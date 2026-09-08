#!/usr/bin/env Rscript
if (!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse")
library(argparse)

# ====================== 参数解析 ======================
parser <- ArgumentParser(
  description = 'RRBS cell proportion estimation (supports EpiDISH + devtEp for PBMC)',
  epilog = paste(
    'Examples:\n',
    '  # CordBlood + uniLIFE7 (hg19)\n',
    '  Rscript script.R -i input_hg19.txt -c cell_ratios.txt -t CordBlood --ref uniLIFE7\n\n',
    '  # CordBlood + FlowSorted (hg38)\n',
    '  Rscript script.R -i input_hg38.txt -c cell_ratios.txt -t CordBlood --ref FlowSorted.CordBlood.450k --hg38ToHg19 /path/to/hg38ToHg19.over.chain\n\n',
    '  # PBMC + cent12CT (default)\n',
    '  Rscript script.R -i input.txt -c cell_ratios.txt -t PBMC\n\n',
    '  # PBMC + devtEp (custom RRBS-native reference)\n',
    '  Rscript script.R -i input.txt -c cell_ratios.txt -t PBMC --ref devtEp --reffile ref.file'
  )
)
parser$add_argument('-i', '--inputf', required = TRUE,
                    help = 'Input methylation ratio file (chrom, pos, samples...)')
parser$add_argument('-c', '--cellratiof', required = TRUE,
                    help = 'Output cell proportion file')
parser$add_argument('--hg38ToHg19', default = NULL,
                    help = 'Path to hg38ToHg19.over.chain. If provided → liftOver; if omitted → assume hg19')
parser$add_argument('-t', '--type', required = TRUE,
                    help = 'Sample type: PBMC or CordBlood')
parser$add_argument('--ref', default = "uniLIFE7",
                    help = paste(
                      'Reference choice:\n',
                      '  PBMC: cent12CT (default) | devtEp\n',
                      '  CordBlood: uniLIFE7 | uniLIFE | FlowSorted.CordBlood.450k\n',
                      'Default: uniLIFE7'
                    ))
parser$add_argument('--reffile', default = NULL,
                    help = 'Custom reference file for --ref devtEp (columns: chr start end celltype.sample...). Required when --ref devtEp.')
args <- parser$parse_args()

inputf         <- args$inputf
cellratiof     <- args$cellratiof
hg38ToHg19file <- args$hg38ToHg19
sample_type    <- args$type
ref_choice     <- args$ref
reffile        <- args$reffile

# ====================== 安装与加载包 ======================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

bioc_pkgs <- c(
  "EpiDISH",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "GenomicRanges",
  "rtracklayer",
  "liftOver",
  "minfi",
  "FlowSorted.CordBlood.450k"
)
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

# 安装 devtEp（仅在需要时，保留包以便兼容其他功能；本脚本在 custom 模式下不再依赖 demo 数据）
if (sample_type == "PBMC" && ref_choice == "devtEp") {
  if (!requireNamespace("devtEp", quietly = TRUE)) {
    cat("正在安装 CORRS-LAB/devtEp ...\n")
    devtools::install_github("CORRS-LAB/devtEp", upgrade = "never")
  }
  library(devtEp)
}

library(EpiDISH)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(GenomicRanges)
library(rtracklayer)
library(minfi)

# ====================== 文件与参数检查 ======================
if (!file.exists(inputf)) stop("Input file does not exist: ", inputf)
do_liftover <- !is.null(hg38ToHg19file)
if (do_liftover && !file.exists(hg38ToHg19file)) {
  stop("Provided --hg38ToHg19 chain file does not exist: ", hg38ToHg19file)
}

use_devtep <- (sample_type == "PBMC" && ref_choice == "devtEp")
if (use_devtep) {
  if (is.null(reffile) || !file.exists(reffile)) {
    stop("When --ref devtEp, --reffile must be provided and exist. Example: --reffile ref.file")
  }
}

cat("Input file      :", inputf, "\n")
cat("Cell ratio file :", cellratiof, "\n")
cat("Sample type     :", sample_type, "\n")
cat("Reference choice:", ref_choice, "\n")
if (use_devtep) cat("Custom reffile  :", reffile, "\n")
cat("Do liftOver     :", do_liftover, ifelse(do_liftover, paste0(" (", hg38ToHg19file, ")"), " (input assumed hg19)"), "\n")

# ====================== 读取数据 ======================
df <- read.table(inputf, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
df$chrom <- as.character(df$chrom)
df$pos   <- as.numeric(as.character(df$pos))

# ====================== 根据参考类型处理数据 ======================
if (use_devtep) {
  # -------------------- custom ref.file 模式（RRBS 原生区域匹配） --------------------
  cat("\n===== 使用自定义 ref.file 构建 RRBS 参考 signature =====\n")

  ref.dat <- read.table(reffile, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(ref.dat) < 4) stop("ref.file 至少需要 4 列: chr start end + 至少一个细胞类型列")

  # 前三列：chr / start / end；染色体加 "chr" 前缀
  chr_raw   <- as.character(ref.dat[[1]])
  chr_raw   <- sub("^chr", "", chr_raw, ignore.case = TRUE)
  chr_vec   <- paste0("chr", chr_raw)
  start_vec <- as.numeric(as.character(ref.dat[[2]]))
  end_vec   <- as.numeric(as.character(ref.dat[[3]]))

  valid <- !is.na(start_vec) & !is.na(end_vec) & (end_vec >= start_vec)
  if (sum(!valid) > 0) {
    cat("警告：过滤掉", sum(!valid), "个无效参考区间\n")
  }
  ref.dat   <- ref.dat[valid, , drop = FALSE]
  chr_vec   <- chr_vec[valid]
  start_vec <- start_vec[valid]
  end_vec   <- end_vec[valid]

  # 提取细胞类型部分并按细胞类型取平均 → 得到区域级参考矩阵
  ref.dt <- ref.dat[, -c(1:3), drop = FALSE]
  celltypes <- colnames(ref.dt) <- as.character(read.table(text = colnames(ref.dt), sep = ".")[, 2])

  ref_mat_region <- matrix(NA, ncol = length(unique(celltypes)), nrow = nrow(ref.dt))
  colnames(ref_mat_region) <- unique(celltypes)
  for (i in colnames(ref_mat_region)) {
    ref_mat_region[, i] <- apply(ref.dt[, celltypes == i, drop = FALSE], 1, mean, na.rm = TRUE)
  }

  cat("devtEp 区域级参考矩阵维度:", paste(dim(ref_mat_region), collapse = " x "), "\n")
  cat("细胞类型:", paste(colnames(ref_mat_region), collapse = ", "), "\n")

  gr_ref <- GRanges(
    seqnames = chr_vec,
    ranges   = IRanges(start = start_vec, end = end_vec)
  )

  gr_sample <- GRanges(
    seqnames = df$chrom,
    ranges   = IRanges(start = df$pos, end = df$pos),
    orig_id  = seq_len(nrow(df))
  )

  if (do_liftover) {
    cat("注意：devtEp 模式建议输入已是与参考相同的基因组版本。当前仍尝试 liftOver...\n")
    chain <- import.chain(hg38ToHg19file)
    gr_sample <- unlist(liftOver(gr_sample, chain))
  }

  overlaps <- findOverlaps(gr_sample, gr_ref)
  if (length(overlaps) == 0) {
    stop("devtEp 参考区域与输入 CpG 没有重叠，请检查基因组版本或坐标。")
  }

  # ---------- 关键修改：按位点展开，而不是按区域平均 ----------
  # queryHits  → 输入位点（df）的索引
  # subjectHits → 参考区域（ref_mat_region）的索引
  site_idx   <- gr_sample$orig_id[queryHits(overlaps)]   # 对应 df 的行号
  region_idx <- subjectHits(overlaps)                    # 对应 ref 区域的行号

  # 1. 构建 beta_mat：每一行是一个重叠的 CpG 位点
  beta_mat <- as.matrix(df[site_idx, 3:ncol(df), drop = FALSE])
  beta_mat <- apply(beta_mat, 2, as.numeric)
  colnames(beta_mat) <- colnames(df)[3:ncol(df)]

  # 2. 构建 ref_mat：把对应区域的参考值复制到每一个重叠位点
  ref_mat <- ref_mat_region[region_idx, , drop = FALSE]

  # 3. 统一行名（epidish 要求两个矩阵行名完全一致）
  common_rownames <- paste0("R", seq_len(nrow(beta_mat)))
  rownames(beta_mat) <- common_rownames
  rownames(ref_mat)  <- common_rownames

  cat("成功匹配到的 CpG 位点数:", nrow(beta_mat), "\n")
  cat("最终 beta_mat / ref_mat 维度:", paste(dim(beta_mat), collapse = " x "), "\n")

} else {
  # -------------------- 原有 450k 匹配流程 --------------------
  anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  gr_450k <- GRanges(
    seqnames = anno450k$chr,
    ranges   = IRanges(start = anno450k$pos, end = anno450k$pos),
    cguid    = rownames(anno450k)
  )

  if (do_liftover) {
    gr_hg38 <- GRanges(
      seqnames = df$chrom,
      ranges   = IRanges(start = pmax(1, df$pos - 1), end = df$pos),
      orig_id  = seq_len(nrow(df))
    )
    chain <- import.chain(hg38ToHg19file)
    gr_hg19 <- unlist(liftOver(gr_hg38, chain))
    cat("liftOver 成功转换位点数:", length(gr_hg19), " / ", nrow(df), "\n")
    overlaps     <- findOverlaps(gr_hg19, gr_450k)
    orig_indices <- gr_hg19$orig_id[queryHits(overlaps)]
    matched_cg   <- gr_450k$cguid[subjectHits(overlaps)]
  } else {
    cat("未提供 --hg38ToHg19，跳过 liftOver，直接使用输入坐标匹配 450k 注释...\n")
    gr_input <- GRanges(
      seqnames = df$chrom,
      ranges   = IRanges(start = df$pos, end = df$pos),
      orig_id  = seq_len(nrow(df))
    )
    overlaps     <- findOverlaps(gr_input, gr_450k)
    orig_indices <- gr_input$orig_id[queryHits(overlaps)]
    matched_cg   <- gr_450k$cguid[subjectHits(overlaps)]
  }

  raw_matrix <- as.matrix(df[orig_indices, 3:ncol(df), drop = FALSE])
  raw_matrix <- apply(raw_matrix, 2, as.numeric)
  rownames(raw_matrix) <- matched_cg

  beta_mat <- aggregate(raw_matrix, by = list(cg = rownames(raw_matrix)), FUN = mean, na.rm = TRUE)
  rownames(beta_mat) <- beta_mat$cg
  beta_mat <- as.matrix(beta_mat[, -1, drop = FALSE])
  cat("成功匹配到的 450k 位点数:", nrow(beta_mat), "\n")

  # -------------------- 选择参考矩阵 --------------------
  if (sample_type == "PBMC") {
    if (ref_choice == "cent12CT") {
      data(cent12CT.m)
      ref_mat <- cent12CT.m
      cat("使用参考: cent12CT.m (PBMC/成人血)\n")
    } else {
      stop("--ref 对 PBMC 只支持: cent12CT | devtEp")
    }
  } else if (sample_type == "CordBlood") {
    if (ref_choice == "uniLIFE7") {
      data(centUniLIFE.m)
      cord_cells <- c("B", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")
      available  <- intersect(cord_cells, colnames(centUniLIFE.m))
      if (length(available) < 5) {
        stop("参考矩阵中找不到足够的脐血细胞类型！实际列名: ",
             paste(colnames(centUniLIFE.m), collapse = ", "))
      }
      ref_mat <- centUniLIFE.m[, available, drop = FALSE]
      cat("使用参考: centUniLIFE.m（仅保留脐血 7 种细胞）\n")
      cat("实际使用的细胞:", paste(available, collapse = ", "), "\n")
    } else if (ref_choice == "uniLIFE") {
      data(centUniLIFE.m)
      ref_mat <- centUniLIFE.m
      cat("使用参考: centUniLIFE.m（完整 19 种细胞）\n")
    } else if (ref_choice == "FlowSorted.CordBlood.450k") {
      cat("正在加载 FlowSorted.CordBlood.450k.ModelPars（官方 700 探针）...\n")
      data(FlowSorted.CordBlood.450k.ModelPars, package = "FlowSorted.CordBlood.450k")
      ref_mat <- as.matrix(FlowSorted.CordBlood.450k.ModelPars)

      colnames(ref_mat) <- gsub("CD4T|CD4", "CD4T", colnames(ref_mat), ignore.case = TRUE)
      colnames(ref_mat) <- gsub("CD8T|CD8", "CD8T", colnames(ref_mat), ignore.case = TRUE)
      colnames(ref_mat) <- gsub("Bcell|B cell|B-cell|^B$", "B", colnames(ref_mat), ignore.case = TRUE)
      colnames(ref_mat) <- gsub("Gran|Granulocyte", "Gran", colnames(ref_mat), ignore.case = TRUE)
      colnames(ref_mat) <- gsub("Mono|Monocyte", "Mono", colnames(ref_mat), ignore.case = TRUE)
      colnames(ref_mat) <- gsub("NK|NKcell", "NK", colnames(ref_mat), ignore.case = TRUE)
      colnames(ref_mat) <- gsub("nRBC|NRBC", "nRBC", colnames(ref_mat), ignore.case = TRUE)

      desired <- c("B", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")
      keep <- intersect(desired, colnames(ref_mat))
      if (length(keep) < 5) {
        stop("ModelPars 中可用细胞类型不足。实际列名: ", paste(colnames(ref_mat), collapse = ", "))
      }
      ref_mat <- ref_mat[, keep, drop = FALSE]
      cat("使用参考: FlowSorted.CordBlood.450k.ModelPars\n")
      cat("实际使用的细胞:", paste(colnames(ref_mat), collapse = ", "), "\n")
    } else {
      stop("--ref 对 CordBlood 只支持: uniLIFE7, uniLIFE, FlowSorted.CordBlood.450k")
    }
  } else {
    stop("sample_type 必须是 'PBMC' 或 'CordBlood'")
  }
}

cat("最终参考矩阵维度:", paste(dim(ref_mat), collapse = " x "), "\n")

# ====================== EpiDISH 去卷积 ======================
out.epidish <- epidish(beta.m = beta_mat, ref.m = ref_mat, method = "RPC")
cat("common features number:", nrow(out.epidish$ref), "\n")
cellFrac <- out.epidish$estF
print(head(cellFrac))
cat("\n各细胞平均比例:\n")
print(round(colMeans(cellFrac), 4))

# ====================== 输出 ======================
write.table(
  cbind(sample_name = rownames(cellFrac), cellFrac),
  file = cellratiof,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n完成！\n")
cat("细胞比例文件:", cellratiof, "\n")