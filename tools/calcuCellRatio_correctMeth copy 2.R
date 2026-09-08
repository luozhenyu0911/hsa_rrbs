#!/usr/bin/env Rscript
if (!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse")
library(argparse)
# ====================== 参数解析 ======================
parser <- ArgumentParser(description = 'RRBS cell proportion estimation + correction (optional hg38 → hg19)\n',
                      epilog = 'This script estimates cell proportions from RRBS methylation data and optionally corrects for differences between hg38 and hg19 coordinates.\n',
                      usage = 'Rscript script.R -i input_hg19.txt -c cell_ratios.txt -t CordBlood --ref uniLIFE7\n, 
                      or \n    Rscript script.R -i input_hg38.txt -c cell_ratios.txt -t CordBlood --ref FlowSorted.CordBlood --hg38ToHg19 /path/to/hg38ToHg19.over.chain')
parser$add_argument('-i', '--inputf', required = TRUE,
                    help = 'Input methylation ratio file (chrom, pos, samples...)')
# parser$add_argument('-o', '--outputf', required = TRUE,
#                     help = 'Output corrected methylation ratio file (currently only cell ratios are written)')
parser$add_argument('-c', '--cellratiof', required = TRUE,
                    help = 'Output cell proportion file')
parser$add_argument('--hg38ToHg19', default = NULL,
                    help = 'Path to hg38ToHg19.over.chain file. If provided, perform liftOver (hg38→hg19); if omitted, assume input is already hg19.')
parser$add_argument('-t', '--type', required = TRUE,
                    help = 'Sample type: PBMC or CordBlood')
parser$add_argument('--ref', default = "uniLIFE7",
                    help = paste('[PBMC, uniLIFE7, uniLIFE, FlowSorted.CordBlood.450k]',
                      'Reference matrix choice:\n',
                      'PBMC: only cent12CT.m is used.\n',
                      '(for CordBlood) uniLIFE7: (7 cord cells from centUniLIFE.m);',
                      'uniLIFE (full 19 cells from centUniLIFE.m);',
                      'or FlowSorted.CordBlood.450k: (Bakulski 7 cord cells).\n',
                      'Default: uniLIFE7'
                    ))
args <- parser$parse_args()

inputf         <- args$inputf
# outputf        <- args$outputf
cellratiof     <- args$cellratiof
hg38ToHg19file <- args$hg38ToHg19
sample_type    <- args$type
ref_choice     <- args$ref

# ====================== 安装与加载包 ======================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
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
library(EpiDISH)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(GenomicRanges)
library(rtracklayer)
library(minfi)

# ====================== 文件与参数检查 ======================
if (!file.exists(inputf)) stop("Input file does not exist: ", inputf)

do_liftover <- !is.null(hg38ToHg19file)
if (do_liftover) {
  if (!file.exists(hg38ToHg19file)) {
    stop("Provided --hg38ToHg19 chain file does not exist: ", hg38ToHg19file)
  }
}

cat("Input file      :", inputf, "\n")
# cat("Output file     :", outputf, "\n")
cat("Cell ratio file :", cellratiof, "\n")
cat("Sample type     :", sample_type, "\n")
cat("Reference choice:", ref_choice, "\n")
cat("Do liftOver     :", do_liftover, ifelse(do_liftover, paste0(" (", hg38ToHg19file, ")"), " (input assumed hg19)"), "\n")

# ====================== 读取数据 ======================
df <- read.table(inputf, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
df$chrom <- as.character(df$chrom)
df$pos   <- as.numeric(as.character(df$pos))

# 450k 注释（hg19）
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
gr_450k <- GRanges(
  seqnames = anno450k$chr,
  ranges   = IRanges(start = anno450k$pos, end = anno450k$pos),
  cguid    = rownames(anno450k)
)

if (do_liftover) {
  # ====================== liftOver (hg38 → hg19) ======================
  # 很多 RRBS 是 1-based，用 pos-1 更稳妥
  gr_hg38 <- GRanges(
    seqnames = df$chrom,
    ranges   = IRanges(start = pmax(1, df$pos - 1), end = df$pos),
    orig_id  = seq_len(nrow(df))
  )
  chain <- import.chain(hg38ToHg19file)
  gr_hg19_list <- liftOver(gr_hg38, chain)
  gr_hg19 <- unlist(gr_hg19_list)
  cat("liftOver 成功转换位点数:", length(gr_hg19), " / ", nrow(df), "\n")

  overlaps     <- findOverlaps(gr_hg19, gr_450k)
  orig_indices <- gr_hg19$orig_id[queryHits(overlaps)]
  matched_cg   <- gr_450k$cguid[subjectHits(overlaps)]
} else {
  # 直接按坐标匹配（输入已是 hg19）
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

# 提取甲基化矩阵
raw_matrix <- as.matrix(df[orig_indices, 3:ncol(df)])
raw_matrix <- apply(raw_matrix, 2, as.numeric)
rownames(raw_matrix) <- matched_cg

# 去重取平均
beta_mat <- aggregate(raw_matrix, by = list(cg = rownames(raw_matrix)), FUN = mean, na.rm = TRUE)
rownames(beta_mat) <- beta_mat$cg
beta_mat <- as.matrix(beta_mat[, -1])
cat("成功匹配到的 450k 位点数:", nrow(beta_mat), "\n")

# ====================== 选择 / 构建参考矩阵 ======================
if (sample_type == "PBMC") {
  data(cent12CT.m)
  ref_mat <- cent12CT.m
  cat("使用参考: cent12CT.m (PBMC/成人血)\n")

} else if (sample_type == "CordBlood") {

  if (ref_choice == "uniLIFE7") {
    data(centUniLIFE.m)
    cord_cells <- c("B", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")
    available  <- intersect(cord_cells, colnames(centUniLIFE.m))
    if (length(available) < 5) {
      stop("参考矩阵中找不到足够的脐血细胞类型，请检查列名！\n实际列名: ",
           paste(colnames(centUniLIFE.m), collapse = ", "))
    }
    ref_mat <- centUniLIFE.m[, available, drop = FALSE]
    cat("使用参考: centUniLIFE.m（仅保留脐血 7 种细胞）\n")
    cat("实际使用的细胞:", paste(available, collapse = ", "), "\n")

  } else if (ref_choice == "uniLIFE") {
    data(centUniLIFE.m)
    ref_mat <- centUniLIFE.m
    cat("使用参考: centUniLIFE.m（完整 19 种细胞）\n")

  # } else if (ref_choice == "FlowSorted.CordBlood.450k") {
  #   library(FlowSorted.CordBlood.450k)
  #   # 从 FlowSorted.CordBlood.450k 构建均值参考矩阵
  #   cat("正在从 FlowSorted.CordBlood.450k 构建参考矩阵...\n")
  #   data(FlowSorted.CordBlood.450k)
  #   Mset  <- preprocessRaw(FlowSorted.CordBlood.450k)
  #   beta  <- getBeta(Mset)
  #   cell_types <- pData(FlowSorted.CordBlood.450k)$CellType

  #   # 标准化常见命名
  #   cell_types <- gsub("CD4T", "CD4T", cell_types)
  #   cell_types <- gsub("CD8T", "CD8T", cell_types)
  #   cell_types <- gsub("Bcell|B cell|B-cell", "B", cell_types, ignore.case = TRUE)
  #   cell_types <- gsub("Gran|Granulocyte", "Gran", cell_types, ignore.case = TRUE)
  #   cell_types <- gsub("Mono|Monocyte", "Mono", cell_types, ignore.case = TRUE)
  #   cell_types <- gsub("NK|NKcell", "NK", cell_types, ignore.case = TRUE)
  #   cell_types <- gsub("nRBC|NRBC|nRBC", "nRBC", cell_types, ignore.case = TRUE)

  #   unique_cts <- unique(cell_types)
  #   cat("检测到的细胞类型:", paste(unique_cts, collapse = ", "), "\n")

  #   ref_list <- lapply(unique_cts, function(ct) {
  #     idx <- which(cell_types == ct)
  #     if (length(idx) == 0) return(NULL)
  #     rowMeans(beta[, idx, drop = FALSE], na.rm = TRUE)
  #   })
  #   names(ref_list) <- unique_cts
  #   ref_list <- ref_list[!sapply(ref_list, is.null)]
  #   ref_mat  <- do.call(cbind, ref_list)

  #   # 只保留标准 7 种脐血细胞
  #   desired <- c("B", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")
  #   keep    <- intersect(desired, colnames(ref_mat))
  #   if (length(keep) < 5) {
  #     stop("FlowSorted.CordBlood.450k 中可用脐血细胞类型不足。实际列名: ",
  #          paste(colnames(ref_mat), collapse = ", "))
  #   }
  #   ref_mat <- ref_mat[, keep, drop = FALSE]
  #   cat("使用参考: FlowSorted.CordBlood.450k（构建的 7 种脐血细胞均值矩阵）\n")
  #   cat("实际使用的细胞:", paste(colnames(ref_mat), collapse = ", "), "\n")
  } else if (ref_choice == "FlowSorted.CordBlood.450k") {
    cat("正在加载 FlowSorted.CordBlood.450k.ModelPars（官方 700 探针参考矩阵）...\n")
    # 确保包已加载
    if (!requireNamespace("FlowSorted.CordBlood.450k", quietly = TRUE)) {
      stop("请先安装 FlowSorted.CordBlood.450k 包")
    }
    data(FlowSorted.CordBlood.450k.ModelPars, package = "FlowSorted.CordBlood.450k")
    
    ref_mat <- as.matrix(FlowSorted.CordBlood.450k.ModelPars)
    
    # 查看并标准化列名
    cat("原始列名:", paste(colnames(ref_mat), collapse = ", "), "\n")
    
    # 常见列名映射（根据实际输出可能需要微调）
    colnames(ref_mat) <- gsub("CD4T|CD4", "CD4T", colnames(ref_mat), ignore.case = TRUE)
    colnames(ref_mat) <- gsub("CD8T|CD8", "CD8T", colnames(ref_mat), ignore.case = TRUE)
    colnames(ref_mat) <- gsub("Bcell|B cell|B-cell|B", "B", colnames(ref_mat), ignore.case = TRUE)
    colnames(ref_mat) <- gsub("Gran|Granulocyte", "Gran", colnames(ref_mat), ignore.case = TRUE)
    colnames(ref_mat) <- gsub("Mono|Monocyte", "Mono", colnames(ref_mat), ignore.case = TRUE)
    colnames(ref_mat) <- gsub("NK|NKcell", "NK", colnames(ref_mat), ignore.case = TRUE)
    colnames(ref_mat) <- gsub("nRBC|NRBC|nRBC", "nRBC", colnames(ref_mat), ignore.case = TRUE)
    
    # 只保留标准 7 种细胞
    desired <- c("B", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")
    keep <- intersect(desired, colnames(ref_mat))
  if (length(keep) < 5) {
    stop("ModelPars 中可用细胞类型不足。实际列名: ", paste(colnames(ref_mat), collapse = ", "))
  }
  ref_mat <- ref_mat[, keep, drop = FALSE]
  
  cat("使用参考: FlowSorted.CordBlood.450k.ModelPars（官方 700 探针）\n")
  cat("实际使用的细胞:", paste(colnames(ref_mat), collapse = ", "), "\n")
  cat("参考矩阵维度:", paste(dim(ref_mat), collapse = " x "), "\n")

  } else {
    stop("--ref 对 CordBlood 只支持: uniLIFE7, uniLIFE, FlowSorted.CordBlood.450k")
  }

} else {
  stop("sample_type 必须是 'PBMC' 或 'CordBlood'")
}

cat("参考矩阵维度:", paste(dim(ref_mat), collapse = " x "), "\n")

# ====================== EpiDISH 去卷积 ======================
out.epidish <- epidish(beta.m = beta_mat, ref.m = ref_mat, method = "RPC")
cat("common CpG number:", nrow(out.epidish$ref), "\n")
cellFrac <- out.epidish$estF
print(head(cellFrac))
cat("\n各细胞平均比例:\n")
print(round(colMeans(cellFrac), 4))

# ====================== 输出细胞比例 ======================
write.table(
  cbind(sample_name = rownames(cellFrac), cellFrac),
  file = cellratiof,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n完成！\n")
cat("细胞比例文件:", cellratiof, "\n")
# 注意：原脚本中甲基化校正部分仍被注释，因此 outputf 目前不会写入内容。