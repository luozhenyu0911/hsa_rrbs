#!/usr/bin/env Rscript
###############################################
# 差异甲基化分析完整流程（limma + DMRcate）
# 支持命令行参数
###############################################

suppressPackageStartupMessages({
  library(argparse)
  library(limma)
  library(DMRcate)
  library(GenomicRanges)
  library(dplyr)
})

###############################################
# 命令行参数解析
###############################################
parser <- ArgumentParser(description = "Differential methylation analysis (limma + DMRcate)")

parser$add_argument("-i", "--input", required = TRUE,
                    help = "Input beta matrix file (chr POS sample1 sample2 ...), tab-separated")
parser$add_argument("-g", "--groups", required = TRUE,
                    help = "Comma-separated group labels matching sample order, e.g. Control,Control,Case,Case")
parser$add_argument("-o", "--outdir", default = ".",
                    help = "Output directory [default: current directory]")
parser$add_argument("--prefix", default = "result",
                    help = "Output file prefix [default: result]")

# DMP 筛选阈值
parser$add_argument("--dmp_fdr", type = "double", default = 0.05,
                    help = "FDR threshold for significant DMPs [default: 0.05]")
parser$add_argument("--dmp_delta", type = "double", default = 0.1,
                    help = "|deltaBeta| threshold for significant DMPs [default: 0.1]")

# DMR 相关参数
parser$add_argument("--lambda", type = "integer", default = 1000,
                    help = "Bandwidth for dmrcate (lambda) [default: 1000]")
parser$add_argument("--C", type = "double", default = 2,
                    help = "Scaling factor for dmrcate (C) [default: 2]")
parser$add_argument("--pcutoff", type = "double", default = 0.05,
                    help = "Individual site FDR cutoff passed to dmrcate [default: 0.05]")
parser$add_argument("--dmr_fdr", type = "double", default = 0.05,
                    help = "minfdr / Stouffer threshold for significant DMRs [default: 0.05]")
parser$add_argument("--dmr_meandiff", type = "double", default = 0.1,
                    help = "|meanbetafc| threshold for significant DMRs [default: 0.1]")

# 其他
parser$add_argument("--min_coverage", type = "double", default = 0.8,
                    help = "Minimum proportion of non-NA samples to keep a CpG [default: 0.8]. Set to 0 to disable filtering.")
parser$add_argument("--coverage_by", default = "global",
                    choices = c("global", "group"),
                    help = "Coverage calculation mode for filtering: global (overall non-NA proportion) or group (require every group to meet min_coverage) [default: global]")
parser$add_argument("--impute", default = "global",
                    choices = c("global", "group"),
                    help = "NA imputation method: global (row mean) or group (within-group mean) [default: global]")

args <- parser$parse_args()

# 创建输出目录
if (!dir.exists(args$outdir)) {
  dir.create(args$outdir, recursive = TRUE)
}

cat("===== Parameters =====\n")
print(args)
cat("======================\n\n")

###############################################
# 1. 读取数据
###############################################
dat <- read.table(args$input, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)

coords <- dat[, 1:2]
colnames(coords) <- c("chr", "pos")
beta  <- as.matrix(dat[, -c(1:2)])
rownames(beta) <- paste(coords$chr, coords$pos, sep = "_")

# 解析分组
group_labels <- strsplit(args$groups, ",")[[1]]
if (length(group_labels) != ncol(beta)) {
  stop("Number of group labels (", length(group_labels),
       ") does not match number of samples (", ncol(beta), ")")
}
group <- factor(group_labels)
cat("Group counts:\n")
print(table(group))

if (nlevels(group) != 2) {
  stop("Currently only supports two-group comparison.")
}
group_levels <- levels(group)
contrast_name <- paste0(group_levels[2], "_vs_", group_levels[1])

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

###############################################
# 2. NA 处理
###############################################
cat("\n原始位点数:", nrow(beta), "\n")
cat("缺失比例概况:\n")
print(summary(rowMeans(is.na(beta))))

# 过滤高缺失位点
if (args$min_coverage > 0) {
  if (args$coverage_by == "global") {
    # 整体非缺失比例
    keep <- rowMeans(!is.na(beta)) >= args$min_coverage
  } else {
    # 组间：要求每一个组都满足 min_coverage
    keep <- rep(TRUE, nrow(beta))
    for (g in levels(group)) {
      idx <- which(group == g)
      cov_g <- rowMeans(!is.na(beta[, idx, drop = FALSE]))
      keep <- keep & (cov_g >= args$min_coverage)
    }
  }
  beta   <- beta[keep, , drop = FALSE]
  coords <- coords[keep, ]
  cat("过滤后剩余位点数 (coverage >=", args$min_coverage,
      ", by", args$coverage_by, "):", nrow(beta), "\n")
}

# 填充 NA
if (any(is.na(beta))) {
  if (args$impute == "global") {
    beta[is.na(beta)] <- rowMeans(beta, na.rm = TRUE)[row(beta)[is.na(beta)]]
  } else {
    # 按组内均值填充
    for (g in levels(group)) {
      idx <- which(group == g)
      mat <- beta[, idx, drop = FALSE]
      na_idx <- which(is.na(mat), arr.ind = TRUE)
      if (nrow(na_idx) > 0) {
        mat[na_idx] <- rowMeans(mat, na.rm = TRUE)[na_idx[, 1]]
        beta[, idx] <- mat
      }
    }
  }
}

cat("填充后是否还有 NA:", any(is.na(beta)), "\n")

###############################################
# 3. Beta 转 M 值
###############################################
beta[beta <= 0] <- 1e-6
beta[beta >= 1] <- 1 - 1e-6
M <- log2(beta / (1 - beta))

###############################################
# 4. 差异甲基化位点分析（DMPs）—— limma
###############################################
fit <- lmFit(M, design)

contrast.matrix <- makeContrasts(
  contrasts = paste0(group_levels[2], " - ", group_levels[1]),
  levels = design
)
colnames(contrast.matrix) <- contrast_name

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 提取全部结果
DMP <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")

# 添加坐标
DMP$chr <- coords$chr[match(rownames(DMP), rownames(beta))]
DMP$pos <- coords$pos[match(rownames(DMP), rownames(beta))]

# 计算 deltaBeta
mean_g1 <- rowMeans(beta[, group == group_levels[1], drop = FALSE])
mean_g2 <- rowMeans(beta[, group == group_levels[2], drop = FALSE])
DMP$deltaBeta <- mean_g2[match(rownames(DMP), names(mean_g2))] -
                 mean_g1[match(rownames(DMP), names(mean_g1))]

# 筛选显著 DMP
sig_DMP <- DMP[DMP$adj.P.Val < args$dmp_fdr & abs(DMP$deltaBeta) > args$dmp_delta, ]

# 保存
write.table(DMP,
            file.path(args$outdir, paste0(args$prefix, "_DMP_all_results.txt")),
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(sig_DMP,
            file.path(args$outdir, paste0(args$prefix, "_DMP_significant.txt")),
            sep = "\t", quote = FALSE, row.names = TRUE)

cat("\n显著 DMP 数量 (FDR <", args$dmp_fdr, "& |Δβ| >", args$dmp_delta, "):",
    nrow(sig_DMP), "\n")

###############################################
# 5. 差异甲基化区域分析（DMRs）
###############################################
tt <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
stopifnot(all(rownames(tt) == rownames(beta)))

gr <- GRanges(
  seqnames = coords$chr,
  ranges   = IRanges(start = as.integer(coords$pos), width = 1)
)

mcols(gr) <- DataFrame(
  stat     = tt$t,
  diff     = DMP$deltaBeta[match(rownames(tt), rownames(DMP))],
  rawpval  = tt$P.Value,
  ind.fdr  = tt$adj.P.Val,
  is.sig   = tt$adj.P.Val < args$pcutoff
)
names(gr) <- rownames(tt)

gr <- sort(sortSeqlevels(gr))

myannotation <- new("CpGannotated", ranges = gr)

dmrcoutput <- dmrcate(myannotation,
                      lambda  = args$lambda,
                      C       = args$C,
                      pcutoff = args$pcutoff)

# results.ranges <- extractRanges(dmrcoutput)
# DMR_df <- as.data.frame(results.ranges)

###############################################
# 离线提取 DMR（不依赖 ExperimentHub）
###############################################
# 从 DMResults 对象中直接提取坐标和统计量
coords_df <- DMRcate:::extractCoords(dmrcoutput@coord)

DMR_df <- data.frame(
  seqnames          = coords_df$chrom,
  start             = as.integer(coords_df$chromStart),
  end               = as.integer(coords_df$chromEnd),
  width             = as.integer(coords_df$chromEnd) - as.integer(coords_df$chromStart) + 1,
  no.cpgs           = dmrcoutput@no.cpgs,
  min_smoothed_fdr  = dmrcoutput@min_smoothed_fdr,   # 或 minfdr
  Stouffer          = dmrcoutput@Stouffer,
  HMFDR             = dmrcoutput@HMFDR,
  Fisher            = dmrcoutput@Fisher,
  maxdiff           = dmrcoutput@maxdiff,            # 或 maxbetafc
  meandiff          = dmrcoutput@meandiff,           # 或 meanbetafc
  stringsAsFactors  = FALSE
)

# 兼容不同版本的列名
colnames(DMR_df) <- gsub("min_smoothed_fdr", "minfdr", colnames(DMR_df))
colnames(DMR_df) <- gsub("maxdiff", "maxbetafc", colnames(DMR_df))
colnames(DMR_df) <- gsub("meandiff", "meanbetafc", colnames(DMR_df))

# 保存全部 DMR
write.table(DMR_df,
            file.path(args$outdir, paste0(args$prefix, "_DMR_all_results.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("检测到的 DMR 总数:", nrow(DMR_df), "\n")

# 筛选显著 DMR
fdr_col <- if ("minfdr" %in% colnames(DMR_df)) "minfdr" else "Stouffer"
meandiff_col <- if ("meanbetafc" %in% colnames(DMR_df)) "meanbetafc" else "maxbetafc"

sig_DMR <- DMR_df[
  DMR_df[[fdr_col]] < args$dmr_fdr &
  abs(DMR_df[[meandiff_col]]) > args$dmr_meandiff,
]

write.table(sig_DMR,
            file.path(args$outdir, paste0(args$prefix, "_DMR_significant.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("显著 DMR 数量 (", fdr_col, "<", args$dmr_fdr,
    "& |", meandiff_col, "| >", args$dmr_meandiff, "):",
    nrow(sig_DMR), "\n")

cat("\n分析完成！结果已保存到:", args$outdir, "\n")
