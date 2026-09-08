#!/usr/bin/env Rscript

if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}

# 1. 安装必要的包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioc_pkgs <- c(
  "rtracklayer", "EpiDISH",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "liftOver", "FlowSorted.CordBlood.450k"
)
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
}

# 加载 argparse
library(argparse)

# 创建解析器
parser <- ArgumentParser(description = 'calculate cell proportion and correct methylation rates by cell ratio')

parser$add_argument('-i', '--inputf', required = TRUE, 
                    help = 'the methylation ratio file (tab-separated)')
parser$add_argument('-o', '--outputf', required = TRUE, 
                    help = 'the corrected methylation ratio output file')
parser$add_argument('-c', '--cellratiof', required = TRUE, 
                    help = 'the cell ratio output file')
parser$add_argument('--hg38ToHg19', required = TRUE, 
                    help = 'path to hg38ToHg19.over.chain file;wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz')
parser$add_argument('-t', '--type', required = TRUE, 
                    help = 'sample type: PBMC or CordBlood')

args <- parser$parse_args()

# 获取参数
inputf        <- args$inputf
outputf       <- args$outputf
cellratiof    <- args$cellratiof
hg38ToHg19file <- args$hg38ToHg19
sample_type   <- args$type


# 2. 检查文件是否存在
if (!file.exists(inputf)) {
  stop("Input file does not exist!")
}
if (!file.exists(hg38ToHg19file)) {
  stop("hg38ToHg19.chain file does not exist!")
}


cat("Input file:", inputf, "\n")
cat("Output file:", outputf, "\n")
cat("Cell ratio file:", cellratiof, "\n")

library(rtracklayer)
library(EpiDISH)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# 1. 确保读取时数据类型正确
df <- read.table(inputf, header = TRUE, check.names = FALSE)

# 强制转换坐标和数据类型 (防止出现字符串)
df$chrom <- as.character(df$chrom)
df$pos   <- as.numeric(as.character(df$pos))

# 2. 仅使用坐标构建 GRanges，并记录原始行号 ID
gr_hg38 <- GRanges(
  seqnames = df$chrom,
  ranges   = IRanges(start = df$pos-1, end = df$pos),
  orig_id  = 1:nrow(df)  # 用行号记录对应关系
)

# 1. 定义本地链文件名
# # 2. 如果本地不存在，则自动从 UCSC 官网下载并解压
# wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
# gunzip hg38ToHg19.over.chain.gz

# 3. 导入 chain 文件并转换坐标
chain <- import.chain(hg38ToHg19file)
gr_hg19_list <- liftOver(gr_hg38, chain)
gr_hg19 <- unlist(gr_hg19_list)

# 4. 获取 450k 芯片坐标
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
gr_450k <- GRanges(
  seqnames = anno450k$chr,
  ranges   = IRanges(start = anno450k$pos, end = anno450k$pos),
  cguid    = rownames(anno450k)
)

# 5. 求交集并匹配 cg 号
overlaps <- findOverlaps(gr_hg19, gr_450k)

# 提取重合位点在原始 df 中的行号和对应的 cg 号
orig_indices <- gr_hg19$orig_id[queryHits(overlaps)]
matched_cg   <- gr_450k$cguid[subjectHits(overlaps)]

# 6. 提取对应的甲基化矩阵 (转换为纯 numeric 矩阵)
raw_matrix <- as.matrix(df[orig_indices, 3:ncol(df)])
# 确保矩阵内部全部为 numeric 类型
raw_matrix <- apply(raw_matrix, 2, as.numeric)
rownames(raw_matrix) <- matched_cg

# 7. 去重 (如果同一个 cg 位点匹配到了多个坐标，取平均值)
beta_mat <- aggregate(raw_matrix, by = list(rownames(raw_matrix)), FUN = mean, na.rm = TRUE)
rownames(beta_mat) <- beta_mat$Group.1
beta_mat <- as.matrix(beta_mat[, -1])

cat("成功转换并构建 Beta 矩阵，匹配到的 450k 位点数:", nrow(beta_mat), "\n")


library(EpiDISH)
# 设置样本类型："PBMC" = 外周血/成人血, "CordBlood" = 脐带血
if (sample_type == "PBMC") {
  # ---------- 外周血 / 成人血：EpiDISH 内置参考 ----------
  library(EpiDISH)
  data(cent12CT.m)        # EPIC 阵列用 cent12CT.m；若是 450k 阵列改为 data(cent12CT450k.m)
  ref_mat <- cent12CT.m
  # 若是 450k 阵列：data(cent12CT450k.m); ref_mat <- cent12CT450k.m

} else if (sample_type == "CordBlood") {
  # ---------- 脐带血：FlowSorted.CordBlood.450k ----------
  library(FlowSorted.CordBlood.450k)
  library(minfi)
  data("FlowSorted.CordBlood.450k")

  mSet      <- preprocessRaw(FlowSorted.CordBlood.450k)
  all_beta  <- getBeta(mSet)
  cell_types <- pData(FlowSorted.CordBlood.450k)$CellType

  valid_idx  <- !is.na(cell_types) & cell_types != "WholeBlood"
  sub_beta   <- all_beta[, valid_idx, drop = FALSE]
  sub_types  <- as.character(cell_types[valid_idx])

  ref_mat <- do.call(cbind, lapply(split(seq_along(sub_types), sub_types), function(idx) {
    rowMeans(sub_beta[, idx, drop = FALSE], na.rm = TRUE)
  }))

} else {
  stop("sample_type 必须是 'PBMC' 或 'CordBlood'")
}

cat("成功加载参考甲基化矩阵，样本类型:", sample_type, "\n")
dim(ref_mat)

# 2.2 运行 EpiDISHDeconvolution (RPC 方法)
out.epidish <- epidish(beta.m = beta_mat, ref.m = ref_mat, method = 'RPC')

cat("common CpG number:", nrow(out.epidish$ref), "\n")

# 提取估计出的细胞比例矩阵
cellFrac <- out.epidish$estF
print(head(cellFrac))


# 3.1 提取原始全部位点的 Beta 矩阵 (包括无法映射到 450k 的所有 RRBS 位点)
raw_beta <- as.matrix(df[, 3:ncol(df)])

# 确保样本顺序完全一致
cellFrac <- cellFrac[colnames(raw_beta), ]

# 3.2 自定义细胞比例线性校正函数 (基于 ChAMP champ.refbase 原理)
correct_cell_proportion <- function(beta, cellFrac) {
  # 剔除平均比例最低的一种细胞类型，防止多重共线性 (Collinearity)
  min_col <- which.min(colMeans(cellFrac))
  predictors <- cellFrac[, -min_col]
  
  # 对每个 CpG 位点拟合线性模型: t(beta) ~ cellFrac
  lm.o <- lm(t(beta) ~ predictors)
  
  # 提取残差 (Residuals) 并加回全组位点的 Mean Beta 值
  tmp.m <- t(lm.o$residuals) + rowMeans(beta, na.rm = TRUE)
  
  # 截断极值，确保甲基化 Beta 值在 [0, 1] 范围内
  tmp.m[tmp.m <= 0] <- min(tmp.m[tmp.m > 0], na.rm = TRUE)
  tmp.m[tmp.m >= 1] <- max(tmp.m[tmp.m < 1], na.rm = TRUE)
  
  return(tmp.m)
}

# 3.3 执行校正
corrected_beta <- correct_cell_proportion(raw_beta, cellFrac)

# 3.4 拼回原有的 chrom 和 pos 信息并导出结果
final_df <- cbind(df[, 1:2], corrected_beta)

# 导出校正后的甲基化率文件
write.table(final_df, outputf, sep = "\t", quote = FALSE, row.names = FALSE)
# 导出估计的细胞比例结果
write.table(
  cbind(sample_name = rownames(cellFrac), cellFrac),
  file = cellratiof,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("细胞比例估计与甲基化率校正完成！\n")
