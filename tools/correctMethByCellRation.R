#!/usr/bin/env Rscript

if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}

# 加载 argparse
library(argparse)

# 创建解析器
parser <- ArgumentParser(description = 'Correct methylation rates by cell ratio')

parser$add_argument('-m', '--methylationfile', required = TRUE, 
                    help = 'path to methylation ratio file (tab-separated)')
parser$add_argument('-c', '--cellratiofile', required = TRUE, 
                    help = 'path to the corrected methylation ratio file (tab-separated)')
parser$add_argument('-o', '--outputfile', required = TRUE, 
                    help = 'path to the corrected methylation rate file (tab-separated)')

args <- parser$parse_args()

# 获取参数
methylationfile <- args$methylationfile
cellratiofile   <- args$cellratiofile
outputfile      <- args$outputfile


# 2. 检查文件是否存在
if (!file.exists(methylationfile)) {
  stop("Input file does not exist!")
}
if (!file.exists(cellratiofile)) {
  stop("hg38ToHg19.chain file does not exist!")
}

cat("methylation file:",  methylationfile, "\n")
cat("cellratio file:", cellratiofile, "\n")

# ==============================================================================
# 细胞比例矫正脚本 (ChAMP Linear Regression Principle)
# ==============================================================================

# 1. 读取输入文件
# ------------------------------------------------------------------------------
# (1) 读取甲基化率文件
# 请将 "methylation_rates.txt" 替换为你本地的甲基化率文件名
meth_df <- read.table(methylationfile, header = TRUE, sep = "\t", check.names = FALSE)

# (2) 读取细胞比例文件
# 请将 "cell_proportions.txt" 替换为你本地的细胞比例文件名
cell_df <- read.table(cellratiofile, header = TRUE, sep = "\t", check.names = FALSE)


# 2. 数据预处理与样本对齐
# ------------------------------------------------------------------------------
# 提取第一、二列的基因组位置信息 (chrom 和 pos)
loc_info <- meth_df[, 1:2]

# 提取甲基化率数值矩阵 (从第三列开始)
beta_mat <- as.matrix(meth_df[, 3:ncol(meth_df)])

# 整理细胞比例矩阵 (将第一列 bigcsID 设为行名)
cell_mat <- as.matrix(cell_df[, -1])
rownames(cell_mat) <- cell_df[, 1]

# 确保甲基化矩阵与细胞比例矩阵中的样本顺序完全一致
common_samples <- intersect(colnames(beta_mat), rownames(cell_mat))

if (length(common_samples) < length(colnames(beta_mat))) {
  cat("警告：有部分样本在细胞比例文件中未找到，仅对交集样本进行校正。\n")
}

beta_mat <- beta_mat[, common_samples]
cell_mat <- cell_mat[common_samples, ]


# 3. 线性回归校正函数 (依据 ChAMP champ.refbase 原理)
# ------------------------------------------------------------------------------
correct_cell_proportion <- function(beta, cellFrac) {
  # 剔除平均比例最低的一种细胞类型，防止多元线性回归时的共线性 (Collinearity) 报错
  min_col <- which.min(colMeans(cellFrac, na.rm = TRUE))
  predictors <- cellFrac[, -min_col]
  
  cat("拟合线性模型中... (剔除的最小比例细胞类为:", colnames(cellFrac)[min_col], ")\n")
  
  # 对每一个 CpG 位点拟合线性模型: t(beta) ~ predictors
  lm.o <- lm(t(beta) ~ predictors)
  
  # 提取残差 (Residuals) 并加回该位点在所有样本中的原始均值
  tmp.m <- t(lm.o$residuals) + rowMeans(beta, na.rm = TRUE)
  
  # 截断极值，确保校正后的甲基化 Beta 值严格限定在 [0, 1] 物理区间内
  tmp.m[tmp.m <= 0] <- min(tmp.m[tmp.m > 0], na.rm = TRUE)
  tmp.m[tmp.m >= 1] <- max(tmp.m[tmp.m < 1], na.rm = TRUE)
  
  return(tmp.m)
}


# 4. 执行校正与导出
# ------------------------------------------------------------------------------
cat("开始执行细胞比例矫正...\n")
corrected_beta <- correct_cell_proportion(beta_mat, cell_mat)

# 拼回原始的 chrom 和 pos 两列
final_df <- cbind(loc_info, corrected_beta)

# 保存矫正后的甲基化率文件
write.table(final_df, file = outputfile, 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("矫正完成！结果已成功导出至 ", outputfile, " 文件中。\n")