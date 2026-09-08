#!/usr/bin/env Rscript
###############################################
# Differential methylation analysis pipeline (limma + DMRcate)
# Supports command-line arguments
# Supports optional cell-type proportion correction
###############################################

suppressPackageStartupMessages({
  library(argparse)
  library(limma)
  library(DMRcate)
  library(GenomicRanges)
  library(dplyr)
})

###############################################
# Command-line argument parsing
###############################################
parser <- ArgumentParser(description = "Differential methylation analysis (limma + DMRcate). example: Rscript script.R -i methylation_matrix.txt -g Case,Case,Control,Control -o result_dir")

parser$add_argument("-i", "--input", required = TRUE,
                    help = "Input beta matrix file (chr POS sample1 sample2 ...), tab-separated")
parser$add_argument("-g", "--groups", required = TRUE,
                    help = "Comma-separated group labels matching sample order, e.g. Case,Case,Control,Control")
parser$add_argument("-o", "--outdir", default = ".",
                    help = "Output directory [default: current directory]")
parser$add_argument("--prefix", default = "result",
                    help = "Output file prefix [default: result]")

# Character(s) representing missing values in the input file
parser$add_argument("--na_string", default = ".",
                    help = "Character(s) representing NA in the input beta matrix. Multiple values can be separated by comma, e.g. '.,NA,NaN' [default: .]")

# DMP filtering thresholds
parser$add_argument("--dmp_fdr", type = "double", default = 0.05,
                    help = "FDR threshold for significant DMPs [default: 0.05]")
parser$add_argument("--dmp_delta", type = "double", default = 0.1,
                    help = "|deltaBeta| threshold for significant DMPs [default: 0.1]")

# DMR-related parameters
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

# NA handling
parser$add_argument("--min_coverage", type = "double", default = 0.8,
                    help = "Minimum proportion of non-NA samples to keep a CpG [default: 0.8]. Set to 0 to disable filtering.")
parser$add_argument("--coverage_by", default = "global",
                    choices = c("global", "group"),
                    help = "Coverage calculation mode for filtering: global (overall non-NA proportion) or group (require every group to meet min_coverage) [default: global]")
parser$add_argument("--impute", default = "global",
                    choices = c("global", "group"),
                    help = "NA imputation method: global (row mean) or group (within-group mean) [default: global]")

# Optional cell-type proportion correction
parser$add_argument("--cell_props", default = NULL,
                    help = "Optional cell type proportion file (rows = sample names matching beta columns, columns = cell types). If provided, used as covariates for DMP/DMR correction.")

args <- parser$parse_args()

# Create output directory
if (!dir.exists(args$outdir)) {
  dir.create(args$outdir, recursive = TRUE)
}

cat("===== Parameters =====\n")
print(args)
cat("======================\n\n")

###############################################
# 1. Read data
###############################################
# Parse na_string (support multiple values separated by comma)
na_strings <- strsplit(args$na_string, ",")[[1]]
na_strings <- trimws(na_strings)   # remove possible spaces

dat <- read.table(args$input,
                  header = TRUE,
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  check.names = FALSE,
                  na.strings = na_strings)

coords <- dat[, 1:2]
colnames(coords) <- c("chr", "pos")
beta  <- as.matrix(dat[, -c(1:2)])
rownames(beta) <- paste(coords$chr, coords$pos, sep = "_")

# Parse group labels
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

# Contrast: group_levels[1] vs group_levels[2]  (level1 - level2)
contrast_name <- paste0(group_levels[1], "_vs_", group_levels[2])

###############################################
# 2. Design matrix (optional cell-type proportion correction)
###############################################
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

if (!is.null(args$cell_props)) {
  cat("\nReading cell proportion file and applying correction...\n")
  cell <- read.table(args$cell_props, header = TRUE, row.names = 1,
                     sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
  # Check sample name matching
  if (!all(colnames(beta) %in% rownames(cell))) {
    stop("Row names of the cell proportion file must exactly match the column names (sample names) of the beta matrix!\n",
         "Missing samples: ", paste(setdiff(colnames(beta), rownames(cell)), collapse = ", "))
  }
  
  # Reorder according to beta column order
  cell <- cell[colnames(beta), , drop = FALSE]
  cell_mat <- as.matrix(cell)
  
  # Check for missing values
  if (any(is.na(cell_mat))) {
    stop("NA values found in the cell proportion file. Please handle them first!")
  }
  
  # Append to design matrix
  design <- cbind(design, cell_mat)
  
  cat("Cell proportion covariates added, total", ncol(cell_mat), "cell types:\n")
  print(colnames(cell_mat))
  cat("Design matrix dimensions:", paste(dim(design), collapse = " x "), "\n")
} else {
  cat("\nNo cell proportion file provided; cell composition correction is not performed.\n")
}

###############################################
# 3. NA handling
###############################################
cat("\nOriginal number of sites:", nrow(beta), "\n")
cat("Missingness summary:\n")
print(summary(rowMeans(is.na(beta))))

# Filter high-missingness sites
if (args$min_coverage > 0) {
  if (args$coverage_by == "global") {
    # Overall non-missing proportion
    keep <- rowMeans(!is.na(beta)) >= args$min_coverage
  } else {
    # By group: every group must meet min_coverage
    keep <- rep(TRUE, nrow(beta))
    for (g in levels(group)) {
      idx <- which(group == g)
      cov_g <- rowMeans(!is.na(beta[, idx, drop = FALSE]))
      keep <- keep & (cov_g >= args$min_coverage)
    }
  }
  beta   <- beta[keep, , drop = FALSE]
  coords <- coords[keep, ]
  cat("Sites remaining after filtering (coverage >=", args$min_coverage,
      ", by", args$coverage_by, "):", nrow(beta), "\n")
}

# Impute NAs
if (any(is.na(beta))) {
  if (args$impute == "global") {
    beta[is.na(beta)] <- rowMeans(beta, na.rm = TRUE)[row(beta)[is.na(beta)]]
  } else {
    # Within-group mean imputation
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

cat("Any NA remaining after imputation:", any(is.na(beta)), "\n")

# Final safety check before M-value conversion
if (!is.numeric(beta)) {
  stop("beta is still non-numeric before M-value conversion. Please check the input file and --na_string parameter!")
}
storage.mode(beta) <- "numeric"

###############################################
# 4. Convert beta to M-values
###############################################
beta[beta <= 0] <- 1e-6
beta[beta >= 1] <- 1 - 1e-6
M <- log2(beta / (1 - beta))

###############################################
# 5. Differentially methylated positions (DMPs) — limma
###############################################
fit <- lmFit(M, design)

# Contrast: group_levels[1] - group_levels[2]
contrast.matrix <- makeContrasts(
  contrasts = paste0(group_levels[1], " - ", group_levels[2]),
  levels = design
)
colnames(contrast.matrix) <- contrast_name

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract full results
DMP <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")

# Add coordinates
DMP$chr <- coords$chr[match(rownames(DMP), rownames(beta))]
DMP$pos <- coords$pos[match(rownames(DMP), rownames(beta))]

# Calculate deltaBeta (level1 - level2; still based on original beta values)
mean_g1 <- rowMeans(beta[, group == group_levels[1], drop = FALSE])
mean_g2 <- rowMeans(beta[, group == group_levels[2], drop = FALSE])
DMP$deltaBeta <- mean_g1[match(rownames(DMP), names(mean_g1))] -
                 mean_g2[match(rownames(DMP), names(mean_g2))]

# Filter significant DMPs
sig_DMP <- DMP[DMP$adj.P.Val < args$dmp_fdr & abs(DMP$deltaBeta) > args$dmp_delta, ]

# Save results
write.table(DMP,
            file.path(args$outdir, paste0(args$prefix, "_DMP_all_results.txt")),
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(sig_DMP,
            file.path(args$outdir, paste0(args$prefix, "_DMP_significant.txt")),
            sep = "\t", quote = FALSE, row.names = TRUE)

cat("\nNumber of significant DMPs (FDR <", args$dmp_fdr, "& |Δβ| >", args$dmp_delta, "):",
    nrow(sig_DMP), "\n")

###############################################
# 6. Differentially methylated regions (DMRs)
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

###############################################
# Offline extraction of DMRs (no ExperimentHub dependency)
###############################################
coords_df <- DMRcate:::extractCoords(dmrcoutput@coord)

DMR_df <- data.frame(
  seqnames          = coords_df$chrom,
  start             = as.integer(coords_df$chromStart),
  end               = as.integer(coords_df$chromEnd),
  width             = as.integer(coords_df$chromEnd) - as.integer(coords_df$chromStart) + 1,
  no.cpgs           = dmrcoutput@no.cpgs,
  min_smoothed_fdr  = dmrcoutput@min_smoothed_fdr,
  Stouffer          = dmrcoutput@Stouffer,
  HMFDR             = dmrcoutput@HMFDR,
  Fisher            = dmrcoutput@Fisher,
  maxdiff           = dmrcoutput@maxdiff,
  meandiff          = dmrcoutput@meandiff,
  stringsAsFactors  = FALSE
)

# Compatibility with different column name versions
colnames(DMR_df) <- gsub("min_smoothed_fdr", "minfdr", colnames(DMR_df))
colnames(DMR_df) <- gsub("maxdiff", "maxbetafc", colnames(DMR_df))
colnames(DMR_df) <- gsub("meandiff", "meanbetafc", colnames(DMR_df))

# Save all DMRs
write.table(DMR_df,
            file.path(args$outdir, paste0(args$prefix, "_DMR_all_results.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Total number of detected DMRs:", nrow(DMR_df), "\n")

# Filter significant DMRs
fdr_col <- if ("minfdr" %in% colnames(DMR_df)) "minfdr" else "Stouffer"
meandiff_col <- if ("meanbetafc" %in% colnames(DMR_df)) "meanbetafc" else "maxbetafc"

sig_DMR <- DMR_df[
  DMR_df[[fdr_col]] < args$dmr_fdr &
  abs(DMR_df[[meandiff_col]]) > args$dmr_meandiff,
]

write.table(sig_DMR,
            file.path(args$outdir, paste0(args$prefix, "_DMR_significant.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Number of significant DMRs (", fdr_col, "<", args$dmr_fdr,
    "& |", meandiff_col, "| >", args$dmr_meandiff, "):",
    nrow(sig_DMR), "\n")

cat("\nAnalysis finished! Results saved to:", args$outdir, "\n")