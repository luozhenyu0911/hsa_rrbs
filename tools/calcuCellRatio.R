#!/usr/bin/env Rscript

# ==============================================================================
# RRBS Cell Proportion Estimation
# Supports: EpiDISH (450k-based) + custom RRBS-native reference (devtEp)
# ==============================================================================

# -------------------------- 依赖安装与加载 ------------------------------------
if (!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse")
library(argparse)

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

# -------------------------- 参数解析 ------------------------------------------
parser <- ArgumentParser(
  description = "RRBS cell proportion estimation (EpiDISH + optional RRBS-native reference)",
  epilog = paste(
    "Examples:\n",
    "  # CordBlood + uniLIFE7 (hg19)\n",
    "  Rscript script.R -i input_hg19.txt -c cell_ratios.txt -t CordBlood --ref uniLIFE7\n\n",
    "  # CordBlood + FlowSorted.CordBlood.450k.ModelPars (hg38)\n",
    "  Rscript script.R -i input_hg38.txt -c cell_ratios.txt -t CordBlood \\\n",
    "    --ref FlowSorted.CordBlood.450k.ModelPars --hg38ToHg19 /path/to/hg38ToHg19.over.chain\n\n",
    "  # CordBlood + FlowSorted.CordBlood.450k (full mean matrix)\n",
    "  Rscript script.R -i input.txt -c cell_ratios.txt -t CordBlood --ref FlowSorted.CordBlood.450k\n\n",
    "  # PBMC + cent12CT (default)\n",
    "  Rscript script.R -i input.txt -c cell_ratios.txt -t PBMC\n\n",
    "  # PBMC + devtEp (custom RRBS-native reference)\n",
    "  Rscript script.R -i input.txt -c cell_ratios.txt -t PBMC --ref devtEp --reffile ref.file"
  )
)

parser$add_argument("-i", "--inputf", required = TRUE,
                    help = "Input methylation ratio file (chrom, pos, samples...)")
parser$add_argument("-c", "--cellratiof", required = TRUE,
                    help = "Output cell proportion file")
parser$add_argument("--hg38ToHg19", default = NULL,
                    help = "Path to hg38ToHg19.over.chain. If provided → liftOver; if omitted → assume hg19")
parser$add_argument("-t", "--type", required = TRUE,
                    help = "Sample type: PBMC or CordBlood")
parser$add_argument("--ref", default = "uniLIFE7",
                    help = paste(
                      "Reference choice:\n",
                      "  PBMC     : cent12CT (default) | devtEp\n",
                      "  CordBlood: uniLIFE7 | uniLIFE | FlowSorted.CordBlood.450k.ModelPars | FlowSorted.CordBlood.450k\n",
                      "Default: uniLIFE7"
                    ))
parser$add_argument("--reffile", default = NULL,
                    help = "Custom reference file for --ref devtEp (chr start end celltype.sample...). Required when --ref devtEp.")

args <- parser$parse_args()

inputf         <- args$inputf
cellratiof     <- args$cellratiof
hg38ToHg19file <- args$hg38ToHg19
sample_type    <- args$type
ref_choice     <- args$ref
reffile        <- args$reffile

# -------------------------- 条件安装 devtEp -----------------------------------
use_devtep <- (sample_type == "PBMC" && ref_choice == "devtEp")
if (use_devtep) {
  if (!requireNamespace("devtEp", quietly = TRUE)) {
    cat("Installing CORRS-LAB/devtEp ...\n")
    devtools::install_github("CORRS-LAB/devtEp", upgrade = "never")
  }
  library(devtEp)
}

# -------------------------- 加载核心包 ----------------------------------------
library(EpiDISH)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(GenomicRanges)
library(rtracklayer)
library(minfi)

# -------------------------- 参数与文件检查 ------------------------------------
if (!file.exists(inputf)) stop("Input file does not exist: ", inputf)

do_liftover <- !is.null(hg38ToHg19file)
if (do_liftover && !file.exists(hg38ToHg19file)) {
  stop("Provided --hg38ToHg19 chain file does not exist: ", hg38ToHg19file)
}

if (use_devtep) {
  if (is.null(reffile) || !file.exists(reffile)) {
    stop("When --ref devtEp, --reffile must be provided and exist.")
  }
}

cat("Input file      :", inputf, "\n")
cat("Cell ratio file :", cellratiof, "\n")
cat("Sample type     :", sample_type, "\n")
cat("Reference choice:", ref_choice, "\n")
if (use_devtep) cat("Custom reffile  :", reffile, "\n")
cat("Do liftOver     :", do_liftover,
    ifelse(do_liftover, paste0(" (", hg38ToHg19file, ")"), " (input assumed hg19)"), "\n")

# -------------------------- 读取输入数据 --------------------------------------
df <- read.table(inputf, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
df$chrom <- as.character(df$chrom)
df$pos   <- as.numeric(as.character(df$pos))

# ==============================================================================
# 主处理逻辑
# ==============================================================================
if (use_devtep) {
  # --------------------------------------------------------------------------
  # Mode 1: custom RRBS-native reference (devtEp)
  # --------------------------------------------------------------------------
  cat("\n===== Using custom ref.file (RRBS-native region matching) =====\n")

  ref.dat <- read.table(reffile, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(ref.dat) < 4) stop("ref.file must have at least 4 columns: chr start end + ≥1 cell-type column")

  # Standardize chromosome names
  chr_raw   <- sub("^chr", "", as.character(ref.dat[[1]]), ignore.case = TRUE)
  chr_vec   <- paste0("chr", chr_raw)
  start_vec <- as.numeric(as.character(ref.dat[[2]]))
  end_vec   <- as.numeric(as.character(ref.dat[[3]]))

  valid <- !is.na(start_vec) & !is.na(end_vec) & (end_vec >= start_vec)
  if (sum(!valid) > 0) cat("Warning: filtered", sum(!valid), "invalid reference intervals\n")

  ref.dat   <- ref.dat[valid, , drop = FALSE]
  chr_vec   <- chr_vec[valid]
  start_vec <- start_vec[valid]
  end_vec   <- end_vec[valid]

  # Build region-level reference matrix (average replicates per cell type)
  ref.dt <- ref.dat[, -c(1:3), drop = FALSE]
  celltypes <- colnames(ref.dt) <- as.character(read.table(text = colnames(ref.dt), sep = ".")[, 2])

  ref_mat_region <- matrix(NA, nrow = nrow(ref.dt), ncol = length(unique(celltypes)))
  colnames(ref_mat_region) <- unique(celltypes)
  for (ct in colnames(ref_mat_region)) {
    ref_mat_region[, ct] <- rowMeans(ref.dt[, celltypes == ct, drop = FALSE], na.rm = TRUE)
  }

  cat("Region-level reference matrix:", paste(dim(ref_mat_region), collapse = " x "), "\n")
  cat("Cell types:", paste(colnames(ref_mat_region), collapse = ", "), "\n")

  # GRanges objects
  gr_ref <- GRanges(seqnames = chr_vec,
                    ranges   = IRanges(start = start_vec, end = end_vec))

  gr_sample <- GRanges(seqnames = df$chrom,
                       ranges   = IRanges(start = df$pos, end = df$pos),
                       orig_id  = seq_len(nrow(df)))

  if (do_liftover) {
    cat("Note: for devtEp it is recommended that input and reference share the same genome build.\n")
    chain     <- import.chain(hg38ToHg19file)
    gr_sample <- unlist(liftOver(gr_sample, chain))
  }

  overlaps <- findOverlaps(gr_sample, gr_ref)
  if (length(overlaps) == 0) {
    stop("No overlapping CpG sites found between input and reference regions. Check genome build / coordinates.")
  }

  # Expand to site level (one row per overlapping CpG)
  site_idx   <- gr_sample$orig_id[queryHits(overlaps)]   # rows in df
  region_idx <- subjectHits(overlaps)                    # rows in ref_mat_region

  # beta_mat: one row per overlapping CpG site
  beta_mat <- as.matrix(df[site_idx, 3:ncol(df), drop = FALSE])
  beta_mat <- apply(beta_mat, 2, as.numeric)
  colnames(beta_mat) <- colnames(df)[3:ncol(df)]

  # ref_mat: replicate the corresponding region values for each site
  ref_mat <- ref_mat_region[region_idx, , drop = FALSE]

  # Identical rownames required by epidish
  common_rownames <- paste0("R", seq_len(nrow(beta_mat)))
  rownames(beta_mat) <- common_rownames
  rownames(ref_mat)  <- common_rownames

  cat("Successfully matched CpG sites:", nrow(beta_mat), "\n")
  cat("Final beta_mat / ref_mat dimensions:", paste(dim(beta_mat), collapse = " x "), "\n")

} else {
  # --------------------------------------------------------------------------
  # Mode 2: classic 450k-based matching (EpiDISH references)
  # --------------------------------------------------------------------------
  anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  gr_450k  <- GRanges(seqnames = anno450k$chr,
                      ranges   = IRanges(start = anno450k$pos, end = anno450k$pos),
                      cguid    = rownames(anno450k))

  if (do_liftover) {
    gr_hg38 <- GRanges(seqnames = df$chrom,
                       ranges   = IRanges(start = pmax(1, df$pos - 1), end = df$pos),
                       orig_id  = seq_len(nrow(df)))
    chain   <- import.chain(hg38ToHg19file)
    gr_hg19 <- unlist(liftOver(gr_hg38, chain))
    cat("liftOver successful sites:", length(gr_hg19), "/", nrow(df), "\n")

    overlaps     <- findOverlaps(gr_hg19, gr_450k)
    orig_indices <- gr_hg19$orig_id[queryHits(overlaps)]
    matched_cg   <- gr_450k$cguid[subjectHits(overlaps)]
  } else {
    cat("No --hg38ToHg19 provided → assuming input coordinates are already hg19\n")
    gr_input <- GRanges(seqnames = df$chrom,
                        ranges   = IRanges(start = df$pos, end = df$pos),
                        orig_id  = seq_len(nrow(df)))
    overlaps     <- findOverlaps(gr_input, gr_450k)
    orig_indices <- gr_input$orig_id[queryHits(overlaps)]
    matched_cg   <- gr_450k$cguid[subjectHits(overlaps)]
  }

  raw_matrix <- as.matrix(df[orig_indices, 3:ncol(df), drop = FALSE])
  raw_matrix <- apply(raw_matrix, 2, as.numeric)
  rownames(raw_matrix) <- matched_cg

  # Aggregate multiple hits to the same CpG
  beta_mat <- aggregate(raw_matrix, by = list(cg = rownames(raw_matrix)), FUN = mean, na.rm = TRUE)
  rownames(beta_mat) <- beta_mat$cg
  beta_mat <- as.matrix(beta_mat[, -1, drop = FALSE])
  cat("Successfully matched 450k CpGs:", nrow(beta_mat), "\n")

  # -------------------- Select reference matrix --------------------
  if (sample_type == "PBMC") {
    if (ref_choice == "cent12CT") {
      data(cent12CT.m)
      ref_mat <- cent12CT.m
      cat("Using reference: cent12CT.m (PBMC / adult blood)\n")
    } else {
      stop("--ref for PBMC only supports: cent12CT | devtEp")
    }

  } else if (sample_type == "CordBlood") {

    if (ref_choice == "uniLIFE7") {
      data(centUniLIFE.m)
      cord_cells <- c("B", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")
      available  <- intersect(cord_cells, colnames(centUniLIFE.m))
      if (length(available) < 5) {
        stop("Insufficient cord-blood cell types in reference. Columns: ",
             paste(colnames(centUniLIFE.m), collapse = ", "))
      }
      ref_mat <- centUniLIFE.m[, available, drop = FALSE]
      cat("Using reference: centUniLIFE.m (7 cord-blood cell types)\n")
      cat("Cells used:", paste(available, collapse = ", "), "\n")

    } else if (ref_choice == "uniLIFE") {
      data(centUniLIFE.m)
      ref_mat <- centUniLIFE.m
      cat("Using reference: centUniLIFE.m (full 19 cell types)\n")

    } else if (ref_choice == "FlowSorted.CordBlood.450k.ModelPars") {
      library(FlowSorted.CordBlood.450k)
      # Official pre-computed ModelPars (700 probes)
      cat("Loading FlowSorted.CordBlood.450k.ModelPars (official 700 probes)...\n")
      data(FlowSorted.CordBlood.450k.ModelPars, package = "FlowSorted.CordBlood.450k")
      ref_mat <- as.matrix(FlowSorted.CordBlood.450k.ModelPars)

      # Harmonize cell-type names
      colnames(ref_mat) <- gsub("CD4T|CD4", "CD4T", colnames(ref_mat), ignore.case = TRUE)
      colnames(ref_mat) <- gsub("CD8T|CD8", "CD8T", colnames(ref_mat), ignore.case = TRUE)
      colnames(ref_mat) <- gsub("Bcell|B cell|B-cell|^B$", "B", colnames(ref_mat), ignore.case = TRUE)
      colnames(ref_mat) <- gsub("Gran|Granulocyte", "Gran", colnames(ref_mat), ignore.case = TRUE)
      colnames(ref_mat) <- gsub("Mono|Monocyte", "Mono", colnames(ref_mat), ignore.case = TRUE)
      colnames(ref_mat) <- gsub("NK|NKcell", "NK", colnames(ref_mat), ignore.case = TRUE)
      colnames(ref_mat) <- gsub("nRBC|NRBC", "nRBC", colnames(ref_mat), ignore.case = TRUE)

      desired <- c("B", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")
      keep    <- intersect(desired, colnames(ref_mat))
      if (length(keep) < 5) {
        stop("Insufficient cell types in ModelPars. Columns: ",
             paste(colnames(ref_mat), collapse = ", "))
      }
      ref_mat <- ref_mat[, keep, drop = FALSE]
      cat("Using reference: FlowSorted.CordBlood.450k.ModelPars\n")
      cat("Cells used:", paste(colnames(ref_mat), collapse = ", "), "\n")

    } else if (ref_choice == "FlowSorted.CordBlood.450k") {
      library(FlowSorted.CordBlood.450k)
      # Build mean reference matrix from the full FlowSorted.CordBlood.450k dataset
      cat("Building reference matrix from FlowSorted.CordBlood.450k ...\n")
      data(FlowSorted.CordBlood.450k)
      Mset  <- preprocessRaw(FlowSorted.CordBlood.450k)
      beta  <- getBeta(Mset)
      cell_types <- pData(FlowSorted.CordBlood.450k)$CellType

      # Standardize common cell-type names
      cell_types <- gsub("CD4T", "CD4T", cell_types)
      cell_types <- gsub("CD8T", "CD8T", cell_types)
      cell_types <- gsub("Bcell|B cell|B-cell", "B", cell_types, ignore.case = TRUE)
      cell_types <- gsub("Gran|Granulocyte", "Gran", cell_types, ignore.case = TRUE)
      cell_types <- gsub("Mono|Monocyte", "Mono", cell_types, ignore.case = TRUE)
      cell_types <- gsub("NK|NKcell", "NK", cell_types, ignore.case = TRUE)
      cell_types <- gsub("nRBC|NRBC|nRBC", "nRBC", cell_types, ignore.case = TRUE)

      unique_cts <- unique(cell_types)
      cat("Detected cell types:", paste(unique_cts, collapse = ", "), "\n")

      ref_list <- lapply(unique_cts, function(ct) {
        idx <- which(cell_types == ct)
        if (length(idx) == 0) return(NULL)
        rowMeans(beta[, idx, drop = FALSE], na.rm = TRUE)
      })
      names(ref_list) <- unique_cts
      ref_list <- ref_list[!sapply(ref_list, is.null)]
      ref_mat  <- do.call(cbind, ref_list)

      # Keep only the standard 7 cord-blood cell types
      desired <- c("B", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")
      keep    <- intersect(desired, colnames(ref_mat))
      if (length(keep) < 5) {
        stop("Insufficient cord-blood cell types in FlowSorted.CordBlood.450k. Actual columns: ",
             paste(colnames(ref_mat), collapse = ", "))
      }
      ref_mat <- ref_mat[, keep, drop = FALSE]
      cat("Using reference: FlowSorted.CordBlood.450k (mean matrix of 7 cord-blood cell types)\n")
      cat("Cells used:", paste(colnames(ref_mat), collapse = ", "), "\n")

    } else {
      stop("--ref for CordBlood only supports: uniLIFE7 | uniLIFE | FlowSorted.CordBlood.450k.ModelPars | FlowSorted.CordBlood.450k")
    }

  } else {
    stop("sample_type must be 'PBMC' or 'CordBlood'")
  }
}

cat("Final reference matrix dimensions:", paste(dim(ref_mat), collapse = " x "), "\n")

# -------------------------- EpiDISH deconvolution -----------------------------
out.epidish <- epidish(beta.m = beta_mat, ref.m = ref_mat, method = "RPC")
cat("Number of common features:", nrow(out.epidish$ref), "\n")

cellFrac <- out.epidish$estF
print(head(cellFrac))
cat("\nMean cell proportions:\n")
print(round(colMeans(cellFrac), 4))

# -------------------------- 输出结果 ------------------------------------------
write.table(
  cbind(sample_name = rownames(cellFrac), cellFrac),
  file      = cellratiof,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("\nFinished!\n")
cat("Cell proportion file written to:", cellratiof, "\n")
