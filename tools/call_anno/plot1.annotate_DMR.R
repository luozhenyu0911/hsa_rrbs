# setwd('E:\\tmp_other\\04_lhh_meth\\pheno')

# BiocManager::install('ChIPseeker')
# BiocManager::install('ChIPpeakAnno')
# BiocManager::install('org.Hs.eg.db')
# BiocManager::install('HDO.db')
# ### https://www.jianshu.com/p/e0515071f175
# BiocManager::install("TxDb.Drerio.UCSC.danRer11.refGene")
# #BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# #BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# BiocManager::install("clusterProfiler")

# options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))
# install.packages("magick")
# install.packages("ggimage")
# install.packages('ggupset')

print("###################################  Processing first step: annotate DMR/DMCs with ChIPpeakAnno  ##################################")
print("Loading required packages...")

.libPaths(c("/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta2/lib/R/library"))

library(ggplot2)
library('GO.db')
library(ChIPseeker)
library(org.Hs.eg.db)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggupset)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
args=commandArgs(T)

prefix <- args[1]
inBed <- paste0("01.DMR_DMC_calling/",prefix,".DMR_qval.0.05.out")

directory_path <- "02.anno_plot/Out_anno"

if (dir.exists(directory_path)) {
  cat("The directory '02.anno_plot/Out_anno' exists.\n")
} else {
  dir.create ("02.anno_plot/Out_anno", recursive = TRUE)
}

# prefix <- sub("\\.[^.]+$", "", inBed)

outPeakFig <- paste0("02.anno_plot/Out_anno/outPeakFig.",prefix,".pdf")
outResult <- paste0("02.anno_plot/Out_anno/outAnno.",prefix,".csv")
outR <- paste0("02.anno_plot/Out_anno/out.",prefix,".Rdata")

dat <- readPeakFile(inBed, header=F)
dat_peakAnno <- annotatePeak(dat, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb="org.Hs.eg.db")

# plot figures
pdf(outPeakFig, height = 5, width = 11, onefile = T)
#plotAnnoBar(dat_peakAnno)
plotAnnoPie(dat_peakAnno)
#vennpie(dat_peakAnno)
upsetplot(dat_peakAnno, vennpie=T)
plotDistToTSS(dat_peakAnno,title="Distribution of transcription methylation loci relative to TSS")
dev.off()

# save to csv
dat_Anno = as.data.frame(dat_peakAnno)

colnames(dat_Anno) <- c("seqnames","start","end","width","strand","qvalue",
    "mean_met_diff","CpGs","mean_group1",
    "mean_group2","annotation","geneChr","geneStart","geneEnd","geneLength",
    "geneStrand","geneId","transcriptId","distanceToTSS","ENSEMBL","SYMBOL","GENENAME")
write.csv(dat_Anno, file=outResult, row.names = FALSE)
save(dat_peakAnno, dat_Anno, file = outR)

print("###################################  Completed first step: annotate DMR with ChIPpeakAnno  ##################################")
