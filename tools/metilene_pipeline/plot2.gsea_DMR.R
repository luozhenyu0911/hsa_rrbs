print("################  Processing second step: Function enrichment for genes annotated with DMRs/DMCs ##################################")

.libPaths(c("/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta2/lib/R/library"))

library(ggplot2)
library(ggimage)
library(ChIPseeker)
library(org.Hs.eg.db)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(patchwork)
library(enrichplot)
library(DOSE)
library(pathview)


# remotes::install_github("YuLab-SMU/createKEGGdb")
# library(createKEGGdb)
# 
# species <-c("hsa", "mmu", "rno") 
# create_kegg_db(species)
# install.packages("./KEGG.db_1.0.tar.gz", repos=NULL)
# 确保使用双反斜杠并且路径正确
#install.packages("D:\\A Tutu\\01 博士\\01 博士课题\\001 结果-方法\\3-甲基化-LHH\\db_KEGG\\KEGG.db_hsa_1.0.tar.gz", repos=NULL, type="source")
library(KEGG.db)

args=commandArgs(T)

# inBed <- args[1]
# prefix <- sub("\\.[^.]+$", "", inBed)

prefix <- args[1]
inBed <- paste0("01.DMR_DMC_calling/",prefix,".DMRs_qval.0.05.out")

################ GSEA MsigDB ########################
## 只保留基因及基因临近区域: Promoter (<=1kb), gene, Downstream (<=300bp)
#dat_Anno1 <- subset(dat_Anno, annotation!="Distal Intergenic" & annotation!="Promoter (1-2kb)" & annotation!="Promoter (2-3kb)")
#dat_Anno1 <- filter(dat_Anno, !grepl("Distal Intergenic|Promoter|Downstream",annotation))

## 聚合同一个基因的甲基化水平，用mean_met_diff和CpGs做加权平均
#dat_Anno1$met_diff <- dat_Anno1$mean_met_diff * dat_Anno1$CpGs
#x1 <- aggregate(dat_Anno1$met_diff, by = list(dat_Anno1$geneId), FUN=sum)
#x2 <- aggregate(dat_Anno1$CpGs, by = list(dat_Anno1$geneId), FUN=sum)
#glist <- x1$x / x2$x
#names(glist) <- x1$Group.1
#glist <- sort(glist,decreasing = T)
#length(glist)

# dir.create ("02.anno_plot/Out_anno")
# dir.create ("02.anno_plot/Out_GO_KEGG")

# Define the names of the folders to create
folders <- c("02.anno_plot/Out_anno", "02.anno_plot/Out_GO_KEGG")

# Loop through each folder name and create the directory if it doesn't exist
for (folder in folders) {
    if (!dir.exists(folder)) {
        dir.create(folder, recursive = TRUE)
        cat("Directory", folder, "created.\n")
    } else {
        cat("Directory", folder, "already exists.\n")
    }
}


analysisGSEA <- function(varname, qcut=0.01){
  
  inDat <- paste0("02.anno_plot/Out_anno/out.",varname,".Rdata")
  load(inDat)
  
  ### differential methylation genes enrich analysis
  gene = unique(subset(dat_Anno, qvalue<qcut)$geneId)
  print(length(gene))
  # Run GO enrichment analysis 
  ego <- enrichGO(gene = gene, 
                  keyType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "ALL",  # BP, MF, CC, ALL
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)
  # 进行KEGG富集
  ekegg <- enrichKEGG(gene = gene,
                      organism = "hsa",
                      keyType = "kegg",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      use_internal_data=T)
  # plot
  fig.ego <- dotplot(ego, showCategory=6, split="ONTOLOGY",
                     color="pvalue", orderBy ="pvalue",decreasing=FALSE,
                     label_format = 35,font.size=10) + 
    facet_grid(ONTOLOGY~.,scale="free") + 
    ggtitle("Ehrich GO result of \ndifferential methylation genes (qvalue<0.1)")
  
  fig.ekegg <- dotplot(ekegg, showCategory=20, color="pvalue", 
                       orderBy ="pvalue",decreasing=FALSE,label_format = 35,font.size=10) +
    ggtitle("Ehrich KEGG pathway result of \ndifferential methylation genes (qvalue<0.1)")
  
  x <- pairwise_termsim(ekegg)
  fig.emap <- emapplot(x,layout="fr",color="pvalue")
  
  write.csv(ego@result, file = paste0("02.anno_plot/Out_GO_KEGG/dmg_",varname,"_go.csv"))
  write.csv(ekegg@result, file = paste0("02.anno_plot/Out_GO_KEGG/dmg_",varname,"_kegg.csv"))
  try(ggsave(fig.ego, filename = paste0("02.anno_plot/Out_GO_KEGG/dmg_",varname,"_go.pdf"), 
             width = 10, height = 8), silent = TRUE)
  try(ggsave(fig.ekegg, filename = paste0("02.anno_plot/Out_GO_KEGG/dmg_",varname,"_kegg.pdf"), 
             width = 10, height = 8), silent = TRUE)
  try(ggsave(fig.emap, filename = paste0("02.anno_plot/Out_GO_KEGG/dmg_",varname,"_emap.pdf"),
             width = 10, height = 8), silent = TRUE)
  
  #if (FALSE){
  ## 只保留基因及基因临近区域: Promoter (<=1kb), gene, Downstream (<=300bp)
  dat_Anno1 <- subset(dat_Anno, annotation!="Distal Intergenic" & 
                        annotation!="Promoter (1-2kb)" & 
                        annotation!="Promoter (2-3kb)")
  # dat_Anno1 <- filter(dat_Anno, !grepl("Distal Intergenic|Promoter|Downstream",annotation))
  
  ## 聚合同一个基因的甲基化水平，用mean_met_diff和CpGs做加权平均
  dat_Anno1$met_diff <- dat_Anno1$mean_met_diff * dat_Anno1$CpGs
  x1 <- aggregate(dat_Anno1$met_diff, by = list(dat_Anno1$geneId), FUN=sum)
  x2 <- aggregate(dat_Anno1$CpGs, by = list(dat_Anno1$geneId), FUN=sum)
  glist <- x1$x / x2$x
  names(glist) <- x1$Group.1
  glist <- sort(glist,decreasing = T)
  length(glist)
  
  dir <- "/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/database/msigdb_v7.4_GMTs/"
  gmt1 <- read.gmt(paste0(dir,"c2.cp.kegg.v7.4.entrez.gmt"))
  gmt2 <- read.gmt(paste0(dir,"c5.go.v7.4.entrez.gmt"))
  
  gsea_kegg <- GSEA(glist, TERM2GENE = gmt1, pvalueCutoff=1, seed=T)
  gsea_go <- GSEA(glist, TERM2GENE = gmt2, pvalueCutoff=1, seed=T)
  
  write.csv(gsea_kegg@result, file = paste0("02.anno_plot/Out_GO_KEGG/gsea_",varname,"_kegg.csv"))
  write.csv(gsea_go@result, file = paste0("02.anno_plot/Out_GO_KEGG/gsea_",varname,"_go.csv"))
  fig1 <- dotplot(gsea_kegg,split=".sign",
                  color="pvalue",
                  orderBy ="pvalue",decreasing=FALSE,
                  label_format = 35,
                  font.size=10)+facet_grid(~.sign) #出点图，并且分面激活和抑制
  fig2 <- dotplot(gsea_go,split=".sign",
                  color="pvalue",
                  orderBy ="pvalue",decreasing=FALSE,
                  label_format = 35,
                  font.size=10)+facet_grid(~.sign) #出点图，并且分面激活和抑制
  
  try(ggsave(fig1, filename = paste0("02.anno_plot/Out_GO_KEGG/gsea_",varname,"_kegg.pdf"),
             width = 10, height = 8),silent = TRUE)
  try(ggsave(fig2, filename = paste0("02.anno_plot/Out_GO_KEGG/gsea_",varname,"_go.pdf"), 
             width = 10, height = 8),silent = TRUE)
  
  result <- list(enrich_go=ego,gsea_go=gsea_go,gsea_kegg=gsea_kegg,glist=glist)#,enrich_kegg=ekegg
  return(result)
}

analysis_1 <- analysisGSEA(prefix,qcut = 0.01)
# analysis_2 <- analysisGSEA(varname = vars[2],qcut = 0.01)

analysis_1$glist[analysis_1$glist>0]
analysis_1$glist[analysis_1$glist<0]
##############view pathway#############################
# 查看*_kegg.csv文件的ID列，去掉hsa前缀就是pathway.id
pathview(gene.data = analysis_1$glist, pathway.id = "05130", species = "hsa",
         out.suffix = prefix, kegg.native = T)

print("################  Completed second step: Function enrichment for genes annotated with DMRs/DMCs ##################################")
