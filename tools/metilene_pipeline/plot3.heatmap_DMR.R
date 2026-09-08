print("################  Processing third step: third step: Plot the heatmap to show the differential methylation of DMRs/DMCs ##################################")

.libPaths(c("/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta2/lib/R/library"))

library(reshape2)
library(ggplot2)
library(dplyr)
require(cowplot)
library(ComplexHeatmap)

args=commandArgs(T)

varname <- args[1] 
pheno <- args[2]  

folders <- c("02.anno_plot/Out_heatmap")

# Loop through each folder name and create the directory if it doesn't exist
for (folder in folders) {
    if (!dir.exists(folder)) {
        dir.create(folder, recursive = TRUE)
        cat("Directory", folder, "created.\n")
    } else {
        cat("Directory", folder, "already exists.\n")
    }
}

###############读取文件合并成绘图文件
anno = read.csv(paste0("02.anno_plot/Out_anno/","outAnno.",varname,".csv"),header = 1)
met = read.csv(paste0("01.DMR_DMC_calling/",varname,".DMR_qval.0.05.sampleMet.bed"),header = F,sep = "\t")
# met = met[,-71] ## 删除最后一列NA
###读取跑甲基化之前的phe文件,注意列名要与varname一致
phe= read.csv(pheno,header = T)
t<-data.frame(c("seqnames_1","start_1","end_1","qvalue_1","mean_met_diff_1",
                "CpGs_1","mean_group1_1","mean_group2_1"))
colnames(t)<-strsplit(varname, "\\.")[[1]][1]
vs_name<-strsplit(varname, "\\.")[[1]][1]
colname<-rbind(t,phe[vs_name])
colnames(met) <- as.list(colname[[1]])#$之后加变量名
mydata<-met

# mydata<-cbind(anno,met)
# ##############检查2个文件列是否一致
# all(mydata$seqnames == mydata$seqnames_1)
# all(mydata$start == mydata$start_1+1)
# all(mydata$end == mydata$end_1)

# dat = subset(mydata,qvalue_1<0.0001)
# #dat <- dplyr::arrange(mydata,qvalue)[1:500,]
# 
# dim(dat)
# 
# tmp1 <- dat[,grepl("g1", colnames(dat))]
# tmp2 <- dat[,grepl("g2", colnames(dat))]
# group1 <- c(rep("g1",ncol(tmp1)),rep("g2",ncol(tmp2)))
# df2 <- cbind(tmp1,tmp2)
# 
# #df2 <- dat[,25:84]
# df2 <- as.matrix(df2)
# sampleGroup <- factor(group1, levels = c("g1","g2"))
# annot_df2 <- data.frame(Group=sampleGroup)
# colo = list(Group = c("g1" = '#66C2A5', "g2" = '#377EB8'))
# .violin = anno_density(df2, type = "violin", gp = gpar(fill = "lightblue"))
# ha_mix_top = HeatmapAnnotation(df = annot_df2, col=colo, violin = .violin,height = unit(2, "cm"))
# 
# ### 按行中心标准化，减均值除方差，是Z-score转换
# df2_scaled = apply(df2, 1, scale)
# # 继续原数据表列名
# rownames(df2_scaled) = colnames(df2)
# # 转置才与原方向一致
# df2_scaled = t(df2_scaled)
# P1 <- Heatmap(df2_scaled,name = "rel.methylation",
#               show_row_names = FALSE,
#               show_column_names = FALSE,
#               #show_row_dend = FALSE,
#               cluster_rows = TRUE,
#               ##row_split = 4, ##行分为4组
#               #cluster_columns = FALSE,
#               top_annotation = ha_mix_top, na_col="white") 
# ##获取差异甲基化聚类后的分组信息
# #site_group <- row_order(P1)
# #site_group1 <- mydata[df2_scaled[1],]
# #site_group2 <- mydata[df2_scaled[2],]
# #site_group3 <- mydata[df2_scaled[3],]
# #site_group4 <- mydata[df2_scaled[4],]
# outPDF <- paste0("02.anno_plot/Out_heatmap/Fig.",varname,".heatmap.pdf")
# pdf(file = outPDF, height = 8, width = 8)
# P1
# dev.off()
# 
# ############### no scaled ##############
# color = colorRampPalette(colors = c ("#7e52fd","white","#ff5636")) (100)
# P2 <- Heatmap(df2,name = "rel.methylation",
#               #col = color,
#               show_row_names = FALSE,
#               show_column_names = FALSE,
#               show_row_dend = FALSE,
#               cluster_rows = TRUE,
#               #cluster_columns = FALSE,
#               top_annotation = ha_mix_top, na_col="white") 
# 
# outPDF2 <- paste0("02.anno_plot/Out_heatmap/Fig.",varname,".heatmap2.pdf")
# pdf(file = outPDF2, height = 8, width = 8)
# P2
# dev.off()


############ 小提琴图 ########################
#varname <- "WBC3" # "IFNg", "IL1b", "IL6", "TNFa", "WBC1", "WBC2", "WBC3"
#varname <- "IL10" # "IL10","IP10"

outPDF3 <- paste0("02.anno_plot/Out_heatmap/Fig.",varname,".violin.pdf")

###############
dat_up = subset(mydata,qvalue_1<0.05 & mean_met_diff_1>(0.1))
dat_down = subset(mydata,qvalue_1<0.05 & mean_met_diff_1<(-0.1))

dat1control <- dat_up[,grepl("g2|seqnames", colnames(dat_up))]
dat1case <- dat_up[,grepl("g1|seqnames", colnames(dat_up))]
dat2control <- dat_down[,grepl("g2|seqnames", colnames(dat_down))]
dat2case <- dat_down[,grepl("g1|seqnames", colnames(dat_down))]


hyperN <- paste0("hyper-dmrs (",nrow(dat_up),")")
hypoN <- paste0("hypo-dmrs (",nrow(dat_down),")")

dat1case_long <- melt(dat1case,value.name=c("methylation"))
dat1case_long["Group1"] <- rep(hyperN,nrow(dat1case))
dat1case_long["Group"] <- rep("High",nrow(dat1case))
dat1control_long <- melt(dat1control,value.name=c("methylation"))
dat1control_long["Group1"] <- rep(hyperN,nrow(dat1control))
dat1control_long["Group"] <- rep("Low",nrow(dat1control))
dat2case_long <- melt(dat2case,value.name=c("methylation"))
dat2case_long["Group1"] <- rep(hypoN,nrow(dat2case))
dat2case_long["Group"] <- rep("High",nrow(dat2case))
dat2control_long <- melt(dat2control,value.name=c("methylation"))
dat2control_long["Group1"] <- rep(hypoN,nrow(dat2control))
dat2control_long["Group"] <- rep("Low",nrow(dat2control))
dat_long <- rbind(dat1case_long,dat1control_long,dat2case_long,dat2control_long)
dat_long$Group <- factor(dat_long$Group,levels = c("High","Low"))
dat_long3 <- subset(dat_long,methylation>0)
P2 <- ggplot(data = dat_long3,aes(x=Group,y=methylation,fill=Group)) + geom_violin() + 
  geom_boxplot(width=0.1,color="black",notch=TRUE,size=0.8) + facet_grid(~Group1)+theme_bw()+
  ggtitle(paste0("High vs Low ",varname))+theme(plot.title = element_text(hjust = 0.5,size=18))

pdf(file = outPDF3, height = 5, width = 7)
P2
dev.off()

print("################  Completed third step: Plot the heatmap to show the differential methylation of DMRs/DMCs ##################################")
