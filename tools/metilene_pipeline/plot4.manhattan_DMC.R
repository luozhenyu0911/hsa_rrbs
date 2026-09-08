print("################  Processing forth step: Plot the mahattan plot for DMC sites ##################################")
library(CMplot)
args=commandArgs(T)


var <- args[1]
filef <- paste0("01.DMR_DMC_calling/",var,".DMC_pval.0.05.out")

data<- read.table(filef,sep = '\t')
# colnames(data)<-c('Chr', 'start', 'End', 'Diff', 
#                   'qval', 'count','Amean', 'Bmean')
head(data)

folders <- c("02.anno_plot/Out_Manhtn")

# Loop through each folder name and create the directory if it doesn't exist
for (folder in folders) {
    if (!dir.exists(folder)) {
        dir.create(folder, recursive = TRUE)
        cat("Directory", folder, "created.\n")
    } else {
        cat("Directory", folder, "already exists.\n")
    }
}

data$name <- paste0(data$V1,data$V2)

data2<-data[,c(8,1,2,3)]

CMplot(data2,plot.type="d",bin.size=1e6,col=c("darkgreen", "yellow", "red"),
       file="pdf",file.name='plot2',main = "DMC distribution diagram",dpi=300,file.output=TRUE, verbose=TRUE)


CMplot(data2,plot.type="m",LOG10=TRUE,threshold=NULL,chr.den.col=NULL,
       file="pdf",file.name='plot3',main = "DMC distribution diagram",dpi=300,file.output=TRUE, verbose=TRUE)

print("################  Processing forth step: Plot the mahattan plot for DMC sites ##################################")
