library(getopt)

spec <- matrix(
  c("sample", "s", 1, "character", "name,expression matrix",
    "meta", "m", 1, "character", "meta.information",
    "help", "h", 0, "logical", "This is the help!"),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$meta) || is.null(opt$sample)){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}

library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)
utils::methods(class = 'Seurat')
library(tools)

data<-read.table(file_path_as_absolute(opt$sample),header=T,sep="\t")
meta<-read.table(file_path_as_absolute(opt$meta),header=T,sep="\t")

cell_stat<-as.data.frame(matrix(nrow=0,ncol=2))
merge<-merge(data,meta,by="sample")

nrows<-nrow(data)
cols<-colnames(merge)

for(subrow in 1:nrows){

    sam_name<-merge[subrow,1]
    sam_path<-file_path_as_absolute(merge[subrow,2])
    
    if(grepl("h5$",sam_path)){
       sam_data <-Read10X_h5(sam_path)
    }else{
        sam_data<-Read10X(data.dir=sam_path)
    }  
    seuobj<-CreateSeuratObject(counts = sam_data , project = sam_name, min.cells = 3, min.features = 200)

    seuobj[["percent.mt"]] <- PercentageFeatureSet(seuobj, pattern = "^MT")
    VlnPlot(seuobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(paste("step1_",sam_name,"_raw_Count_MT_VlnPlot.pdf",sep=""),width=8,height=6,dpi=600)
    ggsave(paste("step1_",sam_name,"_raw_Count_MT_VlnPlot.png",sep=""),width=8,height=6,dpi=600)

    plot1 <- FeatureScatter(seuobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(seuobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    plot1 + plot2
    ggsave(paste("step1_",sam_name,"_raw_Count_MT_Scatter.pdf",sep=""),width=8,height=6,dpi=600)
    ggsave(paste("step1_",sam_name,"_raw_Count_MT_Scatter.png",sep=""),width=8,height=6,dpi=600)

    raw_cell <- nrow(seuobj@meta.data)
    for(order in 3:length(cols)){
      seuobj[[cols[order]]]<-merge[subrow,order]
    }
    seuobj$sample<-sam_name

    tmp_cell_num<-as.data.frame(matrix(c(sam_name,raw_cell),nrow=1,ncol=2))
    cell_stat<-rbind(cell_stat,tmp_cell_num)

    saveRDS(seuobj,paste("step1_",sam_name,"_raw.Rds",sep=""))
    
}
colnames(cell_stat)<-c("sample","raw_cell")
write.table(cell_stat,paste("step1_samples_raw_cell_num.csv",sep=""),sep=",",row.names=FALSE,quote=FALSE)
