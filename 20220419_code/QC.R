library(getopt)

spec <- matrix(
  c("sample", "s", 1, "character", "name,rds,min_nFeature_RNA,max_nFeature_RNA,min_nCount_RNA,max_nCount_RNA,percent.mt!",
    "help",  "h", 0, "logical", "This is the help!"),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$sample)){
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
cell_stat<-as.data.frame(matrix(nrow=0,ncol=2))

nrows<-nrow(data)

for(subrow in 1:nrows){

    sam_name<-data[subrow,1]
    sam_rds<-file_path_as_absolute(data[subrow,2])
    min_nFeature<-data[subrow,3]
    max_nFeature<-data[subrow,4]
    min_nCount<-data[subrow,5]
    max_nCount<-data[subrow,6]
    max_mt<-data[subrow,7]
    
    seuobj<-readRDS(sam_rds)
    seuobj <- subset(seuobj, subset = nFeature_RNA>min_nFeature & nFeature_RNA<max_nFeature & percent.mt<max_mt & nCount_RNA>min_nCount & nCount_RNA<max_nCount)
    
    seuobj[["percent.mt"]] <- PercentageFeatureSet(seuobj, pattern = "^MT")
    VlnPlot(seuobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(paste("step2_",sam_name,"_filter_Count_MT_VlnPlot.pdf",sep=""),width=8,height=6,dpi=600)
    ggsave(paste("step2_",sam_name,"_filter_Count_MT_VlnPlot.png",sep=""),width=8,height=6,dpi=600)

    plot1 <- FeatureScatter(seuobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(seuobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    plot1 + plot2
    ggsave(paste("step2_",sam_name,"_filter_Count_MT_Scatter.pdf",sep=""),width=8,height=6,dpi=600)
    ggsave(paste("step2_",sam_name,"_filter_Count_MT_Scatter.png",sep=""),width=8,height=6,dpi=600)

    clean_cell <- nrow(seuobj@meta.data)
    cell_num <- as.data.frame(matrix(c(sam_name,clean_cell),nrow=1,ncol=2))
    cell_stat <- rbind(cell_stat,cell_num)

    seuobj <- NormalizeData(seuobj, normalization.method = "LogNormalize", scale.factor = 10000)

    seuobj <- FindVariableFeatures(seuobj, selection.method = "vst", nfeatures = 2000)
    seuobj <- RenameCells(object = seuobj, add.cell.id = sam_name)
    saveRDS(seuobj, file = paste("step2_",sam_name,"_filter.Rds",sep=""))
}

colnames(cell_stat)<-c("sample","clean_cell")
cell_stat$clean_cell<-as.numeric(cell_stat$clean_cell)
write.table(cell_stat,"step2_samples_clean_cell_num.csv",sep=",",quote=FALSE,col.names=TRUE,row.names=FALSE)
