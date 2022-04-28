library(getopt)

spec <- matrix(
  c("rds", "r", 1, "character", "10x RDS file!",
    "infor","i", 1, "character", "information of cell",
    "resolution", "s", 2, "numeric", "resolution!",
    "dims", "d", 2, "numeric", "the number of dimension of PCA/harmony!",   
    "help",  "h", 0, "logical", "This is the help!"),
  byrow=TRUE, ncol=5)

opt<-getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$rds) || is.null(opt$infor)){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}

library(dplyr)
library(Seurat)
library(cowplot)
utils::methods(class = 'Seurat')
library(ggplot2)
library(tools)

data <- readRDS(file_path_as_absolute(opt$rds))

res<-"0.8"
if(!is.null(opt$resolution)){
        res<-opt$resolution
}

res<-as.numeric(res)

ndim<-"20"
if(!is.null(opt$dims)){
        ndim<-opt$dims
}
ndim<-as.numeric(ndim)

if(!is.null (opt$infor)){
        infor<-opt$infor

        infors<-unlist(strsplit(infor,split=";"))

        lens<-length(infors)
        for(order in 1:lens){
                subinfor<-unlist(strsplit(infors[order],split=":"))
                ssinfor<-unlist(strsplit(subinfor[2],split=","))
                if(subinfor[1] == "idents"){
                        data<-subset(data,idents = ssinfor)
                }else{
                        rnames<-rownames(data@meta.data[which(data@meta.data[[subinfor[1]]] %in% ssinfor),])
                        data<-subset(data,cells = rnames)
                }
        }
}

DefaultAssay(data) <- "RNA"
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

DefaultAssay(data) <- "integrated"
data <- ScaleData(data, verbose = FALSE)

data <- RunPCA(data, npcs = 30, verbose = FALSE)

## Determine the ‘dimensionality’ of the dataset ##
data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)
JackStrawPlot(data, dims = 1:20)
ggsave("recluster.JackStrawPlot.pdf",width=9,height=6,dpi=600)
ggsave("recluster.JackStrawPlot.png",width=9,height=6,dpi=600)
ElbowPlot(data)
ggsave("recluster.Elbowplot.pdf",width=9,height=6,dpi=600)
ggsave("recluster.Elbowplot.png",width=9,height=6,dpi=600)

###### Cluster the cells ######
if("harmony" %in% names(data@reductions)){
	data <- FindNeighbors(data, reduction = "harmony", dims = 1:ndim) %>% FindClusters()
}else{
	data <- FindNeighbors(data, dims = 1:ndim)
	data <- FindClusters(data, resolution = res)
}
data$celltype <- Idents(data)

##### Run non-linear dimensional reduction (UMAP/tSNE) ######
data <- RunUMAP(data, dims = 1:20)
DimPlot(data, reduction = "umap", label = TRUE)
ggsave("recluster.raw_umap.pdf",width=9,height=6,dpi=600)
ggsave("recluster.raw_umap.png",width=9,height=6,dpi=600)

p2 <- DimPlot(data, reduction = "umap", group.by = "sample")
p1 <- DimPlot(data, reduction = "umap", label = TRUE)
p1 +  p2
ggsave("recluster.sample.raw_umap.pdf",width=8,height=6,dpi=600)
ggsave("recluster.sample.raw_umap.png",width=8,height=6,dpi=600)

sams<-length(unique(data@meta.data$sample))

numwidth<-""
if(sams >4){
    numwidth<-as.numeric("12")
    
}else{
    numwidth<-sams*3
}
numheight<-""
numcols<-ceiling(sams/4)
numheight<-as.numeric(numcols*3)

DimPlot(data, reduction = "umap", split.by = "sample",ncol=4)
ggsave("recluster.raw_umap.sample.pdf",width=numwidth,height=numheight,dpi=600)

class.num <- table(data$sample,data@meta.data$seurat_clusters)
write.table("sample.recluster.num.xls",col.names =NA,quote =FALSE,sep ="\t")
saveRDS(data,"recluster.dimreduce.Rds")

data.markers <- FindAllMarkers(data, test.use = "MAST",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(data.markers,"recluster.exp_gene.diff.pos.xls",sep="\t",row.names =FALSE, col.names =TRUE,  quote =FALSE)

DefaultAssay(data) <- "integrated"
top20 <- data.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
p <- DoHeatmap(data, features = top20$gene,label=T,size=5,angle=90,group.bar.height = 0.05,group.by ="celltype") + NoLegend()
p<-p+theme(axis.text.y=element_blank())
p<-p+theme(plot.margin = unit(c(4,1,1,1),"cm"))
ggsave("recluster.diff_gene.top20.png",width=12,height=9)
ggsave("recluster.diff_gene.top20.pdf",width=12,height=9)
