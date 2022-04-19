library(getopt)

spec <- matrix(
  c("sample", "s", 1, "character", "name,rds",
    "batch", "b", 1,"character", "remove batch effect by information!",
    "method", "m", 2, "character", "the method of integrating the data!",
    "resolution", "r", 2, "numeric", "resolution!",
    "dims", "d", 2, "numeric", "the number of dimension of PCA/harmony!",
    "help",  "h", 0, "logical", "This is the help!"),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$sample) || is.null(opt$batch)){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}

library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)
utils::methods(class = 'Seurat')
library(tools)

batch<-opt$batch

res<-"0.8"
if(!is.null(opt$resolution)){
	res<-opt$resolution	
}

method<-"integrated"
if(!is.null(opt$method)){
	method<-opt$method	
}
ndim<-"20"
if(!is.null(opt$dims)){
	ndim<-opt$dims
}

data<-read.table(file_path_as_absolute(opt$sample),header=T,sep="\t")
cell_stat<-as.data.frame(matrix(nrow=0,ncol=2))

nrows<-nrow(data)
obj.list<-list()
sumcell<-"0"
sumcell<-as.numeric(sumcell)

combined<-""

for(subrow in 1:nrows){

    sam_name<-data[subrow,1]
    sam_rds<-file_path_as_absolute(data[subrow,2])
    
    seuobj<-readRDS(sam_rds)
  
    tmpcell<-nrow(seuobj@meta.data)
    sumcell<-tmpcell+sumcell
    #obj.list<-append(obj.list,seuobj)
    if(subrow ==1){
    	combined<-seuobj
     }else{
	combined<-merge(combined,seuobj,merge.data=TRUE)
    }
}

obj.list<-SplitObject(combined, split.by = batch)

features <- SelectIntegrationFeatures(object.list = obj.list)

###### IntegrateData ######

data.anchors<-""
data.combined<-""
dimmax<-"30"
if(sumcell >200000 & nrows>2){
	dimmax<-"50"
	dimmax <- as.numeric(dimmax)
	obj.list <- lapply(X = obj.list, FUN = function(x) {
		x <- ScaleData(x, features = features, verbose = FALSE)
		x <- RunPCA(x, features = features, verbose = FALSE)
	})
	data.anchors <- FindIntegrationAnchors(object.list = obj.list, reference = c(1,2), reduction = "rpca",dims = 1:50)
     
}else if(sumcell >200000 & nrows<=2){
    	dimmax<-"50"
	dimmax <- as.numeric(dimmax)
	obj.list <- lapply(X = obj.list, FUN = function(x) {
                x <- ScaleData(x, features = features, verbose = FALSE)
                x <- RunPCA(x, features = features, verbose = FALSE)
	})
	data.anchors <- FindIntegrationAnchors(object.list = obj.list, reference = c(1), reduction = "rpca",dims = 1:50)
}else{
	data.anchors<-FindIntegrationAnchors(object.list = obj.list,anchor.features = features)
}

if(method == "integrated"){
	data.combined <- IntegrateData(anchorset = data.anchors,dim=1:dimmax)
	
	DefaultAssay(data.combined) <- "integrated"
	data.combined <- ScaleData(data.combined, verbose = FALSE)
	data.combined <- RunPCA(data.combined, npcs = 30 , verbose = FALSE)

	data.combined <- JackStraw(data.combined, num.replicate = 100,dims = 30)
	data.combined <- ScoreJackStraw(data.combined, dims = 1:30)
	JackStrawPlot(data.combined, dims = 1:30)
	ggsave("step3_combined_JackStrawPlot.pdf",width=8,height=6,dpi=600)
	ggsave("step3_combined_JackStrawPlot.png",width=8,height=6,dpi=600)
	ElbowPlot(data.combined,ndims = 30)
	ggsave("step3_combined_Elbowplot.pdf",width=8,height=6,dpi=600)
	ggsave("step3_combined_Elbowplot.png",width=8,height=6,dpi=600)

	###### Cluster the cells ######
	data.combined <- FindNeighbors(data.combined, dims = 1:ndim)
	data.combined <- FindClusters(data.combined, resolution = res)
	
	##### Run non-linear dimensional reduction (UMAP/tSNE) ######
	data.combined <- RunUMAP(data.combined, dims = 1:30)
	data.combined$celltype <- Idents(data.combined)
	saveRDS(data.combined, file = "step3_comb_dimreduce.Rds")

}else if(method == "harmony"){

	data.combined <- NormalizeData(combined)
	data.combined <- FindVariableFeatures(data.combined)
	data.combined <- ScaleData(data.combined, verbose = FALSE)
	data.combined <- RunPCA(data.combined, npcs = 30 , verbose = FALSE)
		
	data.combined <- RunHarmony(data.combined, group.by.vars = batch)
	data.combined <- RunUMAP(data.combined, reduction = "harmony", dims = 1:30)
	data.combined <- FindNeighbors(data.combined, reduction = "harmony", dims = 1:ndim) %>% FindClusters()
	data.combined$celltype <- Idents(data.combined)
	saveRDS(data.combined, file = "step3_comb_dimreduce.Rds")
	
}
	
DimPlot(data.combined, reduction = "umap", label = TRUE)
ggsave("step3_raw_umap.pdf",width=8,height=6,dpi=600)
ggsave("step3_raw_umap.png",width=8,height=6,dpi=600)

p1 <- DimPlot(data.combined, reduction = "umap", group.by = "group")
p2 <- DimPlot(data.combined, reduction = "umap", group.by = "sample")
p1 + p2
ggsave("step3_umap_group_sample.pdf",width=12,height=6,dpi=600)
ggsave("step3_umap_group_sample.png",width=12,height=6,dpi=600)

numwidth<-""
if(nrows >4){
	numwidth<-"12"
    	numwidth<-as.numeric(numwidth)
}else{
    	numwidth<-nrows*3
}
numheight<-""
numcols<-ceiling(nrows/4)
numheight<-numcols*3

numheight<-as.numeric(numheight)

DimPlot(data.combined, reduction = "umap", split.by = "sample",ncol=4)
ggsave("step3_umap_split_sample.pdf",width=numwidth,height=numheight)
ggsave("step3_umap_split_sample.png",width=numwidth,height=numheight)

class.num<-table(data.combined$sample,data.combined@meta.data$seurat_clusters)
write.table(class.num,"step3_sample_cluster_num.csv",col.names =NA,quote =FALSE,sep =",")

DefaultAssay(data.combined) <- "RNA"

data.combined.markers <- FindAllMarkers(data.combined, test.use = "MAST",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(data.combined.markers,"step3_exp_gene_diff_pos.csv",sep=",",row.names =FALSE,col.names =TRUE, quote =FALSE)

if(method == "integrated"){
	DefaultAssay(data.combined) <- "integrated"
	top20 <- data.combined.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
	p<-DoHeatmap(data.combined, features = top20$gene,label=T,size=5,angle=90,group.bar.height = 0.05,group.by ="celltype") + NoLegend()
	p<-p+theme(axis.text.y=element_blank())
	p<-p+theme(plot.margin = unit(c(4,1,1,1),"cm"))
	ggsave("step3_top20_gene_Heatmap.pdf",width=9,height=12,dpi=600)
	ggsave("step3_top20_gene_Heatmap.png",width=9,height=12,dpi=600)

}else if(method == "harmony"){
	top20 <- data.combined.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
	p<-DoHeatmap(data.combined, features = top20$gene,label=T,size=5,angle=90,group.bar.height = 0.05,group.by ="celltype") + NoLegend()
	p<-p+theme(axis.text.y=element_blank())
	p<-p+theme(plot.margin = unit(c(4,1,1,1),"cm"))
	ggsave("step3_top20_gene_Heatmap.pdf",width=9,height=12,dpi=600)
	ggsave("step3_top20_gene_Heatmap.png",width=9,height=12,dpi=600)
}
