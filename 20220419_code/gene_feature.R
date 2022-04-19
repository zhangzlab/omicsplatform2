library(getopt)

spec <- matrix(
  c("rds",  "r", 1, "character", "10x RDS file!",
    "gene", "g", 1, "character", "gene list file!",
    "infor","i", 2, "complex", "information of cell",
    "help",  "h", 0, "logical", "help document!"),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)

if( !is.null(opt$help) || (is.null(opt$rds) && is.null(opt$gene))){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}


library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)
utils::methods(class = 'Seurat')
library(tools)

data <- readRDS(file_path_as_absolute(opt$rds))
genes <-read.table(file_path_as_absolute(opt$gene),header=T,sep="\t")

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

cell_num<-length(unique(data@meta.data$celltype))
ngenes<-nrow(genes)
nwidths<-cell_num*0.8+2
nheights<-ngenes*0.8+1

for(num in 1:ngenes){

	gene_name<-genes[num,1]
	FeaturePlot(data,features = genes[num,1],pt.size=0.6, raster=FALSE)
	ggsave(paste("step4_",gene_name,"_exp_features.pdf",sep=""),width=8,height=6)
	ggsave(paste("step4_",gene_name,"_exp_features.png",sep=""),width=8,height=6)

	VlnPlot(data, features =genes[num,1],group.by='celltype')
	ggsave(paste("step4_",gene_name,"_exp_vlnplot.pdf",sep=""),width=nwidths,height=5)
	ggsave(paste("step4_",gene_name,"_exp_vlnplot.png",sep=""),width=nwidths,height=5)
}

p<-DotPlot(data, features = genes[,1])+coord_flip()
p<-p+theme(axis.text.x=element_text(size=14,angle=45,hjust=1,vjust=1),axis.text.y=element_text(size=14),axis.title.x=element_text(),axis.title.y=element_text())
ggsave("step4_all_genes_exp_point_features.pdf",width=nwidths,height=nheights)
ggsave("step4_all_genes_exp_point_features.png",width=nwidths,height=nheights)
