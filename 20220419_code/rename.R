library(getopt)

spec <- matrix(
  c("rds",  "r", 2, "character", "Support 10x RDS file!",
    "name", "n", 1, "character", "the rename information including two column: cluster and celltype!",
    "help",  "h", 0, "logical", "This is the help!"),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$rds) || is.null(opt$name)){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}

library(dplyr)
library(plyr)
library(Seurat)
library(cowplot)
utils::methods(class = 'Seurat')
library(ggplot2)
library(tools)
library(reshape2)
library(ggpubr)

data<-readRDS(file_path_as_absolute(opt$rds))
infor<-read.table(file_path_as_absolute(opt$name),header=T,sep="\t")
nrows<-nrow(infor)

infor[,1]<-as.character(infor[,1])
infor[,2]<-as.character(infor[,2])

for(num in 1:nrows){
    cells.use <- WhichCells(data, idents = infor[num,1])
    data <- SetIdent(data, cells = cells.use, value = infor[num,2])
    #data <- RenameIdents(data,infor[num,1] = infor[num,2])
}

data$celltype <- Idents(data)

p <- DimPlot(data,reduction = "umap", label = TRUE, pt.size = 1,label.size=6) + NoLegend()
p <- p+theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30),axis.title.x=element_text(size=30),axis.title.y=element_text(size=30))
ggsave("step5_rename_umap.png",width=8,height=6)
ggsave("step5_rename_umap.pdf",width=8,height=6)
saveRDS(data, file = "step5_rename_dim.Rds")

p1 <- DimPlot(data, reduction = "umap",group.by = "group")
p2 <- DimPlot(data, reduction = "umap", group.by = "sample")
p1 + p2
ggsave("step5_rename_umap_group_sample.pdf",width=12,height=6,dpi=600)
ggsave("step5_rename_umap_group_sample.png",width=12,height=6,dpi=600)

grp<-unique(data@meta.data$group)
sams<-unique(data@meta.data$sample)
ngrp<-length(grp)
nsams<-length(sams)

numwidth<-""
if(nsams >4){
    numwidth<-"12"
    numwidth<-as.numeric(numwidth)
}else{
    numwidth<-nsams*3
}
numheight<-""
numcols<-ceiling(nsams/4)
numheight<-numcols*3

numheight<-as.numeric(numheight)

DimPlot(data, reduction = "umap", split.by = "sample",ncol=4)
ggsave("step5_rename_umap_split_sample.pdf",width=numwidth,height=numheight)

class.num<-table(data$sample,data@meta.data$celltype)
write.table(class.num,"step5_rename_sample_cluster_num.csv",col.names =NA,quote =FALSE,sep =",")
group.num<-table(data$group,data@meta.data$celltype)
write.table(group.num,"step5_rename_group_cluster_num.csv",col.names =NA,quote =FALSE,sep =",")

data.markers <- FindAllMarkers(data, test.use = "MAST",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(data.markers,"step5_rename_positive_diff.csv",sep=",",row.names =FALSE, col.names =TRUE,  quote =FALSE)
top20 <- data.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
p<-DoHeatmap(data, features = top20$gene,label=T,size=4,angle=90,group.bar.height = 0.05,group.by ="celltype")+ NoLegend()
p<-p+theme(axis.text.y=element_blank())
p<-p+theme(plot.margin = unit(c(4,1,1,1),"cm"))
ggsave("step5_rename_top20_gene.pdf",width=8,height=6)
ggsave("step5_rename_top20_gene.png",width=8,height=6)

#####group ratio ######
maxwidth<-1*ngrp+4
group.num<-as.data.frame(group.num)
colnames(group.num)<-c("group","celltype","num")

p <- ggplot(group.num, aes(x =group,y = num,fill=celltype))+geom_bar(stat = 'identity', position = "fill",width=0.8)
p<-p+labs(x = "",y = "Ratio", title = "")
p<-p+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),plot.margin = unit(c(2,2,2,2), "lines"),panel.border = element_blank(),axis.line=element_line(size=1))
p<-p+theme(axis.ticks =element_line(size=0.5),axis.title.y = element_text(size=20),axis.text.x = element_text(size=20,angle=270,vjust=0.5,hjust=0),axis.text.y =element_text(size=20),plot.title=element_text(size=20,hjust = 0.5))
p<-p+ theme(strip.text = element_text(face = 'bold', size = 10), strip.background = element_rect(size = 1))
p<-p+guides(fill = guide_legend(reverse = F))
p<-p+theme(legend.title = element_blank(),legend.text = element_text(size = 14), legend.key.size=unit(0.5,'cm'),legend.position = "right",legend.direction="vertical",legend.spacing.x = unit(0.5, 'cm'),legend.key.width = unit(0.5, "cm"))
ggsave("step5_rename_group_cell_ratio.pdf",width=maxwidth,height=8)
ggsave("step5_rename_group_cell_ratio.png",width=maxwidth,height=8)

######sample ratio ######
maxwidth<-1*nsams+4
class.num<-as.data.frame(class.num)
colnames(class.num)<-c("sample","celltype","num")
p<-ggplot(class.num, aes(x =sample,y = num,fill=celltype))+geom_bar(stat = 'identity', position = "fill",width=0.8)
p<-p+labs(x = "",y = "Ratio", title = "")
p<-p+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),plot.margin = unit(c(2,2,2,2), "lines"),panel.border = element_blank(),axis.line=element_line(size=1))
p<-p+theme(axis.ticks =element_line(size=0.5),axis.title.y = element_text(size=20),axis.text.x = element_text(size=20,angle=270,vjust=0.5,hjust=0),axis.text.y =element_text(size=20),plot.title=element_text(size=20,hjust = 0.5))
p<-p+theme(strip.text = element_text(face = 'bold', size = 10), strip.background = element_rect(size = 1))
p<-p+guides(fill = guide_legend(reverse = F))
p<-p+theme(legend.title = element_blank(),legend.text = element_text(size = 14), legend.key.size=unit(0.5,'cm'),legend.position = "right",legend.direction="vertical",legend.spacing.x = unit(0.5, 'cm'),legend.key.width = unit(0.5, "cm"))
ggsave("step5_rename_sample_cell_ratio.pdf",width=maxwidth,height=8)
ggsave("step5_rename_sample_cell_ratio.png",width=maxwidth,height=8)

###### sample cell ######
metas<-cbind(as.data.frame(rownames(data@meta.data)),data@meta.data$sample,data@meta.data$group,data@meta.data$celltype)

colnames(metas)<-c("cell","sample","group","celltype")
cellsum<-ddply(metas,.(sample,celltype),function(subss){data.frame(num=nrow(subss))})
write.table(cellsum,"step5_rename_sample_celltype_number.csv",sep=",",quote=FALSE)

samsum<-ddply(metas,.(sample),function(subst){data.frame(sum=nrow(subst))})
write.table(cellsum,"step5_rename_sample_celltype_number.csv",sep=",",quote=FALSE)

samgr<-cbind(as.data.frame(data@meta.data$sample),as.data.frame(data@meta.data$group))
colnames(samgr)<-c("sample","group")
samgr<-unique(samgr)
comb<-merge(cellsum,samsum,by="sample")
comb<-merge(comb,samgr,by="sample")
comb$ratio<-round(comb$num/comb$sum,4)
write.table(comb,"step5_rename_sample_cell_ratio.csv",sep=",",quote=FALSE)

two_com<-t(combn(grp,2))
ntwo_com<-nrow(two_com)

my_comparisons<-split(t(combn(grp,2)), 1:nrow(t(combn(grp,2))))

metass<-unique(comb$celltype)

for( subcell in metass){

    subdata<-comb[which(comb$celltype==subcell),]
    maxratio<-max(subdata$ratio)
    #maxwidth<-1*ngrp+4

    yvalue<-c(maxratio+0.025)
    for(order in 2:ntwo_com ){
        yvalue<-append(yvalue,maxratio+0.025*order)
    }
    p <- ggplot(subdata, aes(x =group,y = ratio,fill=group))+geom_boxplot(outlier.size=NA, outlier.shape = NA,width=0.7,lwd=0.9)
    p<-p+geom_jitter(width = 0.3, size=2)
    p<-p+stat_compare_means(comparisons = my_comparisons,bracket.size=0.5,tip.length=0.015,label="p.signif",label.y.npc=1,label.y =yvalue,hide.ns = TRUE)
    p<-p+theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),plot.margin = unit(c(2,2,2,2), "lines"),panel.border = element_blank(),axis.line=element_line(size=1))
    p<-p+theme(axis.ticks =element_line(size=0.5),axis.title.y = element_text(size=14),axis.text.x = element_text(size=14,angle=270,vjust=0.5,hjust=0),axis.text.y =element_text(size=14),plot.title=element_text(size=14,hjust = 0.5))
    #p<-p+theme(strip.text = element_text(face = 'bold', size = 10), strip.background = element_rect(size = 0.8))
    p<-p+labs(x="",y="",title=subcell)
    p<-p+guides(fill = FALSE)

    tmp_cell<-gsub(" ","_",subcell)
    ggsave(paste("step5_rename_",tmp_cell,"_cell_ratio_stat.pdf",sep=""),width=5,height=5)
    ggsave(paste("step5_rename_",tmp_cell,"_cell_ratio_stat.png",sep=""),width=5,height=5)
}
