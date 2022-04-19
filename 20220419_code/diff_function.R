library(getopt)

spec <- matrix(
  c("diff",  "f", 1, "character", "differential expression genes, one column should be named \"gene\".",
    "help",  "h", 0, "logical", "help document."),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$diff)){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}

library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tools)

data<-read.table(file_path_as_absolute(opt$diff),sep="\t",header=T)

data$gene<-toupper(data$gene)
gene<- bitr(data$gene, fromType = "SYMBOL",toType = c("ENTREZID", "ENSEMBL"),OrgDb = org.Hs.eg.db)

ALL <- enrichGO(gene$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',pvalueCutoff  = 1,pAdjustMethod = "BH", qvalueCutoff  = 1, readable=T)
BP <- enrichGO(gene$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff  = 1,pAdjustMethod = "BH",qvalueCutoff  = 1, readable=T)
MF <- enrichGO(gene$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff  = 1,pAdjustMethod = "BH",qvalueCutoff  = 1, readable=T)
CC <- enrichGO(gene$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff  = 1,pAdjustMethod = "BH",qvalueCutoff  = 1, readable=T)

write.table(as.data.frame(ALL@result), file="step7_all_GO.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=FALSE)

CC_simp <- simplify(CC, cutoff=0.7, by="p.adjust",select_fun=min)
BP_simp <- simplify(BP, cutoff=0.7, by="p.adjust",select_fun=min)
MF_simp <- simplify(MF, cutoff=0.7, by="p.adjust",select_fun=min)

CC_simps <- as.data.frame(CC_simp@result)
BP_simps <- as.data.frame(BP_simp@result)
MF_simps <- as.data.frame(MF_simp@result)

CC_sign<-CC_simps[which(CC_simps$p.adjust<=0.05),]
CC_sign$p.adjust<-as.numeric(CC_sign$p.adjust)
maxvalue<-max(-1*log10(CC_sign$p.adjust))

num_cc<-nrow(CC_sign)
height_cc<-0.34*num_cc+1
CC_sign$Description<-factor(CC_sign$Description,levels=CC_sign$Description)
p <- ggplot(CC_sign,aes(x=Description,y=-1*log10(CC_sign$p.adjust),fill="#FF9912"))+ geom_bar(stat="identity",width = 0.8,fill="#FF9912")
p <- p + coord_flip()
p <- p + labs(y=expression(-log[10](p.adjust)),x="",title="CC")
p <- p + theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),plot.margin = unit(c(2,3,2,1), "lines"))
p <- p + theme(panel.background=element_rect(fill="white",colour="black",size=1))
p <- p + theme(axis.ticks =element_line(size=1),axis.title.x = element_text(size=16),axis.title.y = element_text(size=12),axis.text.y = element_text(size=10),axis.text.x =element_text(size=12,angle=45,hjust=1,vjust=1),plot.title=element_text(size=18,hjust = 0.5))
p <- p + scale_y_continuous(limits =c(0,maxvalue+1),breaks=seq(0,maxvalue+1,2),labels=seq(0,maxvalue+1,2))
p <- p + guides(fill=FALSE)
ggsave("step7_CC_function.pdf",width=15,height=height_cc,useDingbats=FALSE)
ggsave("step7_CC_function.png",width=15,height=height_cc)

BP_sign<-BP_simps[which(BP_simps$p.adjust<=0.05),]
BP_sign$p.adjust<-as.numeric(BP_sign$p.adjust)
num_bp<-nrow(BP_sign)
height_bp<-0.34*num_bp+1
maxvalue<-max(-1*log10(BP_sign$p.adjust))
BP_sign$Description<-factor(BP_sign$Description,levels=BP_sign$Description)
p <- ggplot(BP_sign,aes(x=Description,y=-1*log10(BP_sign$p.adjust),fill="#1E90FF"))+ geom_bar(stat="identity",width = 0.8,fill="#1E90FF")
p <- p + coord_flip()
p <- p + labs(y=expression(-log[10](p.adjust)),x="",title="BP")
p <- p + theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),plot.margin = unit(c(2,3,2,1), "lines"))
p <- p + theme(panel.background=element_rect(fill="white",colour="black",size=1))
#p <- p + theme(axis.ticks =element_line(size=1),axis.title.x = element_text(size=16),axis.title.y = element_text(size=12),axis.text.x = element_text(size=14),axis.text.y =element_text(size=10,angle=45,hjust=1,vjust=1),plot.title=element_text(size=18,hjust = 0.5))
p <- p + theme(axis.ticks =element_line(size=1),axis.title.x = element_text(size=16),axis.title.y = element_text(size=12),axis.text.y = element_text(size=10),axis.text.x =element_text(size=12,angle=45,hjust=1,vjust=1),plot.title=element_text(size=18,hjust = 0.5))
p <- p + scale_y_continuous(limits =c(0,maxvalue+1),breaks=seq(0,maxvalue+1,2),labels=seq(0,maxvalue+1,2))
p <- p + guides(fill=FALSE)
ggsave("step7_BP_function.pdf",width=15,height=height_bp,useDingbats=FALSE)
ggsave("step7_BP_function.png",width=15,height=height_bp)

MF_sign<-MF_simps[which(MF_simps$p.adjust<=0.05),]
MF_sign$p.adjust<-as.numeric(MF_sign$p.adjust)
num_mf<-nrow(MF_sign)
height_mf<-0.34*num_mf+1
maxvalue<-max(-1*log10(MF_sign$p.adjust))

MF_sign$Description<-factor(MF_sign$Description,levels=MF_sign$Description)
p <- ggplot(MF_sign,aes(x=Description,y=-1*log10(MF_sign$p.adjust),fill="#FF6347"))+ geom_bar(stat="identity",width = 0.8,fill="#FF6347")
p <- p + coord_flip()
p <- p + labs(y=expression(-log[10](p.adjust)),x="",title="MF")
p <- p + theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),plot.margin = unit(c(2,3,2,1), "lines"))
p <- p + theme(panel.background=element_rect(fill="white",colour="black",size=1))
#p <- p + theme(axis.ticks =element_line(size=1),axis.title.x = element_text(size=16),axis.title.y = element_text(size=12),axis.text.x = element_text(size=14),axis.text.y =element_text(size=10,angle=45,hjust=1,vjust=1),plot.title=element_text(size=18,hjust = 0.5))
p <- p + theme(axis.ticks =element_line(size=1),axis.title.x = element_text(size=16),axis.title.y = element_text(size=12),axis.text.y = element_text(size=10),axis.text.x =element_text(size=12,angle=45,hjust=1,vjust=1),plot.title=element_text(size=18,hjust = 0.5))
p <- p + scale_y_continuous(limits =c(0,maxvalue+1),breaks=seq(0,maxvalue+1,2),labels=seq(0,maxvalue+1,2))
p <- p + guides(fill=FALSE)
ggsave("step7_MF_function.pdf",width=15,height=height_mf,useDingbats=FALSE)
ggsave("step7_MF_function.png",width=15,height=height_mf)

write.table(as.data.frame(CC_simp@result),file="step7_CC_GO_simp.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=FALSE)
write.table(as.data.frame(BP_simp@result),file="step7_BP_GO_simp.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=FALSE)
write.table(as.data.frame(MF_simp@result),file="step7_MF_GO_simp.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=FALSE)

keggs<-enrichKEGG(gene$ENTREZID, organism ="hsa",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 1,use_internal_data =FALSE)

kegg<-as.data.frame(keggs@result)
kegg_sign<-kegg[which(kegg$p.adjust<=0.05),]
kegg_sign$p.adjust<-as.numeric(kegg_sign$p.adjust)
num_kegg<-nrow(kegg_sign)
height_kegg<-0.34*num_kegg+1
maxvalue<-max(-1*log10(kegg_sign$p.adjust))

kegg_sign$Description<-factor(kegg_sign$Description,levels=kegg_sign$Description)
p <- ggplot(kegg_sign,aes(x=Description,y=-1*log10(kegg_sign$p.adjust),fill="#FF6103"))+ geom_bar(stat="identity",width = 0.8,fill="#FF6103")
p <- p + coord_flip()
p <- p + labs(y=expression(-log[10](p.adjust)),x="",title="kegg")
p <- p + theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),plot.margin = unit(c(2,3,2,1), "lines"))
p <- p + theme(panel.background=element_rect(fill="white",colour="black",size=1))
#p <- p + theme(axis.ticks =element_line(size=1),axis.title.x = element_text(size=16),axis.title.y = element_text(size=12),axis.text.x = element_text(size=14),axis.text.y =element_text(size=10,angle=45,hjust=1,vjust=1),plot.title=element_text(size=18,hjust = 0.5))
p <- p + theme(axis.ticks =element_line(size=1),axis.title.x = element_text(size=16),axis.title.y = element_text(size=12),axis.text.y = element_text(size=10),axis.text.x =element_text(size=12,angle=45,hjust=1,vjust=1),plot.title=element_text(size=18,hjust = 0.5))
p <- p + scale_y_continuous(limits =c(0,maxvalue+1),breaks=seq(0,maxvalue+1,2),labels=seq(0,maxvalue+1,2))
p <- p + guides(fill=FALSE)
ggsave("step7_kegg_function.pdf",width=15,height=height_kegg,useDingbats=FALSE)
ggsave("step7_kegg_function.png",width=15,height=height_kegg)

write.table(as.data.frame(keggs@result), file="step7_KEGG.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=FALSE)
