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

data<-read.table(file_path_as_absolute(opt$diff),sep=",",header=T)

data$gene<-toupper(data$gene)
gene<- bitr(data$gene, fromType = "SYMBOL",toType = c("ENTREZID", "ENSEMBL"),OrgDb = org.Hs.eg.db)

ALL <- enrichGO(gene$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',pvalueCutoff  = 1,pAdjustMethod = "BH", qvalueCutoff  = 1, readable=T)
BP <- enrichGO(gene$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff  = 1,pAdjustMethod = "BH",qvalueCutoff  = 1, readable=T)
MF <- enrichGO(gene$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff  = 1,pAdjustMethod = "BH",qvalueCutoff  = 1, readable=T)
CC <- enrichGO(gene$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff  = 1,pAdjustMethod = "BH",qvalueCutoff  = 1, readable=T)

write.table(as.data.frame(ALL@result), file="step7_all_GO.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=FALSE)

CC_simp <- simplify(CC, cutoff=0.7, by="p.adjust",select_fun=min)
BP_simp <- simplify(BP, cutoff=0.7, by="p.adjust",select_fun=min)
MF_simp<- simplify(MF, cutoff=0.7, by="p.adjust",select_fun=min)

write.table(as.data.frame(CC_simp@result),file="step7_CC_GO_simp.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=FALSE)
write.table(as.data.frame(BP_simp@result),file="step7_BP_GO_simp.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=FALSE)
write.table(as.data.frame(MF_simp@result),file="step7_MF_GO_simp.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=FALSE)

kegg<-enrichKEGG(gene$ENTREZID, organism ="hsa",pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 1,use_internal_data =FALSE)

write.table(as.data.frame(kegg@result), file="step7_KEGG.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=FALSE)
