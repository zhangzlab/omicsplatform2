library(getopt)

spec <- matrix(
  c("rds",  "r", 1, "character", "10x RDS file!",
    "comp", "c", 1, "character", "case and control of comparison!",
    "db", "d", 1, "character", "the used database!",
    "help",  "h", 0, "logical", "help document!"),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)

if( !is.null(opt$help) || (is.null(opt$rds) && is.null(opt$comp))  || is.null(opt$db)){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}

library(dplyr)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tools)
library(Seurat)
library(cowplot)
library(ggplot2)
utils::methods(class = 'Seurat')
library(stringr)

db<-"/code/c2.cp.kegg.v7.5.1.symbols.gmt"
if(opt$db == "GO"){
	db<-"/code/c5.go.v7.5.1.symbols.gmt"
}else if(opt$db == "KEGG"){
	db<-"/code/c2.cp.kegg.v7.5.1.symbols.gmt"
}else if(opt$db == "HK"){
	db<-"/code/h.all.v7.5.1.symbols.gmt"
}

all <- readRDS(file_path_as_absolute(opt$rds))
comps<-read.table(file_path_as_absolute(opt$comp),header=T,sep="\t")
nrows<-nrow(comps)

for(num in 1:nrows){
        data<-all
        ctl_data<-all
        case_data<-all

        ctl_sam<-comps[num,1]
        case_sam<-comps[num,2]

        ctl_infor<-unlist(strsplit(ctl_sam,split=";"))
        case_infor<-unlist(strsplit(case_sam,split=";"))

        ctl_len<-length(ctl_infor)
        case_len<-length(case_infor)

        tmp_ctl_name<-NULL
        tmp_case_name<-NULL

        for(order in 1:ctl_len){
                subinfor<-unlist(strsplit(ctl_infor[order],split=":"))
                ssinfor<-unlist(strsplit(subinfor[2],split=","))
                tmp_ctl_name<-append(tmp_ctl_name,ssinfor)

                if(subinfor[1] == "idents"){
                        ctl_data<-subset(ctl_data,idents = ssinfor)
                }else{
                        rnames<-rownames(ctl_data@meta.data[which(ctl_data@meta.data[[subinfor[1]]] %in% ssinfor),])
                        ctl_data<-subset(ctl_data,cells = rnames)
                }
        }
        ctl_cells<-rownames(ctl_data@meta.data)

        for(order in 1:ctl_len){
                subinfor<-unlist(strsplit(case_infor[order],split=":"))
                ssinfor<-unlist(strsplit(subinfor[2],split=","))
                tmp_case_name<-append(tmp_case_name,ssinfor)

                 if(subinfor[1] == "idents"){
                        case_data<-subset(case_data,idents = ssinfor)
                }else{
                        rnames<-rownames(case_data@meta.data[which(case_data@meta.data[[subinfor[1]]] %in% ssinfor),])
                        case_data<-subset(case_data,cells = rnames)
                }
        }
        case_cells<-rownames(case_data@meta.data)
	
	tmp_ctl_name<-str_c(tmp_ctl_name,collapse='_')
        ctl_name<-gsub(" ","_",tmp_ctl_name)
        ctl_name<-gsub("/","_",ctl_name)

        tmp_case_name<-str_c(tmp_case_name,collapse='_')
        case_name<-gsub(" ","_",tmp_case_name)
        case_name<-gsub("/","_",case_name)

        Idents(data,cells = ctl_cells) <- ctl_name
        Idents(data,cells = case_cells) <- case_name

	case_sub<-subset(data,idents=case_name)
	case_exp<-case_sub[['RNA']]@data
	case_cellexp<-as.data.frame(as.matrix(case_exp))
	
	#case_cellexp$cell<-rownames(case_cellexp)
	ctl_sub<-subset(data,idents=ctl_name)
	ctl_exp<-ctl_sub[['RNA']]@data
	ctl_cellexp<-as.data.frame(as.matrix(ctl_exp))
	#ctl_cellexp$cell<-rownames(ctl_cellexp)

	casemean<-apply(case_cellexp,1,mean)
	casemean<-as.data.frame(casemean)
	colnames(casemean)<-"mean_case"
	casemean$gene<-rownames(casemean)

	ctlmean<-apply(ctl_cellexp,1,mean)
	ctlmean<-as.data.frame(ctlmean)
	colnames(ctlmean)<-"mean_ctl"
	ctlmean$gene<-rownames(ctlmean)

	commean<-merge(casemean,ctlmean,by="gene",all=T)

	commean<-commean[which(commean$mean_case >0 |commean$mean_ctl>0),]
	commean[is.na(commean)] <- 0
	commean$log2<-round(log((commean$mean_case+1)/(commean$mean_ctl+1),2),4)
    
	#exp_file<-paste(,case_name,"-vs-",ctl_name,".gene_exp.xls",sep="")
	#write.table(commean,exp_file,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
	dbs <- read.gmt(db)
	glist <- commean$log2
	names(glist) <- toupper(as.character(commean$gene))
	glist <- sort(glist,decreasing = T)
	name<- bitr(commean$gene, fromType = "SYMBOL",toType = c("ENTREZID", "ENSEMBL"),OrgDb = org.Hs.eg.db)
	gsea_out<- GSEA(glist, TERM2GENE=dbs, verbose=FALSE, pvalueCutoff =1)
        #gsea_out$Description <-tolower(gsea_out$Description)
	num_gsea<-nrow(gsea_out)
	for(order in 1:num_gsea){
		 gsea_title<-gsea_out[order,2]
		 preout<-paste("step8_",case_name,"-vs-",ctl_name,"_",opt$db,"_",gsea_title,sep="")
		 p<-gseaplot2(gsea_out,title = gsea_title,order, color="red",base_size = 16, subplots = 1:2, pvalue_table = T)
		 ggsave(paste(preout,".pdf",sep=""),width=8,height=6)
		 ggsave(paste(preout,".png",sep=""),width=8,height=6)
	}
		 
	out_file<-paste("step8_",case_name,"-vs-",ctl_name,"_",opt$db,"_gsea.csv",sep="")
 	write.table(gsea_out,out_file,quote=FALSE,sep=",",col.names=TRUE,row.names=FALSE)
	
}
#},silent=T)
