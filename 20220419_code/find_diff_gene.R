library(getopt)

spec <- matrix(
  c("rds",  "r", 1, "character", "10x RDS file!",
    "comp", "c", 1, "character", "comparsion information",
    "help",  "h", 0, "logical", "help document!"),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$rds) || is.null(opt$comp)){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}

library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)
utils::methods(class = 'Seurat')
library(tools)
library(stringr)

all <- readRDS(file_path_as_absolute(opt$rds))

comps<-read.table(file_path_as_absolute(opt$comp),header=T,sep="\t")
nrows<-nrow(comps)

for(num in 1:nrows){
	data<-all
	ctl_data<-all
	case_data<-all

	ctl_sam<-comps[num,2]
	case_sam<-comps[num,1]

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
	
	data_diff<-FindMarkers(data,test.use = "MAST",ident.1 = case_name, ident.2 = ctl_name, verbose = FALSE, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.585)
    	data_diff<-data_diff[which(data_diff$p_val<=0.05),]
    	name<-data.frame(gene=rownames(data_diff))
    	data_diff<-cbind(name,data_diff)

    	out_name<-paste("step6_",case_name,"-vs-",ctl_name,"_diff_sign_gene.csv",sep="")
    	write.table(data_diff,file=out_name,sep=",",col.names=TRUE,quote=FALSE,row.names=FALSE)
	
}
