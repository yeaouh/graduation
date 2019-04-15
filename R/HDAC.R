options(stringsAsFactors = FALSE)
library(PharmacoGx)
library(VennDiagram)
library(RColorBrewer)
library(R.utils)
library(downloader)
library(gdata)
library(AnnotationDbi)
library(hgu133a.db)
library(gridExtra)
library(biomaRt)
nbcore <- 4

#download and process HDAC signature
mydir <- "1132939s"
if(!file.exists(mydir)){
  downloader::download(paste("http://sciencemag.org/content/suppl/2006/09/29/313.5795.1929.DC1/",
                             mydir,".zip",sep=""),destfile=paste(mydir,".zip",sep=""))
  unzip(paste(mydir,".zip",sep=""))
}
perl<-"C:\\Users\\yuki\\Downloads\\mingw\\bin\\Perl\\bin\\perl.exe"
HDAC_up<-read.xls(paste(mydir,paste(mydir,"sigS1.xls",sep="_"),sep="/"),sheet = 1,header=FALSE,as.is=TRUE,perl = perl)
HDAC_down<-read.xls(paste(mydir,paste(mydir,"sigS1.xls",sep="_"),sep="/"),sheet = 2,header=FALSE,as.is=TRUE,perl = perl)
HDAC<-as.data.frame(matrix(NA,nrow = nrow(HDAC_down)+nrow(HDAC_up),ncol = 2))
annot<-AnnotationDbi::select(hgu133a.db,keys = c(HDAC_up[[1]],HDAC_down[[1]]),columns = c("ENSEMBL"),keytype = "PROBEID")
gene_up<-unique(annot[match(HDAC_up[[1]],annot[,1]),2])
gene_down<-na.omit(unique(annot[match(HDAC_down[[1]],annot[,1]),2]))
HDAC_genes<-as.data.frame(matrix(NA,nrow = length(gene_down)+length(gene_up),ncol = 2))
#ensembl<-useMart("ensembl",dataset = "hsapiens_gene_ensembl")
#entrzID=c(gene_up,gene_down)
#genesymbol<-getBM(attributes = c("entrezgene","hgnc_symbol","ensembl_gene_id"),filters = "entrezgene",values = entrzID,mart = ensembl)
genesymbol<-c(gene_up,gene_down)
HDAC_genes[,1]<-paste(genesymbol,"at",sep = "_")
HDAC_genes[,2]<-c(rep(1,times=length(gene_up)),rep(-1,times=length(gene_down)))
rownames(HDAC_genes)<-HDAC_genes[ ,1]
HDAC<-HDAC_genes[ ,2]
names(HDAC)<-rownames(HDAC_genes)

##download/compute drug perturbation signatures from CMAP
drug.perturbation<-PharmacoGx::downloadPertSig("CMAP")
myfn<-"CMAP_genes_HDAC_connectivity.RData"
if(!file.exists(myfn)){
  message("Be aware that computing sensitivity will take some time...")
  cl<-parallel::makeCluster(nbcore)
  res<-parApply(drug.perturbation[ , ,c("tstat","fdr")],2,function(x,HDAC){return(PharmacoGx::connectivityScore(x=x,y=HDAC,method = "gsea",nperm = 100))},cl=cl,HDAC=HDAC)
  stopCluster(cl)
  rownames(res)<-c("Connectivity","P_Value")
  res<-t(res)
  save(res,file=myfn)
}else{
  load(myfn)
}

HDAC_inhibitors<-c("vorinostat","trichostatin A","HC toxin","valproic acid")
res<-res[order(res[,1],decreasing = T),]
my_ranks<-which(rownames(res) %in% HDAC_inhibitors)
total<-dim(drug.perturbation)[2]
ranks<-my_ranks
positions<-total-ranks+1

pdf("CMAP_Case_Study_all_drugs.pdf")
plot.new()
rect(0,0,.2,1,col = "lightgrey",border = FALSE)
for (position in positions) {
  rect(0,1/total*(position),.2,1/total*(position-1),col = "red",border = FALSE)
}
position=total-50
rect(0,1/total*(position),.2,1/total*(position),col = "black")
dev.off()

total<-50
ranks<-my_ranks
positions<-total-ranks+1

pdf("CMAP_Case_Study_all_drugs_top_50.pdf")
plot.new()
rect(0,0,1,1,col = "lightgrey",border = FALSE)
for (position in positions) {
  rect(0,1/total*(position),1,1/total*(position-1),col = "red",border = FALSE)
}
dev.off()

drugs_v1<-unique(read.xls(paste(mydir,paste(mydir,"tableS1.xls",sep = "_"),sep = "/"),sheet = 1,stringsAsFactors=FALSE)[,"cmap_name"])
drugs_v1<-drugs_v1[!drugs_v1%in%c("",NA)]
int_drugs<-intersect(drugs_v1,rownames(res))
res_v1<-res[int_drugs,]
res_v1<-res_v1[order(res_v1[,1],decreasing = T),]
colnames(res_v1)<-c("Connectivity","P_Value")
my_ranks<-which(rownames(res_v1)%in%HDAC_inhibitors)
total<-nrow(res_v1)
ranks<-my_ranks
positions<-total-ranks+1

pdf("CMAP_Case_Study_v1_drugs.pdf")
plot.new()
rect(0,0,.2,1,col = "lightgrey",border = FALSE)
for (position in positions) {
  rect(0,1/total*(position),.2,1/total*(position-1),col = "red",border = FALSE)
}
dev.off()

pdf("HDAC_res_table.pdf")
grid.table(res[1:20, ])
dev.off()

pdf("HDAC_res_v1_table.pdf")
grid.table(res_v1[1:20, ])
dev.off()
