

nbcore<-48
ket<-drug.perturbation[,"ketoconazole",c("tstat","pvalue")]
myfn1<-"CMAP_genes_ketoconazole_connectivity.RData"
if(!file.exists(myfn1)){
  message("Be aware that computing sensitivity will take some time...")
  cl<-parallel::makeCluster(nbcore)
  res1<-parApply(drug.perturbation[ , ,c("tstat","pvalue")],2,function(x,ket){return(PharmacoGx::connectivityScore(x=x,y=ket,method = "gwc",nperm = 100))},cl=cl,ket=ket)
  stopCluster(cl)
  rownames(res1)<-c("Connectivity","P_Value")
  res1<-t(res1)
  save(res1,file=myfn1)
}else{
  load(myfn1)
}
res1<-res1[order(res1[,1],decreasing = T),]

mic<-drug.perturbation[,"miconazole",c("tstat","pvalue")]
myfn2<-"CMAP_genes_miconazole_connectivity.RData"
if(!file.exists(myfn2)){
  message("Be aware that computing sensitivity will take some time...")
  cl<-parallel::makeCluster(nbcore)
  res2<-parApply(drug.perturbation[ , ,c("tstat","pvalue")],2,function(x,mic){return(PharmacoGx::connectivityScore(x=x,y=mic,method = "gwc",nperm = 100))},cl=cl,mic=mic)
  stopCluster(cl)
  rownames(res2)<-c("Connectivity","P_Value")
  res2<-t(res2)
  save(res2,file=myfn2)
}else{
  load(myfn2)
}
res2<-res2[order(res2[,1],decreasing = T),]

amB<-drug.perturbation[,"amphotericin B",c("tstat","pvalue")]
myfn3<-"CMAP_amB_connectivity.RData"
if(!file.exists(myfn3)){
  message("Be aware that computing sensitivity will take some time...")
  cl<-parallel::makeCluster(nbcore)
  res3<-parApply(drug.perturbation[ , ,c("tstat","pvalue")],2,function(x,amB){return(PharmacoGx::connectivityScore(x=x,y=amB,method = "gwc",nperm = 100))},cl=cl,amB=amB)
  stopCluster(cl)
  rownames(res3)<-c("Connectivity","P_Value")
  res3<-t(res3)
  save(res3,file=myfn3)
}else{
  load(myfn3)
}
res3<-res3[order(res3[,1],decreasing = T),]

nys<-drug.perturbation[,"nystatin",c("tstat","pvalue")]
myfn4<-"CMAP_nystatin_connectivity.RData"
if(!file.exists(myfn4)){
  message("Be aware that computing sensitivity will take some time...")
  cl<-parallel::makeCluster(nbcore)
  res4<-parApply(drug.perturbation[ , ,c("tstat","pvalue")],2,function(x,nys){return(PharmacoGx::connectivityScore(x=x,y=nys,method = "gwc",nperm = 100))},cl=cl,nys=nys)
  stopCluster(cl)
  rownames(res4)<-c("Connectivity","P_Value")
  res4<-t(res4)
  save(res4,file=myfn4)
}else{
  load(myfn4)
}
res4<-res4[order(res4[,1],decreasing = T),]

cop<-drug.perturbation[,"copper sulfate",c("tstat","pvalue")]
myfn5<-"CMAP_cop_connectivity.RData"
if(!file.exists(myfn5)){
  message("Be aware that computing sensitivity will take some time...")
  cl<-parallel::makeCluster(nbcore)
  res5<-parApply(drug.perturbation[ , ,c("tstat","pvalue")],2,function(x,cop){return(PharmacoGx::connectivityScore(x=x,y=cop,method = "gwc",nperm = 100))},cl=cl,cop=cop)
  stopCluster(cl)
  rownames(res5)<-c("Connectivity","P_Value")
  res5<-t(res5)
  save(res5,file=myfn5)
}else{
  load(myfn5)
}
res5<-res5[order(res5[,1],decreasing = T),]

pdf("ketoconazole_table.pdf")
grid.table(res1[1:20, ])
dev.off()

pdf("miconazole_table.pdf")
grid.table(res2[1:20, ])
dev.off()

pdf("amphotericin B_table.pdf")
grid.table(res3[1:20, ])
dev.off()

pdf("nystatin_table.pdf")
grid.table(res4[1:20, ])
dev.off()

pdf("copper_table.pdf")
grid.table(res5[1:20, ])
dev.off()