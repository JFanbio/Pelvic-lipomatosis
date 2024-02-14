library(Seurat)
library(dplyr)
library(cowplot)
library(sctransform)
setwd("/gpfs/share/home/2116392061/01Pelvic_Lipomatosis/04Singlecelltype_all/01Endo")
load("/gpfs/share/home/2116392061/01Pelvic_Lipomatosis/03Identification_all/alldata_rename.Rdata")
DefaultAssay(alldata)<-"integrated"
alldata$ident<-alldata@active.ident
Endo=subset(alldata,subset = ident=="Endothelial")

Endo <- RunPCA(Endo,npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.4, verbose = FALSE)
p1 <- DimPlot(Endo, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(Endo, reduction = "umap", label = TRUE)
p1 | p2


Endo$ident<-Endo@active.ident
DefaultAssay(Endo) <- "RNA"
for(i in unique(Endo$ident)){
  markers <- FindMarkers(Endo, ident.1 = i, min.pct = 0.25)
  markers<-cbind(rownames(markers),markers)
  colnames(markers)[1]<-"gene"
  write.table(markers,paste0("cluster",i,".markers.txt"),sep = "\t",row.names=FALSE)
}


for(i in 0:(length(unique(Endo$ident))-1)){
  markers <- FindMarkers(Endo, ident.1 = "i", min.pct = 0.25)
  markers<-cbind(rownames(markers),markers)
  colnames(markers)[1]<-"gene"
  write.table(markers,paste0("cluster",i,".markers.txt"),sep = "\t",row.names=FALSE)
}

DefaultAssay(Endo) <- "RNA"
pdf("PPAR.pdf",width =5,height = 15)
plots <-VlnPlot(Endo, features = c("CD36","AQP7","DBI","FABP4","LPL","PLIN2","PPARG")
                ,pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()


save(Endo,file = "Endo.Rdata")



FeaturePlot(Endo, features = c("HSPA1A"),cols =colorRampPalette(c("lightgrey","yellow","red","red"))(100)) 


####  vein
pdf("ident_vein1.pdf",width =10,height = 9)
FeaturePlot(Endo, features = c("ACKR1","SELP","IL1R1"),cols =colorRampPalette(c("lightgrey","yellow","red","red"))(100),max.cutoff = c(50,5,10)) 
dev.off()
pdf("ident_vein2.pdf",width =10,height = 15)
plots <-VlnPlot(Endo, features = c("ACKR1","SELP","IL1R1")
                ,pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()

####  artery
pdf("ident_artery1.pdf",width =10,height = 9)
FeaturePlot(Endo, features = c("GJA5","ARL15","IGFBP3","HEY1"),cols =colorRampPalette(c("lightgrey","yellow","red","red"))(100),max.cutoff = c(10,5,25,10)) 
dev.off()
pdf("ident_artery2.pdf",width =10,height = 20)
plots <-VlnPlot(Endo, features = c("GJA5","ARL15","IGFBP3","HEY1")
                ,pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()

### capillary
pdf("ident_capillary1.pdf",width =10,height = 9)
FeaturePlot(Endo, features = c("RGCC","LPL","CD300LG","CA4"),cols =colorRampPalette(c("lightgrey","yellow","red","red"))(100),max.cutoff = c(20,50,20,30)) 
dev.off()
pdf("ident_capillary2.pdf",width =10,height = 20)
plots <-VlnPlot(Endo, features = c("RGCC","LPL","CD300LG","CA4")
                ,pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()




#### cell percentage
library(Seurat)
library(dplyr)
library(cowplot)
library(sctransform)
setwd("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/01Endo/")
load("Endo_rename.Rdata")
DefaultAssay(Endo)<-"RNA"
a <- subset(Endo, orig.ident =="CWL_Case")
b=as.matrix(table(a@active.ident))
a <- subset(Endo, orig.ident =="CWL_Control")
b=cbind(b,as.matrix(table(a@active.ident)))
a <- subset(Endo, orig.ident =="HB_Case")
b=cbind(b,as.matrix(table(a@active.ident)))
a <- subset(Endo, orig.ident =="LZF_Case")
b=cbind(b,as.matrix(table(a@active.ident)))
a <- subset(Endo, orig.ident =="LZF_Control")
b=cbind(b,as.matrix(table(a@active.ident)))
a <- subset(Endo, orig.ident =="ZX_Case")
b=cbind(b,as.matrix(table(a@active.ident)))
a <- subset(Endo, orig.ident =="ZX_Control")
b=cbind(b,as.matrix(table(a@active.ident)))
a <- subset(Endo, orig.ident =="ZZL_Case")
b=cbind(b,as.matrix(table(a@active.ident)))
a <- subset(Endo, orig.ident =="ZZL_Control")
b=cbind(b,as.matrix(table(a@active.ident)))
colnames(b)<-c("CWL_Case","CWL_Control","HB_Case","LZF_Case","LZF_Control","ZX_Case","ZX_Control","ZZL_Case","ZZL_Control")
b<-t(b)
celltype<-rep(colnames(b),each=9)
samples<-rep(rownames(b),3)
Num<-c(b[,1],b[,2],b[,3])
Individual<-c()
Status<-c()
for(i in samples){
  Individual<-c(Individual,strsplit(i,"_")[[1]][1])
  Status<-c(Status,strsplit(i,"_")[[1]][2])
}
b<-data.frame(celltype,samples,Num,Individual,Status)
setwd("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/01Endo/Cell_percentage")
write.table(b,"Cell_percentage.txt",sep = "\t",row.names = F)
