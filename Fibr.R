library(Seurat)
library(dplyr)
library(cowplot)
library(sctransform)
setwd("/gpfs/share/home/2116392061/01Pelvic_Lipomatosis/04Singlecelltype_all/06Fibr")
load("/gpfs/share/home/2116392061/01Pelvic_Lipomatosis/03Identification_all/alldata_rename.Rdata")
DefaultAssay(alldata)<-"integrated"
alldata$ident<-alldata@active.ident
Fibr=subset(alldata,subset = ident=="Pericyte" |ident=="Smooth muscle")

Fibr <- RunPCA(Fibr,npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.8, verbose = FALSE)

save(Fibr,file = "Fibr.Rdata")
pdf("metadata.pdf",width =10,height = 12)
plots <-VlnPlot(Fibr, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo")
                ,pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()

p1 <- DimPlot(Fibr, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(Fibr, reduction = "umap", label = TRUE)
pdf("UMAP1.pdf",width = 12,height = 5)
p1 | p2
dev.off()

pdf("UMAP2.pdf",width = 15,height = 12)
DimPlot(Fibr, reduction = "umap", split.by = "orig.ident",label = TRUE,ncol = 5)
dev.off()

pdf("UMAP3_status.pdf",width = 12,height = 5)
DimPlot(Fibr, reduction = "umap", split.by = "status",label = TRUE)
dev.off()

pdf("UMAP4_library.pdf",width = 12,height = 5)
DimPlot(Fibr, reduction = "umap", split.by = "library",label = TRUE)
dev.off()

pdf("UMAP5_individual.pdf",width = 12,height = 5)
DimPlot(Fibr, reduction = "umap", split.by = "individual",label = TRUE)
dev.off()


p1 <- DimPlot(Fibr, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(Fibr, reduction = "tsne", label = TRUE)
pdf("TSNE1.pdf",width = 12,height = 5)
p1 | p2
dev.off()

pdf("TSNE2.pdf",width = 15,height = 10)
DimPlot(Fibr, reduction = "tsne", split.by = "orig.ident",label = TRUE,ncol = 6)
dev.off()

pdf("TSNE3_status.pdf",width = 12,height = 5)
DimPlot(Fibr, reduction = "tsne", split.by = "status",label = TRUE)
dev.off()

pdf("TSNE4_library.pdf",width = 12,height = 5)
DimPlot(Fibr, reduction = "tsne", split.by = "library",label = TRUE)
dev.off()

pdf("TSNE5_individual.pdf",width = 12,height = 5)
DimPlot(Fibr, reduction = "tsne", split.by = "individual",label = TRUE)
dev.off()

Fibr=subset(Fibr,subset = seurat_clusters!="11" & seurat_clusters!="13")
save(Fibr,file = "Fibr_new.Rdata")


load("Fibr.Rdata")
DefaultAssay(Fibr) <- "RNA"
pdf("tmp1.pdf",width =6,height = 7)
plots <-VlnPlot(Fibr, features = c("STEAP4","MYOCD","COL1A1","ACTA2","THY1")
                ,pt.size = 0, combine = FALSE,split.by = "status")
CombinePlots(plots = plots, ncol = 1)
dev.off()


Fibr <- RenameIdents(Fibr,
                    '0' = "SM_C1",
                    '1' = "SM_C2",
                    '2' = "Per_C1",
                    '3' = "SM_C3",
                    '4' = "SM_C4",
                    '5' = "Per_C2",
                    '6' = "SM_C5",
                    '7' = "Per_C3",
                    '8' = "SM_C4",
                    '9' = "SM_C6")
Fibr$ident1<-Fibr@active.ident
Fibr@active.ident<-Fibr$seurat_clusters
Fibr <- RenameIdents(Fibr,
                     '0' = "SM",
                     '1' = "SM",
                     '2' = "Per",
                     '3' = "SM",
                     '4' = "SM",
                     '5' = "Per",
                     '6' = "SM",
                     '7' = "Per",
                     '8' = "SM",
                     '9' = "SM")
Fibr$ident2<-Fibr@active.ident
Fibr@active.ident<-Fibr$ident1
save(Fibr,file="Fibr_rename.Rdata")
library(ggsci)
library(ggplot2)
col=colorRampPalette(pal_npg()(7))(9)
Fibr$ident1<-factor(Fibr@active.ident,levels=c("Per_C1","Per_C2","Per_C3","SM_C1","SM_C2","SM_C3","SM_C4","SM_C5","SM_C6"))
Idents(Fibr)<-"ident1"
pdf("Identification_TSNE_Fibr.pdf",width =7,height = 6)
DimPlot(Fibr,cols = col,reduction = "tsne",label = TRUE)+theme_classic()
dev.off()
pdf("Identification_UMAP_Fibr.pdf",width =7,height = 6)
DimPlot(Fibr,cols = col,reduction = "umap",label = TRUE)+theme_classic()
dev.off()
save(Fibr,file="Fibr_rename.Rdata")


Fibr@active.ident<-Fibr$ident2
col=colorRampPalette(pal_npg()(7))(9)
col=c("#739FAD","#E64B35")
Fibr$ident1<-factor(Fibr@active.ident,levels=c("Per","SM"))
Idents(Fibr)<-"ident2"
pdf("Identification_TSNE_Fibr2.pdf",width =7,height = 6)
DimPlot(Fibr,cols = col,reduction = "tsne",label = TRUE)+theme_classic()
dev.off()
pdf("Identification_UMAP_Fibr2.pdf",width =7,height = 6)
DimPlot(Fibr,cols = col,reduction = "umap",label = TRUE)+theme_classic()
dev.off()
save(Fibr,file="Fibr_rename.Rdata")


##### marker genes
DefaultAssay(Fibr)<-"SCT"
library(pheatmap)
library(ggplot2)
features = c("CD36","FABP4","POSTN",#Per_C1
             "CFD","CLSTN2","CTSC",#Per_C2
             "RGS5","ARHGDIB","CPM",#Per_C3
             "ZNF331","USP37","CREM",#SM_C1
             "MYH11","KCNMA1","LTBP1",#SM_C2
             "RCAN2","LAMA3","PPM1L",#SM_C3
             "PLN","ACTG2","CNN1",#SM_C4
             "MT1M","MT1E","MT1X",#SM_C5
             "RYR2","IGF1R","CACNA1C"#SM_C6
             )
markergene<-Fibr[["SCT"]][c(features),]
meta<-Fibr@meta.data
SCT<-c()
for(i in (as.vector(unique(meta$ident1)))){
  print(i)
  tmp=markergene[,meta$ident1==i]
  tmp=apply(tmp,1,mean)
  SCT<-cbind(SCT,tmp)
}
colnames(SCT)<-unique(meta$ident1)
SCT<-SCT[,c("Per_C1","Per_C2","Per_C3","SM_C1","SM_C2","SM_C3","SM_C4","SM_C5","SM_C6")]
#data.1 <- decostand(SCT,"standardize",MARGIN = 1)
pdf("markergene_heatmap1.pdf",width =3,height = 4)
pheatmap(SCT,scale="row",cluster_rows = F,cluster_cols = F,border=F,color =colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0" ,"#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))(100),treeheight_col = 10,gaps_row = c(3,6,9,12,15,18,21,24),gaps_col = c(3))+theme(axis.title.x = element_text(vjust = 1))
dev.off()

#####  find marker
DefaultAssay(Fibr)<-"RNA"
for(i in 0:(length(unique(Fibr$seurat_clusters))-1)){
  markers <- FindMarkers(Fibr, ident.1 = i,logfc.threshold =0.25,test.use = "MAST", min.pct = 0.25)
  markers<-cbind(rownames(markers),markers)
  colnames(markers)[1]<-"gene"
  write.table(markers,paste0("cluster",i,".markers.txt"),sep = "\t",row.names=FALSE)
}
Fibr@active.ident<-Fibr$seurat_clusters
markers <- FindMarkers(Fibr, ident.1 = c(4,8),logfc.threshold =0.25,test.use = "MAST", min.pct = 0.25)
markers<-cbind(rownames(markers),markers)
colnames(markers)[1]<-"gene"
write.table(markers,paste0("cluster_4and8.markers.txt"),sep = "\t",row.names=FALSE)


##### cluster_ident1_casecontrol
DefaultAssay(Fibr)<-"RNA"
Idents(Fibr)<-"status"
setwd("/gpfs/share/home/2116392061/01Pelvic_Lipomatosis/04Singlecelltype_all/06Fibr/cluster_ident1_CaseControl")
for(i in unique(Fibr$ident1)){
  a=subset(Fibr,ident1==i)
  markers <- FindMarkers(a, ident.1 = "Case",ident.2 = "Control",logfc.threshold =0.25,test.use = "MAST", min.pct = 0.25)
  markers<-cbind(rownames(markers),markers)
  colnames(markers)[1]<-"gene"
  write.table(markers,paste0("cluster",i,".markers_CaseControl_ident1.txt"),sep = "\t",row.names=FALSE)
}
##### cluster_ident2_casecontrol
setwd("/gpfs/share/home/2116392061/01Pelvic_Lipomatosis/04Singlecelltype_all/06Fibr/cluster_ident1_CaseControl2")
for(i in unique(Fibr$ident2)){
  a=subset(Fibr,ident2==i)
  markers <- FindMarkers(a, ident.1 = "Case",ident.2 = "Control",logfc.threshold =0.25,test.use = "MAST", min.pct = 0.25)
  markers<-cbind(rownames(markers),markers)
  colnames(markers)[1]<-"gene"
  write.table(markers,paste0("cluster",i,".markers_CaseControl_ident2.txt"),sep = "\t",row.names=FALSE)
}


#### cell percentage
Idents(Fibr)<-Fibr$ident1
for(i in c("CWL_Case","CWL_Control","HB_Case","LZF_Case","LZF_Control","ZX_Case","ZX_Control","ZZL_Case","ZZL_Control")){
  a <- subset(Fibr, orig.ident ==i)
  a=as.matrix(table(a@active.ident))
  a=data.frame(a,name=rownames(a))
  if(i=="CWL_Case"){
    b=a
  }
  else {
    b=merge(b,a,by="name",all=T)
  }
}
colnames(b)<-c("celltype2","CWL_Case","CWL_Control","HB_Case","LZF_Case","LZF_Control","ZX_Case","ZX_Control","ZZL_Case","ZZL_Control")
celltype1<-c(rep("Per",3),rep("SM",6))
b<-cbind(celltype1,b)
b<-t(b)
celltype1<-rep(b[1,],each=9)
celltype2<-rep(b[2,],each=9)
Num<-c()
for(i in 1:dim(b)[2]){Num<-c(Num,b[-c(1,2),i])}
samples<-rep(c("CWL_Case","CWL_Control","HB_Case","LZF_Case","LZF_Control","ZX_Case","ZX_Control","ZZL_Case","ZZL_Control"),dim(b)[2])
Individual<-c()
Status<-c()
for(i in samples){
  Individual<-c(Individual,strsplit(i,"_")[[1]][1])
  Status<-c(Status,strsplit(i,"_")[[1]][2])
}
b<-data.frame(celltype1,celltype2,samples,Num,Individual,Status)
b[is.na.data.frame(b)]=0
setwd("/gpfs/share/home/2116392061/01Pelvic_Lipomatosis/04Singlecelltype_all/06Fibr/Cell_percentage")
write.table(b,"Cell_percentage1.txt",sep = "\t",row.names = F)

Idents(Fibr)<-Fibr$ident2
for(i in c("CWL_Case","CWL_Control","HB_Case","LZF_Case","LZF_Control","ZX_Case","ZX_Control","ZZL_Case","ZZL_Control")){
  a <- subset(Fibr, orig.ident ==i)
  a=as.matrix(table(a@active.ident))
  a=data.frame(a,name=rownames(a))
  if(i=="CWL_Case"){
    b=a
  }
  else {
    b=merge(b,a,by="name",all=T)
  }
}
colnames(b)<-c("celltype1","CWL_Case","CWL_Control","HB_Case","LZF_Case","LZF_Control","ZX_Case","ZX_Control","ZZL_Case","ZZL_Control")
b<-t(b)
celltype1<-rep(b[1,],each=9)
Num<-c()
for(i in 1:dim(b)[2]){Num<-c(Num,b[-c(1),i])}
samples<-rep(c("CWL_Case","CWL_Control","HB_Case","LZF_Case","LZF_Control","ZX_Case","ZX_Control","ZZL_Case","ZZL_Control"),dim(b)[2])
Individual<-c()
Status<-c()
for(i in samples){
  Individual<-c(Individual,strsplit(i,"_")[[1]][1])
  Status<-c(Status,strsplit(i,"_")[[1]][2])
}
b<-data.frame(celltype1,samples,Num,Individual,Status)
b[is.na.data.frame(b)]=0
setwd("/gpfs/share/home/2116392061/01Pelvic_Lipomatosis/04Singlecelltype_all/06Fibr/Cell_percentage")
write.table(b,"Cell_percentage2.txt",sep = "\t",row.names = F)
##### cell percentage end



setwd("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/06Fibr")
Idents(Fibr)<-"ident2"
pdf("VEGFA.pdf",width =4,height = 3)
VlnPlot(Fibr,features = "VEGFA",pt.size = 0,split.by = "status",adjust = 2,cols = c("#EFC4CE","#2775AB"))
dev.off()
pdf("VEGFB.pdf",width =4,height = 3)
VlnPlot(Fibr,features = "VEGFB",pt.size = 0,split.by = "status",adjust = 2,cols = c("#EFC4CE","#2775AB"))
dev.off()
pdf("EPAS1.pdf",width =4,height = 3)
VlnPlot(Fibr,features = "EPAS1",pt.size = 0,split.by = "status",adjust = 2,cols = c("#EFC4CE","#2775AB"))
dev.off()
pdf("HIF1A.pdf",width =4,height = 3)
VlnPlot(Fibr,features = "HIF1A",pt.size = 0,split.by = "status",adjust = 2,y.max = 15,cols = c("#EFC4CE","#2775AB"))
dev.off()
Idents(Fibr)<-"ident1"
pdf("VEGFA2.pdf",width =4,height = 3)
VlnPlot(Fibr,features = "VEGFA",pt.size = 0,split.by = "status",adjust = 2,cols = c("#EFC4CE","#2775AB"))
dev.off()
pdf("VEGFB2.pdf",width =4,height = 3)
VlnPlot(Fibr,features = "VEGFB",pt.size = 0,split.by = "status",adjust = 2,cols = c("#EFC4CE","#2775AB"))
dev.off()
pdf("EPAS12.pdf",width =4,height = 3)
VlnPlot(Fibr,features = "EPAS1",pt.size = 0,split.by = "status",adjust = 2,cols = c("#EFC4CE","#2775AB"))
dev.off()
pdf("HIF1A2.pdf",width =4,height = 3)
VlnPlot(Fibr,features = "HIF1A",pt.size = 0,split.by = "status",adjust = 2,y.max = 15,cols = c("#EFC4CE","#2775AB"))
dev.off()