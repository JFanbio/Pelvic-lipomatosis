library(Seurat)
library(dplyr)
library(cowplot)
library(sctransform)
setwd("/gpfs/share/home/2116392061/01Pelvic_Lipomatosis/04Singlecelltype_all/03Adi")
load("/gpfs/share/home/2116392061/01Pelvic_Lipomatosis/03Identification_all/alldata_rename.Rdata")
DefaultAssay(alldata)<-"integrated"
alldata$ident<-alldata@active.ident
Adi=subset(alldata,subset = ident=="Adipocyte")

Adi <- RunPCA(Adi,npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

# save(Adi,file = "Adi.Rdata")
# pdf("metadata.pdf",width =10,height = 12)
# plots <-VlnPlot(Adi, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo")
#                 ,pt.size = 0, combine = FALSE)
# CombinePlots(plots = plots, ncol = 1)
# dev.off()
# 
# p1 <- DimPlot(Adi, reduction = "umap", group.by = "orig.ident")
# p2 <- DimPlot(Adi, reduction = "umap", label = TRUE)
# pdf("UMAP1.pdf",width = 12,height = 5)
# p1 | p2
# dev.off()
# 
# pdf("UMAP2.pdf",width = 15,height = 12)
# DimPlot(Adi, reduction = "umap", split.by = "orig.ident",label = TRUE,ncol = 5)
# dev.off()
# 
# pdf("UMAP3_status.pdf",width = 12,height = 5)
# DimPlot(Adi, reduction = "umap", split.by = "status",label = TRUE)
# dev.off()
# 
# pdf("UMAP4_library.pdf",width = 12,height = 5)
# DimPlot(Adi, reduction = "umap", split.by = "library",label = TRUE)
# dev.off()
# 
# pdf("UMAP5_individual.pdf",width = 12,height = 5)
# DimPlot(Adi, reduction = "umap", split.by = "individual",label = TRUE)
# dev.off()
# 
# 
# 
# p1 <- DimPlot(Adi, reduction = "tsne", group.by = "orig.ident")
# p2 <- DimPlot(Adi, reduction = "tsne", label = TRUE)
# pdf("TSNE1.pdf",width = 12,height = 5)
# p1 | p2
# dev.off()
# 
# pdf("TSNE2.pdf",width = 15,height = 10)
# DimPlot(Adi, reduction = "tsne", split.by = "orig.ident",label = TRUE,ncol = 6)
# dev.off()
# 
# pdf("TSNE3_status.pdf",width = 12,height = 5)
# DimPlot(Adi, reduction = "tsne", split.by = "status",label = TRUE)
# dev.off()
# 
# pdf("TSNE4_library.pdf",width = 12,height = 5)
# DimPlot(Adi, reduction = "tsne", split.by = "library",label = TRUE)
# dev.off()
# 
# pdf("TSNE5_individual.pdf",width = 12,height = 5)
# DimPlot(Adi, reduction = "tsne", split.by = "individual",label = TRUE)
# dev.off()

for(i in 0:(length(unique(Adi$seurat_clusters))-1)){
  markers <- FindMarkers(Adi, ident.1 = i, min.pct = 0.25)
  markers<-cbind(rownames(markers),markers)
  colnames(markers)[1]<-"gene"
  write.table(markers,paste0("cluster",i,".markers.txt"),sep = "\t",row.names=FALSE)
}

load("Adi.Rdata")
DefaultAssay(Adi) <- "RNA"

pdf("tmp1.pdf",width =6,height = 5)
VlnPlot(Adi, features = c("AR")
        ,pt.size = 0, combine = FALSE,split.by = "status")
dev.off()
pdf("tmp2.pdf",width =10,height = 9)
plots <-VlnPlot(Adi, features = c("PDGFRA","PDGFRB")
                ,pt.size = 0, combine = FALSE,split.by = "status")
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("ASC.pdf",width =10,height = 12)
plots <-VlnPlot(Adi, features = c("CD55","DPP4","MFAP5")
                ,pt.size = 0, combine = FALSE,split.by = "status")
CombinePlots(plots = plots, ncol = 1)
dev.off()


pdf("preA.pdf",width =10,height = 12)
plots <-VlnPlot(Adi, features = c("ICAM1","PPARG","GGT5")
                ,pt.size = 0, combine = FALSE,split.by = "status")
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("preA.pdf",width =10,height = 20)
plots <-VlnPlot(Adi, features = c("ICAM1","PPARG","GGT5","LPL","FABP4","AOC3","MGP","APOD")
                ,pt.size = 0, combine = FALSE,split.by = "status")
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("Areg.pdf",width =10,height = 12)
plots <-VlnPlot(Adi, features = c("F3","CLEC11A","FMO2")
                ,pt.size = 0, combine = FALSE,split.by = "status")
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("AR_Adi.pdf",width =4,height = 4)
VlnPlot(Adi, features = c("AR"),pt.size = 0, combine = FALSE,split.by = "status",group.by = "ident")
dev.off()
Adi <- RenameIdents(Adi,
                        '0' = "Venous ECs",
                        '1' = "Capillary ECs",
                        '2' = "Venous ECs",
                        '3' = "Venous ECs",
                        '4' = "Arterial ECs",
                        '5' = "Arterial ECs",
                        '6' = "Lymphatic ECs")
save(Adi,file="Adi_rename.Rdata")
library(ggsci)
library(ggplot2)
col=colorRampPalette(pal_jama()(4))(4)
col=c("#E64B35","#8A7284","#62ABBE","#16A79D","#227487","#D3988E","#91D1C2","#859AB6")
Adi@active.ident<-factor(Adi@active.ident,levels=c("Adipocytes","Adithelial cells","T/NK cells","Fibroblasts","Neutrophils","Monocytes","Mast cells","B cells"))
DimPlot(Adi,cols = col)+theme_test()

cluster0vs23.markers <- FindMarkers(Adi, ident.1 =0,ident.2 = c(2,3), min.pct = 0.25)
write.table(cluster0vs23.markers,"cluster0vs23.markers.txt",sep="\t")

#### vECs1 0
pdf("ident_vECs1_1.pdf",width =10,height = 5)
FeaturePlot(Adi, features = c("AQP1","POSTN"),cols =colorRampPalette(c("lightgrey","yellow","red","red"))(100),max.cutoff = c(60,15)) 
dev.off()
pdf("ident_vECs1_2.pdf",width =10,height = 10)
plots <-VlnPlot(Adi, features = c("AQP1","POSTN")
                ,pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()

#### vECs2 2
pdf("ident_vECs2_1.pdf",width =10,height = 9)
FeaturePlot(Adi, features = c("HLA-DQA1","PLVAP","SCARB1"),cols =colorRampPalette(c("lightgrey","yellow","red","red"))(100),max.cutoff = c(10,30,10)) 
dev.off()
pdf("ident_vECs2_2.pdf",width =10,height = 15)
plots <-VlnPlot(Adi, features = c("HLA-DQA1","PLVAP","SCARB1")
                ,pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()

#### vECs3 3
pdf("ident_vECs3_1.pdf",width =10,height = 9)
FeaturePlot(Adi, features = c("ICAM1","IL6","SELE","ADAMTS4"),cols =colorRampPalette(c("lightgrey","yellow","red","red"))(100),max.cutoff = c(50,100,100,50)) 
dev.off()
pdf("ident_vECs3_2.pdf",width =10,height = 20)
plots <-VlnPlot(Adi, features = c("ICAM1","IL6","SELE","ADAMTS4")
                ,pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()

#### vACs1 5
pdf("ident_vACs1_1.pdf",width =10,height = 9)
FeaturePlot(Adi, features = c("CXCL12","NEBL","SEMA3G"),cols =colorRampPalette(c("lightgrey","yellow","red","red"))(100),max.cutoff = c(20,5,10)) 
dev.off()
pdf("ident_vACs1_2.pdf",width =10,height =15)
plots <-VlnPlot(Adi, features = c("CXCL12","NEBL","SEMA3G")
                ,pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()

#### vACs2 4
pdf("ident_vACs2_1.pdf",width =10,height = 5)
FeaturePlot(Adi, features = c("TSPAN2","FBLIM1"),cols =colorRampPalette(c("lightgrey","yellow","red","red"))(100),max.cutoff = c(5,5)) 
dev.off()
pdf("ident_vACs2_2.pdf",width =10,height = 10)
plots <-VlnPlot(Adi, features = c("TSPAN2","FBLIM1")
                ,pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()

Adi <- RenameIdents(Adi,
                     '0' = "vECs1",
                     '1' = "cECs",
                     '2' = "vECs2",
                     '3' = "vECs3",
                     '4' = "vACs2",
                     '5' = "vACs1",
                     '6' = "lECs")
save(Adi,file="Adi_rename2.Rdata")
library(ggsci)
library(ggplot2)
col=colorRampPalette(pal_jama()(7))(7)
DimPlot(Adi,cols = col)+theme_test()

a <- subset(Adi, subset = sample =="Control1")
b=as.matrix(table(a@active.ident))
a <- subset(Adi, subset = sample =="Control2")
b=cbind(b,as.matrix(table(a@active.ident)))
a <- subset(Adi, subset = sample =="Control3")
b=cbind(b,as.matrix(table(a@active.ident)))
a <- subset(Adi, subset = sample =="Case1")
b=cbind(b,as.matrix(table(a@active.ident)))
a <- subset(Adi, subset = sample =="Case2")
b=cbind(b,as.matrix(table(a@active.ident)))
a <- subset(Adi, subset = sample =="Case3")
b=cbind(b,as.matrix(table(a@active.ident)))
colnames(b)<-c("Control1","Control2","Control3","Case1","Case2","Case3")
write.table(b,"Cell_percentage.txt",sep = "\t",row.names = T)



##### convert to anndata
SaveH5Seurat(Adi, filename = "Adi_fromseurat.h5Seurat")
Convert("Adi_fromseurat.h5Seurat", dest = "h5ad")











load("Adi_rename.Rdata")
DefaultAssay(Adi)<-"integrated"
Adi$ident<-Adi@active.ident
vECs=subset(Adi,subset = ident=="Venous ECs")
vECs<-SCTransform(vECs,vst.flavor = "v2", verbose = FALSE)


vECs <- RunPCA(vECs,npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

DefaultAssay(vECs)<-"RNA"
#save(Adi,file = "Adi.Rdata")
p1 <- DimPlot(vECs, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(vECs, reduction = "umap", label = TRUE)
p1 | p2





Adi<-Adi[,data.frame(Adi@assays[["RNA"]]@counts)["PECAM1",]==0]
Adi<-Adi[,data.frame(Adi@assays[["RNA"]]@counts)["ACTA2",]==0]
Adi<-Adi[,data.frame(Adi@assays[["RNA"]]@counts)["PDGFRA",]==0]
Adi<- subset(Adi,subset =  percent.mt < 5)
Adi <- RunPCA(Adi,npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)





DimPlot(Adi, reduction = "umap", label = TRUE,split.by = "status")
DefaultAssay(Adi) <- "RNA"
FeaturePlot(Adi, features = c("CD3D","CD4","CD8A","NKG7"),cols =colorRampPalette(c("lightgrey","yellow","red","red"))(100)) 




M=SplitObject(vECs,split.by = "orig.ident")

for (i in 1:length(M)) {
  M[[i]] <- NormalizeData(M[[i]], verbose = FALSE)
  M[[i]] <- FindVariableFeatures(
    M[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE
  )
}


M <- lapply(X = M, FUN = function(x) {
  x <- SCTransform(x,vst.flavor = "v2", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = M, nfeatures = 50)
M <- PrepSCTIntegration(object.list = M, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = M, normalization.method = "SCT",anchor.features = features,dims = 1:15)
#save(anchors,file = "anchors.Rdata")
vECs <- IntegrateData(anchorset = anchors, normalization.method = "SCT")


DimPlot(Adi, reduction = "umap", split.by = "status")







features <- SelectIntegrationFeatures(object.list = M)
M <- lapply(X = M, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

alldata <- merge(LZF_Case,c(ZX_Control,ZZL_Case,ZZL_Control,LZF_Case,LZF_Control), add.cell.ids=c("Case1","Control1","Case2","Control2","Case3","Control3"))

ZX_Control <- SCTransform(ZX_Control, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)
ZX_Case <- SCTransform(ZX_Case, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)
ZZL_Case <- SCTransform(ZZL_Case, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)
ZZL_Control <- SCTransform(ZZL_Control, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)
LZF_Case <- SCTransform(LZF_Case, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)
LZF_Control <- SCTransform(LZF_Control, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)

features <- SelectIntegrationFeatures(object.list = list(ZX_Case,ZX_Control,ZZL_Case,ZZL_Control,LZF_Case,LZF_Control), nfeatures = 3000)
M <- PrepSCTIntegration(object.list = list(ZX_Case,ZX_Control,ZZL_Case,ZZL_Control,LZF_Case,LZF_Control), anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = M, normalization.method = "SCT",
                                  anchor.features = features)
M.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
M.combined <- RunPCA(M.combined, verbose = FALSE)
M.combined <- RunUMAP(M.combined, reduction = "pca", dims = 1:30, verbose = FALSE)
M.combined <- FindNeighbors(M.combined, reduction = "pca", dims = 1:30)
M.combined <- FindClusters(M.combined, resolution = 0.7)

ZX_Control <- RunPCA(ZX_Control, verbose = FALSE)
ZX_Case <- RunPCA(ZX_Case, verbose = FALSE)
ZZL_Case <- RunPCA(ZZL_Case, verbose = FALSE)
ZZL_Control <- RunPCA(ZZL_Control, verbose = FALSE)
LZF_Case <- RunPCA(LZF_Case, verbose = FALSE)
LZF_Control <- RunPCA(LZF_Control, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = list(ZX_Case,ZX_Control,ZZL_Case,ZZL_Control,LZF_Case,LZF_Control), reduction = "rpca", dims = 1:50)

DefaultAssay(M.combined) <- "integrated"
M.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
M.integrated <- ScaleData(M.integrated, verbose = FALSE)
M.integrated <- RunPCA(M.integrated, verbose = FALSE)
M.integrated <- RunUMAP(M.integrated, dims = 1:50)
M.integrated <- FindNeighbors(M.integrated,reduction = "pca", dims = 1:30, verbose = FALSE) 
M.integrated <- FindClusters(M.integrated,resolution = 0.7, verbose = FALSE)


M <- FindIntegrationAnchors(object.list = M, dims = 1:30)
Adi <- IntegrateData(anchorset = M, dims = 1:30)


Adi <- SCTransform(Adi, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

p1 <- DimPlot(Adi, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(Adi, reduction = "umap", label = TRUE)
p1 | p2




