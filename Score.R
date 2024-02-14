library(Seurat)
library(dplyr)
library(cowplot)
library(ggpubr)
library(ggplot2)
library(ggsci)
setwd("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/00alldata/")
load("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/00alldata/alldata_rename.Rdata")
alldata@active.ident<-factor(alldata@active.ident,levels=c("Adipocyte","ASPC","Endothelial","Lymphatic endo","Pericyte","Smooth muscle","Macrophage","Monocyte","Dendritic cell","Mast cell","Neutrophil","B cell","NK cell","T cell"))
alldata$ident<-alldata@active.ident
Idents(alldata)<-"ident"
col=colorRampPalette(pal_npg()(7))(14)


#### 50 hallmarker
geneList=readLines("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/01Endo/Score/h.all.v2022.1.Hs.symbols.gmt")
geneList <- strsplit(geneList, "\t")
compaired=list(c("Case","Control"))
for(i in 1:50){
  fun_name<-geneList[[i]][1]
  genelist<-geneList[[i]][-c(1,2)]
  genelist<-list(genelist)
  alldata=AddModuleScore(alldata,features = genelist,name = fun_name,assay="RNA")
  p1=ggecdf(alldata@meta.data, x =paste0(fun_name,"1"),color = "ident",size = 1,palette = col)+theme_test()
  p2=ggecdf(alldata@meta.data, x =paste0(fun_name,"1"),color = "status",size = 1,palette = c("#EFC4CE","#2775AB"))+theme_test()+facet_wrap(.~ident)
  p3=VlnPlot(alldata,features = paste0(fun_name,"1"),pt.size = 0,cols = col,y.max = max(alldata@meta.data[,paste0(fun_name,"1")])*1.1)+NoLegend()+geom_boxplot(width=0.2,outlier.colour = NA,fill="white")+stat_compare_means(comparisons = compaired,method="wilcox.test")+ggtitle(fun_name)+theme(plot.title = element_text(size = 8))
  ggsave((p1|p3),filename =paste0(fun_name,"_Score1.pdf"),device = "pdf",width = 7,height = 4)
  ggsave(p2,filename =paste0(fun_name,"_Score2.pdf"),device = "pdf",width = 9,height = 8)
}
save(alldata,file="alldata_rename.Rdata")



#### hallmarker hypoxia
fun_name<-"HALLMARK_HYPOXIA"
p1=ggecdf(alldata@meta.data, x =paste0(fun_name,"1"),color = "ident",size = 1,palette = col,xlim=c(-1,4))+theme_test()
p2=ggecdf(alldata@meta.data, x =paste0(fun_name,"1"),color = "status",size = 1,palette = c("#EFC4CE","#2775AB"),xlim=c(-1,4))+theme_test()+facet_wrap(.~ident)
p3=VlnPlot(alldata,features = paste0(fun_name,"1"),pt.size = 0,cols = col,y.max = max(alldata@meta.data[,paste0(fun_name,"1")])*1.1)+NoLegend()+geom_boxplot(width=0.2,outlier.colour = NA,fill="white")+stat_compare_means(comparisons = compaired,method="wilcox.test")+ggtitle(fun_name)+theme(plot.title = element_text(size = 8))
ggsave((p1|p3),filename =paste0(fun_name,"_Score1.pdf"),device = "pdf",width = 7,height = 4)
ggsave(p2,filename =paste0(fun_name,"_Score2.pdf"),device = "pdf",width = 9,height = 8)

#### hallmarker TGF_BETA_SIGNALING
fun_name<-"HALLMARK_TGF_BETA_SIGNALING"
p1=ggecdf(alldata@meta.data, x =paste0(fun_name,"1"),color = "ident",size = 1,palette = col,xlim=c(-1,2))+theme_test()
p2=ggecdf(alldata@meta.data, x =paste0(fun_name,"1"),color = "status",size = 1,palette = c("#EFC4CE","#2775AB"),xlim=c(-1,2))+theme_test()+facet_wrap(.~ident)
p3=VlnPlot(alldata,features = paste0(fun_name,"1"),pt.size = 0,cols = col,y.max = max(alldata@meta.data[,paste0(fun_name,"1")])*1.1)+NoLegend()+geom_boxplot(width=0.2,outlier.colour = NA,fill="white")+stat_compare_means(comparisons = compaired,method="wilcox.test")+ggtitle(fun_name)+theme(plot.title = element_text(size = 8))
ggsave((p1|p3),filename =paste0(fun_name,"_Score1.pdf"),device = "pdf",width = 7,height = 4)
ggsave(p2,filename =paste0(fun_name,"_Score2.pdf"),device = "pdf",width = 9,height = 8)



####  LAM_Score (Lipid-Associated Macrophages)
genelist=list(c("LIPA","CTSB","CTSL","FABP4","FABP5","LGALS3","CD9","CD36"))
fun_name<-"LAM_Score"
alldata=AddModuleScore(alldata,features = genelist,name = fun_name,assay="RNA")
p1=ggecdf(alldata@meta.data,x =paste0(fun_name,"1"),color = "status",size = 1,palette=c("#EFC4CE","#2775AB"),xlim=c(0,100))+theme_test()+facet_wrap(.~ident1)
p2=ggecdf(alldata@meta.data,x =paste0(fun_name,"1"),color = "ident1",size = 1,palette = col,xlim=c(0,100))+theme_test()
ggsave(p1+p2,device = "pdf",width = 10,height = 4,filename = "LAM_Score.pdf")

p1=VlnPlot(subset(alldata,ident2=="Mac"), features = c("TREM2"),pt.size = 0,split.by = "status",cols = c("#EFC4CE","#2775AB"),adjust = 2)
p2=VlnPlot(subset(alldata,ident2=="Mac"), features = c("LAM_Score1"),pt.size = 0,split.by = "status",cols = c("#EFC4CE","#2775AB"))
ggsave(p1| p2,device = "pdf",width = 8,height = 3,filename = "TREM2_vlnplot_LAM.pdf")

#### ANGIOGENESIS score
genelist=readLines("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/01Endo/Score/HALLMARK_ANGIOGENESIS.v2022.1.Hs.gmt")
genelist <- strsplit(genelist, "\t")
fun_name<-genelist[[1]][1]
genelist<-genelist[[1]][-c(1,2)]
genelist<-list(genelist)
alldata=AddModuleScore(alldata,features = genelist,name = fun_name,assay="RNA")
p1=ggecdf(alldata@meta.data,x =paste0(fun_name,"1"),color = "ident1",size = 1,palette = col,xlim=c(0,5))+theme_test()
p2=ggecdf(alldata@meta.data,x =paste0(fun_name,"1"),color = "ident",size = 1,palette = col)+theme_test()
p3=ggecdf(subset(alldata,ident1=="Mono_C1")@meta.data,x =paste0(fun_name,"1"),color = "status",size = 1,palette = c("#EFC4CE","#2775AB"),xlim=c(-2,10))+theme_test()
p4=VlnPlot(subset(alldata,ident2=="Mono"), features = c("ANGIOGENESIS1"),pt.size = 0,split.by = "status",cols = c("#EFC4CE","#2775AB"))
ggsave((p1|p2)/(p3|p4) ,device = "pdf",width = 8,height = 6.5,filename = "ANGIOGENESIS_Score.pdf")
p1=VlnPlot(subset(alldata,ident2=="Mono"), features = c("VEGFA"),pt.size = 0,split.by = "status",cols = c("#EFC4CE","#2775AB"))
p2=VlnPlot(Endo, features = c("KDR"),pt.size = 0,split.by = "status",cols = c("#EFC4CE","#2775AB"))
ggsave((p1|p2) ,device = "pdf",width = 8,height = 3,filename = "VEGF.pdf")



#### HALLMARK_HYPOXIA score
genelist=readLines("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/01Endo/Score/HALLMARK_HYPOXIA.v2022.1.Hs.gmt")
genelist <- strsplit(genelist, "\t")
fun_name<-genelist[[1]][1]
genelist<-genelist[[1]][-c(1,2)]
genelist<-list(genelist)
alldata=AddModuleScore(alldata,features = genelist,name = fun_name,assay="RNA")
p1=ggecdf(alldata@meta.data,x =paste0(fun_name,"1"),color = "ident1",size = 1,palette = col,xlim=c(-1,4))+theme_test()
p2=ggecdf(alldata@meta.data,x =paste0(fun_name,"1"),color = "ident",size = 1,palette = col,xlim=c(-1,4))+theme_test()
p3=ggecdf(subset(alldata,ident1=="Mono_C1")@meta.data,x =paste0(fun_name,"1"),color = "status",size = 1,palette = c("#EFC4CE","#2775AB"),xlim=c(-1,4))+theme_test()
p4=VlnPlot(subset(alldata,ident2=="Mono"), features = c("HALLMARK_HYPOXIA1"),pt.size = 0,split.by = "status",cols = c("#EFC4CE","#2775AB"))
ggsave((p1|p2)/(p3|p4) ,device = "pdf",width = 8,height = 6.5,filename = "HALLMARK_HYPOXIA_Score.pdf")
p1=VlnPlot(subset(alldata,ident2=="Mono"), features = c("HIF1A"),pt.size = 0,split.by = "status",cols = c("#EFC4CE","#2775AB"))
p2=VlnPlot(Endo, features = c("HIF1A"),pt.size = 0,split.by = "status",cols = c("#EFC4CE","#2775AB"))
ggsave((p1|p2) ,device = "pdf",width = 8,height = 3,filename = "HIF1A.pdf")


####  Dendritic_activated_Score
genelist=list(c("CCL22","LAMP3","IDO1","CCR7","CXCL10","BIRC3","TNFAIP6","IL7R","RSAD2","SAMSN1","KYNU","CCL5","EBI3","HLA-DQA1","CCL17","BCL2A1","CCL13","CLIC2","PLA2G7","CD40","SLC15A3","CD80","CD86","CCL4","FPR3","CST7","IFI44L","ACP5","CXCL11","MMP12"))
fun_name<-"Dendritic_activated_Score"
alldata=AddModuleScore(alldata,features = genelist,name = fun_name,assay="RNA")
p1=ggecdf(subset(alldata,ident2=="DC")@meta.data,x =paste0(fun_name,"1"),color = "status",size = 1,palette=c("#EFC4CE","#2775AB"),xlim=c(-1,10))+theme_test()
ggsave(p1,device = "pdf",width = 4,height = 3,filename = "Dendritic_activated_Score.pdf")

####  Dendritic_resting_Score
genelist=list(c("MMP12","CD1B","CD1A","ACP5","NCF2","CD1C","CD1E","RNASE6","HLA-DQA1","TREM2","FPR3","PLA2G7","CLEC7A","CLEC10A","C1orf54","MNDA","EGR2","HCK","CLEC4A","CCL22"))


####   Exhausted_Score
genelist=list(c("LAG3","TIGIT","PDCD1","CTLA4","HAVCR2"))
fun_name<-"Exhausted_Score"
alldata=AddModuleScore(alldata,features = genelist,name = fun_name,assay="RNA")
p1=ggecdf(alldata@meta.data,x =paste0(fun_name,"1"),color = "status",size = 1,palette=c("#EFC4CE","#2775AB"),xlim=c(-1,2))+theme_test()
a=subset(alldata,ident1!="CD4_C7")
p2=ggecdf(alldata@meta.data,x =paste0(fun_name,"1"),color = "ident1",size = 1,palette = col,xlim=c(-1,2))+theme_test()
ggsave(p1+p2,device = "pdf",width = 8,height = 3,filename = "Exhausted_Score.pdf")



library(ggpointdensity)
Case<-subset(alldata,status=="Case")
control<-subset(alldata,status=="Control")
data1 <-cbind(Embeddings(object=Case[['umap']]),FetchData(Case, "ident2"))
data2 <-cbind(Embeddings(object=control[['umap']]),FetchData(control, "ident2"))
p1=ggplot(data = data1, mapping = aes(x =UMAP_1,y = UMAP_2)) +
  geom_pointdensity() + #密度散点图(geom pointdensity)
  scale_color_viridis()+theme_classic()
p2=ggplot(data = data2, mapping = aes(x =UMAP_1,y = UMAP_2)) +
  geom_pointdensity() + #密度散点图(geom pointdensity)
  scale_color_viridis()+theme_classic()

DefaultAssay(alldata)<-"integrated"
setwd("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/05alldata/Score")

DoHeatmap(alldata, features = c("ACKR1","IL1R1","GJA5","HEY1","RGCC","CD300LG"), cells = 1:500, size = 4,
          angle = 90) + NoLegend()
DoHeatmap(alldata)

genelist=readLines("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/01Endo/Score/ANGIOGENESIS.v2022.1.Hs.gmt")
genelist=readLines("HALLMARK_ADIPOGENESIS.v2022.1.Hs.gmt")
genelist=readLines("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/01Endo/Score/HALLMARK_HYPOXIA.v2022.1.Hs.gmt")

genelist=readLines("h.all.v2022.1.Hs.symbols.gmt")
genelist <- strsplit(genelist, "\t")
fun_name<-genelist[[1]][1]
genelist<-genelist[[1]][-c(1,2)]
genelist<-list(genelist)


Capillary<-subset(alldata,ident=="Capillary")
Arterial<-subset(alldata,ident=="Arterial")
Venous<-subset(alldata,ident=="Venous")
Idents(Capillary)<-"status"
Idents(Arterial)<-"status"
Idents(Venous)<-"status"
compaired=list(c("Case","Control"))
geneList=readLines("h.all.v2022.1.Hs.symbols.gmt")
geneList <- strsplit(geneList, "\t")

for(i in 1:50){
  fun_name<-geneList[[i]][1]
  genelist<-geneList[[i]][-c(1,2)]
  genelist<-list(genelist)
  alldata=AddModuleScore(alldata,features = genelist,name = fun_name,assay="RNA")
}
for(i in 1:50){
  fun_name<-geneList[[i]][1]
  genelist<-geneList[[i]][-c(1,2)]
  genelist<-list(genelist)
  Capillary=AddModuleScore(Capillary,features = genelist,name = fun_name,assay="RNA")
  Arterial=AddModuleScore(Arterial,features = genelist,name = fun_name,assay="RNA")
  Venous=AddModuleScore(Venous,features = genelist,name = fun_name,assay="RNA")
  p1=ggecdf(Capillary@meta.data, x =paste0(fun_name,"1"),color = "status",size = 1,palette = c("#EFC4CE","#2775AB"))+theme_test()+ggtitle("Capillary")
  p2=ggecdf(Arterial@meta.data, x =paste0(fun_name,"1"),color = "status",size = 1,palette = c("#EFC4CE","#2775AB"))+theme_test()+ggtitle("Arterial")
  p3=ggecdf(Venous@meta.data, x =paste0(fun_name,"1"),color = "status",size = 1,palette = c("#EFC4CE","#2775AB"))+theme_test()+ggtitle("Venous")
  ggsave(p1+p2+p3,filename =paste0("cdf_",fun_name,".pdf"),device = "pdf",width = 11,height = 3)
  p4=VlnPlot(Capillary,features = paste0(fun_name,"1"),pt.size = 0,cols = c("#EFC4CE","#2775AB"),y.max = max(Capillary@meta.data[,paste0(fun_name,"1")])*1.1)+NoLegend()+geom_boxplot(width=0.2,outlier.colour = NA,fill="white")+stat_compare_means(comparisons = compaired,method="wilcox.test")+ggtitle("Capillary")+theme(plot.title = element_text(size = 8))
  p5=VlnPlot(Arterial,features = paste0(fun_name,"1"),pt.size = 0,cols = c("#EFC4CE","#2775AB"),y.max = max(Arterial@meta.data[,paste0(fun_name,"1")])*1.1)+NoLegend()+geom_boxplot(width=0.2,outlier.colour = NA,fill="white")+stat_compare_means(comparisons = compaired,method="wilcox.test")+ggtitle("Arterial")+theme(plot.title = element_text(size = 8))
  p6=VlnPlot(Venous,features = paste0(fun_name,"1"),pt.size = 0,cols = c("#EFC4CE","#2775AB"),y.max = max(Venous@meta.data[,paste0(fun_name,"1")])*1.1)+NoLegend()+geom_boxplot(width=0.2,outlier.colour = NA,fill="white")+stat_compare_means(comparisons = compaired,method="wilcox.test")+ggtitle("Venous")+theme(plot.title = element_text(size = 8))
  ggsave(p4|p5|p6,filename =paste0("Vlnplot_",fun_name,".pdf"),device = "pdf",width = 6,height = 4)
}



#### hypoxia in different alldatathelial cell
genelist=readLines("HALLMARK_HYPOXIA.v2022.1.Hs.gmt")
genelist <- strsplit(genelist, "\t")
fun_name<-genelist[[1]][1]
genelist<-genelist[[1]][-c(1,2)]
genelist<-list(genelist)
alldata=AddModuleScore(alldata,features = genelist,name = fun_name,assay="RNA")
compaired<-list(c("Capillary","Arterial"),c("Capillary","Venous"),c("Venous","Arterial"))
SAT<-subset(alldata,status=="Control")
VAT<-subset(alldata,status=="Case")
p1=VlnPlot(VAT,features = paste0(fun_name,"1"),pt.size = 0,cols = col,y.max = max(VAT@meta.data[,paste0(fun_name,"1")])*1.5)+NoLegend()+geom_boxplot(width=0.2,outlier.colour = NA,fill="white")+stat_compare_means(comparisons = compaired,method="wilcox.test")+ggtitle("VAT_HYPOXIA")+theme(plot.title = element_text(size = 8))
p2=VlnPlot(SAT,features = paste0(fun_name,"1"),pt.size = 0,cols = col,y.max = max(SAT@meta.data[,paste0(fun_name,"1")])*1.5)+NoLegend()+geom_boxplot(width=0.2,outlier.colour = NA,fill="white")+stat_compare_means(comparisons = compaired,method="wilcox.test")+ggtitle("SAT_HYPOXIA")+theme(plot.title = element_text(size = 8))
p3=ggecdf(VAT@meta.data, x =paste0(fun_name,"1"),color = "ident",size = 1,palette = col)+theme_test()+ggtitle("VAT_HYPOXIA")
p4=ggecdf(SAT@meta.data, x =paste0(fun_name,"1"),color = "ident",size = 1,palette = col)+theme_test()+ggtitle("SAT_HYPOXIA")
(p1|p2)/(p3|p4)
ggsave(p1|p2,filename =paste0("hypoxia in different alldatathelial cell1.pdf"),device = "pdf",width = 5,height = 4)
ggsave(p3|p4,filename =paste0("hypoxia in different alldatathelial cell2.pdf"),device = "pdf",width = 8,height = 3)


##### different gene asscociated with hypoxia
VlnPlot(alldata,features ="HIF1A",fill.by="ident",split.by = "status",cols = c("#EFC4CE","#2775AB"),pt.size=0,stack = F,flip=F,adjust = 2,y.max = 20)
VlnPlot(alldata,features ="EPAS1",fill.by="ident",split.by = "status",cols = c("#EFC4CE","#2775AB"),pt.size=0,stack = F,flip=F,adjust = 2,y.max = 50)
hypoxia_feature<-c("HSP90AA1","HSP90AB1")
VlnPlot(alldata,features = hypoxia_feature,fill.by="ident",split.by = "status",cols = c("#EFC4CE","#2775AB"),pt.size=0,stack = F,flip=F,adjust = 2,y.max = 100)
hypoxia_feature<-c("SLC25A4","SLC25A5","SLC25A6")
VlnPlot(alldata,features = hypoxia_feature,fill.by="ident",split.by = "status",cols = c("#EFC4CE","#2775AB"),pt.size=0,stack = F,flip=F,adjust = 2)
hypoxia_feature<-c("UCP1","UCP2","UCP3")
VlnPlot(alldata,features = hypoxia_feature,fill.by="ident",split.by = "status",cols = c("#EFC4CE","#2775AB"),pt.size=0,stack = F,flip=F,adjust = 2)

hypoxia_feature<-c("DTNA","PPP1R15A","IER3","JUN","FOS","ETS1","DUSP1","PNRC1","CCN1","LDHA","BHLHE40","WSB1","PLIN2","PDGFB","SLC2A3","MAFF")
hypoxia_feature<-c("HIF1A","EPAS1","ARNT")
DefaultAssay(alldata)<-"SCT"
VlnPlot(alldata,features = hypoxia_feature,fill.by="ident",split.by = "status",cols = c("#EFC4CE","#2775AB"),pt.size=0,stack = T,flip=F,adjust = 2)+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+
  xlab("")+ylab("")
