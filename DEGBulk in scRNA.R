library(Seurat)
library(dplyr)
library(cowplot)
library(sctransform)
library(plyr)
library(ggplot2)
setwd("E:/00Scientific_Research/Pelvic_Lipomatosis_new/06BulkRNA/01alldata")
filelist<-list.files("E:/00Scientific_Research/Pelvic_Lipomatosis_new/06BulkRNA/01alldata")[-c(4,5)]
data<-c()
for (i in filelist) {
  a=read.table(i,head=T)
  celltype<-gsub("_withBulk.txt","",i)
  ident=rep(celltype,dim(a)[1])
  a=cbind(ident,a)
  data<-rbind(data,a)
}
congruous<-ifelse(data$avg_log2FC*data$log2FoldChange>0,"congruous","incongruous")
direction<-ifelse(data$avg_log2FC>0,"Up","Down")
data<-cbind(direction,data)
data<-cbind(congruous,data)
data$p_val_adj=ifelse(data$p_val_adj==0,min(data$p_val_adj[data$p_val_adj!=0]),data$p_val_adj)
data$ident<-factor(data$ident,levels=c("Adipocyte","ASPC","Endothelial","Lymphatic endo","Pericyte","Smooth muscle" ,"T cell","NK cell","B cell","Macrophage","Monocyte","Neutrophil","Dendritic cell","Mast cell"))
 
ggplot(data[data$congruous=="congruous",],aes(x=ident,y=gene,size=abs(avg_log2FC),color=direction))+
  #geom_point()+
  geom_point(aes(alpha=-log(p_val_adj)))+
  scale_color_manual(values=c("navy","firebrick3"))+
  scale_alpha_continuous(range = c(0.4, 1))+
  theme_bw()+
  theme(text = element_text(colour = "black"),axis.text.x = element_text(angle = 60,hjust = 1))+
  coord_flip()+coord_polar()
## 6*11

ggplot(data[data$congruous=="congruous",],aes(x=gene,y=ident,size=abs(avg_log2FC),color=direction))+
  #geom_point()+
  geom_point(aes(alpha=-log(p_val_adj)))+
  scale_color_manual(values=c("navy","firebrick3"))+
  scale_alpha_continuous(range = c(0.4, 1))+
  theme_bw()+
  theme(text = element_text(colour = "black"),axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
