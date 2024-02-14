library(Seurat)
library(dplyr)
library(cowplot)
library(CellChat)
load("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/06Fibr/Fibr_rename.Rdata")
load("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/01Endo/Endo_rename.Rdata")
Fibr$ident<-Fibr$ident
Endo$ident<-Endo@active.ident
EndoFibr<-merge(Endo,Fibr)
EndoFibr$ident<-factor(EndoFibr$ident,levels = c("Venous","Arterial","Capillary","Per","SM"))
EndoFibr@active.ident<-EndoFibr$ident
setwd("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/06Fibr/Endofibr_CellChat/")
CellChatDB <- CellChatDB.human
#showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

col= c("#374E55","#DF8F44","#00A1D5","#E64B35","#739FAD")
for(i in c("Control","Case")){
  data<-subset(EndoFibr,status==i)
  data.input <- GetAssayData(data,assay = "RNA",slot = "data")
  data$ident<-data@active.ident
  Idents(data)<-data@active.ident
  labels <- Idents(data)
  identity <- data.frame(group = labels, row.names = names(labels))
  # create a dataframe of the cell
  cellchat <- createCellChat(object = data.input, meta = identity, group.by = "group")
  ####  Preprocessing the expression data for cell-cell communication analysis
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)# This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net<-subsetCommunication(cellchat,slot.name = "net")
  write.csv(df.net,paste0("net_lr_",i,".csv"))
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp<-subsetCommunication(cellchat,slot.name = "netP")
  write.csv(df.netp,paste0("net_pathway_",i,".csv"))
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  ### û??counts????
  pdf(paste0("Number of interactions of ",i,".pdf"),width = 7, height = 6)
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = paste0("Number of interactions of ",i),color.use = col)
  #ggsave(p,device = "pdf",width = 15,height = 5,filename = paste0("Number of interactions of ",i))
  dev.off()
  ### ??counts????
  pdf(paste0("Number of interactions of ",i,"_count.pdf"),width = 7, height = 6)
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = paste0("Number of interactions of ",i),color.use = col)
  dev.off()
  ######  ????ϸ??????
  mat <- cellchat@net$count
  pdf(paste0("Number of interactions of ",i,"_singlecelltype.pdf"),width = 8, height = 5)
  par(mfrow = c(2,3), xpd=TRUE)
  for (j in 1:5) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[j, ] <- mat[j, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[j],color.use = col)
  }
  dev.off()
  pdf(paste0("Number of interactions of ",i,"_singlecelltype_count.pdf"),width = 8, height = 5)
  par(mfrow = c(2,3), xpd=TRUE)
  for (j in 1:5) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[j, ] <- mat[j, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[j],label.edge = T,color.use = col)
  }
  dev.off()
  ########  weight  #############
  ### û??counts????
  pdf(paste0("weight of interactions of ",i,".pdf"),width = 7, height = 6)
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = paste0("weight of interactions of ",i),color.use = col)
  dev.off()
  ######  ????ϸ??????
  mat <- cellchat@net$weight
  pdf(paste0("weight of interactions of ",i,"_singlecelltype.pdf"),width = 8, height = 5)
  par(mfrow = c(2,3), xpd=TRUE)
  for (j in 1:5) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[j, ] <- mat[j, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[j],color.use = col)
  }
  dev.off()
  save(cellchat,file = paste0("cellchat_",i,".Rdata"))
}

la
#### comparison
library(CellChat)
library(Seurat)
library(dplyr)
library(cowplot)
library(sctransform)
setwd("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/06Fibr/Endofibr_CellChat/diff")
load("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/06Fibr/Endofibr_CellChat/cellchat_Case.Rdata")
Case<-cellchat
load("E:/00Scientific_Research/Pelvic_Lipomatosis_new/04Singlecelltype_all/06Fibr/Endofibr_CellChat/cellchat_Control.Rdata")
Control<-cellchat

status.list<-list(Case=Case,Control=Control)
cellchat<-mergeCellChat(status.list,add.names =names(status.list),cell.prefix = TRUE)
####ͨѶ??��??ǿ?ȶԱ?
gg1<-compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = "count",,color.use = c("#EFC4CE","#2775AB"))
gg2<-compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = "weight",,color.use = c("#EFC4CE","#2775AB"))
ggsave("Overview_number_strength.pdf",gg1+gg2,width = 6,height = 4)
####  ??��??ǿ?Ȳ???????ͼ
pdf("number and weight of different interactions.pdf",width = 14, height = 6)
par(mfrow=c(1,2), xpd=T)
netVisual_diffInteraction(cellchat,weight.scale = T,comparison = c(2,1))
netVisual_diffInteraction(cellchat,weight.scale = T,comparison = c(2,1),measure = "weight")
dev.off()

####  ??��??ǿ?Ȳ?????ͼ
pdf("number and weight of different interactions_heamap.pdf",width = 12, height = 6)
par(mfrow=c(1,2), xpd=F)
h1=netVisual_heatmap(cellchat, comparison = c(2,1),color.use=col)
h2=netVisual_heatmap(cellchat,measure = "weight",comparison = c(2,1),color.use=col)
h1+h2
dev.off()
gg1<-rankNet(cellchat,mode = "comparison",stacked = T,do.stat = TRUE,color.use = c("#EFC4CE","#2775AB"))
gg2<-rankNet(cellchat,mode = "comparison",stacked = F,do.stat = TRUE,color.use = c("#EFC4CE","#2775AB"))
ggsave("Compare_pathway_strengh.pdf",gg1+gg2,width = 9,height = 5)


num.link <- sapply(status.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(status.list)) {
  gg[[i]] <- netAnalysis_computeCentrality(status.list[[i]])%>%netAnalysis_signalingRole_scatter(title = names(status.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("interaction strength.pdf",width = 10, height = 4.5)
patchwork::wrap_plots(plots = gg)
dev.off()


######   ??????ͼ??very important!!
for(i in 1:5){
  pdf(paste0(i,"_bubble.pdf"),width = 12, height = 6)
  netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:5),  comparison = c(1, 2), angle.x = 45)
  dev.off()
}
pdf("2.pdf",width = 12, height = 6)
netVisual_bubble(cellchat, sources.use =c(1:14) , targets.use = "Endothelial",  comparison = c(1, 2), angle.x = 45)
dev.off()
##### pathways
pathways.show <- c("ANGPTL")
weight.,ax<-getMaxWeight()
levels(cellchat@idents)
vertex.receiver = c(1:5)
par(mfrow=c(2,1))
netVisual_aggregate(cellchat,signaling = pathways.show,vertex.receiver = vertex.receiver,layout = "hierarchy") 
p1=netVisual_aggregate(Case,signaling = "VEGF",color.use = col)
p2=netVisual_aggregate(Control,signaling = "VEGF",color.use = col) 
netVisual_heatmap(cellchat,signaling = pathways.show)

