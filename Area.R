setwd("E:/00Scientific_Research/Pelvic_Lipomatosis_new/细胞面积/Area")
PL=read.table("E:/00Scientific_Research/Pelvic_Lipomatosis_new/细胞面积/Area/PL_Area.txt",header = T,sep = "\t")
top<-colnames(PL)
area_PL=c()
for(i in 1:dim(PL)[2]){
  a=na.omit(PL[,i])
  b=strsplit(top[i],split = "_")[[1]]
  tmp<-cbind(rep(b[1],length(a)),rep(b[2],length(a)),rep("PL",length(a)),a)
  area_PL=rbind(area_PL,tmp)
}
BC=read.table("E:/00Scientific_Research/Pelvic_Lipomatosis_new/细胞面积/Area/BC_Area.txt",header = T,sep = "\t")
top<-colnames(BC)
area_BC=c()
for(i in 1:dim(BC)[2]){
  a=na.omit(BC[,i])
  b=strsplit(top[i],split = "_")[[1]]
  tmp<-cbind(rep(b[1],length(a)),rep(b[2],length(a)),rep("BC",length(a)),a)
  area_BC=rbind(area_BC,tmp)
}
area<-rbind(area_PL,area_BC)
colnames(area)<-c("Individual","Status","dis","Area")
library(ggplot2)
area<-data.frame(area)
area[,4]<-as.numeric(area[,4])
area[,2]<-factor(area[,2],levels = c("SAT","PAT"))
#### single individal
ggplot(area,aes(colour=Status,fill=Status))+geom_density(aes(x=Area,alpha=0.5))+facet_grid(dis~Individual)+
  theme_bw()+theme(axis.text = element_text(colour = "black"),axis.text.x = element_text(angle = 45))+scale_color_manual(values=c("#2775AB","#EFC4CE"))+scale_fill_manual(values=c("#2775AB","#EFC4CE"))
#### BC vs PL
ggplot(area,aes(colour=Status,fill=Status))+geom_density(aes(x=Area,alpha=0.5))+facet_grid(.~dis)+
  theme_bw()+theme(axis.text = element_text(colour = "black"))+scale_color_manual(values=c("#2775AB","#EFC4CE"))+scale_fill_manual(values=c("#2775AB","#EFC4CE"))
  
  #### BC vs PL
ggplot(area,aes(colour=dis,fill=dis))+geom_density(aes(x=Area,alpha=0.3))+facet_grid(.~Status)+
  theme_bw()+theme(axis.text = element_text(colour = "black"))+scale_color_manual(values=c("#03776A","#EB4B17"))+scale_fill_manual(values=c("#03776A","#EB4B17"))
 


library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(ggforce)
library(dplyr)
library(ggsci)
library(ggpubr)
compaired=list(c("PAT","SAT"))

PL=apply(PL, 2, function(x){mean(x,na.rm=T)}) 
tmp<-strsplit(names(PL),"_")
status<-c()
individual<-c()
for(i in 1:length(tmp)){
  status<-c(status,tmp[[i]][2])
  individual<-c(individual,tmp[[i]][1])
}
dis<-rep("PL",14)
area_PL=data.frame(individual,status,PL)
colnames(area_PL)[3]<-"Area"

df=area_PL %>% pivot_wider(names_from = status,values_from = Area)
df1=df %>% mutate(PAT1=PAT,SAT1=SAT)
df1=df1[,c("individual","PAT","PAT1","SAT","SAT1")]
df2=df1 %>% pivot_longer(-1,names_to = "status",values_to = "Area") 
area_bezier=df2 %>% mutate(Group=case_when(status=="PAT"~1.0,
                                           status=="PAT1"~1.3,
                                           status=="SAT"~1.7,
                                           status=="SAT1"~2.0))


######----plot
ggplot() +
  geom_boxplot(data=area_PL,aes(x=status,y=Area,fill=status),
               width=0.35,position = position_dodge(0),size=0.1,outlier.size = 0) +
  geom_point(data=area_PL,aes(x=status,y=Area,fill=status),shape=21,colour="black",size=2)+
  geom_bezier(data=area_bezier,aes(x=Group,y=Area, group=individual, linetype = 'cubic'),size=0.25,colour="grey20") + #??Ҫ????ggforce??
  scale_fill_manual(values=c("#EFC4CE","#2775AB"))+
  theme_classic2()+
  theme(axis.text = element_text(colour = "black"))+
  stat_compare_means(comparisons = compaired)+ylim(c(700,1800))



BC=apply(BC, 2, function(x){mean(x,na.rm=T)}) 
tmp<-strsplit(names(BC),"_")
status<-c()
individual<-c()
for(i in 1:length(tmp)){
  status<-c(status,tmp[[i]][2])
  individual<-c(individual,tmp[[i]][1])
}
dis<-rep("BC",14)
area_BC=data.frame(individual,status,BC)
colnames(area_BC)[3]<-"Area"

df=area_BC %>% pivot_wider(names_from = status,values_from = Area)
df1=df %>% mutate(PAT1=PAT,SAT1=SAT)
df1=df1[,c("individual","PAT","PAT1","SAT","SAT1")]
df2=df1 %>% pivot_longer(-1,names_to = "status",values_to = "Area") 
area_bezier=df2 %>% mutate(Group=case_when(status=="PAT"~1.0,
                                           status=="PAT1"~1.3,
                                           status=="SAT"~1.7,
                                           status=="SAT1"~2.0))


######----plot
ggplot() +
  geom_boxplot(data=area_BC,aes(x=status,y=Area,fill=status),
               width=0.35,position = position_dodge(0),size=0.1,outlier.size = 0) +
  geom_point(data=area_BC,aes(x=status,y=Area,fill=status),shape=21,colour="black",size=2)+
  geom_bezier(data=area_bezier,aes(x=Group,y=Area, group=individual, linetype = 'cubic'),size=0.25,colour="grey20") + #??Ҫ????ggforce??
  scale_fill_manual(values=c("#EFC4CE","#2775AB"))+
  theme_classic2()+
  theme(axis.text = element_text(colour = "black"))+
  stat_compare_means(comparisons = compaired)+ylim(c(700,1800))