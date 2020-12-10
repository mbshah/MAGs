# Title     : step8-2
# Objective : Run EcoUtils to identify generalists and specialists
# Created by: hb0358
# Created on: 10/12/20
#install.packages("vegan")
#install.packages("spaa")
#devtools::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)
library(RCurl)
x<-"/home/hb0358/PycharmProjects/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles/Bacteria/l6_species_space_replaced2.csv"
comm.tab<-read.table(file=x,sep="\t",row.names=1,header=TRUE,comment.char="@")
comm.tab<-t(comm.tab)
comm.tab<-comm.tab[,which(colSums(comm.tab)>0)]
res<-spec.gen(comm.tab,n=1000)

comm.tab.bin<-ceiling(comm.tab/max(comm.tab))
x<- plot(colSums(comm.tab),colSums(comm.tab.bin)/dim(comm.tab.bin)[1],col=res$sign,pch=19,log="x",xlab="Abundance",ylab="Occurrence")
x<-legend("bottomright",levels(res$sign),col=1:3,pch=19,inset=0.01,cex=0.7)