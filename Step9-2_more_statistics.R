# Title     : TODO
# Objective : TODO
# Created by: hb0358
# Created on: 05/02/21

tsv_file<-"/home/hb0358/PycharmProjects/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles2/tax_annot.tsv"
require(data.table)

data<-as.data.frame(fread(tsv_file))
dat<-read.table(tsv_file,sep="\t",row.names=1,header=TRUE,comment.char="@")
data.tab2<-subset(dat,"Classification!"="NS")
data.tab<-within(read.table(tsv_file,sep="\t",row.names=1,header=TRUE,comment.char="@"),{
  Classification<-as.factor(Classification)
})
detach(data.tab2)
attach(data.tab2)
kruskal.test(GC,Classification)
summary(aov(GC~Classification))
cor.test(GC, Classification)
Classification
