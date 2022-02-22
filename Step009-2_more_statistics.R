

# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(plyr)
library(ggridges)
library(car)
require(data.table)
library(ggpubr)
library(vegan)
#library(modeest)
#library(EValue)
#tsv_file<-"/home/manan/disk2/UDE/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles2/tax_annot_gs.csv"


#data<-as.data.frame(fread(tsv_file))
#data$perc_coding<-data$total_len_ofCDS*100/data$total_Len_Scaffolds
#data$cds_per_kb<-data$Total_number_of_cds/(data$total_Len_Scaffolds/1000)
#colours<-c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677") #colorblind_pallate from Paul Tol https://davidmathlogic.com/colorblind/
#data$color<-ifelse(data$Classification=="GEN",colours[1],ifelse(data$Classification=="SPEC",colours[3],colours[5]))
#data$cl<-ifelse(data$Classification=="GEN",1,ifelse(data$Classification=="SPEC",2,3))
#gen_gc<-data[data$Classification=="GEN",]$GC
#spec_gc<-data[data$Classification=="SPEC",]$GC
#ns_gc<-data[data$Classification=="NS",]$GC
#evalues.HR(0.86, 0.75, 0.99, rare = FALSE)

#fit<-lm(GC~cds_per_kb,data=data)
#plot(data$cl,data$GC,col=data$color)


# Plot
#data2:entire data set
#data3:top 600 of gen and 600 of spec
#data4: filtered
#gs_class: gens_spec_output of ecolutils
#bac_profile: output profile of kraken
#filtered_profile: bac_profile_filtered as per data3
#t_f_p: taransformed filtered profile
#save.image("/home/hb0358/PycharmProjects/mbs_general/.RDataFiles/725233088.RData")
#load(paste0(str_replace(getwd(),"/MAGs",""),"/MAGs/.RDataFiles/entire.RData"))

#folder_source<-"/home/manan/disk2/UDE/mbs_general/WMG/"
folder_source<-"/home/hb0358/PycharmProjects/mbs_general/WMG/"
t2<- paste0(folder_source, "Kraken2_NT_Scaffolds_profiles2/tax_annot_ko_summary.tsv")
t3<-paste0(folder_source, "Kraken2_NT_Scaffolds_profiles2/Bacteria/Genralists_Specialists/gen_spec1000.tsv")
t4<-paste0(folder_source,"Kraken2_NT_Scaffolds_profiles2/Bacteria/l6_species_normalized_deseq2.tsv")
t5<-paste0(getwd(),"/MAGs/Assembled_Metagenomes/CED91220.tsv")
t6<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/test.tsv")

gs_class<-as.data.frame(fread(t3))
bac_profile<-as.data.frame(fread(t4))
rownames(bac_profile)<-bac_profile$V1
bac_profile$V1<-NULL
meta_data<-as.data.frame(fread(t5))
imp_parameters<-as.data.frame(fread(t6))
colours<-c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677") #colorblind_pallate from Paul Tol https://davidmathlogic.com/colorblind/
data<-as.data.frame(fread(t2))
data2<-merge(data,gs_class,by.x="Inferred_Lineage",by.y="V1")
data2$perc_coding<-data2$total_len_ofCDS*100/data2$total_Len_Scaffolds
data2$perc_noncoding<-(data2$total_Len_Scaffolds-data2$total_len_ofCDS)*100/data2$total_Len_Scaffolds
data2$cds_per_kb<-data2$Total_number_of_cds/(data2$total_Len_Scaffolds/1000)
data2$color<-ifelse(data2$Classification=="GEN",colours[1],ifelse(data2$Classification=="SPEC",colours[3],colours[5]))
data2$cl<-ifelse(data2$Classification=="GEN",1,ifelse(data2$Classification=="SPEC",2,3))
data2$dev<-(data2$mean.simulated-data2$observed)
data2<-data2[order(data2$dev),]
data2$perc_of_scaf_with_cds<-100*data2$Scaffolds_with_CDS/data2$Number_fo_Scaffolds
data2$total_classes<-data2$kc_Metabolism+data2$kc_Cellular_Processes+data2$kc_Environmental_Information_Processing+
  data2$kc_Genetic_Information_Processing+data2$kc_Human_Diseases+data2$kc_Organismal_Systems
data3<-rbind(data2 %>% top_n(data2$dev,n=-600),data2 %>% top_n(data2$dev,n=600))
#write.csv(data3, file=paste0(folder_source,"Kraken2_NT_Scaffolds_profiles2/Bacteria/Genralists_Specialists/tax_annot.csv"),quote = F)
data4<-read.csv(paste0(folder_source,"Kraken2_NT_Scaffolds_profiles2/Bacteria/Genralists_Specialists/tax_annot2.csv"))
data4$Classification<-as.factor(data4$Classification)
data5<-data4
data5[data5==0]<-NA
kegg_pathways<-as.array(names(data4)[grep("^ko",names(data4))])
kc_classes<-as.array(names(data4)[grep("kc_",names(data4))])
data4$no_of_KOs
plot(x=rowSums(data4[,grepl("^K",names(data4))]),y=data4$no_of_KOs,col=data4$color)
data6<- data4[!(rowSums(data4[,grepl("^K",names(data4))])==0),]
#
#
plot(data4$total_Len_Scaffolds,data4$no_of_KOs,cols=data4$Classification)
ggplot(data4, aes(x=kc_Translation,y=Classification, fill=stat(x),color=color)) +
  geom_density_ridges_gradient(rel_min_height=0.001,scale=0.9)+
  scale_fill_viridis_c(name="No of Genes") +
  scale_colour_identity(name="Classification",
                        breaks=levels(as.factor(data4$color)),
                        labels=levels(as.factor(data4$Classification)),guide="legend")+
  geom_vline(data=grpdata, aes(xintercept=grp.mean, color=color),linetype="solid")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  axis.title=element_text(size=18,face="bold"))

aov(kc_Folding_sorting_and_degradation~Classification,data = data4)
plot(data4$cds_per_kb,data4$Total_number_of_cds,col=data4$Classification)
shapiro.test(data4$kc_Development_and_regeneration)


grpdata<-ddply(data4,"color",summarise,
               grp.mode=densMode(kc_Translation),
               grp.mean=mean(kc_Translation),
               grp.median=median(kc_Translation))
ggplot(data4, aes(x=Total_number_of_cds,color=color)) +
  geom_density(alpha=0.8)+
  scale_colour_identity(name="Classification",
                        breaks=levels(as.factor(data4$color)),
                        labels=levels(as.factor(data4$Classification)),guide="legend")+
  #geom_vline(data=grpdata, aes(xintercept=grp.mean, color=color), linetype="dashed")+
  #geom_vline(data=grpdata, aes(xintercept=grp.median, color=color),linetype="solid")+
  #geom_vline(data=grpdata, aes(xintercept=grp.mode, color=color),linetype="dotted")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  axis.title=element_text(size=18,face="bold"))+xlim(c(0,1000))

ggplot(data4, aes(x=log(no_of_KOs), y=log(total_Len_Scaffolds), label=Classification,color=color)) +
  geom_point()  +
  scale_colour_identity(name="Classification",
                        breaks=levels(as.factor(data4$color)),
                        labels=levels(as.factor(data4$Classification)),guide="legend")+
  labs(title="Correlation between no. of scaffolds and predicted genes",x="log KO Counts (unique)",y="Log of Length of Scaffolds") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  axis.title=element_text(size=15),plot.title = element_text(size=18,face="bold",hjust=0.5))+
  coord_flip()


cor.test(data4$kc_Folding_sorting_and_degradation, data4$cl)
cor.test(data4$kc_Development_and_regeneration, data4$cl)

densMode <- function(x){
    td <- density(x)
    maxDens <- which.max(td$y)
    td$x[maxDens]
}

##rarefication curve
bp2<-t(floor(bac_profile))
b_n<-specnumber(bp2)
(raremax<-min(rowSums(bp2)))
Bpr<-rarefy(floor(bp2),raremax)
plot(b_n,Bpr)
abline(0,1)
rarecurve(bp2,step=20,sample=raremax,col="blue",cex=0.6)

#PCA_plot
library(ggfortify)
pca_res<-prcomp(bp2)
autoplot(pca_res)

##data4<-data3[data2$Classification %in% c("GEN","SPEC"),]
##data5<-data4[data4$total_Len_Scaffolds>=100000,]
#library(corrplot)
#corrplot(cor(data3[c("Metabolism","Environmental_Information_Processing","GC")]), method = "circle")
#
#data3 %>%
#  ggplot( aes(x=cl, y=kc_Cell_growth_and_death, fill=Classification)) +
#    geom_boxplot(outlier.shape = NA) +
#    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
#    geom_jitter(color="black", size=0.4, alpha=0.9) +
#    theme_ipsum() +
#    theme(
#      legend.position="none",
#      plot.title = element_text(size=11)
#    ) +
#    ggtitle("GC") +
#    xlab("")
##density plots
#library("ggplot2")
#library("plyr")
#grpdata<-ddply(data3,"Classification",summarise,grp.mean=mean(GC))
#grpdata
#ggplot(data3, aes(x=GC, color=Classification)) +
#  geom_density()
#  #geom_vline(data=grpdata, aes(xintercept=grp.mean, color=Classification),
#  #           linetype="dashed")
#tsv_file
#fit<-lm(mean.simulated~dev,data=data3)
#plot(fit,col=data3$color)
#plot(data$cl,data$GC,col=data$color)
#nrow(gs_class[gs_class$sign=="SPECIALIST",])
#nrow(gs_class[gs_class$sign=="GENERALIST",])
#filtered_profile<-subset(bac_profile,V1 %in% data3$Inferred_Lineage)
#rownames(filtered_profile)<-filtered_profile[,1]
#filtered_profile<-filtered_profile[,-1]
#filtered_profile2<-filtered_profile[order(match(rownames(filtered_profile),data3$Inferred_Lineage)),]
#t_f_p<-as.data.frame(t(filtered_profile))
#colnames(t_f_p)<-1:200
#colnames(t_f_p)
#for_cca<-merge(meta_data[,c("ID","pH1","conductivity","temperature","TP","DP","DOC","DN","DRSi","elevation (m)")],t_f_p,by.x="ID",by.y=0)
#write.csv(for_cca, file="/home/hb0358/PycharmProjects/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles2/Bacteria/Genralists_Specialists/for_cca.csv",quote = F)
#
#BiocManager::install('mixOmics')
#library(mixOmics)
#
#
#
#
#
#
#
#
#
#
#
#
#data.tab<-within(data,{
#  Classification<-as.factor(Classification)
#  cl<-as.factor(cl)
#})
#detach(data.tab)
#attach(data.tab)
#kruskal.test(GC,Classification)
#summary(aov(GC~Classification))
#cor.test(data4$kc_Folding_sorting_and_degradation, data4$cl)
#t.test(Classification,GC,paired=TRUE)
#Classification
#
