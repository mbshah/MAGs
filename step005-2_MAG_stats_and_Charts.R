# Libraries
install.packages("tidyverse")
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(plyr)
library(ggridges)
library(data.table)
library(ggpubr)
library(units)
library(ggplot2)
#save.image(paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/.RDataFiles/entire.RData"))
#load(paste0(str_replace(getwd(),"/MAGs",""),"/MAGs/.RDataFiles/entire.RData"))

#folder_source<-"/home/manan/disk2/UDE/mbs_general/"
folder_source<-"/home/hb0358/PycharmProjects/mbs_general/MAGs"
all_bins_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/summary_all.tsv")
good_bins_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/clusters_dastool_90_10_average/clusters_wphyp_sited_m.tsv")
dist_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/dist_matrix_dastool_90_10.tsv")
metadata_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/Assembled_Metagenomes/CED91220.csv")
all_bins<-as.data.frame(fread(all_bins_file))
all_bins$qual<-ifelse(all_bins$completeness>90 & all_bins$contamination<10,"good","bad")
good_bins<-as.data.frame(fread(good_bins_file))
dist_matrix<-as.data.frame(fread(dist_file))
meta_data<-as.data.frame(fread(metadata_file))
row.names(dist_matrix)<-dist_matrix$bins
dist_matrix$bins<-NULL




ggplot(all_bins,aes(x=completeness,y=contamination,color=qual))+
  geom_point()+geom_jitter()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  axis.title=element_text(size=18),plot.title=element_text(size=18,face="bold",hjust=0.5))+
  labs(title="Post Binning Quality",x="Completeness",y="Contamination",colour= "Quality")+
  xlim(c(80,100))+ylim(c(0,15))


good_bins$Proteins_with_functional_assignments<-as.numeric(good_bins$`Proteins with functional assignments`)
ggplot(good_bins, aes(x=Proteins_with_functional_assignments)) +
  geom_histogram(fill="dark green",binwidth = 10)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  axis.title=element_text(size=18,face="bold"))






#library(maps)
#library(ggmap)
#library(osmdata)
#library(ggplot2)
#library(plyr)
#library(sf)
#library(extraoperators)

#maps packages required and versions used:
#'ggrepel_0.9.1       OpenStreetMap_0.3.4 data.table_1.14.2   forcats_0.5.1       stringr_1.4.0
#'dplyr_1.0.7         purrr_0.3.4         readr_2.1.0         tidyr_1.1.4         tibble_3.1.6
#'ggplot2_3.3.5       tidyverse_1.3.1
library(OpenStreetMap)
library(tidyverse)
library(data.table)
library(ggrepel)
metadata_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/Assembled_Metagenomes/CED91220.csv")
meta_data<-as.data.frame(fread(metadata_file))
meta_data$index<-1:44
meta_data$label<-paste(meta_data$index,meta_data$ID,sep=":")
meta_data$ele_lvl<-ifelse(meta_data$elevation>1500,"high",ifelse(meta_data$elevation<500,"low","med"))
caption_str<-paste("Map of sampling sites selected for this study; For every dot: Size: signifies the pH, colour: shows temperature variance, and label colours: denote elevation variation (Pink: High altitude, Orange for Intermediate altitude, Green for low altitude); sample names:"
  ,paste(meta_data$label, collapse = "; "), sep=" ")
#attach(meta_data)
##MAG->zoom_level_6 type:osm
##gen_spec->zoom:5 type:apple-iphoto
##good resolution 1256x829
mp<-autoplot.OpenStreetMap(openproj(openmap(c(35.00,-15.56),c(69.16,30),type='bing',zoom=5)))+
  geom_point(data = meta_data,aes(x=longitude,y=latitude,col=temperature,size=pH1))
mp<-autoplot.OpenStreetMap(openproj(openmap(c(35.00,-15.56),c(69.16,30),type='apple-iphoto',zoom=5)))+
  geom_point(data = meta_data,aes(x=longitude,y=latitude,col=temperature,size=pH1))
mp<-autoplot.OpenStreetMap(openproj(openmap(c(35.00,-15.56),c(69.16,30),type='osm',zoom=5)))+
  geom_point(data = meta_data,aes(x=longitude,y=latitude,col=temperature,size=pH1))
mp+labs(colour="Temperature",size="pH")+
  geom_label_repel(data=meta_data,aes(x=longitude,y=latitude,label=index),
                   segment.color="black",
                   color=rep(c("deeppink","darkgreen","orangered")[as.factor(meta_data$ele_lvl)]),
                   seed=354636,fill = alpha("white",0.6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title=element_blank(),plot.title =element_text(size=18,face="bold",hjust=0.5),axis.ticks=element_blank(),
        axis.text.x=element_blank(),axis.text.y=element_blank(),)



plot(meta_data$conductivity,meta_data$elevation)
legend(400,2000, c("High","Medium","Low"),col=c("deeppink","orangered","darkgreen"))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("High","Medium","Low"), pch=16, pt.cex=3, cex=1.5, bty='n',
    col = c("deeppink","orangered","darkgreen"))
mtext("Elevation", at=0.0, cex=2)

###Codon usage bias analysis
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("coRdon")
library(coRdon)

#all_fasta_combined into one
genome_fasta<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/clusters_dastool_90_10_average/","merged_fna_all_clusters.fna")
reads<-readSet(file=genome_fasta)
genome_codon_table<-codonTable(reads)
cc<-codonCounts(genome_codon_table)
row.names(cc)<-good_bins$Cluster
write.table(cc,file="test.tsv",sep="\t",quote=F,col.names=NA)
Cr_milc<-MILC(genome_codon_table)
Cr_br<-B(genome_codon_table)
Cr_mcb<-MCB(genome_codon_table)
median(Cr_mcb[,1])
good_bins$CUB_MCB<-Cr_mcb[,1]
tc<-t(cc)
freqs <- scale(tc, center = FALSE,
               scale = colSums(tc))
freqs_out_table<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/clusters_dastool_90_10_average/","Cluster_Codon_Table.tsv")
write.table(t(freqs),file=freqs_out_table,sep="\t",quote=F,col.names=NA)




#median of all individual
cub_calc<-function(col_data){
  ret_mat<-matrix()
  for(magid in col_data){
    genome_fasta<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/clusters_dastool_90_10_average/",magid,"/",magid,"_gms.fna")
    reads<-readSet(file=genome_fasta)
    genome_codon_table<-codonTable(reads)
    Cr_mcb<-MCB(genome_codon_table)
    ret_mat[magid]<-median(Cr_mcb)
  }
  ret_mat
}
new_col<-cub_calc(good_bins$Cluster)
good_bins$CUB_MCB2<-new_col[good_bins$Cluster]
write.table(good_bins,file=good_bins_file,quote = FALSE,sep="\t",row.names = FALSE)

##Pathways_images
BiocManager::install("pathview")
library(pathview)
library(tidyverse)
library(data.table)
my_color_pallet<-c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499","#882255")
pathway_of_interes<-"ko00790"
kosummary_table<-as.data.frame(fread("MAGs/outfile/clusters_dastool_90_10_average/gene_differential_presence_actino.tsv")) ##file created by python script function:"make_ko_category_table"
row.names(kosummary_table)<-kosummary_table$KO
kosummary_table$KO<-NULL
out<-pathview(gene.data=kosummary_table[,1:3],
              pathway.id ="02010",
              out.suffix = "actino",
              limit=list(gene=1,cpd=1), bins=list(gene=18,cpd=10),
              both.dirs = list(gene=F,cpd=FALSE),
              high = list(gene = "#003AFF", cpd = "#003AFF"),
              mid=list(gene = "#909090", cpd = "#909090"),
              kegg.dir="ancilary/kegg/",
              key.pos="topright",
              gene.idtype = "KEGG",
              plot.col.key = TRUE,
              species = "ko")


#BOx_plot_members_table with weights keep as is for the p-values, the plots itself also implemented in python.
library(tidyverse)
library(data.table)
library(ggpubr)
members_table<-as.data.frame(fread("/home/hb0358/PycharmProjects/mbs_general/MAGs/outfile/dRep__gff/members_table.csv"))
members_table["wt_mean"]<-members_table$org_size*members_table$abundance
members_table["tp_lev"]<-ifelse(members_table$site_TP>10,"high","low")
p<-ggboxplot(members_table,x="tp_lev",y="size_rep")
p+stat_compare_means(comparisons = list(c("low","high")), method = "t.test")

data <- read.table("/home/hb0358/PycharmProjects/mbs_general/MAGs/outfile/dRep__gff/members_table.csv", header=T, sep=",")
data["tp_lev"]<-ifelse(data$site_TP>10,"High","Low")
library(ggplot2)
library(ggpubr)
library(grid)

q <- ggplot(data, aes(y=size_rep, x=tp_lev))+
  geom_boxplot()+
  theme_bw()
q

t_val<-wtd.t.test(subset(data,tp_lev=="low")$size_mem, subset(data,tp_lev=="high")$size_mem, weight=subset(data,tp_lev=="low")$abundance_m, weighty = subset(data,tp_lev=="high")$abundance_m, samedata=F)$coefficients
t_val["p.value"]
p <- ggplot(data, aes(y=size_rep, x=tp_lev))+
  geom_boxplot(aes(weight=abundance_m))+
  theme_classic2()+
  annotate("rect", xmin = 1, xmax = 2, ymin = 8.5, ymax =8.5, alpha=1,colour = "black")+
  geom_text(x=1.5, y=8.6,label=paste0("p = ",t_val["p.value"]))+xlab("Total Phosphorus Level")+
  ylab("Genome Size (Mb)")

p

+grid.text("test",x=unit(0.15, "npc"),y=unit(0.15, "npc"))


+stat_compare_means(comparisons = list(c("low","high")), method = "t.test")

test <- split(data, f=data$tp_lev)

library(weights)

wtd.t.test(subset(data,tp_lev=="Low")$size_mem, subset(data,tp_lev=="High")$size_mem, weight=subset(data,tp_lev=="Low")$abundance_m, weighty = subset(data,tp_lev=="High")$abundance_m, samedata=T)
wtd.t.test()

subset(data,tp_lev=="low")$size_rep