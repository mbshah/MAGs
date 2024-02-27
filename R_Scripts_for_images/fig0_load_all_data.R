library(ComplexHeatmap)
library(tidyverse)
library(plyr)
library(data.table)
library(ggplot2)
library(vegan)
library(gridExtra)
library(sjmisc)
library(gridExtra)
library(ggpmisc)
library(epade)
library(weights)
library(ggpubr)
library(gghighlight)
library(RColorBrewer)
library(viridis)

#save.image(paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/.RDataFiles/entire.RData"))
#load(paste0(str_replace(getwd(),"/MAGs",""),"/MAGs/.RDataFiles/entire.RData"))

#### DEFINE and READ TABLES
#folder_source<-"/home/manan/disk2/UDE/mbs_general/"
good_bins_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/dRep__gff/post_cluster_data.csv")
metadata_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/Assembled_Metagenomes/CED91220.csv")
raw_md_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/Assembled_Metagenomes/Raw_core_env_data.csv")
members_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/dRep__gff/members_table_super_populated.csv")
abundance_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/dRep__gff/cluster_abundance_profile_normalized_deseq2.tsv")
metbolCompl_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/dRep__gff/MicrobeAnnotator_out/metabolic_summary__module_completeness.tsv")
temp_prec_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/Assembled_Metagenomes/temp_precip.csv")
metabolic_redundancy_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/dRep__gff/kegg_module_redundancy.csv")
metabolic_compl_recal_file <- paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/dRep__gff/kegg_reCal_module_completeness.csv")
metabolic_corelations_file <- paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/dRep__gff/Module_env_Corr.csv")
kegg_ko_profile_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/dRep__gff/post_cluster_data_ko_summary.tsv")
kegg_pathways_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/dRep__gff/post_cluster_data_pathway_summary.tsv")
kegg_subsystem_file<-paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/dRep__gff/post_cluster_data_subSystem_summary.tsv")
MAGs_table<-as.data.frame(fread(good_bins_file))
meta_data<-as.data.frame(fread(metadata_file))
raw_md<-as.data.frame(fread(raw_md_file))
members_table<-as.data.frame(fread(members_file))
metabolCompl<-as.data.frame(fread(metbolCompl_file))
temp_precip<-as.data.frame(fread(temp_prec_file))
ko_profile<-as.data.frame(fread(kegg_ko_profile_file))%>%remove_rownames %>%column_to_rownames(var = "cluster_ID")
pathway_profile<-as.data.frame(fread(kegg_pathways_file))%>%remove_rownames %>%column_to_rownames(var = "cluster_ID")
subsystemprofile<- as.data.frame(fread(kegg_subsystem_file))%>%remove_rownames %>%column_to_rownames(var = "cluster_ID")
metabolic_redundandancy<-read.csv(metabolic_redundancy_file,row.names=1)
metabolic_completeness_new<-read.csv(metabolic_compl_recal_file,row.names = 1 )

##ADD and MODIFY tables
###Changed two names total in raw files:  O151Bu to O151BU      and S031BU to S301BU   in three files raw_md; meta_data; temp_precip in the folder stated above only not original file
row.names(metabolCompl)<-metabolCompl$module
metabolCompl$module<-NULL
abundance_profile<-as.data.frame(fread(abundance_file))
row.names(abundance_profile)<-abundance_profile$abundance
abundance_profile$abundance<-NULL

MAGs_table["tp_lev"]<-ifelse(MAGs_table$weighted_TP>50,"High","Low")
MAGs_table["tns_lev"]<-ifelse(MAGs_table$weighted_TN>1000,"High","Low")
MAGs_table["doc_lev"]<-ifelse(MAGs_table$weighted_TN>3800,"High","Low")
MAGs_table["size_lev"]<-ifelse(MAGs_table$reassembly_size_mb<2,"Small",
                               ifelse(MAGs_table$reassembly_size_mb>4,"Big","Medium"))
MAGs_table["temp_lev"]<-ifelse(MAGs_table$An_Max_Temp>20,"High","Low")
members_table["tp_lev"]<-ifelse(members_table$TP.y>20,"High","Low")
members_table["tn_lev"]<-ifelse(members_table$TNs>600,"High","Low")
members_table["doc_lev"]<-ifelse(members_table$DOC.y>4000,"High","Low")
members_table["size_lev"]<-ifelse(members_table$size_rep<2,"Small",
                               ifelse(members_table$size_rep>4,"Big","Medium"))
members_table["temp_lev"]<-ifelse(members_table$An_Max_Temp>20,"High","Low")
MAGs_table["faa_ko"]<-paste0(MAGs_table$Representative_Genome, ".faa.ko")
MAGs_table<-MAGs_table[order(MAGs_table[,ncol(MAGs_table)],decreasing = F),]
row.names(MAGs_table)<-MAGs_table$cluster_ID
MAGs_table["KOs_for_ABC_Transporter"]<-pathway_profile$ko02010[match(row.names(MAGs_table), row.names(pathway_profile))]
MAGs_table["kc_Cell_motility"]<-subsystemprofile$kc_Cell_motility[match(row.names(MAGs_table), row.names(subsystemprofile))]



####Adding redundancy and completeness to MAGs_table
metabolic_redundandancy_test <-metabolic_redundandancy %>%
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))
rownames(metabolic_redundandancy_test)[rownames(metabolic_redundandancy_test) == "...395"] <- "Total_Redundancy"
metabolic_redundandancy_totals<-metabolic_redundandancy_test["Total_Redundancy",]
MAGs_table_new<-merge(MAGs_table,t(metabolic_redundandancy_totals), by.x="cluster_ID", by.y="row.names")
MAGs_table_new$c80<-NA
MAGs_table_new$c100<-NA
for (rowID in row.names(MAGs_table_new)){
  ClusterID<-MAGs_table_new[rowID,"cluster_ID"]
  redce_c80<-metabolic_completeness_new[ClusterID]>80
  c80_count_tbl<-count(redce_c80)
  completness_count<-c80_count_tbl[c80_count_tbl[ClusterID]==TRUE,"freq"]
  MAGs_table_new[rowID,"c80"]<-completness_count
  redce_c100<-metabolic_completeness_new[ClusterID]==100
  c100_count_tbl<-count(redce_c100)
  completness_count<-c100_count_tbl[c100_count_tbl[ClusterID]==TRUE,"freq"]
  MAGs_table_new[rowID,"c100"]<-completness_count
  #print(complete_completness)
}
MAGs_table_new["No_of_KOs"]<-ko_profile$no_of_KOs[match(MAGs_table_new$cluster_ID,row.names(ko_profile))]
MAGs_table_new["pc_cell_motility"]<-100*MAGs_table_new$kc_Cell_motility/MAGs_table_new$No_of_KOs
MAGs_table_new["pc_abc_transporters"]<-100*MAGs_table_new$KOs_for_ABC_Transporter/MAGs_table_new$No_of_KOs
MAGs_table_new["pc_sigmafactors"]<-100*MAGs_table_new$sigma_factor/MAGs_table_new$no_of_CDS
members_table["kc_Cell_motility"]<-MAGs_table_new$kc_Cell_motility[match(members_table$cluster, MAGs_table_new$cluster_ID)]
members_table["KOs_for_ABC_Transporter"]<-MAGs_table_new$KOs_for_ABC_Transporter[match(members_table$cluster, MAGs_table_new$cluster_ID)]
members_table["sigma_factor"]<-MAGs_table_new$sigma_factor[match(members_table$cluster, MAGs_table_new$cluster_ID)]
members_table["redundant_genes"]<-MAGs_table_new$Total_Redundancy[match(members_table$cluster, MAGs_table_new$cluster_ID)]
members_table["rep_codingpc"]<-MAGs_table_new$coding_perc[match(members_table$cluster, MAGs_table_new$cluster_ID)]

#write.csv(MAGs_table_new,paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/dRep__gff/post_cluster_data_Rver.csv"))