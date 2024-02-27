###LOAD DATA USING 005-4 script
library(tidyverse)
library(plyr)
library(data.table)
library(ggplot2)
library(vegan)

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
MAGs_table<-as.data.frame(fread(good_bins_file))
meta_data<-as.data.frame(fread(metadata_file))
raw_md<-as.data.frame(fread(raw_md_file))
members_table<-as.data.frame(fread(members_file))
metabolCompl<-as.data.frame(fread(metbolCompl_file))
temp_precip<-as.data.frame(fread(temp_prec_file))
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









#KEGG ENRICHMENT ANALYSIS SIZE
# start with KO Profile
ko_profile_file<-"/home/hb0358/PycharmProjects/mbs_general/MAGs/outfile/dRep__gff/ko_profile.tsv"
koProfile<-as.data.frame(fread(ko_profile_file))
row.names(koProfile)<-koProfile$Data
koProfile$Data<-NULL
ko_Small<-koProfile[row.names(koProfile) %in% MAGs_table[MAGs_table$size_lev=='Small',]$cluster_ID, ]
ko_Big<-koProfile[row.names(koProfile) %in% MAGs_table[MAGs_table$size_lev=='Big',]$cluster_ID, ]
ko_Medium<-koProfile[row.names(koProfile) %in% MAGs_table[MAGs_table$size_lev=='Medium',]$cluster_ID, ]
ko_differential_abundance<-data.frame(Small=double(),Big=double(), Intermediate=double())
for (KID in colnames(koProfile)){
  y<-mean(ko_Big[,KID])
  x<-mean(ko_Small[,KID])
  z<-mean(ko_Medium[,KID])
  ko_differential_abundance[KID,]<-c(x,y,z)
}
KOs<-ko_differential_abundance
KOs$KO<-row.names(KOs)
BiocManager::install("gage")
library(gage)
library(pathview)
library(dplyr)

more<-0
less<-0
KOs <- read.table("/home/hb0358/sciebo/mag+gen_spec/MAGs/MAGs_compiled_Paper/tables/gene_differential_presence.tsv", header=T, sep="\t", as.is = T)
# KO should be present in at least 10% of the genomes
KOs2 <- KOs[apply(KOs[,2:3], 1, function(x) any(x>0.1)),] # 2839
# presence in genomes at least twice at high
small <- KOs2[log2((KOs2$Small+0.0001)/(KOs2$Big+0.0001))>=1,'KO'] # 99
big <- KOs2[log2((KOs2$Small+0.0001)/(KOs2$Big+0.0001))<=(-1),'KO'] # 1617
both <- setdiff(setdiff(KOs2$KO, more), less) # 1123
all <- KOs$KO # 5777

KEGGEnrichment <- function (sign.ko, background.ko, kegg.ko)
{
  annList <- list()
  for (i in 1:length(kegg.ko))
  {
    pathway <- names(kegg.ko)[i]
    pathway.ko <- kegg.ko[[i]]
    annotatedMoleculeList <- intersect(pathway.ko, sign.ko)
    annotatedBackgroundList <- intersect(pathway.ko, background.ko)
    ann <- list(pathwayName = "not known", annMoleculeList = character(), annMoleculeNumber = 0, annBgMoleculeList = character(), annBgNumber = 0, moleculeNumber = 0, bgNumber = 0, pvalue = 1)
    ann$pathwayName <- pathway
    ann$annMoleculeList <- annotatedMoleculeList
    ann$annMoleculeNumber <- length(annotatedMoleculeList)
    ann$annBgMoleculeList <- annotatedBackgroundList
    ann$annBgNumber <- length(annotatedBackgroundList)
    ann$moleculeNumber <- length(sign.ko)
    ann$bgNumber <- length(background.ko)
    ann$pvalue <- 1 - phyper(ann$annMoleculeNumber - 1, ann$annBgNumber, ann$bgNumber - ann$annBgNumber, ann$moleculeNumber)
    annList[[i]] <- ann
  }
  return(annList)
}
signKEGGEnrichment <- function (KEGG.enrichment, threshold)
{
  p.value <- sapply(KEGG.enrichment, function(x) return(x$pvalue))
  pathway.table <- unname(t(as.data.frame(lapply(KEGG.enrichment[p.value <= threshold], function(x) c(x$pathwayName, paste(x$annMoleculeList, collapse = ", "), x$pvalue)))))
  if (dim(pathway.table)[1] != 0)
  {
    colnames(pathway.table) <- c("Pathway", "KOs", "P-value")
    pathway.table <- pathway.table[order(as.numeric(pathway.table[, "P-value"])), ]
    return(pathway.table)
  }else
  {
    return(NULL)
  }
}

background.ko <-all
sign.ko <- big
kegg.ko  <- kegg.gsets(species = "ko", id.type = "kegg")$kg.sets
threshold <- 0.01

KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
sign.KEGG.enrichment.big <- signKEGGEnrichment(KEGG.enrichment, threshold)

sign.ko <- small
KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
sign.KEGG.enrichment.small <- signKEGGEnrichment(KEGG.enrichment, threshold)

sign.ko <- both
KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
sign.KEGG.enrichment.both <- signKEGGEnrichment(KEGG.enrichment, threshold)

write.table(sign.KEGG.enrichment.big, file = "MAGs/outfile/dRep__gff/KEGG_Enrichment_big.tsv", col.names = NA, sep = "\t")
write.table(sign.KEGG.enrichment.small, file = "MAGs/outfile/dRep__gff/KEGG_Enrichment_small.tsv", col.names = NA, sep = "\t")
write.table(sign.KEGG.enrichment.both, file = "MAGs/outfile/dRep__gff/KEGG_Enrichment_both.tsv", col.names = NA, sep = "\t")

write.table(big, file= "../outfile/dRep__gff/KOs_big.tsv")
write.table(small, file= "../outfile/dRep__gff/KOs_small.tsv")
write.table(both, file= "../outfile/dRep__gff/KOs_both.tsv")

###END KEGG ENRICHMENT ANALYSIS####

###LOAD DATA USING 005-4 script

#KEGG ENRICHMENT ANALYSIS SIZE @@@@@@2 factors

KEGGEnrichment <- function (sign.ko, background.ko, kegg.ko)
{
  annList <- list()
  for (i in 1:length(kegg.ko))
  {
    pathway <- names(kegg.ko)[i]
    pathway.ko <- kegg.ko[[i]]
    annotatedMoleculeList <- intersect(pathway.ko, sign.ko)
    annotatedBackgroundList <- intersect(pathway.ko, background.ko)
    ann <- list(pathwayName = "not known", annMoleculeList = character(), annMoleculeNumber = 0, annBgMoleculeList = character(), annBgNumber = 0, moleculeNumber = 0, bgNumber = 0, pvalue = 1)
    ann$pathwayName <- pathway
    ann$annMoleculeList <- annotatedMoleculeList
    ann$annMoleculeNumber <- length(annotatedMoleculeList)
    ann$annBgMoleculeList <- annotatedBackgroundList
    ann$annBgNumber <- length(annotatedBackgroundList)
    ann$moleculeNumber <- length(sign.ko)
    ann$bgNumber <- length(background.ko)
    ann$pvalue <- 1 - phyper(ann$annMoleculeNumber - 1, ann$annBgNumber, ann$bgNumber - ann$annBgNumber, ann$moleculeNumber)
    annList[[i]] <- ann
  }
  return(annList)
}
signKEGGEnrichment <- function (KEGG.enrichment, threshold)
{
  p.value <- sapply(KEGG.enrichment, function(x) return(x$pvalue))
  pathway.table <- unname(t(as.data.frame(lapply(KEGG.enrichment[p.value <= threshold], function(x) c(x$pathwayName, paste(x$annMoleculeList, collapse = ", "), x$pvalue)))))
  if (dim(pathway.table)[1] != 0)
  {
    colnames(pathway.table) <- c("Pathway", "KOs", "P-value")
    pathway.table <- pathway.table[order(as.numeric(pathway.table[, "P-value"])), ]
    return(pathway.table)
  }else
  {
    return(NULL)
  }
}

##loop to run for 5 individual phylas
req_phylas<-c("Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota","Verrucomicrobiota")
req_factor<-"tn"
for (req_phyla in req_phylas){
  members_table_sub<-members_table[members_table$gtdb_phylum==req_phyla,]
  # start with KO Profile
  ko_profile_file<-"/home/hb0358/PycharmProjects/mbs_general/MAGs/outfile/dRep__gff/ko_profile.tsv"
  koProfile<-as.data.frame(fread(ko_profile_file))
  row.names(koProfile)<-koProfile$Data
  koProfile$Data<-NULL
  ko_low<-koProfile[row.names(koProfile) %in% unique(sort(members_table_sub[members_table_sub$tn_lev=="Low",]$cluster)),]   #replace tp/tn/doc here
  ko_high<-koProfile[row.names(koProfile) %in% unique(sort(members_table_sub[members_table_sub$tn_lev=="High",]$cluster)),] #replace tp/tn/doc here
  ko_differential_abundance<-data.frame(Low=double(),High=double())
  for (KID in colnames(koProfile)){
    x<-mean(ko_low[,KID])
    y<-mean(ko_high[,KID])
    ko_differential_abundance[KID,]<-c(x,y)
  }
  KOs<-ko_differential_abundance
  KOs$KO<-row.names(KOs)
  #BiocManager::install("gage")
  library(gage)
  library(pathview)
  library(dplyr)

  more<-0
  less<-0
  #KOs <- read.table("/mnt/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/MAGs_compiled_Paper/tables/gene_differential_presence.tsv", header=T, sep="\t", as.is = T)
  # KO should be present in at least 10% of the genomes
  KOs2 <- KOs[apply(KOs[,-which(names(KOs) == "KO")], 1, function(x) any(x>0.1)),] # 2839 $2815
  # presence in genomes at least twice at high
  small <- KOs2[log2((KOs2$Low+0.0001)/(KOs2$High+0.0001))>=1,'KO'] # 99
  big <- KOs2[log2((KOs2$Low+0.0001)/(KOs2$High+0.0001))<=(-1),'KO'] # 1617
  both <- setdiff(setdiff(KOs2$KO, more), less) # 1123
  all <- KOs$KO # 5777

  background.ko <-all
  sign.ko <- big
  kegg.ko  <- kegg.gsets(species = "ko", id.type = "kegg")$kg.sets
  threshold <- 0.01

  KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
  sign.KEGG.enrichment.big <- signKEGGEnrichment(KEGG.enrichment, threshold)

  sign.ko <- small
  KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
  sign.KEGG.enrichment.small <- signKEGGEnrichment(KEGG.enrichment, threshold)

  sign.ko <- both
  KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
  sign.KEGG.enrichment.both <- signKEGGEnrichment(KEGG.enrichment, threshold)

  write.table(sign.KEGG.enrichment.big, file = paste0("MAGs/outfile/dRep__gff/KEGG_Enrichment_Ananlyses/KEGG_Enrichment_High_",req_factor,"_",req_phyla,".tsv"), col.names = NA, sep = "\t")
  write.table(sign.KEGG.enrichment.small, file = paste0("MAGs/outfile/dRep__gff/KEGG_Enrichment_Ananlyses/KEGG_Enrichment_Low_",req_factor,"_",req_phyla,".tsv"), col.names = NA, sep = "\t")
  write.table(sign.KEGG.enrichment.both, file = paste0("MAGs/outfile/dRep__gff/KEGG_Enrichment_Ananlyses/KEGG_Enrichment_Both_",req_factor,"_",req_phyla,".tsv"), col.names = NA, sep = "\t")

}




#for all phylastogether
members_table_sub<-members_table
# start with KO Profile
ko_profile_file<-"/home/hb0358/PycharmProjects/mbs_general/MAGs/outfile/dRep__gff/ko_profile.tsv"
koProfile<-as.data.frame(fread(ko_profile_file))
row.names(koProfile)<-koProfile$Data
koProfile$Data<-NULL
ko_low<-koProfile[row.names(koProfile) %in% unique(sort(members_table_sub[members_table_sub$doc_lev=="Low",]$cluster)),]
ko_high<-koProfile[row.names(koProfile) %in% unique(sort(members_table_sub[members_table_sub$doc_lev=="High",]$cluster)),]
ko_differential_abundance<-data.frame(Low=double(),High=double())
for (KID in colnames(koProfile)){
  x<-mean(ko_low[,KID])
  y<-mean(ko_high[,KID])
  ko_differential_abundance[KID,]<-c(x,y)
}
KOs<-ko_differential_abundance
KOs$KO<-row.names(KOs)
#BiocManager::install("gage")
library(gage)
library(pathview)
library(dplyr)

more<-0
less<-0
#KOs <- read.table("/mnt/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/MAGs_compiled_Paper/tables/gene_differential_presence.tsv", header=T, sep="\t", as.is = T)
# KO should be present in at least 10% of the genomes
KOs2 <- KOs[apply(KOs[,-which(names(KOs) == "KO")], 1, function(x) any(x>0.1)),] # 2839 $2815
# presence in genomes at least twice at high
small <- KOs2[log2((KOs2$Low+0.0001)/(KOs2$High+0.0001))>=1,'KO'] # 99
big <- KOs2[log2((KOs2$Low+0.0001)/(KOs2$High+0.0001))<=(-1),'KO'] # 1617
both <- setdiff(setdiff(KOs2$KO, more), less) # 1123
all <- KOs$KO # 5777

background.ko <-all
sign.ko <- big
kegg.ko  <- kegg.gsets(species = "ko", id.type = "kegg")$kg.sets
threshold <- 0.01

KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
sign.KEGG.enrichment.big <- signKEGGEnrichment(KEGG.enrichment, threshold)

sign.ko <- small
KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
sign.KEGG.enrichment.small <- signKEGGEnrichment(KEGG.enrichment, threshold)

sign.ko <- both
KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
sign.KEGG.enrichment.both <- signKEGGEnrichment(KEGG.enrichment, threshold)

write.table(sign.KEGG.enrichment.big, file = paste0("MAGs/outfile/dRep__gff/KEGG_Enrichment_Ananlyses/KEGG_Enrichment_High_doc_",req_phyla,".tsv"), col.names = NA, sep = "\t")
write.table(sign.KEGG.enrichment.small, file = paste0("MAGs/outfile/dRep__gff/KEGG_Enrichment_Ananlyses/KEGG_Enrichment_Low_doc_",req_phyla,".tsv"), col.names = NA, sep = "\t")
write.table(sign.KEGG.enrichment.both, file = paste0("MAGs/outfile/dRep__gff/KEGG_Enrichment_Ananlyses/KEGG_Enrichment_Both_doc_",req_phyla,".tsv"), col.names = NA, sep = "\t")


###END KEGG ENRICHMENT ANALYSIS####