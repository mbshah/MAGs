# Libraries
BiocManager::install("ComplexHeatmap")
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
members_table["kc_Cell_motility"]<-MAGs_table$kc_Cell_motility[match(members_table$bin, MAGs_table$Representative_Genome)]
members_table["KOs_for_ABC_Transporter"]<-MAGs_table$KOs_for_ABC_Transporter[match(members_table$bin, MAGs_table$Representative_Genome)]
members_table["sigma_factor"]<-MAGs_table$sigma_factor[match(members_table$bin, MAGs_table$Representative_Genome)]

##USE ONLY ONCE
setDT(meta_data)[,paste0("ID",1:2):=tstrsplit(ID,"_")]
phy_cond<-merge(x=meta_data,y=raw_md, by.x="ID2",by.y="Code")
phy_cond<-merge(x=phy_cond,y=temp_precip,by.x="ID2",by.y="Code")
new_members<-merge(x=members_table,y=phy_cond,by.x="site_list",by.y="ID")
write.csv(new_members,paste0(str_replace(getwd(),"/MAGs/",""),"/MAGs/outfile/dRep__gff/members_table_super_populated.csv"))


#NMDS abundance
abundance_matrix<-as.matrix(abundance_profile)
set.seed(2718)
abundance_nmds<-metaMDS(abundance_matrix, distance="bray")
plot(abundance_nmds)
data.scores <-as.data.frame(scores(abundance_nmds)$sites)
data.scores$weightedTP<-MAGs_table$weighted_TP
data.scores$tp_lev<-MAGs_table$tp_lev
data.scores$size_lev<-MAGs_table$size_lev

ggplot(data.scores, aes(x = NMDS1, y = NMDS2, color=tp_lev,shape=size_lev))+
  geom_point(size=4)


#metabolic_completeness
req_genus<-"UBA953"
req_MAGs<-MAGs_table[MAGs_table$gtdb_genus==req_genus,]$Representative_Genome
req_MAGs_faa_ko<- paste0(req_MAGs, ".faa.ko")
req_MAGS_mod_compl<-metabolCompl[,c(c("name","pathway group"),req_MAGs_faa_ko)]
names(req_MAGS_mod_compl)[names(req_MAGS_mod_compl)=='pathway group']<-"pathway_group"

t_rmmc<-t(req_MAGS_mod_compl)
n80<-data.frame(numeric())
names(n80)<-"n80count"
names(n80)
for(mag in req_MAGs_faa_ko)
{
  n80[mag,"n80count"]<-sum(as.numeric(t_rmmc[mag,])>80)
}
n80$size<-MAGs_table[MAGs_table$gtdb_genus==req_genus,]$reassembly_size_mb
n80$tp_lev<-MAGs_table[MAGs_table$gtdb_genus==req_genus,]$tp_lev

library(reshape2)
melted_rmms<-melt(req_MAGS_mod_compl)
ggplot(melted_rmms,aes(y=name,x=variable,fill=value))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1,
    size = 12, hjust = 1))

library(dendextend)
row_dend<-hclust(dist(as.matrix(req_MAGS_mod_compl[,-c(1,2)])))
col_dend<-as.dendrogram(hclust(dist(t(as.matrix(req_MAGS_mod_compl[,-c(1,2)])))))
n80$size_lev<-ifelse(n80$size<2,"small","big")
n80$names<-row.names(n80)
color_codes<-c(small="green",big="red")
labels_colors(col_dend)<-color_codes[n80$size_lev][order.dendrogram(col_dend)]
plot(col_dend)
heatmap(as.matrix(req_MAGS_mod_compl[,-c(1,2)]),cluster_columns=color_branches(col_dend,k=2))

library(gplots)

heatmap.2(as.matrix(req_MAGS_mod_compl[,-c(1,2)]),
          dendrogram = "col",
          Colv=col_dend
)


ggplot(n80,aes(x=n80count,y=size,color=tp_lev))+
  geom_point(size=4)

core_functions<-c("Arginine and proline metabolism","Aromatic amino acid metabolism","ATP synthesis","Branched-chain amino acid metabolism","Carbon fixation","Central carbohydrate metabolism","Cysteine and methionine metabolism","Fatty acid metabolism","Glycan biosynthesis","Glycosaminoglycan metabolism","Histidine metabolism","Lipid metabolism","Lipopolysaccharide metabolism","Lysine metabolism","Macrolide biosynthesis","Metabolic capacity","Methane metabolism","Nitrogen metabolism","Other amino acid metabolism","Photosynthesis","Polyamine biosynthesis","Polyketide sugar unit biosynthesis","Purine metabolism","Pyrimidine metabolism","Serine and threonine metabolism","Sterol biosynthesis","Symbiosis")
secondary_funcitions<-c("Aromatics degradation","Beta-Lactam biosynthesis","Biosynthesis of other secondary metabolites","Cofactor and vitamin metabolism","Drug resistance","Enediyne biosynthesis","Other carbohydrate metabolism","Other terpenoid biosynthesis","Pathogenicity","Plant pathogenicity","Sulfur metabolism","Terpenoid backbone biosynthesis","Type II polyketide biosynthesis")

for (thisFunc in secondary_funcitions){
  req_genus<-"Polynucleobacter"
  req_funcs<-thisFunc
  print(thisFunc)
  outfilepng<-paste("/mnt/Disk2/Data/MAGs_Paper/src/images_charts/heatmaps/Polynucleobacter/Rp_sec",req_genus,req_funcs,"ht.png",sep="_")
  req_MAGS_table<-MAGs_table[MAGs_table$gtdb_genus==req_genus,]
  req_MAGS_table<-req_MAGS_table[order(req_MAGS_table$size_lev),]
  req_MAGs<-req_MAGS_table$Representative_Genome
  req_MAGs_faa_ko<- paste0(req_MAGs, ".faa.ko")
  req_MAGS_mod_compl<-metabolCompl[,c(c("name","pathway group"),req_MAGs_faa_ko)]
  names(req_MAGS_mod_compl)[names(req_MAGS_mod_compl)=='pathway group']<-"pathway_group"
  req_MAGS_func<-req_MAGS_mod_compl[req_MAGS_mod_compl["pathway_group"]==req_funcs,]
  if(sum(req_MAGS_func[,-c(1,2)])==0){
    print ("skip")
    next
  }
  ht<-Heatmap(t(as.matrix(req_MAGS_func[,-c(1,2)])), split=req_MAGS_table$size_lev)
  png(outfilepng)
  draw(ht, padding= unit(c(2,20,2,2),"mm"))
  dev.off()
}

###HEATMAPS_Diviidng big small and medium genomes for specific function type
library(terra)
req_genus<-"UBA953"
req_funcs<-"Arginine and proline metabolism"
outfilepng<-paste("/mnt/Disk2/Data/MAGs_Paper/src/images_charts/heatmaps/Rp",req_genus,req_funcs,"ht.png",sep="_")
req_MAGS_table<-MAGs_table[MAGs_table$gtdb_genus==req_genus,]
req_MAGS_table<-req_MAGS_table[order(req_MAGS_table$size_lev),]
req_MAGs<-req_MAGS_table$Representative_Genome
req_MAGs_faa_ko<- paste0(req_MAGs, ".faa.ko")
req_MAGS_mod_compl<-metabolCompl[,c(c("name","pathway group"),req_MAGs_faa_ko)]
names(req_MAGS_mod_compl)[names(req_MAGS_mod_compl)=='pathway group']<-"pathway_group"
req_MAGS_func<-req_MAGS_mod_compl[req_MAGS_mod_compl["pathway_group"]==req_funcs,]
ht<-Heatmap(t(as.matrix(req_MAGS_func[,-c(1,2)])), split=req_MAGS_table$size_lev)
png(outfilepng)
draw(ht)
dev.off()
col_dend2<-as.dendrogram(hclust(dist(t(as.matrix(req_MAGS_func[,-c(1,2)])))))
labels_colors(col_dend2)<-color_codes[n80$size_lev][order.dendrogram(col_dend2)]
plot(col_dend2)


#NMDS metabolism
meatbolMAT<-t(metabolic_completeness_new[rowSums(metabolic_completeness_new)>0,]) # for metabolis_completeness
meatbolMAT<-t(metabolic_redundandancy[rowSums(metabolic_redundandancy)>0,]) #for Redundance
#meatbolMAT<-meatbolMAT[rowSums(meatbolMAT[])>0,colnames(meatbolMAT)!="EULd_S271TO_unclassified_48_35.fasta.faa.ko"]
#tMeatbolMAT<-t(meatbolMAT)
MeatbolMAT<-meatbolMAT[order(row.names(meatbolMAT),decreasing = F),]
set.seed(2718)
abundance_nmds<-metaMDS(MeatbolMAT)
plot(abundance_nmds)
data.scores <-as.data.frame(scores(abundance_nmds)$sites)
data.scores$weightedTP<-MAGs_table$weighted_TP
data.scores$tp_lev<-MAGs_table$tp_lev
data.scores$size_lev<-MAGs_table$size_lev
data.scores$phyla<-MAGs_table$gtdb_phyla
data.scores$tn_lev<-MAGs_table$tns_lev

ggplot(data.scores, aes(x = NMDS1, y = NMDS2, color=phyla,))+
  geom_point(size=4)

meatbolMAT<-as.data.frame(meatbolMAT)
meatbolMAT$weighted_tp<-MAGs_table$weighted_TP
meatbolMAT$weighted_tn<-MAGs_table$weighted_TN
meatbolMAT$weighted_doc<-MAGs_table$weighted_DOC
#PCA
meatbolMAT<-t(metabolic_completeness_new[rowSums(metabolic_completeness_new)>0,]) # for metabolis_completeness
meatbolMAT<-t(metabolic_redundandancy[rowSums(metabolic_redundandancy)>0,]) #for Redundance
MeatbolMAT.pca<-prcomp(meatbolMAT, center = T, scale. = T)
#install_github("vqv/ggbiplot")
library(ggbiplot)
summary(MeatbolMAT.pca)

row.names(MAGs_table)<-MAGs_table$cluster_ID
MAGs_table<-MAGs_table[order(MAGs_table$cluster_ID),]
ggbiplot(MeatbolMAT.pca,  groups = log(MAGs_table$weighted_DOC),var.axes=FALSE)+theme_bw()+scale_colour_gradientn(colours = topo.colors(10))
ggbiplot(MeatbolMAT.pca,  groups = log(MAGs_table$weighted_TP),var.axes=FALSE)+theme_bw()+scale_colour_gradientn(colours = topo.colors(10))
ggbiplot(MeatbolMAT.pca,  groups = log(MAGs_table$weighted_TN),var.axes=FALSE)+theme_bw()+scale_colour_gradientn(colours = topo.colors(10))
ggbiplot(MeatbolMAT.pca,  groups = log(MAGs_table$reassembly_size_mb),var.axes=FALSE)+theme_bw()+scale_colour_gradientn(colours = topo.colors(10))
PCAMTAX<-ggbiplot(MeatbolMAT.pca,  groups = MAGs_table$gtdb_phyla, var.axes=FALSE)
out_file_svg_PCAMTAX<-"/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/PCAMTAX.svg"
ggsave(file=out_file_svg_PCAMTAX,device="svg",plot=PCAMTAX,width=900, height=715, units="px",dpi=150)

#heatmap of metabolic completness vs TP.y/ TNs/DOC/GS
meatbolMAT<-t(metabolic_completeness_new[rowSums(metabolic_completeness_new)>0,]) # for metabolis_completeness
meatbolDF<-as.data.frame(meatbolMAT)
meatbolDF$weighted_tp<-MAGs_table$weighted_TP
meatbolDF$weighted_tn<-MAGs_table$weighted_TN
meatbolDF$weighted_doc<-MAGs_table$weighted_DOC
meatbolDF$log.wt_tp<-log(meatbolDF$weighted_tp)
meatbolDF$log.wt_tn<-log(meatbolDF$weighted_tn)
meatbolDF$log.wt_doc<-log(meatbolDF$weighted_doc)
x<-data.frame(matrix(NA, nrow=length(colnames(meatbolMAT)), ncol=3),row.names=colnames(meatbolMAT))
colnames(x)<- c("TP", "TN", "DOC")
for(i in colnames(meatbolDF)){
  if(startsWith(i, "M")){
    col_test<-meatbolDF[,i]
    x[i,"TP"]<-cor.test(col_test,meatbolDF$weighted_tp,method = "pearson")["p.value"]
    x[i,"TN"]<-cor.test(col_test,meatbolDF$weighted_tn,method = "pearson")["p.value"]
    x[i,"DOC"]<-cor.test(col_test,meatbolDF$weighted_doc,method = "pearson")["p.value"]

    #if (cor_test['p.value']<0.1){
    #  print(paste(i,cor_test['estimate'],cor_test['p.value']))
    #}
  }
}
meatbolCorr<-merge(x, metabolCompl[c("pathway group","name")],by.x = "row.names",by.y = "row.names")
write.csv(meatbolCorr,metabolic_corelations_file)
#TP
interesting_modules_list<-c("M00118","M00623","M00373","M00862","M00861","M00104","M00726","M00023","M00047","M00034","M00540","M00344","M00014","M00346","M00832","M00051","M00376","M00013","M00555","M00785","M00120","M00525","M00714","M00548","M00744","M00783","M00799","M00127","M00798","M00117","M00129","M00155","M00880","M00027","M00099","M00016","M00345","M00848","M00539")
#TN
interesting_modules_list<-c("M00800","M00039","M00358","M00736","M00801","M00079","M00761","M00075","M00857","M00799","M00364","M00798","M00010","M00065","M00028","M00652","M00014","M00842","M00843","M00860","M00126","M00785","M00035","M00170","M00093","M00718","M00022","M00741","M00172","M00345","M00055","M00047","M00675","M00866","M00044","M00725","M00730","M00129","M00574","M00784","M00171") #TN
#DOC
interesting_modules_list<-c("M00028","M00358","M00023","M00079","M00857","M00844","M00029","M00019","M00063","M00570","M00845","M00372","M00044","M00878","M00016","M00030","M00126","M00015","M00872","M00104","M00051","M00344","M00064","M00345","M00027","M00862","M00048","M00038","M00726","M00012","M00367","M00849","M00840","M00861","M00866","M00039","M00136","M00375","M00801","M00838","M00714","M00623","M00095","M00009","M00155","M00116","M00846","M00014","M00627","M00365","M00033","M00793","M00631","M00010","M00060","M00120","M00376","M00855","M00101","M00026","M00171","M00540","M00031","M00102","M00127","M00784","M00093","M00531","M00034","M00525","M00803","M00826","M00736","M00596","M00081","M00036","M00018","M00848","M00134","M00053","M00089","M00569","M00876","M00877","M00061","M00579","M00555","M00779")
intersting_module<-"M00862"
interesting_variable<-"log.wt_doc"
for(intersting_module in interesting_modules_list){
  formula_oi<-paste(intersting_module,"~",interesting_variable,sep = "")
  reduced_meatbol_df<-meatbolDF[,c(interesting_variable,intersting_module)]
  outfile<-paste0("/mnt/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/MAGs_Paper_V4/graphics/modules_of_interest/",interesting_variable,"_",intersting_module,".png")
  my_plot<-ggplotRegression(lm(formula_oi,data = reduced_meatbol_df))
  png(outfile)
  print(my_plot)
  dev.off()
  print(outfile)
}
formula_oi<-paste(intersting_module,"~log.",intresting_variable,sep = "")
reduced_meatbol_df<-meatbolDF[,c(interesting_variable,intersting_module)]
reduced_meatbol_df[paste("log.",interesting_variable,sep = "")]<-log(reduced_meatbol_df[intresting_variable])
ggplotRegression(lm(formula_oi,data = reduced_meatbol_df))
png(paste("/mnt/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/MAGs_Paper_V4/graphics/modules_of_interest/",intersting_module,"_",interesting_variable,".png"))
ggplotRegression(lm(formula_oi,data = reduced_meatbol_df))
dev.off()


#### NEW STATISTICS ###
#Phylum

cols_phyla <- c("Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota","Verrucomicrobiota")
reduced_members<-members_table[members_table$gtdb_phylum %in% cols_phyla,]
genus_counts<-table(reduced_members$gtdb_genus)
top_genera<-names(genus_counts[genus_counts > 5])
top_genera<-top_genera[!top_genera==""]
reduced_members$gtdb_phylum<-as.factor(reduced_members$gtdb_phylum)
reduced_members$gtdb_genus2<-ifelse(reduced_members$gtdb_genus %in% top_genera, reduced_members$gtdb_genus, "Other")
xs <- split(reduced_members,f = reduced_members$gtdb_phylum)
p1_gs <- ggplot(xs$Actinomycetota,aes(x = An_Mean_Temp,y = size_rep ,group = 1,colour = gtdb_genus2, size = abundance_val, weight=abundance_val )) +
  geom_point() +                                          #replace:TNs/TP.y/DOC.y  #replace:GC_content/size_rep
  stat_smooth(method="lm", formula=y~x) +                 #An_Mean_Temp;An_Max_Temp;An_Min_Temp;Temperature1
  facet_wrap(~gtdb_phylum)+scale_size(guide = 'none')+
  #stat_poly_eq(aes(weight=abundance_val,
  #                 label =  paste(..adj.rr.label..)),
  #             label.x.npc = "right", label.y.npc = 0.15,
  #             formula = y~x, parse = TRUE, size=7)+
  stat_fit_tidy(method = "lm",
                label.x = "right",
                size=7,
                method.args = list(formula = y ~ x, weights = quote(weight)),
                mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                              after_stat(x_estimate),
                                              after_stat(x_p.value))))+
  theme(legend.position="none")
p2_gs <- p1_gs %+% xs$Pseudomonadota
p3_gs <- p1_gs %+% xs$Cyanobacteriota
p4_gs <- p1_gs %+% xs$Bacteroidota
p5_gs <- p1_gs %+% xs$Verrucomicrobiota

gs_grid<-grid.arrange(p1_gs,p2_gs,p3_gs,p4_gs,p5_gs, ncol=1)
p1_gc <- ggplot(xs$Actinomycetota,aes(x = An_Mean_Temp,y = GC_content ,group = 1,colour = gtdb_genus2, size = abundance_val, weight=abundance_val )) +
  geom_point() +                                          #replace:TNs/TP.y/DOC.y  #replace:GC_content/size_rep
  stat_smooth(method="lm", formula=y~x) +                 #An_Mean_Temp;An_Max_Temp;An_Min_Temp;Temperature1
  facet_wrap(~gtdb_phylum)+scale_size(guide = 'none')+
  #stat_poly_eq(aes(weight=abundance_val,
  #                 label =  paste(..adj.rr.label..)),
  #             label.x.npc = "right", label.y.npc = 0.15,
  #             formula = y~x, parse = TRUE, size=7)+
  stat_fit_tidy(method = "lm",
                label.x = "right",
                size=7,
                method.args = list(formula = y ~ x, weights = quote(weight)),
                mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                              after_stat(x_estimate),
                                              after_stat(x_p.value))))
p2_gc <- p1_gc %+% xs$Pseudomonadota
p3_gc <- p1_gc %+% xs$Cyanobacteriota
p4_gc <- p1_gc %+% xs$Bacteroidota
p5_gc <- p1_gc %+% xs$Verrucomicrobiota

gc_grid<-grid.arrange(p1_gc,p2_gc,p3_gc,p4_gc,p5_gc, ncol=1)
gs_gs_grid<-grid.arrange(gs_grid,gc_grid,nrow=1)


ggplot(members_table,aes(x = log(TNs),y = GC_content ,group = 1, size = abundance_val, weight=abundance_val )) +
  geom_point() +                                          #replace:TNs/TP.y/DOC.y  #replace:GC_content/size_rep
  stat_smooth(method="lm", formula=y~x) +                 #An_Mean_Temp;An_Max_Temp;An_Min_Temp;Temperature1
  scale_size(guide = 'none')+
  stat_poly_eq(aes(weight=abundance_val,
                   label =  paste(..adj.rr.label..)),
               label.x.npc = "right", label.y.npc = 0.15,
               formula = y~x, parse = TRUE, size=7)+
  stat_fit_tidy(method = "lm",
                label.x = "right",
                size=7,
                method.args = list(formula = y ~ x, weights = quote(weight)),
                mapping = aes(label = sprintf("Slope = %.3g\np-value = %.3g",
                                              after_stat(x_estimate),
                                              after_stat(x_p.value))))

#summarize 5-10 major findings
#use cutoffs- 1/3rd of data @done
#number of complete modules
#redundancy at phyla and genus @done
#add colour


#genus
cols_genus<- c("Limnohabitans","Rhodoferax","Rhodoluna","Planktophila","UBA2093","UBA3006","Sediminibacterium","Flavobacterium","SXYR01","Cyanobium","UBA953")
reduced_members<-members_table[members_table$gtdb_genus %in% cols_genus,]
reduced_members$gtdb_genus<-as.factor(reduced_members$gtdb_genus)
summary<-reduced_members %>%
  group_by(gtdb_genus)%>%
  nest()%>%
  mutate(model = map(data,~lm(GC_content~TNs,weights=abundance_val, data = .x)),  #replace:TNs/TP.y/DOC.y  #replace:GC_Content/size_rep
         adj.r.squared = map_dbl(model, ~ summary(.x)$adj.r.squared),                          #An_Mean_Temp;An_Max_Temp;An_Min_Temp;Temperature1
         intercept = map_dbl(model, ~ .x$coef[[1]]),
         slope = map_dbl(model, ~ .x$coef[[2]]),
         pvalue = map_dbl(model, ~ summary(.x)$coef[2,4])
         )
# %>%
#   select(-data,-model) %>%
#   left_join(reduced_members)
#   %>%
# ggplot(reduced_members,aes(y=GC_content, x=An_Mean_Temp,size = abundance_val, weight=abundance_val)) + #replace:TNs/TP.y/DOC.y  #replace:GC_Content/size_rep
#   geom_point()+                                                                                     #An_Mean_Temp;An_Max_Temp;An_Min_Temp;Temperature1
#   stat_smooth(method="lm",formula=y~x)+
#   #geom_vline(xintercept=log(3800), linetype='dotted', col = 'red')+ #replace: tp:30,TNs:1000, DOC.y:3800
#   facet_wrap(~gtdb_genus)+
#   stat_poly_eq(aes(label = paste(..rr.label..)),
#                      label.x.npc = "right", label.y.npc = 0.15,
#                      formula = y~x, parse = TRUE, size = 3)+
#   stat_fit_glance(method = 'lm',
#                        method.args = list(formula = y~x),
#                        geom = 'text',
#                        aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
#                        label.x.npc = 'right', label.y.npc = 0.35,
#                        size = 3)
#
#
# geom_text(inherit.aes = FALSE ,aes(x=7, y=4.8,label=paste("Adj R2 = ", adj.r.squared, "\n",  ##adjust x and y to fit graph GC_content tp:0.5,60; tn:8,50 ; doc:6.5,60
#                                  "Slope =", slope, "\n",                                      ##adjust x and y to fit graph genome size tp:0.5,4.8; tn:5,4.8 ; doc:7,4.8
#                                  "P =", pvalue)),size=5)
#

##weighted Box plots

library(epade)
library(weights)
library(ggpubr)
#TP_GC
wtd.t.test(subset(members_table,tp_lev=="Low")$size_rep, subset(members_table,tp_lev=="High")$size_rep,
           weight=subset(members_table,tp_lev=="Low")$abundance_val,
           weighty = subset(members_table,tp_lev=="High")$abundance_val,
           samedata=T)
ggplot(members_table,aes(x=tp_lev,y=size_rep,weight=abundance_val))+
  geom_boxplot()+
  geom_segment(aes(y=9,x=1,xend=2,yend=9),)+annotate("text",label="Weighted t-test p-value: 0.02",x=1.5,y=9.5)

#TN_GC
wtd.t.test(subset(members_table,tn_lev=="Low")$size_rep, subset(members_table,tn_lev=="High")$size_rep,
           weight=subset(members_table,tn_lev=="Low")$abundance_val,
           weighty = subset(members_table,tn_lev=="High")$abundance_val,
           samedata=T)
ggplot(members_table,aes(x=tn_lev,y=size_rep,weight=abundance_val))+
  geom_boxplot()+
  geom_segment(aes(y=9,x=1,xend=2,yend=9),)+annotate("text",label="Weighted t-test p-value: 2.12e-07",x=1.5,y=9.5)

#DOC_GC
wtd.t.test(subset(members_table,doc_lev=="Low")$size_rep, subset(members_table,doc_lev=="High")$size_rep,
           weight=subset(members_table,doc_lev=="Low")$abundance_val,
           weighty = subset(members_table,doc_lev=="High")$abundance_val,
           samedata=T)
ggplot(members_table,aes(x=doc_lev,y=size_rep,weight=abundance_val))+
  geom_boxplot()+
  geom_segment(aes(y=9,x=1,xend=2,yend=9),)+annotate("text",label="Weighted t-test p-value: 0.138",x=1.5,y=9.5)

#TEMP_GC
wtd.t.test(subset(members_table,temp_lev=="Low")$size_rep, subset(members_table,temp_lev=="High")$size_rep,
           weight=subset(members_table,temp_lev=="Low")$abundance_val,
           weighty = subset(members_table,temp_lev=="High")$abundance_val,
           samedata=T)
ggplot(members_table,aes(x=temp_lev,y=size_rep,weight=abundance_val))+
  geom_boxplot()+
  geom_segment(aes(y=9,x=1,xend=2,yend=9),)+annotate("text",label="Weighted t-test p-value:0.02",x=1.5,y=9.5)



#Add Redundancy to MAGS table
metabolic_redundandancy_test <-metabolic_redundandancy %>%
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))
rownames(metabolic_redundandancy_test)[rownames(metabolic_redundandancy_test) == "...395"] <- "Total_Redundancy"
meatbol_red_2<-as.data.frame(t(metabolic_redundandancy_test))
metabolic_redundandancy_totals<-t(metabolic_redundandancy_test["Total_Redundancy",])
MAGs_table_new<-MAGs_table
MAGs_table_new["Total_Redundancy"]<-meatbol_red_2$Total_Redundancy[match(row.names(MAGs_table), row.names(meatbol_red_2))]
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

ggplot(data=MAGs_table_new, aes(y=Total_Redundancy, x=log(site_TP))) + #replace:Total_Redundancy,c80,c100
  geom_point()+
  geom_smooth(method="lm",formula=y~x)+
  stat_poly_eq(aes(label = paste(..rr.label..)),
                     label.x.npc = "right", label.y.npc = 0.15,
                     formula = y~x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                       method.args = list(formula = y~x),
                       geom = 'text',
                       aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                       label.x.npc = 'right', label.y.npc = 0.35,
                       size = 3)

cols_phyla <- c("Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota","Verrucomicrobiota")
reduced_members<-MAGs_table_new[MAGs_table_new$gtdb_phyla %in% cols_phyla,]
my_plot_C<-ggplot(data=reduced_members, aes(y=c80, x=log(weighted_TN))) + #replace:Total_Redundancy,c80,c100
  geom_point()+labs(y="Redundant genes", x="Log of DOC")+
  geom_smooth(method="lm",formula=y~x)+
  facet_wrap(~gtdb_phyla,ncol=1)+
  stat_poly_eq(aes(label = paste(..rr.label..)),
                     label.x.npc = "right", label.y.npc = 0.15,
                     formula = y~x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                       method.args = list(formula = y~x),
                       geom = 'text',label.x = 7.5,label.y = ,
                       aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                       #label.x.npc = 'right', label.y.npc = 0.35,
                       size = 3)
ggarrange(my_plot_n,my_plot_p,my_plot_C, nrow = 1)
out_file<-"/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/scatterc80.svg"
ggsave(file=out_file,device="svg",plot=my_plot,width=500, height=1500, units="px",dpi=170)
cols_genus<- c("Limnohabitans","Rhodoferax","Rhodoluna","Planktophila","UBA2093","UBA3006","Sediminibacterium","Flavobacterium","SXYR01","Cyanobium","UBA953")

reduced_members<-MAGs_table_new[MAGs_table_new$gtdb_genus %in% cols_genus,]
ggplot(data=reduced_members, aes(y=c80, x=reassembly_size_mb)) + #replace:Total_Redundancy,c80,c100  #replace:GC_Content/size_rep
  geom_point()+
  stat_smooth(method="lm",formula=y~x)+
  facet_wrap(~gtdb_genus)+
  stat_poly_eq(aes(label = paste(..rr.label..)),
                     label.x.npc = "right", label.y.npc = 0.15,
                     formula = y~x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                       method.args = list(formula = y~x),
                       geom = 'text',
                       aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                       label.x = 2.5,label.y = 30,
                       size = 3)


ggplotRegression <- function(fit){
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",summary(fit)$adj.r.squared,
                       "Intercept =",fit$coef[[1]],
                       " Slope =",fit$coef[[2]],
                       " P =",summary(fit)$coef[2,4]))
}
