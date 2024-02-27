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

#c80
whole_data_plot_c80<-ggplot(data=MAGs_table_new, aes(y=c80, x=reassembly_size_mb)) + #replace:Total_Redundancy,c80,c100
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
c80_plot_phyla<-ggplot(data=reduced_members, aes(y=c80, x=reassembly_size_mb)) + #replace:Total_Redundancy,c80,c100
  geom_point()+labs(y="Count 80% Complete Pathways", x="Genome Size (Mb)")+
  geom_smooth(method="lm",formula=y~x)+
  facet_wrap(~gtdb_phyla,ncol=1)+
  stat_poly_eq(aes(label = paste(..rr.label..)),
                     label.x.npc = "right", label.y.npc = 0.15,
                     formula = y~x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                       method.args = list(formula = y~x),
                       geom = 'text',label.x = 2.5,label.y = 40,
                       aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                       #label.x.npc = 'right', label.y.npc = 0.35,
                       size = 3)

cols_genus<- c("Limnohabitans","Rhodoferax","Rhodoluna","Planktophila","UBA2093","UBA3006","Sediminibacterium","Flavobacterium","SXYR01","Cyanobium","UBA953")
reduced_members<-MAGs_table_new[MAGs_table_new$gtdb_genus %in% cols_genus,]
c80_plot_genus<-ggplot(data=reduced_members, aes(y=c80, x=reassembly_size_mb)) + #replace:Total_Redundancy,c80,c100  #replace:GC_Content/size_rep
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

#Redundancy
whole_data_plot_redundunacy<-ggplot(data=MAGs_table_new, aes(y=Total_Redundancy, x=reassembly_size_mb)) + #replace:Total_Redundancy,c80,c100
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
redundunacy_plot_phyla<-ggplot(data=reduced_members, aes(y=Total_Redundancy, x=reassembly_size_mb)) + #replace:Total_Redundancy,c80,c100
  geom_point()+labs(y="Count of\nRedundunct Genes", x="Genome Size (Mb)")+
  geom_smooth(method="lm",formula=y~x)+
  facet_wrap(~gtdb_phyla,nrow=1)+
  stat_poly_eq(aes(label = paste(..rr.label..)),
                     label.x.npc = "right", label.y.npc = 0.15,
                     formula = y~x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                       method.args = list(formula = y~x),
                       geom = 'text',label.x = 4,label.y = 40,
                       aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                       #label.x.npc = 'right', label.y.npc = 0.35,
                       size = 3)

cols_genus<- c("Limnohabitans","Rhodoferax","Rhodoluna","Planktophila","UBA2093","UBA3006","Sediminibacterium","Flavobacterium","SXYR01","Cyanobium","UBA953")
reduced_members<-MAGs_table_new[MAGs_table_new$gtdb_genus %in% cols_genus,]
redundunacy_plot_genus<-ggplot(data=reduced_members, aes(y=Total_Redundancy, x=reassembly_size_mb)) + #replace:Total_Redundancy,c80,c100  #replace:GC_Content/size_rep
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

combined_plot_scattter_kegg<-ggarrange(redundunacy_plot_phyla,c80_plot_phyla)
out_file<-"/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/scatterKEGG_redundancy.svg"
ggsave(file=out_file,device="svg",plot=redundunacy_plot_phyla,width=1615, height=330, units="px",dpi=170)