##factors with LOG, (non log ones below)
req_env2<-c("TNs","TP.y","DOC.y" #replace "x = .data[[env]]"  wtih "x = log(.data[[env]])" ###to preserve this comment please add "(" before the X while replaceing
            #"An_Max_Temp","temperature"
)
req_env3<-c("Natural log of total nitrogen (µg/l)", "Natural log of total phosphorus (µg/l)", "Natural log of dissolved organic carbon (µg/l)"
            #"Annual maximum temperature (C)", "Recorded temperature (C)"
)
req_env4<-setNames(as.list(req_env3),req_env2)
cols_phyla <- c("Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota","Verrucomicrobiota")
reduced_members<-members_table[members_table$gtdb_phylum %in% cols_phyla,]
genus_counts<-table(reduced_members$gtdb_genus)
top_genera<-names(genus_counts[genus_counts > 4])
top_genera<-top_genera[!top_genera==""]
reduced_members$gtdb_phylum<-as.factor(reduced_members$gtdb_phylum)
reduced_members$gtdb_genus2<-ifelse(reduced_members$gtdb_genus %in% top_genera, reduced_members$gtdb_genus, "Other")
reduced_members$gtdb_genus2<-factor(reduced_members$gtdb_genus2,levels=c(top_genera,"Other"))
myColors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(length(top_genera))

cols_genus<- c("Limnohabitans","Rhodoferax","Rhodoluna","Planktophila","UBA2093","UBA3006","Sediminibacterium","Flavobacterium","SXYR01","Cyanobium","UBA953")
reduced_members_genus<-members_table[members_table$gtdb_genus %in% cols_genus,]
reduced_members_genus$gtdb_genus<-as.factor(reduced_members_genus$gtdb_genus)
for (env in req_env2){
  p1_gs <- ggplot(reduced_members,aes(x = log(.data[[env]]),y = size_rep ,group = 1,colour = gtdb_genus2, size = abundance_val, weight=abundance_val )) +
    geom_point() +
    stat_smooth(method="lm", formula=y~x) +
    facet_wrap(~gtdb_phylum,ncol = 1)+scale_size(guide = 'none')+stat_poly_eq()+
    stat_fit_tidy(method = "lm",
                  label.x = "right",
                  method.args = list(formula = y ~ x, weights = quote(weight)),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))+
    scale_color_manual(values=c(myColors,"grey"))+
    labs(y="Genome size\n(Mb)",,x=req_env4[[env]]) +theme_linedraw()+
    theme(legend.position="none",strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(colour = "black"))

  p1_gc <- ggplot(reduced_members,aes(x = log(.data[[env]]),y = GC_content ,group = 1,colour = gtdb_genus2, size = abundance_val, weight=abundance_val )) +
    geom_point() +
    stat_smooth(method="lm", formula=y~x) +
    facet_wrap(~gtdb_phylum, ncol=1)+scale_size(guide = 'none')+stat_poly_eq()+
    stat_fit_tidy(method = "lm",
                  label.x = "right",label.y = "bottom",
                  method.args = list(formula = y ~ x, weights = quote(weight)),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))+
    scale_color_manual(values=c(myColors,"grey"))+
    labs(y="GC content\n(%)",x=req_env4[[env]]) +theme_linedraw()+
    theme(strip.background = element_rect(fill = "white", color = "black"), strip.text = element_text(colour = "black"))+
    guides(colour=guide_legend(ncol=1,title="Genera", label.position = "right", label.hjust = 0))

  gs_gc_grid<-grid.arrange(p1_gs,p1_gc,ncol=2, widths = c(1.5,2))
  out_files<-paste0("/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/env_gs_gc_scatters/phyla/supplementary_scatter_gs_gc_vs_",env,".svg")
  out_filep<-paste0("/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/env_gs_gc_scatters/phyla/png/supplementary_scatter_gs_gc_vs_",env,".png")
  ggsave(file=out_files,device="svg",plot=gs_gc_grid,width=1892, height=1600, units="px",dpi=170)
  ggsave(file=out_filep,device="png",plot=gs_gc_grid,width=1892, height=1600, units="px",dpi=170)

  p2_gs <- ggplot(reduced_members_genus,aes(x = log(.data[[env]]),y = size_rep ,group = 1, size = abundance_val, weight=abundance_val )) +
    geom_point() +                                          #replace:TNs/TP.y/DOC.y  #replace:GC_content/size_rep
    stat_smooth(method="lm", formula=y~x) +                 #An_Mean_Temp;An_Max_Temp;An_Min_Temp;Temperature1
    facet_wrap(~gtdb_genus,ncol = 1)+scale_size(guide = 'none')+stat_poly_eq()+
    stat_fit_tidy(method = "lm",
                  label.x = "right",
                  method.args = list(formula = y ~ x, weights = quote(weight)),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))+
    labs(y="Genome size\n(Mb)",,x=req_env4[[env]]) +theme_linedraw()+
    theme(legend.position="none",strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(colour = "black"))

  p2_gc <- ggplot(reduced_members_genus,aes(x = log(.data[[env]]),y = GC_content ,group = 1, size = abundance_val, weight=abundance_val )) +
    geom_point() +                                          #replace:TNs/TP.y/DOC.y  #replace:GC_content/size_rep
    stat_smooth(method="lm", formula=y~x) +                 #An_Mean_Temp;An_Max_Temp;An_Min_Temp;Temperature1
    facet_wrap(~gtdb_genus, ncol=1)+scale_size(guide = 'none')+stat_poly_eq()+
    stat_fit_tidy(method = "lm",
                  label.x = "right",label.y = "bottom",
                  method.args = list(formula = y ~ x, weights = quote(weight)),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))+
    labs(y="GC content\n(%)",x=req_env4[[env]]) +theme_linedraw()+
    theme(strip.background = element_rect(fill = "white", color = "black"), strip.text = element_text(colour = "black"))

  gs_gs_grid<-grid.arrange(p2_gs,p2_gc,ncol=2)
  out_files<-paste0("/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/env_gs_gc_scatters/Genus/supplementary_scatter_gs_gc_vs_",env,"_Genus.svg")
  out_filep<-paste0("/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/env_gs_gc_scatters/Genus/png/supplementary_scatter_gs_gc_vs_",env,"_Genus.png")
  ggsave(file=out_files,device="svg",plot=gs_gs_grid,width=1892, height=1892, units="px",dpi=150)
  ggsave(file=out_filep,device="png",plot=gs_gs_grid,width=1892, height=1892, units="px",dpi=150)
}


#factors without log
req_env2<-c(#"TNs","TP.y","DOC.y" #replace "x = .data[[env]]"  wtih "x = log(.data[[env]])" ###to preserve this comment please add "(" before the X while replaceing
            "An_Max_Temp","temperature"
)
req_env3<-c(#"Log of total nitrogen (µg/l)", "Log of total phosphorus (µg/l)", "Log of dissolved organic carbon (µg/l)"
            "Annual maximum temperature (C)", "Recorded temperature (C)"
)
req_env4<-setNames(as.list(req_env3),req_env2)
cols_phyla <- c("Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota","Verrucomicrobiota")
reduced_members<-members_table[members_table$gtdb_phylum %in% cols_phyla,]
genus_counts<-table(reduced_members$gtdb_genus)
top_genera<-names(genus_counts[genus_counts > 4])
top_genera<-top_genera[!top_genera==""]
reduced_members$gtdb_phylum<-as.factor(reduced_members$gtdb_phylum)
reduced_members$gtdb_genus2<-ifelse(reduced_members$gtdb_genus %in% top_genera, reduced_members$gtdb_genus, "Other")
reduced_members$gtdb_genus2<-factor(reduced_members$gtdb_genus2,levels=c(top_genera,"Other"))
myColors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(length(top_genera))

cols_genus<- c("Limnohabitans","Rhodoferax","Rhodoluna","Planktophila","UBA2093","UBA3006","Sediminibacterium","Flavobacterium","SXYR01","Cyanobium","UBA953")
reduced_members_genus<-members_table[members_table$gtdb_genus %in% cols_genus,]
reduced_members_genus$gtdb_genus<-as.factor(reduced_members_genus$gtdb_genus)
for (env in req_env2){
  p1_gs <- ggplot(reduced_members,aes(x = .data[[env]],y = size_rep ,group = 1,colour = gtdb_genus2, size = abundance_val, weight=abundance_val )) +
    geom_point() +
    stat_smooth(method="lm", formula=y~x) +
    facet_wrap(~gtdb_phylum,ncol = 1)+scale_size(guide = 'none')+stat_poly_eq()+
    stat_fit_tidy(method = "lm",
                  label.x = "right",
                  method.args = list(formula = y ~ x, weights = quote(weight)),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))+scale_color_manual(values=c(myColors,"grey"))+
    labs(y="Genome size\n(Mb)",,x=req_env4[[env]]) +theme_linedraw()+
    theme(legend.position="none",strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(colour = "black"))

  p1_gc <- ggplot(reduced_members,aes(x = .data[[env]],y = GC_content ,group = 1,colour = gtdb_genus2, size = abundance_val, weight=abundance_val )) +
    geom_point() +
    stat_smooth(method="lm", formula=y~x) +
    facet_wrap(~gtdb_phylum, ncol=1)+scale_size(guide = 'none')+stat_poly_eq()+
    stat_fit_tidy(method = "lm",
                  label.x = "right",label.y = "bottom",
                  method.args = list(formula = y ~ x, weights = quote(weight)),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))+scale_color_manual(values=c(myColors,"grey"))+
    labs(y="GC content\n(%)",x=req_env4[[env]]) +theme_linedraw()+
    theme(strip.background = element_rect(fill = "white", color = "black"), strip.text = element_text(colour = "black"))+
    guides(colour=guide_legend(ncol=1,title="Genera", label.position = "right", label.hjust = 0))

  gs_gc_grid<-grid.arrange(p1_gs,p1_gc,ncol=2, widths = c(1.5,2))
  out_files<-paste0("/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/env_gs_gc_scatters/phyla/supplementary_scatter_gs_gc_vs_",env,".svg")
  out_filep<-paste0("/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/env_gs_gc_scatters/phyla/png/supplementary_scatter_gs_gc_vs_",env,".png")
  ggsave(file=out_files,device="svg",plot=gs_gc_grid,width=1892, height=1600, units="px",dpi=170)
  ggsave(file=out_filep,device="png",plot=gs_gc_grid,width=1892, height=1600, units="px",dpi=170)

  p2_gs <- ggplot(reduced_members_genus,aes(x = .data[[env]],y = size_rep ,group = 1, size = abundance_val, weight=abundance_val )) +
    geom_point() +                                          #replace:TNs/TP.y/DOC.y  #replace:GC_content/size_rep
    stat_smooth(method="lm", formula=y~x) +                 #An_Mean_Temp;An_Max_Temp;An_Min_Temp;Temperature1
    facet_wrap(~gtdb_genus,ncol = 1)+scale_size(guide = 'none')+stat_poly_eq()+
    stat_fit_tidy(method = "lm",
                  label.x = "right",
                  method.args = list(formula = y ~ x, weights = quote(weight)),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))+
    labs(y="Genome size\n(Mb)",,x=req_env4[[env]]) +theme_linedraw()+
    theme(legend.position="none",strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(colour = "black"))

  p2_gc <- ggplot(reduced_members_genus,aes(x = .data[[env]],y = GC_content ,group = 1, size = abundance_val, weight=abundance_val )) +
    geom_point() +                                          #replace:TNs/TP.y/DOC.y  #replace:GC_content/size_rep
    stat_smooth(method="lm", formula=y~x) +                 #An_Mean_Temp;An_Max_Temp;An_Min_Temp;Temperature1
    facet_wrap(~gtdb_genus, ncol=1)+scale_size(guide = 'none')+stat_poly_eq()+
    stat_fit_tidy(method = "lm",
                  label.x = "right",label.y = "bottom",
                  method.args = list(formula = y ~ x, weights = quote(weight)),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))+
    labs(y="GC content\n(%)",x=req_env4[[env]]) +theme_linedraw()+
    theme(strip.background = element_rect(fill = "white", color = "black"), strip.text = element_text(colour = "black"))

  gs_gs_grid<-grid.arrange(p2_gs,p2_gc,ncol=2)
  out_files<-paste0("/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/env_gs_gc_scatters/Genus/supplementary_scatter_gs_gc_vs_",env,"_Genus.svg")
  out_filep<-paste0("/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/env_gs_gc_scatters/Genus/png/supplementary_scatter_gs_gc_vs_",env,"_Genus.png")
  ggsave(file=out_files,device="svg",plot=gs_gs_grid,width=1892, height=1892, units="px",dpi=150)
  ggsave(file=out_filep,device="png",plot=gs_gs_grid,width=1892, height=1892, units="px",dpi=150)
}