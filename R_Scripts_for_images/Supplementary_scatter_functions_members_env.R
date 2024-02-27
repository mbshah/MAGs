####Factors needing Log
req_env2<-c("TNs","TP.y","DOC.y"
            #"An_Max_Temp","temperature"
)
req_env3<-c("Natural log of total nitrogen (µg/l)", "Natural log of total phosphorus (µg/l)", "Natural log of dissolved organic carbon (µg/l)"
            #"Annual maximum temperature (C)", "Recorded temperature (C)"
)
req_env4<-setNames(as.list(req_env3),req_env2)
cols_phyla <- c("Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota","Verrucomicrobiota")
reduced_members<-members_table[members_table$gtdb_phylum %in% cols_phyla,]
reduced_members$gtdb_phylum<-as.factor(reduced_members$gtdb_phylum)
for (env in req_env2){
  #p1<-ggplot(data=reduced_members, aes(y=GC_content, x = .data[[env]])) + #replace:Total_Redundancy,c80,c100
  #  geom_point()+labs(y="\nGC content\n(%)",) +theme_linedraw()+
  #theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
  #      strip.background = element_rect(fill = "white", color = "black"),
  #      strip.text = element_text(colour = "black"))+
  #  geom_smooth(method="lm",formula=y~x)+
  #  facet_wrap(~gtdb_phylum,nrow=1)+
  #  stat_fit_tidy(method = "lm",
  #                label.x = "right",label.y="bottom",
  #                size=3,
  #                method.args = list(formula = y ~ x),
  #                mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
  #                                              after_stat(x_estimate),
  #                                              after_stat(x_p.value))))
#
  p2<-ggplot(data=reduced_members, aes(y=redundant_genes, x = log(.data[[env]]))) + #replace:Total_Redundancy,c80,c100
    geom_point()+labs(y="\nCount of\nredundant genes", ) +theme_linedraw()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),strip.text = element_text(colour = "black"))+
    geom_smooth(method="lm",formula=y~x)+
    facet_wrap(~gtdb_phylum,nrow=1)+
    stat_fit_tidy(method = "lm",
                  label.x = "left",
                  size=3,
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))

  p3<-ggplot(data=reduced_members, aes(y=rep_codingpc, x = log(.data[[env]]))) + #replace:Total_Redundancy,c80,c100
    geom_point()+labs(y="\nCoding regions\n(%)", ) +theme_linedraw()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), strip.text.x=element_blank())+
    geom_smooth(method="lm",formula=y~x)+
    facet_wrap(~gtdb_phylum,nrow=1)+
    stat_fit_tidy(method = "lm",
                  label.x = "right", label.y="bottom",
                  size=3,
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))

  p4 <- ggplot(reduced_members,aes(x = log(.data[[env]]),y = KOs_for_ABC_Transporter ,group = 1)) +
    geom_point() +
    stat_smooth(method="lm", formula=y~x) +labs(y="Number of KOs\nfor ABC\ntransporters",)  +theme_linedraw()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), strip.text.x=element_blank())+
    facet_wrap(~gtdb_phylum,nrow = 1)+scale_size(guide = 'none')+
    #stat_poly_eq(aes(,
    #                 label =  paste(..adj.rr.label..)),
    #             label.x.npc = "right", label.y.npc = 0.15,
    #             formula = y~x, parse = TRUE, size=3)+
    stat_fit_tidy(method = "lm",
                  label.x = "right",
                  size=3,
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))

  p5 <- ggplot(reduced_members,aes(x = log(.data[[env]]),y = sigma_factor ,group = 1)) +
    geom_point() +
    stat_smooth(method="lm", formula=y~x) +labs(y="Number of sigma\nfactor\ngenes", )   +theme_linedraw()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), strip.text.x=element_blank())+
    facet_wrap(~gtdb_phylum,nrow = 1)+scale_size(guide = 'none')+
    stat_fit_tidy(method = "lm",
                  label.x = "right",
                  size=3,
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))

  p6 <- ggplot(reduced_members,aes(x = log(.data[[env]]),y = kc_Cell_motility ,group = 1)) +
    geom_point() +
    stat_smooth(method="lm", formula=y~x) +labs(y="Number of KOs\nfor\ncell motility", x=req_env4[[env]])  +
    facet_wrap(~gtdb_phylum,nrow = 1)+scale_size(guide = 'none')+theme_linedraw()+
  theme(strip.text.x=element_blank())+
    stat_fit_tidy(method = "lm",
                  label.x = "right",
                  size=3,
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))

  combined_plot_scattter_kegg_env<-ggarrange(p2,p3,p4,p5,p6,nrow = 5,labels = c("(A)","(B)","(C)","(D)","(E)"),heights=c(rep(0.9,4),1))
  combined_plot_scattter_kegg_env
  out_files<-paste0("/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/env_genome_properties_scatters/supplementary_scatter_genomic_prop_vs_",env,".svg")
  out_filep<-paste0("/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/env_genome_properties_scatters/png/supplementary_scatter_genomic_prop_vs_",env,".png")
  ggsave(file=out_files,device="svg",plot=combined_plot_scattter_kegg_env,width=1892, height=1580, units="px",dpi=170)
  ggsave(file=out_filep,device="png",plot=combined_plot_scattter_kegg_env,width=1892, height=1580, units="px",dpi=170)
}


######NON LOG factors
req_env2<-c(#"TNs","TP.y","DOC.y"
            "An_Max_Temp","temperature"
)
req_env3<-c(#"Natural log of total nitrogen (µg/l)", "Natural log of total phosphorus (µg/l)", "Natural log of dissolved organic carbon (µg/l)"
            "Annual maximum temperature (C)", "Recorded temperature (C)"
)
req_env4<-setNames(as.list(req_env3),req_env2)
cols_phyla <- c("Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota","Verrucomicrobiota")
reduced_members<-members_table[members_table$gtdb_phylum %in% cols_phyla,]
reduced_members$gtdb_phylum<-as.factor(reduced_members$gtdb_phylum)
for (env in req_env2){
  #p1<-ggplot(data=reduced_members, aes(y=GC_content, x = .data[[env]])) + #replace:Total_Redundancy,c80,c100
  #  geom_point()+labs(y="\nGC content\n(%)",) +theme_linedraw()+
  #theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
  #      strip.background = element_rect(fill = "white", color = "black"),strip.text = element_text(colour = "black")
  #      strip.text = element_text(colour = "black"))+
  #  geom_smooth(method="lm",formula=y~x)+
  #  facet_wrap(~gtdb_phylum,nrow=1)+
  #  stat_fit_tidy(method = "lm",
  #                label.x = "right",label.y="bottom",
  #                size=3,
  #                method.args = list(formula = y ~ x),
  #                mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
  #                                              after_stat(x_estimate),
  #                                              after_stat(x_p.value))))
#
  p2<-ggplot(data=reduced_members, aes(y=redundant_genes, x = .data[[env]])) + #replace:Total_Redundancy,c80,c100
    geom_point()+labs(y="\nCount of\nredundant genes", ) +theme_linedraw()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),strip.text = element_text(colour = "black"))+
    geom_smooth(method="lm",formula=y~x)+
    facet_wrap(~gtdb_phylum,nrow=1)+
    stat_fit_tidy(method = "lm",
                  label.x = "left",
                  size=3,
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))

  p3<-ggplot(data=reduced_members, aes(y=rep_codingpc, x = .data[[env]])) + #replace:Total_Redundancy,c80,c100
    geom_point()+labs(y="\nCoding regions\n(%)", ) +theme_linedraw()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), strip.text.x=element_blank())+
    geom_smooth(method="lm",formula=y~x)+
    facet_wrap(~gtdb_phylum,nrow=1)+
    stat_fit_tidy(method = "lm",
                  label.x = "right", label.y="bottom",
                  size=3,
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))

  p4 <- ggplot(reduced_members,aes(x = .data[[env]],y = KOs_for_ABC_Transporter ,group = 1)) +
    geom_point() +
    stat_smooth(method="lm", formula=y~x) +labs(y="Number of KOs\nfor ABC\ntransporters",)  +theme_linedraw()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), strip.text.x=element_blank())+
    facet_wrap(~gtdb_phylum,nrow = 1)+scale_size(guide = 'none')+
    #stat_poly_eq(aes(,
    #                 label =  paste(..adj.rr.label..)),
    #             label.x.npc = "right", label.y.npc = 0.15,
    #             formula = y~x, parse = TRUE, size=3)+
    stat_fit_tidy(method = "lm",
                  label.x = "right",
                  size=3,
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))

  p5 <- ggplot(reduced_members,aes(x = .data[[env]],y = sigma_factor ,group = 1)) +
    geom_point() +
    stat_smooth(method="lm", formula=y~x) +labs(y="Number of sigma\nfactor\ngenes", )   +theme_linedraw()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), strip.text.x=element_blank())+
    facet_wrap(~gtdb_phylum,nrow = 1)+scale_size(guide = 'none')+
    stat_fit_tidy(method = "lm",
                  label.x = "right",
                  size=3,
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))

  p6 <- ggplot(reduced_members,aes(x = .data[[env]],y = kc_Cell_motility ,group = 1)) +
    geom_point() +
    stat_smooth(method="lm", formula=y~x) +labs(y="Number of KOs\nfor\ncell motility", x=req_env4[[env]])  +
    facet_wrap(~gtdb_phylum,nrow = 1)+scale_size(guide = 'none')+theme_linedraw()+
  theme(strip.text.x=element_blank())+
    stat_fit_tidy(method = "lm",
                  label.x = "right",
                  size=3,
                  method.args = list(formula = y ~ x),
                  mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                                after_stat(x_estimate),
                                                after_stat(x_p.value))))

  combined_plot_scattter_kegg_env<-ggarrange(p2,p3,p4,p5,p6,nrow = 5,labels = c("(A)","(B)","(C)","(D)","(E)"),heights=c(rep(0.9,4),1))
  combined_plot_scattter_kegg_env
  out_files<-paste0("/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/env_genome_properties_scatters/supplementary_scatter_genomic_prop_vs_",env,".svg")
  out_filep<-paste0("/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/env_genome_properties_scatters/png/supplementary_scatter_genomic_prop_vs_",env,".png")
  ggsave(file=out_files,device="svg",plot=combined_plot_scattter_kegg_env,width=1892, height=1580, units="px",dpi=170)
  ggsave(file=out_filep,device="png",plot=combined_plot_scattter_kegg_env,width=1892, height=1580, units="px",dpi=170)
}