cols_phyla <- c("Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota","Verrucomicrobiota")
reduced_members<-MAGs_table_new[MAGs_table_new$gtdb_phyla %in% cols_phyla,]
reduced_members$gtdb_phyla<-as.factor(reduced_members$gtdb_phyla)

p1<-ggplot(data=reduced_members, aes(y=reassembly_gc, x=reassembly_size_mb)) + #replace:Total_Redundancy,c80,c100
  geom_point()+labs(y="\nGC content\n(%)",) +theme_linedraw()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(colour = "black"))+
  geom_smooth(method="lm",formula=y~x)+
  facet_wrap(~gtdb_phyla,nrow=1)+
  stat_fit_tidy(method = "lm",
                label.x = "right",
                size=3,
                method.args = list(formula = y ~ x),
                mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                              after_stat(x_estimate),
                                              after_stat(x_p.value))))

p2<-ggplot(data=reduced_members, aes(y=Total_Redundancy, x=reassembly_size_mb)) + #replace:Total_Redundancy,c80,c100
  geom_point()+labs(y="\nCount of\nredundant genes", ) +theme_linedraw()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), strip.text.x=element_blank())+
  geom_smooth(method="lm",formula=y~x)+
  facet_wrap(~gtdb_phyla,nrow=1)+
  stat_fit_tidy(method = "lm",
                label.x = "right",
                size=3,
                method.args = list(formula = y ~ x),
                mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                              after_stat(x_estimate),
                                              after_stat(x_p.value))))

p3<-ggplot(data=reduced_members, aes(y=coding_perc, x=reassembly_size_mb)) + #replace:Total_Redundancy,c80,c100
  geom_point()+labs(y="\nCoding regions\n(%)", ) +theme_linedraw()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), strip.text.x=element_blank())+
  geom_smooth(method="lm",formula=y~x)+
  facet_wrap(~gtdb_phyla,nrow=1)+
  stat_fit_tidy(method = "lm",
                label.x = "left",label.y="bottom",
                size=3,
                method.args = list(formula = y ~ x),
                mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                              after_stat(x_estimate),
                                              after_stat(x_p.value))))

p4 <- ggplot(reduced_members,aes(x = reassembly_size_mb,y = pc_abc_transporters ,group = 1)) +
  geom_point() +
  stat_smooth(method="lm", formula=y~x) +labs(y="Fraction of KOs\nfor ABC\ntransporters (%)",)  +theme_linedraw()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), strip.text.x=element_blank())+
  facet_wrap(~gtdb_phyla,nrow = 1)+scale_size(guide = 'none')+
  #stat_poly_eq(aes(,
  #                 label =  paste(..adj.rr.label..)),
  #             label.x.npc = "right", label.y.npc = 0.15,
  #             formula = y~x, parse = TRUE, size=3)+
  stat_fit_tidy(method = "lm",
                label.x = "left",
                size=3,
                method.args = list(formula = y ~ x),
                mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                              after_stat(x_estimate),
                                              after_stat(x_p.value))))

p5 <- ggplot(reduced_members,aes(x = reassembly_size_mb,y = pc_sigmafactors ,group = 1)) +
  geom_point() +
  stat_smooth(method="lm", formula=y~x) +labs(y="Fraction of\nsigma factor\ngenes (%)", )   +theme_linedraw()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), strip.text.x=element_blank())+
  facet_wrap(~gtdb_phyla,nrow = 1)+scale_size(guide = 'none')+
  stat_fit_tidy(method = "lm",
                label.x = "left",
                size=3,
                method.args = list(formula = y ~ x),
                mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                              after_stat(x_estimate),
                                              after_stat(x_p.value))))

p6 <- ggplot(reduced_members,aes(x = reassembly_size_mb,y = pc_cell_motility ,group = 1)) +
  geom_point() +
  stat_smooth(method="lm", formula=y~x) +labs(y="Fraction of KOs\nfor cell\nmotility (%)", x="Genome size (Mb)")  +
  facet_wrap(~gtdb_phyla,nrow = 1)+scale_size(guide = 'none')+theme_linedraw()+
  theme(strip.text.x=element_blank())+
  stat_fit_tidy(method = "lm",
                label.x = "left",
                size=3,
                method.args = list(formula = y ~ x),
                mapping = aes(label = sprintf("β = %.3g\np-value = %.3g",
                                              after_stat(x_estimate),
                                              after_stat(x_p.value))))

combined_plot_scattter_kegg<-ggarrange(p1,p2,p3,p4,p5,p6,nrow = 6,labels = c("(A)","(B)","(C)","(D)","(E)","(F)"), heights=c(rep(0.9,5),1))
combined_plot_scattter_kegg
out_files<-"/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/fig5/fig5-2.svg"
ggsave(file=out_files,device="svg",plot=combined_plot_scattter_kegg,width=1892, height=1892, units="px",dpi=170)
out_filep<-"/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/fig5/fig5-2.png"
ggsave(file=out_filep,device="png",plot=combined_plot_scattter_kegg,width=1892, height=1892, units="px",dpi=170)