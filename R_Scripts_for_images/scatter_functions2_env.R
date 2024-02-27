cols_phyla <- c("Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota","Verrucomicrobiota")
reduced_members<-MAGs_table_new[MAGs_table_new$gtdb_phyla %in% cols_phyla,]
my_plot_p<-ggplot(data=reduced_members, aes(y=Total_Redundancy, x=log(weighted_TP))) + #replace:Total_Redundancy,c80,c100
  geom_point()+labs(y="Redundant genes", x="Log of Total Phosphorus")+
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
my_plot_n<-ggplot(data=reduced_members, aes(y=Total_Redundancy, x=log(weighted_TN))) + #replace:Total_Redundancy,c80,c100
  geom_point()+labs(y="Redundant genes", x="Log of Nitrogen")+
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
my_plot_C<-ggplot(data=reduced_members, aes(y=Total_Redundancy, x=log(weighted_DOC))) + #replace:Total_Redundancy,c80,c100
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