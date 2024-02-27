req_env2<-c("TNs","TP.y","DOC.y" #replace "x = .data[[env]]"  wtih "x = log(.data[[env]])" ###to preserve this comment please add "(" before the X while replaceing
            #"An_Max_Temp","temperature"
)

nutri_TN_gs <- ggplot(members_table,aes(x = tn_lev,y = size_rep , weight=abundance_val, fill=tn_lev )) +
  geom_boxplot() +                                          #replace:TNs/TP.y/DOC.y  #replace:GC_content/size_rep
  #stat_smooth(method="lm", formula=y~x) +                 #An_Mean_Temp;An_Max_Temp;An_Min_Temp;Temperature1
  #facet_wrap(~gtdb_phylum)+
  scale_size(guide = 'none')+scale_fill_manual(values=c("#A3DADC", "#62EF94"))+
  #stat_poly_eq(aes(weight=abundance_val,
  #                 label =  paste(..adj.rr.label..)),
  #             label.x.npc = "right", label.y.npc = 0.15,
  #             formula = y~x, parse = TRUE, size=7)+
  theme_bw()+theme(legend.position="none")

nutri_TP_gs <- ggplot(members_table,aes(x = tp_lev,y = size_rep , weight=abundance_val , fill=tp_lev)) +
  geom_boxplot() +                                          #replace:TNs/TP.y/DOC.y  #replace:GC_content/size_rep
  #stat_smooth(method="lm", formula=y~x) +                 #An_Mean_Temp;An_Max_Temp;An_Min_Temp;Temperature1
  #facet_wrap(~gtdb_phylum)+
  scale_size(guide = 'none')+scale_fill_manual(values=c("#A3DADC", "#62EF94"))+
  #stat_poly_eq(aes(weight=abundance_val,
  #                 label =  paste(..adj.rr.label..)),
  #             label.x.npc = "right", label.y.npc = 0.15,
  #             formula = y~x, parse = TRUE, size=7)+
  theme_bw()+theme(legend.position="none")

nutri_DC_gs <- ggplot(members_table,aes(x = doc_lev,y = size_rep , weight=abundance_val, fill=doc_lev )) +
  geom_boxplot() +                                          #replace:TNs/TP.y/DOC.y  #replace:GC_content/size_rep
  #stat_smooth(method="lm", formula=y~x) +                 #An_Mean_Temp;An_Max_Temp;An_Min_Temp;Temperature1
  #facet_wrap(~gtdb_phylum)+
  scale_size(guide = 'none')+scale_fill_manual(values=c("#A3DADC", "#62EF94"))+
  #stat_poly_eq(aes(weight=abundance_val,
  #                 label =  paste(..adj.rr.label..)),
  #             label.x.npc = "right", label.y.npc = 0.15,
  #             formula = y~x, parse = TRUE, size=7)+
  theme_bw()+theme(legend.position="none")

grid_image<-grid.arrange(nutri_TN_gs,nutri_TP_gs,nutri_DC_gs,ncol=3)