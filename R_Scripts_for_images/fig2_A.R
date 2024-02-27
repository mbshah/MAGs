cols_phyla <- c("Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota","Verrucomicrobiota")
cols_genus<- c("Limnohabitans","Rhodoferax","Rhodoluna","Planktophila","UBA2093","UBA3006","Sediminibacterium","Flavobacterium","SXYR01","Cyanobium","UBA953")

for (factor1 in c("GC_content","size_rep")){
  my_dir<-"/mnt/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/MAGs_Paper_V4/graphics/"
  table_name<-paste0(factor1,"_slope_pvalue_df")
  new_df<-data.frame(matrix(vector(),17,17,
                            dimnames = list(c(),
                                            c("TNs.slope","TNs.pvalue","TNs.pstr",
                                              "TP.y.slope","TP.y.pvalue","TP.y.pstr",
                                              "DOC.y.slope","DOC.y.pvalue","DOC.y.pstr",
                                              "An_Max_Temp.slope","An_Max_Temp.pvalue","An_Max_Temp.pstr",
                                              "Temperature1.slope","Temperature1.pvalue","Temperature1.pstr",
                                            "Taxonomy","Level"))),stringsAsFactors=F)
  row.names(new_df)<-c("Bacteria","Pseudomonadota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota","Verrucomicrobiota",
                       "Limnohabitans","Rhodoferax","Rhodoluna","Planktophila","UBA2093","UBA3006","Sediminibacterium","Flavobacterium","SXYR01","Cyanobium","UBA953")
  for (factor2 in c("TNs","TP.y","DOC.y","An_Max_Temp","Temperature1")){
    factor2p<-ifelse(str_contains(factor2,"Temp")==F,paste0("log(",factor2,")"),factor2)
    formula<-paste0(factor1,"~",factor2p)
    regression<-summary(lm(formula,weights=abundance_val, data = members_table))["coefficients"]
    slope<-regression$coefficients[factor2p,"Estimate"]
    pvalue<-signif(regression$coefficients[factor2p,"Pr(>|t|)"])
    str_code<-ifelse(between(pvalue,0,0.001),"***",ifelse(between(pvalue,0.001,0.01),"**",ifelse(between(pvalue,0.01,0.05),"*",ifelse(between(pvalue,0.05,0.1),"",""))))
    new_df["Bacteria",paste0(factor2,".slope")]<-slope
    new_df["Bacteria",paste0(factor2,".pvalue")]<-pvalue
    new_df["Bacteria",paste0(factor2,".pstr")]<-str_code
    new_df["Bacteria","Taxonomy"]<-"Bacteria"
    new_df["Bacteria","Level"]<-"Bacteria"
    for(phyla in cols_phyla){
      sub_data<-members_table[members_table$gtdb_phylum==phyla,]
      print(paste("formula:",formula,"phyla:",phyla))
      regression<-summary(lm(formula,weights=abundance_val, data = sub_data))["coefficients"]
      slope<-regression$coefficients[factor2p,"Estimate"]
      pvalue<-signif(regression$coefficients[factor2p,"Pr(>|t|)"])
      str_code<-ifelse(between(pvalue,0,0.001),"***",ifelse(between(pvalue,0.001,0.01),"**",ifelse(between(pvalue,0.01,0.05),"*",ifelse(between(pvalue,0.05,0.1),"",""))))
      new_df[phyla,paste0(factor2,".slope")]<-slope
      new_df[phyla,paste0(factor2,".pvalue")]<-pvalue
      new_df[phyla,paste0(factor2,".pstr")]<-str_code
      new_df[phyla,"Taxonomy"]<-phyla
      new_df[phyla,"Level"]<-"Phylum"
    }
    for(genus in cols_genus){
      sub_data<-members_table[members_table$gtdb_genus==genus,]
      #print(paste("formula:",formula,"genus:",genus))
      regression<-summary(lm(formula,weights=abundance_val, data = sub_data))["coefficients"]
      slope<-regression$coefficients[factor2p,"Estimate"]
      pvalue<-signif(regression$coefficients[factor2p,"Pr(>|t|)"])
      str_code<-ifelse(between(pvalue,0,0.0001),"****",ifelse(between(pvalue,0.0001,0.001),"***",ifelse(between(pvalue,0.001,0.01),"**",ifelse(between(pvalue,0.01,0.1),"*",ifelse(between(pvalue,0.1,1),"","")))))
      new_df[genus,paste0(factor2,".slope")]<-slope
      new_df[genus,paste0(factor2,".pvalue")]<-pvalue
      new_df[genus,paste0(factor2,".pstr")]<-str_code
      new_df[genus,"Taxonomy"]<-genus
      new_df[genus,"Level"]<-"Genus"
    }

  }
  new_df$DOC.y.pvalue.corrected<-p.adjust(new_df$DOC.y.pvalue,method = "bonferroni")
  new_df$DOC.y.pstr.corrected<-ifelse(between(new_df$DOC.y.pvalue.corrected,0,0.001),"***",
                                      ifelse(between(new_df$DOC.y.pvalue.corrected,0.001,0.01),"**",
                                             ifelse(between(new_df$DOC.y.pvalue.corrected,0.01,0.05),"*",
                                                    ifelse(between(new_df$DOC.y.pvalue.corrected,0.05,0.1),"",""))))
  new_df$TP.y.pvalue.corrected<-p.adjust(new_df$TP.y.pvalue,method = "bonferroni")
  new_df$TP.y.pstr.corrected<-ifelse(between(new_df$TP.y.pvalue.corrected,0,0.001),"***",
                                      ifelse(between(new_df$TP.y.pvalue.corrected,0.001,0.01),"**",
                                             ifelse(between(new_df$TP.y.pvalue.corrected,0.01,0.05),"*",
                                                    ifelse(between(new_df$TP.y.pvalue.corrected,0.05,0.1),"",""))))
  new_df$TNs.pvalue.corrected<-p.adjust(new_df$TNs.pvalue,method = "bonferroni")
  new_df$TNs.pstr.corrected<-ifelse(between(new_df$TNs.pvalue.corrected,0,0.001),"***",
                                      ifelse(between(new_df$TNs.pvalue.corrected,0.001,0.01),"**",
                                             ifelse(between(new_df$TNs.pvalue.corrected,0.01,0.05),"*",
                                                    ifelse(between(new_df$TNs.pvalue.corrected,0.05,0.1),"",""))))
  new_df$An_Max_Temp.pvalue.corrected<-p.adjust(new_df$An_Max_Temp.pvalue,method = "bonferroni")
  new_df$An_Max_Temp.pstr.corrected<-ifelse(between(new_df$An_Max_Temp.pvalue.corrected,0,0.001),"***",
                                      ifelse(between(new_df$An_Max_Temp.pvalue.corrected,0.001,0.01),"**",
                                             ifelse(between(new_df$An_Max_Temp.pvalue.corrected,0.01,0.05),"*",
                                                    ifelse(between(new_df$An_Max_Temp.pvalue.corrected,0.05,0.1),"",""))))
  new_df$Temperature1.pvalue.corrected<-p.adjust(new_df$Temperature1.pvalue,method = "bonferroni")
  new_df$Temperature1.pstr.corrected<-ifelse(between(new_df$Temperature1.pvalue.corrected,0,0.001),"***",
                                      ifelse(between(new_df$Temperature1.pvalue.corrected,0.001,0.01),"**",
                                             ifelse(between(new_df$Temperature1.pvalue.corrected,0.01,0.05),"*",
                                                    ifelse(between(new_df$Temperature1.pvalue.corrected,0.05,0.1),"",""))))
  assign(table_name,new_df)
  write.csv(new_df,paste0(my_dir,table_name,".csv"))
}

library(forcats)

req_order<-c("Bacteria",
             "Pseudomonadota", "Actinomycetota", "Bacteroidota",
             "Cyanobacteriota","Verrucomicrobiota","Limnohabitans","Rhodoferax","Rhodoluna",
             "Planktophila","UBA2093","UBA3006","Sediminibacterium","Flavobacterium","SXYR01",
             "Cyanobium","UBA953")


c1_gc<-GC_content_slope_pvalue_df %>%
ggplot(aes(x=TNs.slope,y=fct_relevel(Taxonomy,rev(req_order)),label=TNs.pstr.corrected,
                                      fill = TNs.pvalue.corrected < 0.05,))+labs(x="Total nitrogen")+
  geom_col()+theme(legend.position = "none",
                   #axis.text.y = element_blank(),
                   #axis.ticks.y = element_blank(),
                   axis.title.y = element_blank(),
                   panel.background = element_blank(),axis.line = element_line(colour = "black")
                    )+
  scale_fill_manual(values = c("#777777", "#44AA99")) +
  geom_text(size = 3, position = position_stack(vjust = 0.5))

c2_gc<-GC_content_slope_pvalue_df %>%
ggplot(aes(x=TP.y.slope,y=fct_relevel(Taxonomy,rev(req_order)),label=TP.y.pstr.corrected,
                                      fill = TP.y.pvalue.corrected < 0.05,))+labs(x="Total phosphorus")+
  geom_col()+theme(legend.position = "none",axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y = element_blank(),
                   panel.background = element_blank(),axis.line.x = element_line(colour = "black"))+
  scale_fill_manual(values = c("#777777", "#44AA99")) +
  geom_text(size = 3, position = position_stack(vjust = 0.5))

c3_gc<-GC_content_slope_pvalue_df %>%
ggplot(aes(x=DOC.y.slope,y=fct_relevel(Taxonomy,rev(req_order)),label=DOC.y.pstr.corrected,
                                      fill = DOC.y.pvalue.corrected < 0.05,))+labs(x="Dissolved organic carbon")+
  geom_col()+theme(legend.position = "none",axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y = element_blank(),
                   panel.background = element_blank(),axis.line.x = element_line(colour = "black"))+
  scale_fill_manual(values = c("#777777", "#44AA99")) +
  geom_text(size = 3, position = position_stack(vjust = 0.5))

c4_gc<-GC_content_slope_pvalue_df %>%
ggplot(aes(x=An_Max_Temp.slope,y=fct_relevel(Taxonomy,rev(req_order)),label=An_Max_Temp.pstr.corrected,
                                      fill = An_Max_Temp.pvalue.corrected < 0.05,))+labs(x="Annual max temperature")+
  geom_col()+theme(legend.position = "none",axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y = element_blank(),
                   panel.background = element_blank(),axis.line.x = element_line(colour = "black"))+
  scale_fill_manual(values = c("#777777", "#44AA99")) +
  geom_text(size = 3, position = position_stack(vjust = 0.5))

c5_gc<-GC_content_slope_pvalue_df %>%
ggplot(aes(x=Temperature1.slope,y=fct_relevel(Taxonomy,rev(req_order)),label=Temperature1.pstr.corrected,
                                      fill = Temperature1.pvalue.corrected < 0.05,))+labs(x="Recorded temperature")+
  geom_col()+theme(legend.position = "none",axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y = element_blank(),
                   panel.background = element_blank(),axis.line.x = element_line(colour = "black"))+
  scale_fill_manual(values = c("#777777", "#44AA99")) +
  geom_text(size = 3, position = position_stack(vjust = 0.5))


my_plot_gc<-ggarrange(c1_gc,c2_gc,c3_gc,c4_gc,c5_gc,nrow = 1, widths=c(1,0.7,0.7,0.7,0.7),label.y=)+
  geom_hline(yintercept = 11.5)
pvalue_gc_plot_combined<-annotate_figure(my_plot_gc, left = textGrob("                 Genus                                              Phylum", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                    top = textGrob("Slopes of regression of GC content against environmental parameters", gp = gpar(cex = 1)))+
  annotation_custom(grid.polygon(x=c(0.005, 0.005, 0.995, 0.995,0.005, 0.005, 0.995, 0.995,0.005, 0.005, 0.995, 0.995,0.005, 0.005, 0.995, 0.995),
                                 y=c(0.64, 0.9, 0.90, 0.64,0.64,0.07,0.07,0.64,0.9, 0.95, 0.95, 0.9,0.95,0.999,0.999,0.95),
                                 id.lengths=c(4,4,4,4), gp=gpar(fill=NA)))


pvalue_gc_plot_combined


#out_file<-"/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/pvalue_GC.svg"
#ggsave(file=out_file,device="svg",plot=pvalue_gc_plot_combined,width=1892, height=901, units="px",dpi=-150)