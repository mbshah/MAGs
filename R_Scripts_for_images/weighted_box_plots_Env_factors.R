
library(epade)
library(weights)
library(ggpubr)

#add required Coloumns
members_table["tp_lev"]<-ifelse(members_table$TP.y>20,"High","Low")
members_table["tn_lev"]<-ifelse(members_table$TNs>600,"High","Low")
members_table["doc_lev"]<-ifelse(members_table$DOC.y>4000,"High","Low")
members_table["size_lev"]<-ifelse(members_table$size_rep<2,"Small",
                               ifelse(members_table$size_rep>4,"Big","Medium"))
members_table["temp_lev"]<-ifelse(members_table$An_Max_Temp>20,"High","Low")

#Genome Sizes  
#TP
t_p_val_tp<-wtd.t.test(subset(members_table,tp_lev=="Low")$size_rep, subset(members_table,tp_lev=="High")$size_rep,
           weight=subset(members_table,tp_lev=="Low")$abundance_val,
           weighty = subset(members_table,tp_lev=="High")$abundance_val,
           samedata=T)$coefficients["p.value"]
c1<-ggplot(members_table,aes(x=tp_lev,y=size_rep,weight=abundance_val))+
  geom_boxplot()+labs(y="Total Phosphorus\n\n\nGenome Size(Mb)", title="Genome Sizes")+
  annotate("text",label=paste("Weighted t-test p-value:",signif(t_p_val_tp, 3)),x=1.5,y=7.5)+
  theme(legend.position = "none",
                   axis.text.x = element_blank(),
                   #axis.ticks.y = element_blank(),
                   axis.title.x = element_blank(),
        plot.title=element_text(hjust=0.5)
                    )

#TN
t_p_val_tn<-wtd.t.test(subset(members_table,tn_lev=="Low")$size_rep, subset(members_table,tn_lev=="High")$size_rep,
           weight=subset(members_table,tn_lev=="Low")$abundance_val,
           weighty = subset(members_table,tn_lev=="High")$abundance_val,
           samedata=T)$coefficients["p.value"]
c2<-ggplot(members_table,aes(x=tn_lev,y=size_rep,weight=abundance_val))+
  geom_boxplot()+labs(y="Total Nitrogen\n\n\nGenome Size(Mb)")+
  annotate("text",label=paste("Weighted t-test p-value:",signif(t_p_val_tn, 3)),x=1.5,y=7.5)+
  theme(legend.position = "none",
                   axis.text.x = element_blank(),
                   #axis.ticks.y = element_blank(),
                   axis.title.x = element_blank()
                    )

#DOC
t_p_val_doc<-wtd.t.test(subset(members_table,doc_lev=="Low")$size_rep, subset(members_table,doc_lev=="High")$size_rep,
           weight=subset(members_table,doc_lev=="Low")$abundance_val,
           weighty = subset(members_table,doc_lev=="High")$abundance_val,
           samedata=T)$coefficients["p.value"]
c3<-ggplot(members_table,aes(x=doc_lev,y=size_rep,weight=abundance_val))+
  geom_boxplot()+labs(y="Dissolved Oxygen\nContent\n\nGenome Size(Mb)")+
  annotate("text",label=paste("Weighted t-test p-value:",signif(t_p_val_doc, 3)),x=1.5,y=7.5)+
  theme(legend.position = "none",
                   axis.text.x = element_blank(),
                   #axis.ticks.y = element_blank(),
                   axis.title.x = element_blank()
                    )

#TEMP
members_table["temp_lev"]<-ifelse(members_table$An_Max_Temp>20,"High","Low")
t_p_val_tempar<-wtd.t.test(subset(members_table,temp_lev=="Low")$size_rep, subset(members_table,temp_lev=="High")$size_rep,
           weight=subset(members_table,temp_lev=="Low")$abundance_val,
           weighty = subset(members_table,temp_lev=="High")$abundance_val,
           samedata=T)$coefficients["p.value"]
c4<-ggplot(members_table,aes(x=temp_lev,y=size_rep,weight=abundance_val))+
  geom_boxplot()+labs(y="Annual Max\nTemperature\n\nGenome Size(Mb)")+
  annotate("text",label=paste("Weighted t-test p-value:",signif(t_p_val_tempar, 3)),x=1.5,y=7.5)+
  theme(legend.position = "none",
                   #axis.text.y = element_blank(),
                   #axis.ticks.y = element_blank(),
                   axis.title.x = element_blank()
                    )


#GC_Content 
#TP
t_p_val_tp<-wtd.t.test(subset(members_table,tp_lev=="Low")$GC_content, subset(members_table,tp_lev=="High")$GC_content,
           weight=subset(members_table,tp_lev=="Low")$abundance_val,
           weighty = subset(members_table,tp_lev=="High")$abundance_val,
           samedata=T)$coefficients["p.value"]
c5<-ggplot(members_table,aes(x=tp_lev,y=GC_content,weight=abundance_val))+
  geom_boxplot()+labs(y="GC Content (%)", title="GC content",)+
  annotate("text",label=paste("Weighted t-test p-value:",signif(t_p_val_tp, 3)),x=1.5,y=71)+
  theme(legend.position = "none",
                   axis.text.x = element_blank(),
                   #axis.ticks.y = element_blank(),
                   axis.title.x = element_blank(),
        plot.title=element_text(hjust=0.5)
                    )

#TN
t_p_val_tn<-wtd.t.test(subset(members_table,tn_lev=="Low")$GC_content, subset(members_table,tn_lev=="High")$GC_content,
           weight=subset(members_table,tn_lev=="Low")$abundance_val,
           weighty = subset(members_table,tn_lev=="High")$abundance_val,
           samedata=T)$coefficients["p.value"]
c6<-ggplot(members_table,aes(x=tn_lev,y=GC_content,weight=abundance_val))+
  geom_boxplot()+labs(y="GC Content (%)")+
  annotate("text",label=paste("Weighted t-test p-value:",signif(t_p_val_tn, 3)),x=1.5,y=71)+
  theme(legend.position = "none",
                   axis.text.x = element_blank(),
                   #axis.ticks.y = element_blank(),
                   axis.title.x = element_blank()
                    )

#DOC
t_p_val_doc<-wtd.t.test(subset(members_table,doc_lev=="Low")$GC_content, subset(members_table,doc_lev=="High")$GC_content,
           weight=subset(members_table,doc_lev=="Low")$abundance_val,
           weighty = subset(members_table,doc_lev=="High")$abundance_val,
           samedata=T)$coefficients["p.value"]
c7<-ggplot(members_table,aes(x=doc_lev,y=GC_content,weight=abundance_val))+
  geom_boxplot()+labs(y="GC Content (%)")+
  annotate("text",label=paste("Weighted t-test p-value:",signif(t_p_val_doc, 3)),x=1.5,y=71)+
  theme(legend.position = "none",
                   axis.text.x = element_blank(),
                   #axis.ticks.y = element_blank(),
                   axis.title.x = element_blank()
                    )

#TEMP
members_table["temp_lev"]<-ifelse(members_table$An_Max_Temp>20,"High","Low")
t_p_val_tempar<-wtd.t.test(subset(members_table,temp_lev=="Low")$GC_content, subset(members_table,temp_lev=="High")$GC_content,
           weight=subset(members_table,temp_lev=="Low")$abundance_val,
           weighty = subset(members_table,temp_lev=="High")$abundance_val,
           samedata=T)$coefficients["p.value"]
c8<-ggplot(members_table,aes(x=temp_lev,y=GC_content,weight=abundance_val))+
  geom_boxplot()+labs(y="GC Content (%)")+
  annotate("text",label=paste("Weighted t-test p-value:",signif(t_p_val_tempar, 3)),x=1.5,y=71)+
  theme(legend.position = "none",
                   #axis.text.y = element_blank(),
                   #axis.ticks.y = element_blank(),
                   axis.title.x = element_blank()
                    )

my_plot<-ggarrange(c1,c5,c2,c6,c3,c7,c4,c8,ncol=2,nrow = 4)
out_file<-"/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/Box_plot_Env_wtd.svg"
ggsave(file=out_file,device="svg",plot=my_plot,width=1800, height=1500, units="px",dpi=150)