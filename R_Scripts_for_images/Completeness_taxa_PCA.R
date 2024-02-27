meatbolMAT<-t(metabolic_completeness_new[rowSums(metabolic_completeness_new)>0,]) # for metabolis_completeness
#meatbolMAT<-t(metabolic_redundandancy[rowSums(metabolic_redundandancy)>0,]) #for Redundancy
MeatbolMAT.pca<-prcomp(meatbolMAT, center = T, scale. = T)
#install_github("vqv/ggbiplot")
library(ggbiplot)
summary(MeatbolMAT.pca)

row.names(MAGs_table)<-MAGs_table$cluster_ID
MAGs_table<-MAGs_table[order(MAGs_table$cluster_ID),]
#ggbiplot(MeatbolMAT.pca,  groups = log(MAGs_table$weighted_DOC),var.axes=FALSE)+theme_bw()+scale_colour_gradientn(colours = topo.colors(10))
#ggbiplot(MeatbolMAT.pca,  groups = log(MAGs_table$weighted_TP),var.axes=FALSE)+theme_bw()+scale_colour_gradientn(colours = topo.colors(10))
#ggbiplot(MeatbolMAT.pca,  groups = log(MAGs_table$weighted_TN),var.axes=FALSE)+theme_bw()+scale_colour_gradientn(colours = topo.colors(10))
#ggbiplot(MeatbolMAT.pca,  groups = log(MAGs_table$reassembly_size_mb),var.axes=FALSE)+theme_bw()+scale_colour_gradientn(colours = topo.colors(10))
pca_plot<-ggbiplot(MeatbolMAT.pca,  groups = MAGs_table$gtdb_phyla, var.axes=FALSE)+theme_bw()+labs(color="Phyla")
out_files<-"/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/PCATAXA.svg"
ggsave(file=out_files,device="svg",plot=pca_plot,width=1000, height=800, units="px",dpi=150)
out_filep<-"/home/hb0358/Disk2/Disk2/Data/nextcloud/Sciebo/mag+gen_spec/MAGs/REDONE_NEW/Final_Set/supplementary/PCATAXA.png"
ggsave(file=out_filep,device="png",plot=pca_plot,width=1000, height=800, units="px",dpi=150)