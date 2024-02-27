tabundance<-t(abundance_profile)
row.names(meta_data)<-meta_data$ID
example_NMDS<-metaMDS(tabundance,k=2)
ordiplot(example_NMDS,type="n")
ordisurf(example_NMDS,meta_data$TP)
orditorp(example_NMDS,display="sites",air=0.1, cex=1)