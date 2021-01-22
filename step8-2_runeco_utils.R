# Title     : step8-2
# Objective : Run Spec-gen and spaa to identify generalists and specialists
# Created by: Manan Shah
# Created on: 10/12/20

###install and load_req_libraries###
#install.packages("vegan")
#install.packages("spaa")
#devtools::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)
library(RCurl)

###load and prepare table###
x<-"/home/hb0358/PycharmProjects/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles/Bacteria/l6_species_normalized_deseq2.tsv"
ofl<-"/home/hb0358/PycharmProjects/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles/Bacteria/Genralists_Specialists/gen_spec1000.tsv"
comm.tab<-read.table(file=x,sep="\t",row.names=1,header=TRUE,comment.char="@")
comm.tab<-t(comm.tab)
comm.tab<-comm.tab[,which(colSums(comm.tab)>0)]
comm.tab2<-comm.tab
ct2<-ceiling(comm.tab2)
comm.tab<-ct2
comm.tab.bin<-ceiling(comm.tab/max(comm.tab)) #table of present or absent


###run_nichwidth###
set.seed(45) #increase reproducibility
##all are differeent bootstrapped values###
res10<-spec.gen(comm.tab,n=10)
res100<-spec.gen(comm.tab,n=100)
res500<-spec.gen(comm.tab,n=500)
res1000<-spec.gen(comm.tab,n=1000)
res5000<-spec.gen(comm.tab,n=5000)




###assign colours and plot
cols<-c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677") #colorblind_pallate from Paul Tol https://davidmathlogic.com/colorblind/
res1000$color<-ifelse(res1000$sign=="GENERALIST",cols[1],ifelse(res1000$sign=="SPECIALIST",cols[3],cols[2]))
plot(colSums(comm.tab),colSums(comm.tab.bin)/dim(comm.tab.bin)[1],col=res1000$color,main = ,pch=19,log="x",xlab="Abundance",ylab="Occurrence")
legend("bottomright",levels(res1000$sign),col=cols,pch=19,inset=0.01,cex=0.7)


###write output###
write.table(res1000, file=ofl , sep="\t", quote=F, col.names=NA)

###reread modifiedfiles
ofl2<-"/home/hb0358/PycharmProjects/mbs_general/WMG/Kraken2_NT_Scaffolds_profiles2/Bacteria/Genralists_Specialists/gen_specret.tsv"
res_reread<-read.table(file=ofl2,sep="\t",row.names=1,header=TRUE,comment.char="@")
levels(res1000$sign)
plot(colSums(comm.tab),colSums(comm.tab.bin)/dim(comm.tab.bin)[1],col=res_reread$Color,main = ,pch=19,log="x",xlab="Abundance",ylab="Occurrence")
legend("bottomright",c("Gen","Spec","NS","Gen_MAG","Spec_MAG","NS_MAG"),col=cols,pch=19,inset=0.01,cex=0.7)

#to determine which bootstrap value will give ideal results###
gen_arr<-c(length(which(res10$sign=="GENERALIST")),length(which(res100$sign=="GENERALIST")),length(which(res500$sign=="GENERALIST")),length(which(res1000$sign=="GENERALIST")),length(which(res5000$sign=="GENERALIST")))
spec_arr<-c(length(which(res10$sign=="SPECIALIST")),length(which(res100$sign=="SPECIALIST")),length(which(res500$sign=="SPECIALIST")),length(which(res1000$sign=="SPECIALIST")),length(which(res5000$sign=="SPECIALIST")))
nd_arr<-c(length(which(res10$sign=="NON SIGNIFICANT")),length(which(res100$sign=="NON SIGNIFICANT")),length(which(res500$sign=="NON SIGNIFICANT")),length(which(res1000$sign=="NON SIGNIFICANT")),length(which(res5000$sign=="NON SIGNIFICANT")))
plot (gen_arr)
plot (spec_arr)
plot (nd_arr)
c(length(which(res1000$sign=="GENERALIST")),length(which(res1000$sign=="SPECIALIST")),length(which(res1000$sign=="NON SIGNIFICANT")))
