if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")
profile_tsv <- "/home/hb0358/PycharmProjects/mbs_general/MAGs/outfile/dRep__gff/cluster_abundance_profile.tsv"
metadata_file <- "/home/hb0358/PycharmProjects/mbs_general/MAGs/Assembled_Metagenomes/CED91220.tsv"
metaData <- read.table(file=metadata_file, sep="\t", header = TRUE)
readCount <- na.omit(t(read.table(file=profile_tsv, sep="\t", header=TRUE, row.names=1,)))
dds<- DESeqDataSetFromMatrix(countData = readCount, colData = metaData, design=~nutrient_level)
View (counts(dds))
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
readCountNormalized <- counts(dds, normalized=TRUE)
write.table(data.frame("abundance"=rownames(readCountNormalized),readCountNormalized), file="/home/hb0358/PycharmProjects/mbs_general/MAGs/outfile/dRep__gff/cluster_abundance_profile_normalized_deseq2.tsv", sep="\t", quote=F,row.names=F)
