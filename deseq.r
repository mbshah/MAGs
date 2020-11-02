#if (!requireNamespace("BiocManager", quietly = TRUE))
#  BiocManager::install(version = "3.11")
#BiocManager::install("DESeq2")
library("DESeq2")
profile_tsv <- "/home/hb0358/PycharmProjects/mbs_general/MAGs/outfile/cluster_abundance_profile_transformed.tsv"
metadata_file <- "/home/hb0358/PycharmProjects/mbs_general/MAGs/outfile/meta.csv"
metaData <- read.table(file=metadata_file, sep="\t", header = TRUE)
readCount <- read.table(file=profile_tsv, sep="\t", header=TRUE, row.names=1)
dds<- DESeqDataSetFromMatrix(countData = readCount, colData = metaData, design=~elab)
View (counts(dds))
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
readCountNormalized <- counts(dds, normalized=TRUE)
write.table(readCountNormalized, file="/home/hb0358/PycharmProjects/mbs_general/MAGs/outfile/cluster_abundance_profile_normalized_deseq2.tsv", sep="\t", quote=F, col.names=NA)
