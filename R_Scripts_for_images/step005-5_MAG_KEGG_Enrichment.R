library(ggplot2)

setwd("~/Documents/Manuscripte/Manan_Metagenome/Streamlining/")
mags <- read.table("member_table.tsv", header=T, sep="\t", as.is = T)

small <- mags[as.numeric(mags$org_size) < 2,]
big <- mags[as.numeric(mags$org_size) > 4,]
medium <- mags[as.numeric(mags$org_size) > 2 & as.numeric(mags$org_size) < 4,]

# percentage with different genome sizes
dim(small)[1] / dim(mags)[1] # 26.8%
dim(medium)[1] / dim(mags)[1] # 58.3%
dim(big)[1] / dim(mags)[1] # 14.6%

# most abundant phyla in these groups
df0 <- data.frame(phylum=sort(table(mags$Phylum), decreasing = T), genome_size="All")
df1 <- data.frame(phylum=sort(table(small$Phylum), decreasing = T), genome_size="Small")
df2 <- data.frame(phylum=sort(table(medium$Phylum), decreasing = T),  genome_size="Medium")
df3 <- data.frame(phylum=sort(table(big$Phylum), decreasing = T),  genome_size="Big")

df <- rbind(df0, df1, df2, df3)

phyla <- c("Actinomycetota", "Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria", "Oligoflexia", "Deltaproteobacteria", "Bacteroidetes", "Firmicutes", "Planctomycetes", "Cyanobacteriota", "Verrucomicrobia", "Gemmatimonadetes", "Deinococcus-Thermus", "Armatimonadetes")
colors <- c("#88ccee", "#44aa99", "#332288", "#aaaa00", "#ee8866", "#77aadd", "#117733", "#ddcc77", "#999933", "#cc6677", "#882255", "#aa4499", "#dddddd", "#bbcc33")
names(colors) <- phyla

f <- ggplot(data = df, aes(x=genome_size, y = phylum.Freq, fill=reorder(phylum.Var1, -phylum.Freq)))+
  geom_bar(position = "fill", stat="identity", )+
  scale_fill_manual(values = colors)+
  theme_bw()+
  xlab("Genome Size")+
  ylab("Proportion of phyla")+
  labs(fill = "Phylum")
f 

ggsave("Composition_phyla.svg", f, width = 5, height = 5)

# Redo this part, might be wrong table:

mags <- read.csv("Final_Paper/tables/clusters_wphyp_sited_m.tsv", header=T, sep="\t", as.is=T)

mags$reassembly_size_mb

small <- mags[as.numeric(mags$reassembly_size_mb) < 2,]
big <- mags[as.numeric(mags$reassembly_size_mb) > 4,]
medium <- mags[as.numeric(mags$reassembly_size_mb) > 2 & as.numeric(mags$reassembly_size_mb) < 4,]

# percentage with different genome sizes
dim(small)[1] / dim(mags)[1] # 28.6%, 90
dim(medium)[1] / dim(mags)[1] # 58.3%, 183
dim(big)[1] / dim(mags)[1] # 12.7%, 40

# most abundant phyla in these groups
df0 <- data.frame(phylum=sort(table(mags$Kraken_Phyla), decreasing = T), genome_size="All")
df1 <- data.frame(phylum=sort(table(small$Kraken_Phyla), decreasing = T), genome_size="Small")
df2 <- data.frame(phylum=sort(table(medium$Kraken_Phyla), decreasing = T),  genome_size="Medium")
df3 <- data.frame(phylum=sort(table(big$Kraken_Phyla), decreasing = T),  genome_size="Big")

df <- rbind(df0, df1, df2, df3)
df$genome_size <- factor(df$genome_size, levels=c("All", "Small", "Medium", "Big"))

phyla <- c("Actinomycetota", "Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria", "Oligoflexia", "Deltaproteobacteria", "Bacteroidetes", "Firmicutes", "Planctomycetes", "Cyanobacteriota", "Verrucomicrobia", "Gemmatimonadetes", "Deinococcus-Thermus", "Armatimonadetes")
colors <- c("#88ccee", "#44aa99", "#332288", "#aaaa00", "#ee8866", "#77aadd", "#117733", "#ddcc77", "#999933", "#cc6677", "#882255", "#aa4499", "#dddddd", "#bbcc33")
names(colors) <- phyla

f <- ggplot(data = df, aes(x=genome_size, y = phylum.Freq, fill=reorder(phylum.Var1, -phylum.Freq)))+
  geom_bar(position = "fill", stat="identity", )+
  scale_fill_manual(values = colors)+
  theme_bw()+
  xlab("Genome Size")+
  ylab("Proportion of phyla")+
  labs(fill = "Phylum")
f 

ggsave("Composition_phyla.svg", f, width = 5, height = 5)

# KOs
BiocManager::install("gage")
library(gage)
library(pathview)
library(dplyr)

more<-0
less<-0
KOs <- read.table("/home/hb0358/sciebo/mag+gen_spec/MAGs/MAGs_compiled_Paper/tables/gene_differential_presence.tsv", header=T, sep="\t", as.is = T)
# KO should be present in at least 10% of the genomes
KOs2 <- KOs[apply(KOs[,2:3], 1, function(x) any(x>0.1)),] # 2839
# presence in genomes at least twice at high
small <- KOs2[log2((KOs2$Small+0.0001)/(KOs2$Big+0.0001))>=1,'KO'] # 99
big <- KOs2[log2((KOs2$Small+0.0001)/(KOs2$Big+0.0001))<=(-1),'KO'] # 1617
both <- setdiff(setdiff(KOs2$KO, more), less) # 1123
all <- KOs$KO # 5777

KEGGEnrichment <- function (sign.ko, background.ko, kegg.ko)
{
  annList <- list()
  for (i in 1:length(kegg.ko))
  {
    pathway <- names(kegg.ko)[i]
    pathway.ko <- kegg.ko[[i]]
    annotatedMoleculeList <- intersect(pathway.ko, sign.ko)
    annotatedBackgroundList <- intersect(pathway.ko, background.ko)
    ann <- list(pathwayName = "not known", annMoleculeList = character(), annMoleculeNumber = 0, annBgMoleculeList = character(), annBgNumber = 0, moleculeNumber = 0, bgNumber = 0, pvalue = 1)
    ann$pathwayName <- pathway
    ann$annMoleculeList <- annotatedMoleculeList
    ann$annMoleculeNumber <- length(annotatedMoleculeList)
    ann$annBgMoleculeList <- annotatedBackgroundList
    ann$annBgNumber <- length(annotatedBackgroundList)
    ann$moleculeNumber <- length(sign.ko)
    ann$bgNumber <- length(background.ko)
    ann$pvalue <- 1 - phyper(ann$annMoleculeNumber - 1, ann$annBgNumber, ann$bgNumber - ann$annBgNumber, ann$moleculeNumber)
    annList[[i]] <- ann
  }
  return(annList)
}
signKEGGEnrichment <- function (KEGG.enrichment, threshold)
{
  p.value <- sapply(KEGG.enrichment, function(x) return(x$pvalue))
  pathway.table <- unname(t(as.data.frame(lapply(KEGG.enrichment[p.value <= threshold], function(x) c(x$pathwayName, paste(x$annMoleculeList, collapse = ", "), x$pvalue)))))
  if (dim(pathway.table)[1] != 0)
  {
    colnames(pathway.table) <- c("Pathway", "KOs", "P-value")
    pathway.table <- pathway.table[order(as.numeric(pathway.table[, "P-value"])), ]
    return(pathway.table)
  }else
  {
    return(NULL)
  }
}

background.ko <-all
sign.ko <- big
#kegg.ko  <- kegg.gsets(species = "ko", id.type = "kegg")$kg.sets
threshold <- 0.01

KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
sign.KEGG.enrichment.big <- signKEGGEnrichment(KEGG.enrichment, threshold)

sign.ko <- small
KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
sign.KEGG.enrichment.small <- signKEGGEnrichment(KEGG.enrichment, threshold)

sign.ko <- both
KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
sign.KEGG.enrichment.both <- signKEGGEnrichment(KEGG.enrichment, threshold)

write.table(sign.KEGG.enrichment.big, file = "MAGs/outfile/dRep__gff/KEGG_Enrichment_big.tsv", col.names = NA, sep = "\t")
write.table(sign.KEGG.enrichment.small, file = "MAGs/outfile/dRep__gff/KEGG_Enrichment_small.tsv", col.names = NA, sep = "\t")
write.table(sign.KEGG.enrichment.both, file = "MAGs/outfile/dRep__gff/KEGG_Enrichment_both.tsv", col.names = NA, sep = "\t")

write.table(big, file= "../outfile/dRep__gff/KOs_big.tsv")
write.table(small, file= "../outfile/dRep__gff/KOs_small.tsv")
write.table(both, file= "../outfile/dRep__gff/KOs_both.tsv")

# MAGs, Fig5

cor.test(mags$reassembly_size_mb,mags$perc_coding_combined)

# look at other annotation, e.g. transporters!