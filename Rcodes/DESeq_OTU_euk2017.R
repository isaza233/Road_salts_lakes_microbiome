## ressources https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf
## ressources https://lashlock.github.io/compbio/R_presentation.html

prok_abondance_taxo <- read.table("data/Euk_samples_PR2_2.csv", sep=";",dec=".", head=T)

chem <- read.table("data/chem.csv", sep=";",dec=".", head=T)

library(DESeq2)

##S?lectionner les stations de surface seulement
euk_surf <- prok_abondance_taxo[,-c(24,27,29,32,34,35,36)]

## select year 2017 for chemistry
chem17 <- subset(chem, subset=chem$Year=="2017", select=colnames(chem))

## mean per triplicates for chemistry
chem_mean <- aggregate(.~ chem17$ID+chem17$Lake+chem17$Month, data=chem17, FUN="mean")

## s'assurer que chem_mean, genus_t_surf
euk_sort <- euk_surf[,match(chem_mean$`chem17$ID`,colnames(euk_surf))]
rownames(euk_sort) <- prok_abondance_taxo$OTU

## Cr?er un vecteur pour conductivity low or high 
cond <- vector()
for (i in 1:nrow(chem_mean))
{
  if (chem_mean$`chem17$Lake`[i]=="Saint-Charles") {cond[i] <- "low"}
  else if (chem_mean$`chem17$Lake`[i]=="Saint-Augustin") {cond[i] <- "high"}
  else if (chem_mean$`chem17$Lake`[i]=="Clement") {cond[i] <- "high"}
  else if (chem_mean$`chem17$Lake`[i]=="Clair") {cond[i] <- "low"}
}

## Cr?er un dataframe de metadata
grouping<- data.frame(chem_mean$`chem17$ID`,chem_mean$`chem17$Lake`,cond, chem_mean$Season)

for(i in 1:nrow(grouping))
{
  if (grouping[i,4]==2) {grouping[i,4] <- "Winter"}
  else if (grouping[i,4]==1) {grouping[i,4] <- "Summer"}
}

colnames(grouping) <- c("ID", "lake","cond", "season")
grouping$season <- as.factor(grouping$season)

## add one to every OTU for every sample to be able to calculate the geometric mean, which can't be done if at least
## one entry is 0. As the log is compute, adding 1 won't have that much effect
## see forum: https://biostar.usegalaxy.org/p/30097/
euk_sort1 <- euk_sort # to store a copy of euk_sort

for (i in 1:nrow(euk_sort))
{
  for (j in 1:ncol(euk_sort))
  {
    euk_sort[i,j] <- euk_sort[i,j]+1
  }
}

## DESeq analysis seasons
ddss <- DESeqDataSetFromMatrix(countData=euk_sort, colData=grouping, design=~cond+season)
ddss <- DESeq(ddss)
ress <- results(ddss)

## results 
seasonres <- matrix(nrow=length(which(ress$padj<0.05)), ncol=5)
seasonres[,1] <- rownames(euk_sort)[which(ress$padj< 0.05)]
for (i in 1:length(which(ress$padj<0.05))){seasonres[i,2] <- as.character(euk_surf$taxo_3.1[which(euk_surf$OTU==seasonres[i,1])])}
for (i in 1:length(which(ress$padj<0.05))){seasonres[i,3] <- as.character(euk_surf$genus.1[which(euk_surf$OTU==seasonres[i,1])])}
for (i in 1:length(which(ress$padj<0.05))){seasonres[i,4] <- as.character(euk_surf$specie.1[which(euk_surf$OTU==seasonres[i,1])])}
seasonres[,5] <- ress$log2FoldChange[which(ress$padj< 0.05)]

seasonres <- as.data.frame(seasonres)
colnames(seasonres) <- c("OTU","Phylum", "Genus","Species", "log2Foldchange")

seasonagg <- aggregate(as.numeric(as.character(seasonres$log2Foldchange)), by=list(seasonres$Phylum, seasonres$Genus, seasonres$Species), 
                       FUN="mean")

## DESeq analysis cond
ddss <- DESeqDataSetFromMatrix(countData=euk_sort, colData=grouping, design=~season+cond)
ddss <- DESeq(ddss)
ress <- results(ddss)

## results 
condres <- matrix(nrow=length(which(ress$padj<0.05)), ncol=5)
condres[,1] <- rownames(euk_sort)[which(ress$padj< 0.05)]
for (i in 1:length(which(ress$padj<0.05))){condres[i,2] <- as.character(euk_surf$taxo_3.1[which(euk_surf$OTU==condres[i,1])])}
for (i in 1:length(which(ress$padj<0.05))){condres[i,3] <- as.character(euk_surf$genus.1[which(euk_surf$OTU==condres[i,1])])}
for (i in 1:length(which(ress$padj<0.05))){condres[i,4] <- as.character(euk_surf$specie.1[which(euk_surf$OTU==condres[i,1])])}
condres[,5] <- ress$log2FoldChange[which(ress$padj< 0.05)]

condres <- as.data.frame(condres)
colnames(condres) <- c("OTU","Phylum", "Genus","Species", "log2Foldchange")

condagg <- aggregate(as.numeric(as.character(condres$log2Foldchange)), by=list(condres$Phylum, condres$Genus, condres$Species), 
                       FUN="mean")