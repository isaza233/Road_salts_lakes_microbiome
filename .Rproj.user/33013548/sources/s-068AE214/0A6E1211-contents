library(vegan)
library(pairwiseAdonis)
euk_abondance_taxo <- read.table("data/Euk_samples_PR2_2.csv", sep=";",dec=".", head=T)
prok_abondance_taxo <- read.table("data/Prok_abondance_taxo.csv", sep=";",dec=".", head=T)
chloroplasts_abondance_taxo <- read.table("data/Chloroplasts.csv", sep=";",dec=".", head=T)
chem <- read.table("data/chem.csv", sep=";",dec=".", head=T)
microscopy2 <- microscopy <- read.table("data/Sommaire2.csv", sep=";",dec=".", head=T)
pigments <- read.table("data/Sommaire_hplc.csv", sep=";",dec=".", head=F, colClasses = "character")

## transform areas in concentrations
ve=2.7
vi=0.04

#filtration volume in Liters from mL
for (i in 5:ncol(pigments))
{
  pigments[7,i] <- as.numeric(pigments[7,i])/1000
}

#actual calculation
for (i in 5:ncol(pigments))
{
  for (j in 9:nrow(pigments))
  {
    pigments[j,i] <- as.numeric(pigments[j,i])*as.numeric(pigments[j,4])*(ve/vi)*(1/as.numeric(pigments[7,i]))
  }
}

## calculate mean per sample 
pigments2 <- t(pigments[c(1,5,8:nrow(pigments)),c(5:ncol(pigments))])
pigments2 <- as.data.frame(pigments2)
colnames(pigments2) <- c("Lake","Season","ID",pigments[c(9:nrow(pigments)),2])

for (j in 4:ncol(pigments2))
{pigments2[,j] <- as.numeric(as.character(pigments2[,j]))}

pigments_mean <- aggregate(.~pigments2$ID, pigments2, FUN=mean)

## Aggregate 18S and 16S par groupes et transposer en sites=rows, groupes=colonnes
pgroups <- aggregate(prok_abondance_taxo[,c(10:52)],list(prok_abondance_taxo$taxo2), sum)
egroups <- aggregate(euk_abondance_taxo[,c(10:52)],list(euk_abondance_taxo$taxo_3), sum)

pgroupst <- as.data.frame(t(pgroups[,c(2:ncol(pgroups))]))
colnames(pgroupst) <- pgroups[,1]
rownames(pgroupst) <- colnames(pgroups)[-1]
egroupst <- as.data.frame(t(egroups[,c(2:ncol(egroups))]))
colnames(egroupst) <- egroups[,1]
rownames(egroupst) <- colnames(egroups)[-1]

## Aggregate 18S and 16S par genres et transposer en sites=rows, genres=colonnes
pgenus <- aggregate(prok_abondance_taxo[,c(10:52)],list(prok_abondance_taxo$taxo6), sum)
egenus <- aggregate(euk_abondance_taxo[,c(10:52)],list(euk_abondance_taxo$genus), sum)

pgenust <- as.data.frame(t(pgenus[,c(2:ncol(pgenus))]))
colnames(pgenust) <- pgenus[,1]
rownames(pgenust) <- colnames(pgenus)[-1]
egenust <- as.data.frame(t(egenus[,c(2:ncol(egenus))]))
colnames(egenust) <- egenus[,1]
rownames(egenust) <- colnames(egenus)[-1]

## transposer en sites=row, OTUs=columns
prok_t <- t(prok_abondance_taxo[,c(10:52)])
prok_t <- as.data.frame(prok_t)
colnames(prok_t) <- prok_abondance_taxo[,1]
rownames(prok_t) <- colnames(prok_abondance_taxo)[10:52]

euk_t <- t(euk_abondance_taxo[,c(10:52)])
euk_t <- as.data.frame(euk_t)
colnames(euk_t) <- euk_abondance_taxo[,1]
rownames(euk_t) <- colnames(euk_abondance_taxo)[10:52]

chloro_t <- t(chloroplasts_abondance_taxo[,c(5:47)])
chloro_t <- as.data.frame(chloro_t)
colnames(chloro_t) <- chloroplasts_abondance_taxo[,1]
rownames(chloro_t) <- colnames(chloroplasts_abondance_taxo)[5:47]

# Calculer les biovolumes 
microscopy1 <- microscopy2 #not to modified the original, which is needed for the counts 

for (j in 7:ncol(microscopy1))
{
  for (i in 1:nrow(microscopy1))
  {
    if (is.na(microscopy1$Volume[i])==TRUE)
    {microscopy1[i,j] <- microscopy1[i,j]}
    else 
    {microscopy1[i,j] <- microscopy1[i,j]*microscopy1$Volume[i]}
  }
}

# transposer mais retirer colonne 13 "LCT_oct_16", 26 "LSA_oct_16" pour retirer les ?chantilons de 2016 
microscopy <- t(microscopy2[,c(6:12,14:25,27:34)])
microscopy <- as.data.frame(microscopy)
colnames(microscopy) <- microscopy2[,2]
rownames(microscopy) <- colnames(microscopy2[c(6:12,14:25,27:34)])

micro_biovol <- t(microscopy1[,c(6:12,14:25,27:34)])
micro_biovol <- as.data.frame(micro_biovol)
colnames(micro_biovol) <- microscopy1[,2]
rownames(micro_biovol) <- colnames(microscopy1[c(6:12,14:25,27:34)])

## average chemistry for three sampling replicates in one row 
chem_mean <- aggregate(.~ chem$ID+chem$Lake+chem$Year+chem$Month, data=chem, FUN="mean")
chem_mean17 <- chem_mean[-c(which(chem_mean$Year=="2016")),]
colnames(chem_mean17)[14] <- "Temperature"
colnames(chem_mean)[11:16] <- c("Column_conductivity", "Conductivity",
                                "Column_temperature", "Temperature", "Water_color", 
                                "Alcalinity")

## s'assurer que chem_mean et euk_t sont dans le m?me ordre 
chloro_OTU_sort <- chloro_t[match(chem_mean17$`chem$ID`,rownames(chloro_t)),]
eOTU_sort <- euk_t[match(chem_mean17$`chem$ID`,rownames(euk_t)),]
egroups_sort <- egroupst[match(chem_mean17$`chem$ID`,rownames(egroupst)),]
egenus_sort <- egenust[match(chem_mean17$`chem$ID`,rownames(egenust)),]
pOTU_sort <- prok_t[match(chem_mean17$`chem$ID`,rownames(prok_t)),]
pgroups_sort <- pgroupst[match(chem_mean17$`chem$ID`,rownames(pgroupst)),]
pgenus_sort <- pgenust[match(chem_mean17$`chem$ID`,rownames(pgenust)),]

#microscopy_sort <- microscopy[match(chem_mean17$`chem$ID`,rownames(microscopy)),] #### microscopy counts if needed
micro_biovol_sort <- micro_biovol[match(chem_mean17$`chem$ID`,rownames(micro_biovol)),]

pigments_mean_sort <- pigments_mean[match(chem_mean17$`chem$ID`,pigments_mean$`pigments2$ID`),]
rownames(pigments_mean_sort) <- pigments_mean_sort[,1]
pigments_mean_sort <- pigments_mean_sort[,-c(1:4)]

##################################### NMDS #################################### 
tiff(filename="images/Figure4_Fournier.tif", width=25, height=18, units="cm", pointsize=10, res=300)
graphics::layout(mat=matrix(c(3,2,1,0,6,5,4,10,6,5,4,11,7,8,9,0), nrow=4, ncol=4, byrow=TRUE), 
                 heights=c(1,0.5,0.5,1))
par(mar=c(4,4,0,0.5))
# Prok OTUs
permp <- adonis(pOTU_sort~chem_mean17$`chem$Lake`*chem_mean17$Season, permutations = 999, method = "bray") #Lake pvalue=0.001, season pvalue=0.001
aaa <- pairwise.adonis(vegdist(pOTU_sort, method="bray"), factors=chem_mean17$`chem$Lake`) #LSC=LCR everything else is different
nmds_p <- metaMDS(pOTU_sort, distance = "bray") #stress 0.11

plot(pch=c(21,22,23,24)[chem_mean17$Lake],x=nmds_p$points[,1], y=nmds_p$points[,2],
     bg=c("red","blue")[chem_mean17$Season], cex=2, ylab="NMDS axis 2", 
     xlab=NA, las=1, cex.axis=1.2, cex.lab=1.4)

mtext(outer=F,side=1,text="NMDS axis 1", line=2.3, cex=0.9)
text(x=((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1], y=((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3], labels=" C. 16S OTU - St. 0.11", cex=1.5, adj=0)

# Prok genus
permpge <- adonis2(pgenus_sort~chem_mean17$`chem$Lake`*chem_mean17$Season, permutations = 999, method = "bray") #Lake pvalue=0.001, season pvalue=0.03
aaa <- pairwise.adonis(vegdist(pgenus_sort, method="bray"), factors=chem_mean17$`chem$Lake`) #LSC=tout everything else is different
nmds_pge <- metaMDS(pgenus_sort, distance = "bray") #stress 0.15

plot(pch=c(21,22,23,24)[chem_mean17$Lake],x=nmds_pge$points[,1], y=nmds_pge$points[,2],
     bg=c("red","blue")[chem_mean17$Season], cex=2, ylab="NMDS axis 2", 
     xlab=NA, las=1, cex.axis=1.2, cex.lab=1.4)

mtext(outer=F,side=1,text="NMDS axis 1", line=2.3, cex=0.9)
text(x=((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1], y=((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3], labels="B. 16S Genus - St. 0.15", cex=1.5, adj=0)

# Prok groups
permpgr <- adonis2(pgroups_sort~chem_mean17$`chem$Lake`*chem_mean17$Season, permutations = 999, method = "bray") #Lake pvalue=0.002, season pvalue=NS
aaa <- pairwise.adonis(vegdist(pgroups_sort, method="bray"), factors=chem_mean17$`chem$Lake`) #LCR and LCL =! LSA
nmds_pgr <- metaMDS(pgroups_sort, distance = "bray") #stress 0.24

plot(pch=c(21,22,23,24)[chem_mean17$Lake],x=-nmds_pgr$points[,1], y=-nmds_pgr$points[,2],
     bg=c("red","blue")[chem_mean17$Season], cex=2, ylab="NMDS axis 2", 
     xlab=NA, las=1, cex.axis=1.2, cex.lab=1.4)

mtext(outer=F,side=1,text="NMDS axis 1", line=2.3, cex=0.9)
text(x=((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1], y=((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3], labels="A. 16S Phylum - St. 0.24", cex=1.5, adj=0)

# Euk OTU
perme <- adonis2(eOTU_sort~chem_mean17$`chem$Lake`*chem_mean17$Season, permutations = 999, method = "bray") #Lake pvalue=0.001, season pvalue=0.001
aaa <- pairwise.adonis(vegdist(eOTU_sort, method="bray"), factors=chem_mean17$`chem$Lake`) #LCR=LSC, LSA=LCL
nmds_e <- metaMDS(eOTU_sort, distance = "bray") #stress 0.13

plot(pch=c(21,22,23,24)[chem_mean17$Lake],x=nmds_e$points[,1], y=nmds_e$points[,2],
     bg=c("red","blue")[chem_mean17$Season], cex=2, ylab="NMDS axis 2", 
     xlab=NA, las=1, cex.axis=1.2, cex.lab=1.4)

mtext(outer=F,side=1,text="NMDS axis 1", line=2.3, cex=0.9)
text(x=((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1], y=((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3], labels="F. 18S OTU - St. 0.13", cex=1.5, adj=0)

# Euk genus
permege <- adonis2(egenus_sort~chem_mean17$`chem$Lake`*chem_mean17$Season, permutations = 999, method = "bray") #Lake pvalue=0.001, season pvalue=0.001
aaa <- pairwise.adonis(vegdist(egenus_sort, method="bray"), factors=chem_mean17$`chem$Lake`) #LCR=! le reste is the same
nmds_ege <- metaMDS(egenus_sort, distance = "bray") #stress 0.16

plot(pch=c(21,22,23,24)[chem_mean17$Lake],x=nmds_ege$points[,1], y=nmds_ege$points[,2],
     bg=c("red","blue")[chem_mean17$Season], cex=2, ylab="NMDS axis 2", 
     xlab=NA, las=1, cex.axis=1.2, cex.lab=1.4)

mtext(outer=F,side=1,text="NMDS axis 1", line=2.3, cex=0.9)
text(x=((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1], y=((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3], labels="E. 18S Genus - St. 0.16", cex=1.5, adj=0)

# Euk groups
permegr <- adonis2(egroups_sort~chem_mean17$`chem$Lake`*chem_mean17$Season, permutations = 999, method = "bray") #Lake pvalue=0.009, season pvalue=0.001
aaa <- pairwise.adonis(vegdist(egroups_sort, method="bray"), factors=chem_mean17$`chem$Lake`) #everything is the same
nmds_egr <- metaMDS(egroups_sort, distance = "bray") #stress 0.18

plot(pch=c(21,22,23,24)[chem_mean17$Lake],x=nmds_egr$points[,1], y=nmds_egr$points[,2],
     bg=c("red","blue")[chem_mean17$Season], cex=2, ylab="NMDS axis 2", 
     xlab=NA, las=1, cex.axis=1.2, cex.lab=1.4)

mtext(outer=F,side=1,text="NMDS axis 1", line=2.3, cex=0.9)
text(x=((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1], y=((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3], labels="D. 18S Phylum - St. 0.18", cex=1.5, adj=0)

# Pigments 
permpig <- adonis2(pigments_mean_sort[,c(5,7:9,11:37,39:42,44:69,71:137,139:ncol(pigments_mean_sort))]~chem_mean17$`chem$Lake`*chem_mean17$Season, permutations = 999, method = "bray") #Lake pvalue=0.001, season pvalue=0.037
aaa <- pairwise.adonis(vegdist(pigments_mean_sort[,c(5,7:9,11:37,39:42,44:69,71:137,139:ncol(pigments_mean_sort))], method="bray"), factors=chem_mean17$`chem$Lake`) #LCR=!LSA
pigments_nmds <- metaMDS(pigments_mean_sort[,c(5,7:9,11:37,39:42,44:69,71:137,139:ncol(pigments_mean_sort))], distance = "bray", noshare=TRUE,
                         try=40)
pigments_nmds2 <- metaMDS(pigments_mean_sort[,c(5,7:9,11:37,39:42,44:69,71:137,139:ncol(pigments_mean_sort))], distance = "bray", noshare=TRUE,
                          try=40, previous.best=pigments_nmds)
pigments_nmds3 <- metaMDS(pigments_mean_sort[,c(5,7:9,11:37,39:42,44:69,71:137,139:ncol(pigments_mean_sort))], distance = "bray", noshare=TRUE,
                          try=40, previous.best=pigments_nmds2) #stress 0.24

plot(pch=c(21,22,23,24)[chem_mean17$Lake],x=pigments_nmds3$points[,1], y=pigments_nmds3$points[,2],
     bg=c("red","blue")[chem_mean17$Season], cex=2, ylab="NMDS axis 2", 
     xlab=NA, las=1, cex.axis=1.2, cex.lab=1.4)

mtext(outer=F,side=1,text="NMDS axis 1", line=2.3, cex=0.9)
text(x=((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1], y=((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3], labels="G. Pigments - St. 0.24", cex=1.5, adj=0)

# Microscopy biovolumes
permm <- adonis2(micro_biovol_sort~chem_mean17$`chem$Lake`*chem_mean17$Season, permutations = 999, method = "bray") #Lake pvalue=0.001, season pvalue=0.001
aaa <- pairwise.adonis(vegdist(micro_biovol_sort, method="bray"), factors=chem_mean17$`chem$Lake`) #everything is the same
nmds_mb <- metaMDS(micro_biovol_sort, distance = "bray") #stress 0.23

plot(pch=c(21,22,23,24)[chem_mean17$Lake],x=nmds_mb$points[,1], y=nmds_mb$points[,2],
     bg=c("red","blue")[chem_mean17$Season], cex=2, ylab="NMDS axis 2", 
     xlab=NA, las=1, cex.axis=1.2, cex.lab=1.4)

mtext(outer=F,side=1,text="NMDS axis 1", line=2.3, cex=0.9)
text(x=((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1], y=((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3], labels="H. Biovolumes - St. 0.23", cex=1.5, adj=0)

# 16S Chloroplasts
permc <- adonis2(chloro_OTU_sort~chem_mean17$`chem$Lake`*chem_mean17$Season, permutations = 999, method = "bray") #Lake pvalue=0.001, season pvalue=0.001
aaa <- pairwise.adonis(vegdist(chloro_OTU_sort, method="bray"), factors=chem_mean17$`chem$Lake`) #LSC and LCR =! LSA
nmds_chloro <- metaMDS(chloro_OTU_sort, distance = "bray") #stress 0.13

plot(pch=c(21,22,23,24)[chem_mean17$Lake],x=nmds_chloro$points[,1], y=-nmds_chloro$points[,2],
     bg=c("red","blue")[chem_mean17$Season], cex=2, ylab="NMDS axis 2", 
     xlab=NA, las=1, cex.axis=1.2, cex.lab=1.4)

mtext(outer=F,side=1,text="NMDS axis 1", line=2.3, cex=0.9)
text(x=((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1], y=((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3], labels="I. Chloro OTU - St. 0.13", cex=1.5, adj=0)

#lake legend
par(mar=c(0,2.5,0,0))
plot(xaxt="n", yaxt="n", x=c(3:8), y=c(3:8), col="white", xlab=NA, ylab=NA, bty="n")
legend(x=3, y=5.5, yjust=0.5, xjust=0, pch=c(21,24,22,23), legend=c("Lake Clair","Lake Saint-Charles","Lake Cl?ment",
                                                                        "Lake Saint-Augustin"), pt.cex=3, cex=2, bty="n")

#season legend
par(mar=c(0,2.5,0,0))
plot(x=c(3:8), y=c(3:8), col="white", xlab=NA, ylab=NA, bty="n", xaxt="n", yaxt="n")
legend(x=3, y=5.5, yjust=0.5, xjust=0, pch=21, legend=c("Ice-cover","Open-water"),
       pt.bg=c("blue","red"), pt.cex=3, cex=2, bty="n")

dev.off()
