#https://jokergoo.github.io/ComplexHeatmap-reference/book/legends.html#heatmap-and-annotation-legends


## RUN DESeq_OTU_euk2017.R first 
source("Rcodes/DESeq_OTU_euk2017.R")

### results data.frame at OTU level ### 
euk_p <- prop.table(as.matrix(t(euk_sort1)),1)
euk_p_chem <- data.frame(cond, chem_mean$Lake, chem_mean$Season, euk_p)
colnames(euk_p_chem) <- c("cond", "lake", "season", rownames(euk_sort))

euk_agg <- aggregate(euk_p_chem[,4:ncol(euk_p_chem)], by=list(euk_p_chem$season, euk_p_chem$lake), FUN="mean")
euk_aggcc <- euk_agg

ab_IC <- vector()
for (i in 1:nrow(seasonres))
{
  ab_IC[i] <- mean(euk_agg[c(2,4,6,8),which(colnames(euk_agg)==seasonres$OTU[i])])
}

ab_OW <- vector()
for (i in 1:nrow(seasonres))
{
  ab_OW[i] <- mean(euk_agg[c(1,3,5,7),which(colnames(euk_agg)==seasonres$OTU[i])])
}

ab_comp <- vector()
dummy <- vector()
for (i in 1:nrow(seasonres))
{
  if (ab_IC[i] < ab_OW[i]) {ab_comp[i] <- 1}
  else if (ab_OW[i] < ab_IC[i]) {ab_comp[i] <- 2}
  else {ab_comp[i] <- 0}
  dummy[i] <- 1
}

seasonres2 <- data.frame(seasonres, ab_IC, ab_OW, ab_comp,dummy)
seasonres2$log2Foldchange <- as.numeric(as.character(seasonres2$log2Foldchange))
colnames(seasonres2) <- c(colnames(seasonres), "ab_IC", "ab_OW", "ab_comp","dummy")

### results data.frame at species level ### 
euk_genus <- aggregate(prok_abondance_taxo[,10:ncol(prok_abondance_taxo)], by=list(prok_abondance_taxo$taxo_3.1, prok_abondance_taxo$genus.1, 
                                                                                   prok_abondance_taxo$specie.1), FUN="sum")
v <- c(1:nrow(euk_genus))
euk_genus <- data.frame(v, euk_genus)

euk_sort <- euk_genus[,match(chem_mean$`chem17$ID`,colnames(euk_genus))]

euk_sort_p <- prop.table(as.matrix(t(euk_sort)),1)
euk_sort_p_chem <- data.frame(cond, chem_mean$Lake, chem_mean$Season, euk_sort_p)
colnames(euk_sort_p_chem) <- c("cond", "lake", "season", v)

euk_agg <- aggregate(euk_sort_p_chem[,4:ncol(euk_sort_p_chem)], by=list(euk_sort_p_chem$season, euk_sort_p_chem$lake), FUN="mean")

ab_low <- vector()
for (i in 1:nrow(seasonagg))
{
  a <- euk_genus$v[which(euk_genus$Group.2==as.character(seasonagg$Group.2[i])&euk_genus$Group.3==as.character(seasonagg$Group.3[i]))]
  ab_low[i] <- mean(euk_agg[c(2,4,6,8),which(colnames(euk_agg)==a)])
}

ab_high <- vector()
for (i in 1:nrow(seasonagg))
{
  a <- euk_genus$v[which(euk_genus$Group.2==as.character(seasonagg$Group.2[i])&euk_genus$Group.3==as.character(seasonagg$Group.3[i]))]
  ab_high[i] <- mean(euk_agg[c(1,3,5,7),which(colnames(euk_agg)==a)])
}

seasonagg2 <- data.frame(seasonagg, ab_low, ab_high)
seasonagg2$x<- as.numeric(as.character(seasonagg$x))
colnames(seasonagg2) <- c(colnames(seasonagg), "ab_low", "ab_high")

#for ice-cover
high_IC <- subset(seasonres2, subset=seasonres2$log2Foldchange>0, select=colnames(seasonres2))
high_ICagg3 <- aggregate(high_IC$log2Foldchange, by=list(high_IC$Phylum), FUN="sum")
pp <- vector()
for (i in 1:nrow(high_ICagg3)){pp[i] <- (high_ICagg3$x[i]/sum(high_ICagg3$x))*100}
high_res <- data.frame(high_ICagg3, pp)

#for open water
high_OW <- subset(seasonres2, subset=seasonres2$log2Foldchange<0, select=colnames(seasonres2))
high_OWagg3 <- aggregate(high_OW$log2Foldchange, by=list(high_OW$Phylum), FUN="sum")
pp <- vector()
for (i in 1:nrow(high_OWagg3)){pp[i] <- (high_OWagg3$x[i]/sum(high_OWagg3$x))*100}
high_res <- data.frame(high_OWagg3, pp)

#Phylum for Ic OW L2FC
#("Cercozoa", "Chlorophytes", "Choanoflagellates", "Chrysophytes", "Ciliates", "Cryptophytes", "Diatoms", "Dinoflagellates", "Euglenophytes", 
# "Haptophytes", "Others")
a <- rev(c(11.4,6.4,0,34.2,13.2,16.6,0.7,2.4,0,0,15.1)) #IC
b <- rev(c(5.9,2.58,0.88,32.9,45.6,2.07,0,3.19,0,1.09,5.79)) #OW
col=rev(c("purple","dark green","yellow", "green","aquamarine2","navy","darkseagreen","orange","tan4","red","gray65"))

png(file="images/Rplot002.png", width=5, height=6, units="cm", pointsize=10, res=500)
par(mar=c(2,2,2,2))
barplot(as.matrix(data.frame(a,b)), col=col, las=1, names.arg=c("", ""),
        cex.axis=0.6, cex.names=0.6)
mtext(side=1, line=0, text="IC", at=0.7, cex=0.6)
mtext(side=1, line=0, text="OW", at=1.917, cex=0.6)
dev.off()

### heatmap ### 
library(graphicsutils)

## select all SIGNIFICANT DESeq
chosenOTU <- vector()
chosenPhylum <- vector()

for (i in 1:nrow(seasonres))
{
  chosenOTU[i] <- which(colnames(euk_aggcc)==seasonres[i,1])
  if (seasonres$Phylum[i]=="Bicoecea"){chosenPhylum[i] <- "Others"}
  else if (seasonres$Phylum[i]=="Bolidophyceae"){chosenPhylum[i] <- "Others"}
  else if (seasonres$Phylum[i]==" Bacillariophyceae"){chosenPhylum[i] <- "Diatoms"}
  else if (seasonres$Phylum[i]=="Centroheliozoa"){chosenPhylum[i] <- "Others"}
  else if (seasonres$Phylum[i]=="Alveolata"){chosenPhylum[i] <- "Others"}
  else if (seasonres$Phylum[i]=="Apicomplexa"){chosenPhylum[i] <- "Others"}
  else if (seasonres$Phylum[i]=="Dictyochophyceae"){chosenPhylum[i] <- "Others"}
  else if (seasonres$Phylum[i]=="Fungi"){chosenPhylum[i] <- "Others"}
  else if (seasonres$Phylum[i]=="Ochrophyta_uncl"){chosenPhylum[i] <- "Others"}
  else if (seasonres$Phylum[i]=="Stramenopiles_uncl"){chosenPhylum[i] <- "Others"}
  else if (seasonres$Phylum[i]=="Eukaryota_uncl"){chosenPhylum[i] <- "Others"}
  else if (seasonres$Phylum[i]=="Telonemia"){chosenPhylum[i] <- "Others"}
  else if (seasonres$Phylum[i]=="Perkinsea"){chosenPhylum[i] <- "Others"}
  else if (seasonres$Phylum[i]=="Katablepharidophyta"){chosenPhylum[i] <- "Others"}
  else if (seasonres$Phylum[i]=="Unclassified"){chosenPhylum[i] <- "Others"}
  else {chosenPhylum[i] <- as.character(seasonres$Phylum[i])}
}

## colors from phylum 
love <- as.matrix(euk_aggcc[c(1,2,7,8,3,4,5,6),c(chosenOTU)])
rownames(love) <- c("LCR-OW","LCR-IC","LSC-OW","LSC-IC","LSA-OW","LSA-IC","LCL-OW","LCL-IC")

################## Figure ##################
library(ComplexHeatmap)
library(vegan)
library(circlize)
library(dendsort)
col_fun = colorRamp2(c(0, 0.005, 0.01, 0.05,0.1), c("bisque", "orange","chocolate","dark red","black"))

seasons <- rowAnnotation(Period=c("open","ice","open","ice","open","ice","open","ice"),
                         col=list(Period=c("ice"="blue", "open"="red")), 
                         border=TRUE, annotation_legend_param=list(border="black"),
                         show_annotation_name=FALSE, show_legend=FALSE)

phyl <- HeatmapAnnotation(Phylum=chosenPhylum, col=list(Phylum=c("Cercozoa"="purple",
                                                                 "Chlorophyta"="dark green",
                                                                 "Choanoflagellida"="yellow",
                                                                 "Chrysophyta"="green",
                                                                 "Ciliophora"="aquamarine2",
                                                                 "Cryptophyta"="navy",
                                                                 "Diatoms"="darkseagreen",
                                                                 "Dinoflagellata"="orange",
                                                                 "Euglenoids"="tan4",
                                                                 "Haptophyta"="dark red",
                                                                 "Others"="gray65")),
                          show_annotation_name=FALSE, show_legend=FALSE)
reads <- Legend(title="Reads (%)", labels=rev(c("0","0.5","1","5","10")), at=c(0,0.005,0.01,0.05,0.1),
                border="gray12",legend_gp=gpar(fill=col_fun(rev(c(0,0.005,0.01,0.05,0.1)))))
seasons_leg <- Legend(title="Period", labels=c("Ice-cover","Open-water"), at=c(1,2),
                      border="gray12",legend_gp=gpar(fill=c("blue","red")))
phyl_leg <- Legend(title="Phylum", labels=c("Cercozoa","Chlorophyta","Choanoflagellida", "Chrysophyta","Ciliophora",
                                            "Cryptophyta","Diatoms","Dinoflagellata","Euglenoids","Haptophta","Others"),
                   at=c(1:8),border="gray12", legend_gp=gpar(fill=c("purple","dark green","yellow", "green",
                                                                    "aquamarine2","navy","darkseagreen","orange","tan4","dark red","gray65")))

dendcol = dendsort(hclust(vegdist(t(love)), method="ward.D2"))
dendrow = dendsort(hclust(vegdist(love), method="ward.D2"))

ht <-Heatmap(love, column_dend_side="top",column_labels=rep("", ncol(love)), col=col_fun,
             row_title=NULL, column_title=NULL, row_split=2, column_split=3,
             show_heatmap_legend=FALSE, top_annotation=phyl, left_annotation=seasons,
             cluster_rows=dendrow, cluster_columns=dendcol)

Objet <- packLegend(phyl_leg,reads,seasons_leg,row_gap = unit(0.6,"cm"))

png(file="images/18Sseas.png", width=17, height=12, units="cm", pointsize=11, res=500)

draw(ht, annotation_legend_list=Objet)

dev.off()