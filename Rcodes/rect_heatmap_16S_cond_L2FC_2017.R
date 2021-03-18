#https://jokergoo.github.io/ComplexHeatmap-reference/book/legends.html#heatmap-and-annotation-legends


## RUN DESeq_OTU_prok.R first 
source("Rcodes/DESeq_OTU_prok2017.R")

### results data.frame at OTU level ### 
euk_p <- prop.table(as.matrix(t(euk_sort1)),1)
euk_p_chem <- data.frame(cond, chem_mean$`chem17$Lake`, chem_mean$Season, euk_p)
colnames(euk_p_chem) <- c("cond", "lake", "season", rownames(euk_sort))

euk_agg <- aggregate(euk_p_chem[,4:ncol(euk_p_chem)], by=list(euk_p_chem$season, euk_p_chem$lake), FUN="mean")
euk_aggaa <- euk_agg

ab_IC <- vector()
for (i in 1:nrow(condres))
{
  ab_IC[i] <- mean(euk_agg[c(1,2,7,8),which(colnames(euk_agg)==condres$OTU[i])])
}

ab_OW <- vector()
for (i in 1:nrow(condres))
{
  ab_OW[i] <- mean(euk_agg[c(2,4,5,6),which(colnames(euk_agg)==condres$OTU[i])])
}

ab_comp <- vector()
dummy <- vector()
for (i in 1:nrow(condres))
{
  if (ab_IC[i] < ab_OW[i]) {ab_comp[i] <- 1}
  else if (ab_OW[i] < ab_IC[i]) {ab_comp[i] <- 2}
  else {ab_comp[i] <- 0}
  dummy[i] <- 1
}

condres2 <- data.frame(condres, ab_IC, ab_OW, ab_comp,dummy)
condres2$log2Foldchange <- as.numeric(as.character(condres2$log2Foldchange))
colnames(condres2) <- c(colnames(condres), "ab_IC", "ab_OW", "ab_comp","dummy")

### results data.frame at species level ### 
euk_genus <- aggregate(prok_abondance_taxo[,10:ncol(prok_abondance_taxo)], by=list(prok_abondance_taxo$taxo3, prok_abondance_taxo$taxo6, 
                                                                                   prok_abondance_taxo$taxo7), FUN="sum")
v <- c(1:nrow(euk_genus))
euk_genus <- data.frame(v, euk_genus)

euk_sort <- euk_genus[,match(chem_mean$`chem17$ID`,colnames(euk_genus))]

euk_sort_p <- prop.table(as.matrix(t(euk_sort)),1)
euk_sort_p_chem <- data.frame(cond, chem_mean$Lake, chem_mean$Season, euk_sort_p)
colnames(euk_sort_p_chem) <- c("cond", "lake", "season", v)

euk_agg <- aggregate(euk_sort_p_chem[,4:ncol(euk_sort_p_chem)], by=list(euk_sort_p_chem$season, euk_sort_p_chem$lake), FUN="mean")

ab_low <- vector()
for (i in 1:nrow(condagg))
{
  a <- euk_genus$v[which(euk_genus$Group.1==as.character(condagg$Group.1[i])&euk_genus$Group.2==as.character(condagg$Group.2[i])&euk_genus$Group.3==as.character(condagg$Group.3[i]))]
  ab_low[i] <- mean(euk_agg[c(1,2,7,8),which(colnames(euk_agg)==a)])
}

ab_high <- vector()
for (i in 1:nrow(condagg))
{
  a <- euk_genus$v[which(euk_genus$Group.1==as.character(condagg$Group.1[i])&euk_genus$Group.2==as.character(condagg$Group.2[i])&euk_genus$Group.3==as.character(condagg$Group.3[i]))]
  ab_high[i] <- mean(euk_agg[c(3,4,5,6),which(colnames(euk_agg)==a)])
}

condagg2 <- data.frame(condagg, ab_low, ab_high)
condagg2$x<- as.numeric(as.character(condagg$x))
colnames(condagg2) <- c(colnames(condagg), "ab_low", "ab_high")

#for low cond
high_IC <- subset(condres2, subset=condres2$log2Foldchange>0, select=colnames(condres2))
high_ICagg3 <- aggregate(high_IC$log2Foldchange, by=list(high_IC$Phylum, high_IC$Class), FUN="sum")
pp <- vector()
for (i in 1:nrow(high_ICagg3)){pp[i] <- (high_ICagg3$x[i]/sum(high_ICagg3$x))*100}
high_res <- data.frame(high_ICagg3, pp)

#for high cond
high_OW <- subset(condres2, subset=condres2$log2Foldchange<0, select=colnames(condres2))
high_OWagg3 <- aggregate(high_OW$log2Foldchange, by=list(high_OW$Phylum, high_OW$Class), FUN="sum")
pp <- vector()
for (i in 1:nrow(high_OWagg3)){pp[i] <- (high_OWagg3$x[i]/sum(high_OWagg3$x))*100}
high_res <- data.frame(high_OWagg3, pp)

#Phylum for low high L2FC
#alphaproteo, Bacteroidetes, betaproteo, cyano, deltaproteo, gammaproteo, Planctomycetes, verrucomicrobia, others
a <- c(23.2,18.3,11.4,1.2,1.5,2.7,5.6,24.2,11.9) #low
b <- c(12.5,25,9.9,5.96,6.7,8.5,12.3,5.9,13.24) #high

col <- rev(c("purple","dark green","aquamarine2","navy", "light blue", "orange", "red","black","gray65"))
png(file="images/Rplot002.png", width=5, height=6, units="cm", pointsize=10, res=500)
par(mar=c(2,2,2,2))
barplot(as.matrix(data.frame(a,b)), col=col, las=1, names.arg=c("", ""),
        cex.axis=0.6, cex.names=0.6)
mtext(side=1, line=0, text="Low", at=0.7, cex=0.6)
mtext(side=1, line=0, text="High", at=1.917, cex=0.6)
dev.off()

### heatmap ### 
library(graphicsutils)

## select all SIGNIFICANT DESeq
chosenOTU <- vector()
chosenPhylum <- vector()

for (i in 1:nrow(condres))
{
  chosenOTU[i] <- which(colnames(euk_aggaa)==condres[i,1])
  if (condres$Phylum[i]=="Spirochaetes"){chosenPhylum[i] <- "Others"}
  else if (condres$Phylum[i]=="Acidobacteria"){chosenPhylum[i] <- "Others"}
  else if (condres$Phylum[i]=="Armatimonadetes"){chosenPhylum[i] <- "Others"}
  else if (condres$Phylum[i]=="Actinobacteria"){chosenPhylum[i] <- "Others"}
  else if (condres$Phylum[i]=="Gemmatimonadetes"){chosenPhylum[i] <- "Others"}
  else if (condres$Phylum[i]=="Margulisbacteria"){chosenPhylum[i] <- "Others"}
  else if (condres$Phylum[i]=="Dependentiae"){chosenPhylum[i] <- "Others"}
  else if (condres$Phylum[i]=="Chloroflexi"){chosenPhylum[i] <- "Others"}
  else if (condres$Phylum[i]=="Chlamydiae"){chosenPhylum[i] <- "Others"}
  else if (condres$Phylum[i]=="Unclassified"){chosenPhylum[i] <- "Others"}
  else if (condres$Phylum[i]=="Firmicutes"){chosenPhylum[i] <- "Others"}
  else if (condres$Phylum[i]=="Epsilonbacteraeota"){chosenPhylum[i] <- "Others"}
  else {chosenPhylum[i] <- as.character(condres$Phylum[i])}
}

for (i in 1:length(chosenPhylum)) ## this cause "unclassified" to appear in the vector
{
  if (chosenPhylum[i]=="Proteobacteria") {chosenPhylum[i] <- as.character(condres$Class[i])}
  else {chosenPhylum[i] <- chosenPhylum[i] }
}

## colors from phylum 
love <- as.matrix(euk_aggaa[c(1,2,7,8,3,4,5,6),c(chosenOTU)])
rownames(love) <- c("LCR-OW","LCR-IC","LSC-OW","LSC-IC","LSA-OW","LSA-IC","LCL-OW","LCL-IC")

################## Figure ##################
library(ComplexHeatmap)
library(vegan)
library(circlize)
library(dendsort)
col_fun = colorRamp2(c(0, 0.005, 0.01, 0.05,0.1), c("bisque", "orange","chocolate","dark red","black"))

seasons <- rowAnnotation(Period=c("low","low","low","low","high","high","high","high"),
                         col=list(Period=c("low"="white", "high"="gray12")), 
                         border=TRUE, annotation_legend_param=list(border="black"),
                         show_annotation_name=FALSE, show_legend=FALSE)

phyl <- HeatmapAnnotation(Phylum=chosenPhylum, col=list(Phylum=c("Alphaproteobacteria"="purple",
                                                                 "Bacteroidetes"="dark green",
                                                                 "Betaproteobacteria"="aquamarine2",
                                                                 "Cyanobacteria"="navy",
                                                                 "Deltaproteobacteria"="light blue",
                                                                 "Gammaproteobacteria"="orange",
                                                                 "Unclassified"="gray65", ## this is the unclassified proteobacteria
                                                                 "Others"="gray65",
                                                                 "Planctomycetes"="red",
                                                                 "Verrucomicrobia"="black")),
                          show_annotation_name=FALSE, show_legend=FALSE)

reads <- Legend(title="Reads (%)", labels=rev(c("0","0.5","1","5","10")), at=c(0,0.005,0.01,0.05,0.1),
                border="gray12",legend_gp=gpar(fill=col_fun(rev(c(0,0.005,0.01,0.05,0.1)))))
seasons_leg <- Legend(title="Salinization", labels=c("Low","High"), at=c(1,2),
                      border="gray12",legend_gp=gpar(fill=c("white","gray12")))
phyl_leg <- Legend(title="Phylum", labels=c("Alphaproteobacteria","Bacteroidetes",
                                            "Betaproteobacteria","Cyanobacteria","Deltaproteobacteria","Gammaproteobacteria",
                                            "Planctomycetes","Verrucomicrobia","Others"),
                   at=c(1:8),border="gray12", legend_gp=gpar(fill=c("purple", "dark green", "aquamarine2","navy","light blue",
                                                                    "orange", "red", "black", "gray65")))
dendcol = dendsort(hclust(vegdist(t(love)), method="ward.D2"))
dendrow = dendsort(hclust(vegdist(love), method="ward.D2"))

ht <-Heatmap(love, column_dend_side="top",column_labels=rep("", ncol(love)), col=col_fun,
             row_title=NULL, column_title=NULL, row_split=2, column_split=3,
             show_heatmap_legend=FALSE, top_annotation=phyl, left_annotation=seasons,
             cluster_rows=dendrow, cluster_columns=dendcol)

Objet <- packLegend(phyl_leg,reads,seasons_leg,row_gap = unit(0.6,"cm"))

png(file="images/16Scond.png", width=17, height=12, units="cm", pointsize=11, res=500)

draw(ht, annotation_legend_list=Objet)

dev.off()