library(dunn.test)

chem <- read.table("data/chem.csv", sep=";",dec=".", head=T)

chem17 <- subset(chem, subset=chem$Year=="2017", select=colnames(chem))

chem_mean <- aggregate(.~ chem17$ID+chem17$Lake+chem17$Year+chem17$Month, data=chem17, FUN="mean")
rownames(chem_mean) <- chem_mean$`chem17$ID`

## correlation heatmap 
library(corrplot)
a <- cor(chem_mean[,c(12,14:ncol(chem_mean))])
corrplot.mixed(a,lower.col="black", number.cex=0.7, sig.level=0.05, order="hclust", hclust.method="ward.D")

## NMDS on chem alone 
library(vegan)
##
xx <- metaMDS(chem_mean[,c(12,15:33)], distance="bray")
bbb <- adonis2(chem_mean[,c(12,15:33)]~chem_mean$`chem17$Lake`*chem_mean$Season, dist="Bray")
bbb #lakes p=0.001 season NS
rownames(xx$species)[3] <- "Alkalinity"
rownames(xx$species)[2] <- "CDOM"

# figure 
tiff(filename="images/Figure2_Fournier.tif", width=8, height=7, units="cm", pointsize=10, res=300)
par(mar=c(3, 3.5, 0, 0))

plot(scores(xx, display = "sites")[,1], scores(xx, display = "sites")[,2], col="gray12", pch=c(21,22,23,24)[chem_mean$`chem17$Lake`], bg=c("red","blue")[chem_mean$Season],
     ylim=c(-0.35, 0.3), xlim=c(-0.3, 0.4), cex=1.2, xlab=NA, ylab=NA, las=1, yaxt="n", xaxt="n")
axis(side=2, at=c(-0.2,0,0.2), labels=NA, las=1, cex.axis=0.5, tick=TRUE)
mtext(side=2, at=c(-0.2,0,0.2), text=c("-0.2","0.0","0.2"), las=1, cex=0.8, line=1)
axis(side=1, at=c(-0.2,0,0.2,0.4), labels=NA, las=1, cex.axis=0.5, tick=TRUE)
mtext(side=1, at=c(-0.2,0,0.2,0.4), text=c("-0.2","0.0","0.2","0.4"), las=1, cex=0.8, line=0.5)

mtext(side=1, at=0.05, adj=0.5, text=c("NMDS Axis 1"), las=1, cex=0.7, line=1.5)
mtext(side=2, at=0, adj=0.5, text=c("NMDS Axis 2") , cex=0.7, line=2.2)

#points(x=xx$species[c(2,3,9,11,13,19),1], y=xx$species[c(2,3,9,11,13,19),2],col="blue")
text(x=xx$species[c(2,3,9,10,11,13,14,17,19),1], y=xx$species[c(2,3,9,10,11,13,14,17,19),2],col="black",
     labels=rownames(xx$species)[c(2,3,9,10,11,13,14,17,19)])

for(i in c(2,3,9,10,11,13,14,17,19))
{
lines(c(scores(xx, display = "species")[i,1],0), c(scores(xx, display = "species")[i,2],0), col="black",
      lwd=0.6)
}

legend(x=0.33, y=-0.342, yjust=0.5, xjust=0.5, pch=21, legend=c("Ice-cover","Open-water"),
       pt.bg=c("blue","red"), pt.cex=1.2, cex=0.6, bty="n")

legend(x=-0.03, y=-0.35, yjust=0.5, xjust=0.5, pch=c(21,24,22,23), legend=c("LCR","LSC", "LCL","LSA"), pt.cex=1.2, cex=0.6, bty="n", horiz=TRUE,
     text.width=c(0.1, 0.1, 0.1, 0.1))

dev.off()