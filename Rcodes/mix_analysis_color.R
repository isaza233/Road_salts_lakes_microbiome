######################################## Microscopie ###############################################

## calcul du biovolume
microscopy <- read.table("data/Sommaire2.csv", sep=";",dec=".", head=T)
mean_biovolume <- vector()
total_biovolume <- vector()
total_cell <- vector()

for(i in 7:ncol(microscopy))
{
  mean_biovolume[i-6] <- weighted.mean(microscopy$Volume[-c(which(is.na(microscopy$Volume)==TRUE))],
                                       microscopy[,i][-c(which(is.na(microscopy$Volume)==TRUE))])
  total_cell[i-6] <-  sum(microscopy[,i][-c(which(is.na(microscopy$Volume)==TRUE))])
  a <- vector()
  for (j in 1:nrow(microscopy))
  {
    a[j] <- microscopy[j,i]*microscopy$Volume[j]
  }
  total_biovolume[i-6] <- sum(a, na.rm=TRUE)
}

## calcul du carbone par ?chantillon (biomass)
ID <- c("Unknown_small", "Unknown_big", "Diatom_small", "Diatom_big", "Chlorophytes", "Chrysophytes", "Dinoflagellates", "Ciliates")
loga <- c(-0.583, -0.665, -0.541, -0.933, -1.026, -1.694, -0.353, -0.639)
b <- c(0.860, 0.939, 0.811, 0.881, 1.088, 1.218, 0.864, 0.984)
pgC <- vector()

for(j in 1:73) #Calcul du nombre de pg par cellule avec les ?quations de: Menden-Deuer et Lessard (2000)
{
  if (is.na(microscopy$Volume_code[j])==TRUE){pgC[j] <- NA}
  else {
    pgC[j] <- loga[which(ID==microscopy$Volume_code[j])]+b[which(ID==microscopy$Volume_code[j])]*log10(microscopy$Volume[j])
  }
}

microscopypgCL <- microscopy #Calcul du nombre de pg par litre d'eau
for (i in 7:34)
{
  for(j in 1:73)
  {
    if (is.na(microscopy$Volume_code[j])==TRUE){microscopypgCL[j,i] <- 0}
    else{
      microscopypgCL[j,i] <- microscopy[j,i]*pgC[j]
    }
  }
}

biomass <- vector() ## calcul de la biomasse totale 

for (i in 7:ncol(microscopypgCL))
{
  biomass[i-6] <- sum(microscopypgCL[,i])
}


microscopy <- microscopy[,-c(13,26,34)]
microscopypgCL <- microscopypgCL[,-c(13,26,34)]
biomass <- biomass[-c(7,20,28)]
mean_biovolume <- mean_biovolume[-c(7,20,28)]
total_biovolume <- total_biovolume[-c(7,20,28)]
total_cell <- total_cell[-c(7,20,28)]


## Cr?ation des histogrammes de biovolumes par saison
micro_ice <- microscopy[,c(6,7,8,9,13,14,15,19,20,21,25,26,27)]
micro_open <- microscopy[,c(6,10,11,12,17,17,18,22,23,24,28,29,30,31)]

a <- vector()

for (i in 2:ncol(micro_ice))
{
  for (j in 1:nrow(micro_ice))
  {
    if (is.na(micro_ice$Volume[j])==TRUE){}
    else{a <- c(a,rep(micro_ice$Volume[j],micro_ice[j,i]))}
  }
}

b <- vector()

for (i in 2:ncol(micro_open))
{
  for (j in 1:nrow(micro_open))
  {
    if (is.na(micro_open$Volume[j])==TRUE){}
    else{b <- c(b,rep(micro_open$Volume[j],micro_open[j,i]))}
  }
}

aa <- hist(a, freq=T, breaks=c(0,100,250,1000,35000))
bb <- hist(b, freq=T, breaks=c(0,100,250,1000,35000))

aap <- vector()
bbp <- vector()

for (i in 1:length(aa$counts))
{
  aap[i] <- (aa$counts[i]/sum(aa$counts))*100
  bbp[i] <- (bb$counts[i]/sum(bb$counts))*100
}

######################################## Pigments ###############################################

pigments <- read.table("data/Sommaire_hplc.csv", sep=";",dec=".", head=F, colClasses = "character")

## transform areas in concentrations
ve=2.7
vi=0.04

#filtration volume in Liters from mL
for (i in 6:ncol(pigments))
{
  pigments[7,i] <- as.numeric(pigments[7,i])/1000
}

#actual calculation from area to ug/L
for (i in 6:ncol(pigments))
{
  for (j in 9:nrow(pigments))
  {
    pigments[j,i] <- as.numeric(pigments[j,i])*as.numeric(pigments[j,4])*(ve/vi)*(1/as.numeric(pigments[7,i]))
  }
}

#actual calculation from ug/L to mol/L
for (i in 6:ncol(pigments))
{
  for (j in 9:nrow(pigments))
  {
    if (is.na(pigments[j,5])==TRUE){pigments[j,i]}
    else{pigments[j,i] <- as.numeric(pigments[j,i])/(10^6*as.numeric(pigments[j,5]))}
  }
}

## Calcul de diff?rentes variables 
PPC <- vector() ## Calculate Photoprotective carotenoids (PPC)
PSC <- vector() ## Calculate Photosynthetic carotenoids (PSC)
TChla <- vector() ## Calculate total chlorophyll a as sum of Chla+allomers of chla (TChla)
Chlacc <- vector() ## Calculate Accessory chlorophylls (Chl acc)
TKCaro <- vector() ## Calculate Total known carotenoids (TKCaro)
TKChl <- vector() ## Calculate Total known chlorophylls (TKChl)
TUCaro <- vector() ## Calculate Total known active carotenoids (TUCaro)
TUChl <- vector() ## Calculate Total known chlorophylls (TUChl)
Chldeg <- vector() ## Calculate Total degraded chlorophylls (Chldeg)

for (i in 6:ncol(pigments))
{
  PPC[i] <- sum(as.numeric(pigments[c(104,111,127,131,137,159,160,130,102),i]))
  PSC[i] <- sum(as.numeric(pigments[c(123,96,86),i]))
  TChla[i] <- sum(as.numeric(pigments[c(66,63,62,61,60,58,33),i]))
  Chlacc[i] <- sum(as.numeric(pigments[c(56,57,26,25,17,23,14,27,24),i]))
  TKCaro[i] <- sum(as.numeric(pigments[c(100,123,110,117,98,104,166,165,137,156,111,127,159,96,131,106,108,86,101,102,130),i]))
  TUCaro[i] <- sum(as.numeric(pigments[c(setdiff(which(pigments$V3=="pda"),c(100,123,110,117,98,104,166,165,137,156,111,127,159,96,131,106,108,86,101,102,130))),i])) #setdiff is for knowing which elements of vector 1 are not in vector 2 setdiff(1,2)
  TKChl[i] <- sum(as.numeric(pigments[c(33,38,58,60,61,62,63,66,56,57,26,25,17,23,14,27,24),i]))
  TUChl[i] <- sum(as.numeric(pigments[c(setdiff(which(pigments$V3=="fluo"),c(33,38,58,60,61,62,63,66,56,57,26,25,17,23,14,27,24))),i]))
  Chldeg[i] <- sum(as.numeric(pigments[c(34,35,36,41,69,70,72,73),i]))
}

PPC[2] <- "PPC"
PSC[2] <- "PSC"
Chlacc[2] <- "Chlacc"
TChla[2] <- "TChla"
TKCaro[2] <- "TKCaro"
TUCaro[2] <- "TUCaro"
TUChl[2] <- "TUChl"
TKChl[2] <- "TKChl"
Chldeg[2] <- "Chldeg"
pigments <- rbind(pigments, "PPC" = PPC)
pigments <- rbind(pigments, "PSC" = PSC)
pigments <- rbind(pigments, "Chlacc" = Chlacc)
pigments <- rbind(pigments, "TChla" = TChla)
pigments <- rbind(pigments, "TKCaro" = TKCaro)
pigments <- rbind(pigments, "TUCaro" = TUCaro)
pigments <- rbind(pigments, "TKChl" = TKChl)
pigments <- rbind(pigments, "TUChl" = TUChl)
pigments <- rbind(pigments, "Chldeg" = Chldeg)


######################################## Pigments preparation ###############################################

## Calculate mean per sample 
pigments2 <- t(pigments[c(1,5,8:nrow(pigments)),c(6:ncol(pigments))])
pigments2 <- as.data.frame(pigments2)
colnames(pigments2) <- c("Lake","Season","ID",pigments[c(9:nrow(pigments)),2])

for (j in 4:ncol(pigments2)) {pigments2[,j] <- as.numeric(as.character(pigments2[,j]))}

pigments_mean <- aggregate(pigments2[,c(4:ncol(pigments2))], by=list(pigments2$ID), FUN=mean)

pigments_mean2 <- as.data.frame(t(pigments_mean))

colnames(pigments_mean2) <- pigments_mean[,1]
pigments_mean2 <- pigments_mean2[-1,]

## rearrange pigments_mean to fit microscopy order (parce qu'il y a moins d'?chantillons dans microscopie
## puisqu'il y a seulement l'ann?e 2017)

pigments_sort <- pigments_mean2[,match(colnames(microscopy)[7:34],colnames(pigments_mean2))[c(which(is.na(match(colnames(microscopy)[7:34],colnames(pigments_mean2)))==FALSE))]]

pigments_sort <- as.data.frame(t(pigments_sort))

######################################## Chemistry ###############################################

chem <- read.table("data/chem.csv", sep=";",dec=".", head=T)

chem_mean <- aggregate(.~ chem$ID, data=chem, FUN="mean")

## rearrange chem_mean to fit microscopy order (parce qu'il y a moins d'?chantillons dans microscopie
## puisqu'il y a seulement l'ann?e 2017)
rownames(chem_mean) <- chem_mean[,1]
chem_mean <- chem_mean[,-c(1,2)]

chem_sort <- chem_mean[match(colnames(microscopy)[7:31],rownames(chem_mean)),]

######################################## Mix analyses ###############################################

lakes <- c("LCT","LCT","LCT","LCT","LCT","LCT","LCL","LCL","LCL","LCL","LCL","LCL","LSA","LSA","LSA","LSA"
           ,"LSA","LSA","LSC","LSC","LSC","LSC","LSC","LSC","LSC")
seasons <- c("w","w","w","s","s","s","w","w","w","s","s","s","w","w","w","s","s","s","w","w","w",
             "s","s","s","s")

tiff(filename="images/Figure3_Fournier.tif", width=20, height=15, units="cm", pointsize=10, res=300)
graphics::layout(mat=matrix(c(4,1,2,5,3,6,5,3,6), nrow=3, ncol=3, byrow=TRUE), 
                 heights=c(1,0.5,0.5))

## 1. Total caro x total chl 
f <- (as.numeric(as.character(pigments_sort$TKCaro))+as.numeric(as.character(pigments_sort$TUCaro)))*10^9
g <- (as.numeric(as.character(pigments_sort$TUChl))+as.numeric(as.character(pigments_sort$TKChl)))*10^9
c <- lm(log10(f)~log10(g))

par(mar=c(4,4.8,0.8,0.5))
plot(y=f, x=g, pch=c(22,21,23,24)[as.factor(lakes)],
     bg=c("red","blue")[as.factor(seasons)], cex=2.5, cex.lab=1.4, cex.axis=1.2, las=1, log="xy", 
     ylab=expression("Total carotenoids"~(nmol~L^-1)), 
     xlab=expression("Total chlorophylls"~(nmol~L^-1)), xaxt="n")
axis(side=1, at=c(0.5, 1, 5, 10, 50), labels=c(0.5, 1, 5, 10, 50), cex.axis=1.2)
abline(c, col="red")
text(x=30, y=1.8, adj=0.5, labels=expression(R^2==0.95))
text(x=30, y=3, adj=0.5, labels=expression("y~0.98x-0.06"))
text(x=10^(((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1]), y=10^(((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3]), labels="B.", cex=1.5, adj=0)

## 2. Total P x Total chl
par(mar=c(4,4.8,0.8,0.5))
f=as.numeric(as.character(pigments_sort$TKChl))*10^9+as.numeric(as.character(pigments_sort$TUChl))*10^9
g <- chem_sort$TP
c <- lm(log10(f)~log10(g))

plot(y=f, x=chem_sort$TP, 
     pch=c(22,21,23,24)[as.factor(lakes)],
     bg=c("red","blue")[as.factor(seasons)], cex=2, las=1, log="xy", 
     ylab=expression("Total chlorophylls"~(nmol~L^-1)), 
     xlab=expression("Total phosphorus"~(mu*g~L^-1)), yaxt="n", cex.lab=1.4, cex.axis=1.2)
axis(side=2, at=c(0.2, 0.5, 1, 5, 10, 50,100), labels=c(0.2, 0.5, 1, 5, 10, 50,100), las=1, cex.axis=1.2)
abline(c, col="red")
text(x=60, y=1.8, adj=0.5, labels=expression(R^2==0.46))
text(x=60, y=3, adj=0.5, labels=expression("y~1.04x-0.06"))
text(x=10^(((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1]), y=10^(((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3]), labels="C.", cex=1.5, adj=0)

## 3. Total chlorophylls/total biovolume x total biovolume
f <- ((as.numeric(as.character(pigments_sort$TKChl))*10^9+as.numeric(as.character(pigments_sort$TUChl))*10^9)/total_biovolume)*10^9
g <- total_biovolume/10^9
c <- lm(log10(f)~log10(g))

par(mar=c(4,4.8,0.8,0.5))
plot(y=f, x=g, pch=c(22,21,23,24)[as.factor(lakes)],
     bg=c("red","blue")[as.factor(seasons)], cex=2, las=1, log="xy", 
     ylab=expression("Total chl/total biovolume"~(10^9~nmol~mm^3^-1)),
     xlab=expression("Total biovolume"~(mm^3~L^-1)), xaxt="n", cex.lab=1.4, cex.axis=1.2)
axis(side=1, at=c(0.05, 0.1, 0.5, 1, 5, 10), labels=c(0.05, 0.1, 0.5, 1, 5, 10), cex.axis=1.2)
abline(c, col="red")
text(x=0.15, y=3.5, adj=0.5, labels=expression(R^2==0.35))
text(x=0.15, y=4.8, adj=0.5, labels=expression("y~-0.51x+0.67"))
text(x=10^(((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1]), y=10^(((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3]), labels="E.", cex=1.5, adj=0)

## 4. Total cells x total biovolume
f <- total_cell/10^6
g <- total_biovolume/10^9
c <- lm(log10(f)~log10(g))

par(mar=c(4,4.8,0.8,0.5))
plot(y=total_cell/10^6, x=total_biovolume/10^9, pch=c(22,21,23,24)[as.factor(lakes)],
     bg=c("red","blue")[as.factor(seasons)], cex=2, las=1, log="xy", 
     ylab=expression("Cells concentration"~~(10^6~cells~L^-1)),
     xlab=expression("Total biovolume"~(mm^3~L^-1)), xaxt="n", cex.lab=1.4, cex.axis=1.2)
axis(side=1, at=c(0.05, 0.1, 0.5, 1, 5, 10), labels=c(0.05, 0.1, 0.5, 1, 5, 10), cex.axis=1.2)
abline(c, col="red")
text(x=4.8, y=0.7, adj=0.5, labels=expression(R^2==0.84))
text(x=4.8, y=1, adj=0.5, labels=expression("y~0.8x+0.23"))
text(x=10^(((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1]), y=10^(((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3]), labels="A.", cex=1.5, adj=0)

## 5. graphique de la taille des cellules
barplot(t(matrix(c(aap,bbp), nrow=4, ncol=2)), beside=TRUE, las=1, ylim=c(0,60), 
        col=c("blue", "red"), ylab="Proportion of cells (%)", 
        names.arg=c("0-100", "100-250", "250-1000", ">1000"),
        xlab=expression("Cell volume"~(mu*m^3)), cex.lab=1.4, cex.axis=1.2)
lines(y=rep(0,12), x=c(1:12))
text(adj=0.5, x=c(1.5,2.5,4.5,5.5,7.5,8.5,10.5,11.5), y=c(31.89,17.06,56.66,39.27,9.58,37.24,9.86,14.43), 
     labels=c("29","15","56","37","7","35","8","13"), cex=0.9)
text(x=((par("usr")[2]-par("usr")[1])*0.05)+par("usr")[1], y=((par("usr")[4]-par("usr")[3])*0.94)+par("usr")[3], labels="D.", cex=1.5, adj=0)

## 6. legends
#lake legend
par(mar=c(0,2.5,0,0))
plot(xaxt="n", yaxt="n", x=c(3:8), y=c(3:8), col="white", xlab=NA, ylab=NA, bty="n")
legend(x=3, y=5.5, yjust=0.5, xjust=0, pch=c(21,24,22,23,NA,21,21), legend=c("Lake Clair","Lake Saint-Charles","Lake Clément",
      "Lake Saint-Augustin","","Ice-cover","Open-water"), col="black", pt.bg=c("white","white","white","white",NA,"blue","red"), pt.cex=3, cex=2, bty="n")

dev.off()
