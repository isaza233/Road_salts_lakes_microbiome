chem <- read.table("data/chem.csv", sep=";",dec=".", head=T)

## calculer le score sur l'axe 1 de chaque site de 2017 sur la NMDS
chem17 <- subset(chem, subset=chem$Year=="2017", select=colnames(chem))
chem_mean <- aggregate(.~ chem17$ID+chem17$Lake+chem17$Year+chem17$Month, data=chem17, FUN="mean")
rownames(chem_mean) <- chem_mean$`chem17$ID`
library(vegan)
xx <- metaMDS(chem_mean[,c(12,15:33)], distance="bray")

## ajouter le score dans le data.frame chem_mean et subdiviser par p?riode
chem1 <- data.frame(scores(xx, display = "sites")[,1], chem_mean)
colnames(chem1)[1] <- "scores"
chem_meanIC <- subset(chem1, subset=chem1$Season=="2", select=colnames(chem1)) ## winter
chem_meanOW <- subset(chem1, subset=chem1$Season=="1", select=colnames(chem1)) ## summer

################################### For phototrophs only
######################### Microscopy biovolumes
microscopy <- read.table("data/Sommaire2.csv", sep=";",dec=".", head=T)

for (j in 7:ncol(microscopy))
{
  for (i in 1:nrow(microscopy))
  {
    if (is.na(microscopy$Volume[i])==TRUE)
    {microscopy[i,j] <- microscopy[i,j]}
    else 
    {microscopy[i,j] <- microscopy[i,j]*microscopy$Volume[i]}
  }
}

micro_vol <- microscopy[c(which(is.na(microscopy$Volume)==FALSE)),]

rr <- vector()
for (i in 1:nrow(micro_vol))
{
  if (micro_vol$Family[i]=="Ciliates"){rr <- c(rr, i)}
  else if (micro_vol$Family[i]=="Flagellates"){rr <- c(rr, i)}
  else if (micro_vol$Family[i]=="Unclassified"){rr <- c(rr, i)}
  else if (micro_vol$Family[i]=="Cyanobacteria"){rr <- c(rr, i)}
  else {rr <- rr}
}

micro_vol <- micro_vol[-rr,]
v <- c(1:nrow(micro_vol))
micro_vol1 <- data.frame(v, micro_vol[,c(7:ncol(micro_vol))])

## transposer en sites=row, genus=columns
micro_t <- t(micro_vol1)
micro_t <- as.data.frame(micro_t)
colnames(micro_t) <- micro_vol1[,1]
micro_t <- micro_t[-1,]

## s'assurer que chem_mean, genus_t_surf
micro_sortIC <- micro_t[match(chem_meanIC$chem17.ID,rownames(micro_t)),]
micro_sortOW <- micro_t[match(chem_meanOW$chem17.ID,rownames(micro_t)),]

## RDA
library(vegan)

dd <- list(chem_meanIC, chem_meanOW, micro_sortIC, micro_sortOW)

for (i in 1:2)
{
  chem_mean <- dd[[i]]
  micro_sort <- dd[[i+2]]
  
  p_micro <- vector() #DOC, scores, Si, TP, TN, Lake, Al
  
  a_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+Condition(chem_mean$DOC)))
  p_micro[1] <- (a_micro$partial.chi/a_micro$tot.chi)*100
  
  b_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$DOC+Condition(chem_mean$scores)))
  p_micro[2] <- (b_micro$partial.chi/b_micro$tot.chi)*100
  
  c_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Si)))
  p_micro[3] <- (c_micro$partial.chi/c_micro$tot.chi)*100
  
  d_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TP)))
  p_micro[4] <- (d_micro$partial.chi/d_micro$tot.chi)*100
  
  e_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TN)))
  p_micro[5] <- (e_micro$partial.chi/e_micro$tot.chi)*100
  
  g_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Lake)))
  p_micro[6] <- (g_micro$partial.chi/g_micro$tot.chi)*100
  
  h_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Lake+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Al)))
  p_micro[7] <- (h_micro$partial.chi/h_micro$tot.chi)*100
  
  if (i==1){p_microIC <- p_micro}
  else{p_microOW <- p_micro}
}

barplot(as.matrix(data.frame(p_microIC, p_microOW)), col=rainbow(7))

######################### EUK 16s chloroplasts
chloroplasts <- read.table("data/Chloroplasts.csv", sep=";",dec=".", head=T)

## transposer en sites=row, genus=columns
chloro_t <- t(chloroplasts)
chloro_t <- as.data.frame(chloro_t)
colnames(chloro_t) <- chloroplasts[,1]
chloro_t <- chloro_t[-c(1:4),]

##S?lectionner les stations de surface seulement
chloro_t_surf <- chloro_t[-c(15,18,20,23,25,26,27), ]

## s'assurer que chem_mean, genus_t_surf
chloro_sortIC <- chloro_t_surf[match(chem_meanIC$chem17.ID,rownames(chloro_t_surf)),]
for (i in 1:ncol(chloro_sortIC))
{
  chloro_sortIC[,i] <- as.numeric(as.character(chloro_sortIC[,i]))
}

chloro_sortOW <- chloro_t_surf[match(chem_meanOW$chem17.ID,rownames(chloro_t_surf)),]
for (i in 1:ncol(chloro_sortOW))
{
  chloro_sortOW[,i] <- as.numeric(as.character(chloro_sortOW[,i]))
}

# RDA
dd <- list(chem_meanIC, chem_meanOW, chloro_sortIC, chloro_sortOW)

for (i in 1:2)
{
  chem_mean <- dd[[i]]
  chloro_sort <- dd[[i+2]]
  
  p_chloro <- vector() #DOC, scores, Si, TP, TN, Lake, Al
  
  a_chloro <- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+Condition(chem_mean$DOC)))
  p_chloro[1] <- (a_chloro$partial.chi/a_chloro$tot.chi)*100
  
  b_chloro<- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$DOC+Condition(chem_mean$scores)))
  p_chloro[2] <- (b_chloro$partial.chi/b_chloro$tot.chi)*100
  
  c_chloro <- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Si)))
  p_chloro[3] <- (c_chloro$partial.chi/c_chloro$tot.chi)*100
  
  d_chloro <- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TP)))
  p_chloro[4] <- (d_chloro$partial.chi/d_chloro$tot.chi)*100
  
  e_chloro <- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TN)))
  p_chloro[5] <- (e_chloro $partial.chi/e_chloro $tot.chi)*100
  
  g_chloro <- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Lake)))
  p_chloro[6] <- (g_chloro$partial.chi/g_chloro$tot.chi)*100
  
  h_chloro <- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Lake+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Al)))
  p_chloro[7] <- (h_chloro$partial.chi/h_chloro$tot.chi)*100
  
  if (i==1){p_chloroIC <- p_chloro}
  else{p_chloroOW <- p_chloro}
}

barplot(as.matrix(data.frame(p_chloroIC, p_chloroOW)), col=rainbow(7))

######################### EUK 18S
euk_abondance_taxo <- read.table("data/Euk_samples_PR2_2.csv", sep=";",dec=".", head=T)

ss <- vector()
for (i in 1:nrow(euk_abondance_taxo))
{
  if (euk_abondance_taxo$taxo_3.1[i]=="Alveolata"){ss <- c(ss, i)}
  else if (euk_abondance_taxo$taxo_3.1[i]=="Apicomplexa"){ss <- c(ss, i)}
  else if (euk_abondance_taxo$taxo_3.1[i]=="Bicoecea"){ss <- c(ss, i)}
  else if (euk_abondance_taxo$taxo_3.1[i]=="Bolidophyceae"){ss <- c(ss, i)}
  else if (euk_abondance_taxo$taxo_3.1[i]=="Centroheliozoa"){ss <- c(ss, i)}
  else if (euk_abondance_taxo$taxo_3.1[i]=="Choanoflagellida"){ss <- c(ss, i)}
  else if (euk_abondance_taxo$taxo_3.1[i]=="Ciliophora"){ss <- c(ss, i)}
  else if (euk_abondance_taxo$taxo_3.1[i]=="Dictyochophyceae"){ss <- c(ss, i)}
  else if (euk_abondance_taxo$taxo_3.1[i]=="Eukaryota_uncl"){ss <- c(ss, i)}
  else if (euk_abondance_taxo$taxo_3.1[i]=="Fungi"){ss <- c(ss, i)}
  else if (euk_abondance_taxo$taxo_3.1[i]=="Katablepharidophyta"){ss <- c(ss, i)}
  else if (euk_abondance_taxo$taxo_3.1[i]=="Perkinsea"){ss <- c(ss, i)}
  else if (euk_abondance_taxo$taxo_3.1[i]=="Telonemia"){ss <- c(ss, i)}
  else {ss <- ss}
}

euk_abondance_taxo <- euk_abondance_taxo[-ss,]

## transposer en sites=row, genus=columns
euk_t <- t(euk_abondance_taxo[,10:52])
euk_t <- as.data.frame(euk_t)
colnames(euk_t) <- euk_abondance_taxo[,1]

##S?lectionner les stations de surface seulement
euk_t_surf <- euk_t[-c(15,18,20,23,25,26,27), ]

## s'assurer que chem_mean, genus_t_surf
euk_sortIC <- euk_t_surf[match(chem_meanIC$chem17.ID,rownames(euk_t_surf)),]
euk_sortOW <- euk_t_surf[match(chem_meanOW$chem17.ID,rownames(euk_t_surf)),]

# RDA
dd <- list(chem_meanIC, chem_meanOW, euk_sortIC, euk_sortOW)

for (i in 1:2)
{
  chem_mean <- dd[[i]]
  euk_sort <- dd[[i+2]]
  
  p_euk <- vector() #DOC, scores, Si, TP, TN, Lake, Al
  
  a_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+Condition(chem_mean$DOC)))
  p_euk[1] <- (a_euk$partial.chi/a_euk$tot.chi)*100
  
  b_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$DOC+Condition(chem_mean$scores)))
  p_euk[2] <- (b_euk$partial.chi/b_euk$tot.chi)*100
  
  c_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Si)))
  p_euk[3] <- (c_euk$partial.chi/c_euk$tot.chi)*100
  
  d_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TP)))
  p_euk[4] <- (d_euk$partial.chi/d_euk$tot.chi)*100
  
  e_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TN)))
  p_euk[5] <- (e_euk$partial.chi/e_euk$tot.chi)*100
  
  g_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Lake)))
  p_euk[6] <- (g_euk$partial.chi/g_euk$tot.chi)*100
  
  h_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Lake+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Al)))
  p_euk[7] <- (h_euk$partial.chi/h_euk$tot.chi)*100
  
  if (i==1){p_eukIC <- p_euk}
  else{p_eukOW <- p_euk}
}

barplot(as.matrix(data.frame(p_eukIC, p_eukOW)), col=rainbow(7))

######################### PROK 16S
prok_abondance_taxo <- read.table("data/Prok_abondance_taxo.csv", sep=";",dec=".", head=T)

## transposer en sites=row, genus=columns
prok_t <- t(prok_abondance_taxo[,10:52])
prok_t <- as.data.frame(prok_t)
colnames(prok_t) <- prok_abondance_taxo[,1]

##S?lectionner les stations de surface seulement
prok_t_surf <- prok_t[-c(15,18,20,23,25,26,27), ]

## s'assurer que chem_mean, genus_t_surf
prok_sortIC <- prok_t_surf[match(chem_meanIC$chem17.ID,rownames(prok_t_surf)),]
prok_sortOW <- prok_t_surf[match(chem_meanOW$chem17.ID,rownames(prok_t_surf)),]

# RDA
dd <- list(chem_meanIC, chem_meanOW, prok_sortIC, prok_sortOW)

for (i in 1:2)
{
  chem_mean <- dd[[i]]
  prok_sort <- dd[[i+2]]
  
  p_prok <- vector() #DOC, scores, Si, TP, TN, Lake, Al
  
  a_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+Condition(chem_mean$DOC)))
  p_prok[1] <- (a_prok$partial.chi/a_prok$tot.chi)*100
  
  b_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$DOC+Condition(chem_mean$scores)))
  p_prok[2] <- (b_prok$partial.chi/b_prok$tot.chi)*100
  
  c_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Si)))
  p_prok[3] <- (c_prok$partial.chi/c_prok$tot.chi)*100
  
  d_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TP)))
  p_prok[4] <- (d_prok$partial.chi/d_prok$tot.chi)*100
  
  e_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TN)))
  p_prok[5] <- (e_prok$partial.chi/e_prok$tot.chi)*100
  
  g_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Lake)))
  p_prok[6] <- (g_prok$partial.chi/g_prok$tot.chi)*100
  
  h_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Lake+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Al)))
  p_prok[7] <- (h_prok$partial.chi/h_prok$tot.chi)*100
  
  if (i==1){p_prokIC <- p_prok}
  else{p_prokOW <- p_prok}
}

barplot(as.matrix(data.frame(p_prokIC, p_prokOW)), col=rainbow(7))

#### barplot of variance partitioning
bbb2<- matrix(c(p_microIC,p_microOW,p_chloroIC,p_chloroOW,p_eukIC,p_eukOW), nrow=7, ncol=6)
bbb2 <- bbb2[c(7,2,1,6,3,5,4),]
sums <- c(round(sum(bbb2[,1])),round(sum(bbb2[,2])),round(sum(bbb2[,3])),round(sum(bbb2[,4])),round(sum(bbb2[,5])),round(sum(bbb2[,6])))
bbb <- bbb2/max(sums)*100

###################################### For everything 
######################### Microscopy biovolumes
microscopy <- read.table("data/Sommaire2.csv", sep=";",dec=".", head=T)

for (j in 7:ncol(microscopy))
{
  for (i in 1:nrow(microscopy))
  {
    if (is.na(microscopy$Volume[i])==TRUE)
    {microscopy[i,j] <- microscopy[i,j]}
    else 
    {microscopy[i,j] <- microscopy[i,j]*microscopy$Volume[i]}
  }
}

micro_vol <- microscopy[c(which(is.na(microscopy$Volume)==FALSE)),]

v <- c(1:nrow(micro_vol))
micro_vol1 <- data.frame(v, micro_vol[,c(7:ncol(micro_vol))])

## transposer en sites=row, genus=columns
micro_t <- t(micro_vol1)
micro_t <- as.data.frame(micro_t)
colnames(micro_t) <- micro_vol1[,1]
micro_t <- micro_t[-1,]

## s'assurer que chem_mean, genus_t_surf
micro_sortIC <- micro_t[match(chem_meanIC$chem17.ID,rownames(micro_t)),]
micro_sortOW <- micro_t[match(chem_meanOW$chem17.ID,rownames(micro_t)),]

## RDA
library(vegan)

dd <- list(chem_meanIC, chem_meanOW, micro_sortIC, micro_sortOW)

for (i in 1:2)
{
  chem_mean <- dd[[i]]
  micro_sort <- dd[[i+2]]
  
  p_micro <- vector() #DOC, scores, Si, TP, TN, Lake, Al
  
  a_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+Condition(chem_mean$DOC)))
  p_micro[1] <- (a_micro$partial.chi/a_micro$tot.chi)*100
  
  b_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$DOC+Condition(chem_mean$scores)))
  p_micro[2] <- (b_micro$partial.chi/b_micro$tot.chi)*100
  
  c_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Si)))
  p_micro[3] <- (c_micro$partial.chi/c_micro$tot.chi)*100
  
  d_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TP)))
  p_micro[4] <- (d_micro$partial.chi/d_micro$tot.chi)*100
  
  e_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TN)))
  p_micro[5] <- (e_micro$partial.chi/e_micro$tot.chi)*100
  
  g_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Lake)))
  p_micro[6] <- (g_micro$partial.chi/g_micro$tot.chi)*100
  
  h_micro <- summary(capscale(vegdist(micro_sort, method="bray")~chem_mean$Lake+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Al)))
  p_micro[7] <- (h_micro$partial.chi/h_micro$tot.chi)*100
  
  if (i==1){p_microIC <- p_micro}
  else{p_microOW <- p_micro}
}

barplot(as.matrix(data.frame(p_microIC, p_microOW)), col=rainbow(7))

######################### EUK 16s chloroplasts
chloroplasts <- read.table("data/Chloroplasts.csv", sep=";",dec=".", head=T)

## transposer en sites=row, genus=columns
chloro_t <- t(chloroplasts)
chloro_t <- as.data.frame(chloro_t)
colnames(chloro_t) <- chloroplasts[,1]
chloro_t <- chloro_t[-c(1:4),]

##S?lectionner les stations de surface seulement
chloro_t_surf <- chloro_t[-c(15,18,20,23,25,26,27), ]

## s'assurer que chem_mean, genus_t_surf
chloro_sortIC <- chloro_t_surf[match(chem_meanIC$chem17.ID,rownames(chloro_t_surf)),]
for (i in 1:ncol(chloro_sortIC))
{
  chloro_sortIC[,i] <- as.numeric(as.character(chloro_sortIC[,i]))
}

chloro_sortOW <- chloro_t_surf[match(chem_meanOW$chem17.ID,rownames(chloro_t_surf)),]
for (i in 1:ncol(chloro_sortOW))
{
  chloro_sortOW[,i] <- as.numeric(as.character(chloro_sortOW[,i]))
}

# RDA
dd <- list(chem_meanIC, chem_meanOW, chloro_sortIC, chloro_sortOW)

for (i in 1:2)
{
  chem_mean <- dd[[i]]
  chloro_sort <- dd[[i+2]]
  
  p_chloro <- vector() #DOC, scores, Si, TP, TN, Lake, Al
  
  a_chloro <- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+Condition(chem_mean$DOC)))
  p_chloro[1] <- (a_chloro$partial.chi/a_chloro$tot.chi)*100
  
  b_chloro<- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$DOC+Condition(chem_mean$scores)))
  p_chloro[2] <- (b_chloro$partial.chi/b_chloro$tot.chi)*100
  
  c_chloro <- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Si)))
  p_chloro[3] <- (c_chloro$partial.chi/c_chloro$tot.chi)*100
  
  d_chloro <- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TP)))
  p_chloro[4] <- (d_chloro$partial.chi/d_chloro$tot.chi)*100
  
  e_chloro <- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TN)))
  p_chloro[5] <- (e_chloro $partial.chi/e_chloro $tot.chi)*100
  
  g_chloro <- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Lake)))
  p_chloro[6] <- (g_chloro$partial.chi/g_chloro$tot.chi)*100
  
  h_chloro <- summary(capscale(vegdist(chloro_sort, method="bray")~chem_mean$Lake+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Al)))
  p_chloro[7] <- (h_chloro$partial.chi/h_chloro$tot.chi)*100
  
  if (i==1){p_chloroIC <- p_chloro}
  else{p_chloroOW <- p_chloro}
}

barplot(as.matrix(data.frame(p_chloroIC, p_chloroOW)), col=rainbow(7))

######################### EUK 18S
euk_abondance_taxo <- read.table("data/Euk_samples_PR2_2.csv", sep=";",dec=".", head=T)

## transposer en sites=row, genus=columns
euk_t <- t(euk_abondance_taxo[,10:52])
euk_t <- as.data.frame(euk_t)
colnames(euk_t) <- euk_abondance_taxo[,1]

##S?lectionner les stations de surface seulement
euk_t_surf <- euk_t[-c(15,18,20,23,25,26,27), ]

## s'assurer que chem_mean, genus_t_surf
euk_sortIC <- euk_t_surf[match(chem_meanIC$chem17.ID,rownames(euk_t_surf)),]
euk_sortOW <- euk_t_surf[match(chem_meanOW$chem17.ID,rownames(euk_t_surf)),]

# RDA
dd <- list(chem_meanIC, chem_meanOW, euk_sortIC, euk_sortOW)

for (i in 1:2)
{
  chem_mean <- dd[[i]]
  euk_sort <- dd[[i+2]]
  
  p_euk <- vector() #DOC, scores, Si, TP, TN, Lake, Al
  
  a_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+Condition(chem_mean$DOC)))
  p_euk[1] <- (a_euk$partial.chi/a_euk$tot.chi)*100
  
  b_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$DOC+Condition(chem_mean$scores)))
  p_euk[2] <- (b_euk$partial.chi/b_euk$tot.chi)*100
  
  c_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Si)))
  p_euk[3] <- (c_euk$partial.chi/c_euk$tot.chi)*100
  
  d_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TP)))
  p_euk[4] <- (d_euk$partial.chi/d_euk$tot.chi)*100
  
  e_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TN)))
  p_euk[5] <- (e_euk$partial.chi/e_euk$tot.chi)*100
  
  g_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Lake)))
  p_euk[6] <- (g_euk$partial.chi/g_euk$tot.chi)*100
  
  h_euk <- summary(capscale(vegdist(euk_sort, method="bray")~chem_mean$Lake+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Al)))
  p_euk[7] <- (h_euk$partial.chi/h_euk$tot.chi)*100
  
  if (i==1){p_eukIC <- p_euk}
  else{p_eukOW <- p_euk}
}

barplot(as.matrix(data.frame(p_eukIC, p_eukOW)), col=rainbow(7))

######################### PROK 16S
prok_abondance_taxo <- read.table("data/Prok_abondance_taxo.csv", sep=";",dec=".", head=T)

## transposer en sites=row, genus=columns
prok_t <- t(prok_abondance_taxo[,10:52])
prok_t <- as.data.frame(prok_t)
colnames(prok_t) <- prok_abondance_taxo[,1]

##S?lectionner les stations de surface seulement
prok_t_surf <- prok_t[-c(15,18,20,23,25,26,27), ]

## s'assurer que chem_mean, genus_t_surf
prok_sortIC <- prok_t_surf[match(chem_meanIC$chem17.ID,rownames(prok_t_surf)),]
prok_sortOW <- prok_t_surf[match(chem_meanOW$chem17.ID,rownames(prok_t_surf)),]

# RDA
dd <- list(chem_meanIC, chem_meanOW, prok_sortIC, prok_sortOW)

for (i in 1:2)
{
  chem_mean <- dd[[i]]
  prok_sort <- dd[[i+2]]
  
  p_prok <- vector() #DOC, scores, Si, TP, TN, Lake, Al
  
  a_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+Condition(chem_mean$DOC)))
  p_prok[1] <- (a_prok$partial.chi/a_prok$tot.chi)*100
  
  b_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$DOC+Condition(chem_mean$scores)))
  p_prok[2] <- (b_prok$partial.chi/b_prok$tot.chi)*100
  
  c_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Si)))
  p_prok[3] <- (c_prok$partial.chi/c_prok$tot.chi)*100
  
  d_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TN+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TP)))
  p_prok[4] <- (d_prok$partial.chi/d_prok$tot.chi)*100
  
  e_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Lake+chem_mean$Al+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$TN)))
  p_prok[5] <- (e_prok$partial.chi/e_prok$tot.chi)*100
  
  g_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Al+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Lake)))
  p_prok[6] <- (g_prok$partial.chi/g_prok$tot.chi)*100
  
  h_prok <- summary(capscale(vegdist(prok_sort, method="bray")~chem_mean$Lake+chem_mean$TN+chem_mean$TP+chem_mean$Si+chem_mean$scores+chem_mean$DOC+Condition(chem_mean$Al)))
  p_prok[7] <- (h_prok$partial.chi/h_prok$tot.chi)*100
  
  if (i==1){p_prokIC <- p_prok}
  else{p_prokOW <- p_prok}
}

barplot(as.matrix(data.frame(p_prokIC, p_prokOW)), col=rainbow(7))

#### barplot of variance partitioning
aaa2<- matrix(c(p_microIC,p_microOW,p_chloroIC,p_chloroOW,p_eukIC,p_eukOW,p_prokIC,p_prokOW), nrow=7, ncol=8)
aaa2 <- aaa2[c(7,2,1,6,3,5,4),]
sums <- c(round(sum(aaa2[,1])),round(sum(aaa2[,2])),round(sum(aaa2[,3])),round(sum(aaa2[,4])),round(sum(aaa2[,5])),round(sum(aaa2[,6])),round(sum(aaa2[,7])),round(sum(aaa2[,8])))
aaa <- aaa2/max(sums)*100

## Figure 
png(file="images/partial_RDA.png", width=14, height=18, units="cm", pointsize=10, res=300)
graphics::layout(mat=matrix(c(1,0,2,0,2,5,3,5,4,5,4,0), nrow=6, ncol=2, byrow=TRUE), widths=c(1,0.35), heights=c(0.05,0.5,0.5,0.05,0.5,0.5))

#1
par(mar=c(0,0,0,0))
plot(x=c(1:10), y=c(1:10), ylab=NA, xlab=NA, col="white", xaxt="n", yaxt="n", bty="n")
text(x=1, y=5, adj=0, labels="A. Heterotrophic and phototrophic prokaryotes and eukaryotes", cex=1.3)

#2
place <- c(0.9,1.95,3.35,4.4,5.8,6.85,8.25,9.3)
par(mar=c(3,4,0,0))
barplot(aaa[c(7:1),], col=rainbow(7), las=1, names.arg=c("","","","","","","",""),
        ylab="Relative variance partitioning (%)", ylim=c(0,105), bty="n",space=c(0.4,0.05,0.4,0.05,0.4,0.05,0.4,0.05))
lines(x=c(1:9), y=rep(0,9), col="black")
text(x=place, y=c(round(sum(aaa[,1])),round(sum(aaa[,2])),round(sum(aaa[,3])),round(sum(aaa[,4])),round(sum(aaa[,5])),round(sum(aaa[,6])),round(sum(aaa[,7])),round(sum(aaa[,8])))+3, labels=c(round(sum(aaa[,1])),round(sum(aaa[,2])),round(sum(aaa[,3])),round(sum(aaa[,4])),round(sum(aaa[,5])),round(sum(aaa[,6])),round(sum(aaa[,7])),round(sum(aaa[,8]))), cex=0.9)

for(i in 1:8)
{
  if (which(aaa[,i]==sort(aaa[,i])[7])==7)
  {text(x=place[i], y=sort(aaa[,i])[7]/2, labels="*")}
  else {text(x=place[i], y=sum(aaa[c((which(aaa[,i]==sort(aaa[,i])[7])+1):7),i])+sort(aaa[,i])[7]/2, labels="*")}
  if (which(aaa[,i]==sort(aaa[,i])[6])==7)
  {text(x=place[i], y=sort(aaa[,i])[6]/2, labels="* *")}
  else {text(x=place[i], y=sum(aaa[c((which(aaa[,i]==sort(aaa[,i])[6])+1):7),i])+sort(aaa[,i])[6]/2, labels="* *")}
  if (which(aaa[,i]==sort(aaa[,i])[5])==7)
  {text(x=place[i], y=sort(aaa[,i])[5]/2, labels="* * *")}
  else {text(x=place[i], y=sum(aaa[c((which(aaa[,i]==sort(aaa[,i])[5])+1):7),i])+sort(aaa[,i])[5]/2, labels="* * *")}
}

mtext(side=1, oma=TRUE, at=c(1.425,3.875,6.3225,8.7725), text=c("Micro","Chloro","18S","16S"), line=1.75, cex=0.8)
mtext(side=1, oma=TRUE, at=place, text=c("Ice","Open"), line=0, cex=0.75)
mtext(side=1, oma=TRUE, at=place, text=c("cover","water"), line=0.75, cex=0.75)

#3
par(mar=c(0,0,0,0))
plot(x=c(1:10), y=c(1:10), ylab=NA, xlab=NA, col="white", xaxt="n", yaxt="n", bty="n")
text(x=1, y=5, adj=0, labels="B. Phototrophic eukaryotes only", cex=1.3)

#4
place <- c(0.9,1.95,3.35,4.4,5.8,6.85)
par(mar=c(3,4,0,0))
barplot(bbb[c(7:1),], col=rainbow(7), las=1, names.arg=c("","","","","",""),
        ylab="Relative variance partitioning (%)", ylim=c(0,105), bty="n",space=c(0.4,0.05,0.4,0.05,0.4,0.05))
lines(x=c(1:9), y=rep(0,9), col="black")
text(x=place, y=c(round(sum(bbb[,1])),round(sum(bbb[,2])),round(sum(bbb[,3])),round(sum(bbb[,4])),round(sum(bbb[,5])),round(sum(bbb[,6])))+3, labels=c(round(sum(bbb[,1])),round(sum(bbb[,2])),round(sum(bbb[,3])),round(sum(bbb[,4])),round(sum(bbb[,5])),round(sum(bbb[,6]))), cex=0.9)


for(i in 1:6)
{
  text(x=place[i], y=sum(bbb[c((which(bbb[,i]==sort(bbb[,i])[7])+1):7),i])+sort(bbb[,i])[7]/2, labels="*")
  if (which(bbb[,i]==sort(bbb[,i])[6])==7)
  {text(x=place[i], y=sort(bbb[,i])[6]/2, labels="* *")}
  else {text(x=place[i], y=sum(bbb[c((which(bbb[,i]==sort(bbb[,i])[6])+1):7),i])+sort(bbb[,i])[6]/2, labels="* *")}
  if (which(bbb[,i]==sort(bbb[,i])[5])==7)
  {text(x=place[i], y=sort(bbb[,i])[5]/2, labels="* * *")}
  else {text(x=place[i], y=sum(bbb[c((which(bbb[,i]==sort(bbb[,i])[5])+1):7),i])+sort(bbb[,i])[5]/2, labels="* * *")}
}

mtext(side=1, oma=TRUE, at=c(1.425,3.875,6.3225), text=c("Micro","Chloro","18S"), line=1.75, cex=0.8)
mtext(side=1, oma=TRUE, at=place, text=c("Ice","Open"), line=0, cex=0.75)
mtext(side=1, oma=TRUE, at=place, text=c("cover","water"), line=0.75, cex=0.75)

#5
par(mar=c(0,0,0,0))
plot(x=c(1:10), y=c(1:10), ylab=NA, xlab=NA, col="white", xaxt="n", yaxt="n", bty="n")
legend(x=5, xjust=0.5, yjust=0.5, y=5, horiz=F, bty="n", legend=c("Al","NMDS Axis 1","DOC","Lake","Si","TN","TP"),
       fill=rev(rainbow(7)), xpd=T, cex=1.5)

dev.off()

