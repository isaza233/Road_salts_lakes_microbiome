source("Rcodes/rect_heatmap_16S_seas_L2FC_2017.R")
source("Rcodes/rect_heatmap_16S_cond_L2FC_2017.R")
source("Rcodes/rect_heatmap_18S_seas_L2FC_2017.R")
source("Rcodes/rect_heatmap_18S_cond_L2FC_2017.R")

library(magick)
#ressources: https://cran.r-project.org/web/packages/magick/vignettes/intro.html#Read_and_write
prokseas <- image_read("images/16Sseas.png")
eukseas <- image_read("images/18Sseas.png")
eukseas <-image_annotate(eukseas, "A. 18S", size=80, gravity="northwest", location="+110+100", color="black")
prokseas <-image_annotate(prokseas, "B. 16S", size=80, gravity="northwest", location="+110+100", color="black")

img <- c(eukseas,prokseas)
img <- image_scale(img, "3000x1333")
tiger <- image_append(img)
image_write(tiger, path="images/Figure7_Fournier.jpeg", format="jpeg")

#vertical
#img <- c(eukseas,prokseas)
#img <- image_scale(img, "3000x1333")
#tiger <- image_append(img, stack=TRUE)
#image_write(tiger, path="C:/Users/Isaza/Documents/Rplot002.png", format="png")

prokcond <- image_read("images/16Scond.png")
eukcond <- image_read("images/18Scond.png")
eukcond  <-image_annotate(eukcond, "A. 18S", size=80, gravity="northwest", location="+110+100", color="black")
prokcond<-image_annotate(prokcond, "B. 16S", size=80, gravity="northwest", location="+110+100", color="black")
img <- c(eukcond,prokcond)
img <- image_scale(img, "3000x1333")
tiger <- image_append(img)
image_write(tiger, path="images/Figure6_Fournier.jpeg", format="jpeg")