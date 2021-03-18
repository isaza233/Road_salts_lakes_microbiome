library(mapedit)
library(sf)

north.arrow = function(x, y, h,lab="North",lab.pos="below") {
  polygon(c(x, x, x + h/0.95), c(y - (0.9*h), y, y - (1 + sqrt(3)/2) * h), col = "black", border = NA)
  polygon(c(x, x + h/0.95, x, x - h/0.95), c(y - (0.9*h), y - (1 + sqrt(3)/2) * h, y, y - (1 + sqrt(3)/2) * h))
  if(lab.pos=="below") text(x, y-(2.5*h), lab, adj = c(0.5, 0), cex = 1)
  else text(x, y+(0.25*h), lab, adj = c(0.5, 0), cex = 1.5)
}

## how to make new polygons/lines
#new<- editMap()
#save(new, file="map/new.RData")
#then save new.RData under another name

## load existing polygons/lines
load("map/LCT.RData")
load("map/LSA.RData")
load("map/stlawrence.RData")
load("map/LCL.RData")
load("map/LSC.RData")
load("map/quebec.RData")
load("map/desmaures.RData")
load("map/highway1.RData")
load("map/highway2.RData")
load("map/highway3.RData")
load("map/highway4.RData")
load("map/highway5.RData")
load("map/highway6.RData")
load("map/highway7.RData")


## figure starts here
library(graphicsutils)

tiff(width=12, height=12, units="cm", res=300, pointsize=9, filename="images/Figure1_Fournier.tif")
par(mar=c(3,3,1,1))
plot(x=c(1:10), y=c(1:10), col="black",xlim=c(288.27,288.92), ylim=c(46.68,47.015), yaxs="i",xaxs="i",
     xlab=NA, ylab=NA, xaxt="n", yaxt="n")
axis(side=2, at=seq(46.7, 47, 0.1), labels=c("46.7", "46.8", "46.9", "47.0"), las=2)
axis(side=1, at=seq(288.3, 288.9, 0.2), labels=seq(-71.70, -71.10, 0.2))
polygon(x=LCT$geometry[[1]][[1]][,1], y=LCT$geometry[[1]][[1]][,2], col="light blue", border="black")
polygon(x=LSA$geometry[[1]][[1]][,1], y=LSA$geometry[[1]][[1]][,2], col="light blue", border="black")
polygon(x=LCL$geometry[[1]][[1]][,1], y=LCL$geometry[[1]][[1]][,2], col="light blue", border="black")
polygon(x=LSC$geometry[[1]][[1]][,1], y=LSC$geometry[[1]][[1]][,2], col="light blue", border="black")
polygon(x=stlawrence$geometry[[1]][[1]][,1], y=stlawrence$geometry[[1]][[1]][,2], col="light blue", 
        border="black")
lines(x=st_coordinates(highway1$geometry)[,1], y=st_coordinates(highway1$geometry)[,2], col="red", lty=1, lwd=1.5)
lines(x=st_coordinates(highway2$geometry)[,1], y=st_coordinates(highway2$geometry)[,2], col="red", lty=1, lwd=1.5)
lines(x=st_coordinates(highway3$geometry)[,1], y=st_coordinates(highway3$geometry)[,2], col="red", lty=1, lwd=1.5)
lines(x=st_coordinates(highway4$geometry)[,1], y=st_coordinates(highway4$geometry)[,2], col="red", lty=1, lwd=1.5)
lines(x=st_coordinates(highway5$geometry)[,1], y=st_coordinates(highway5$geometry)[,2], col="red", lty=1, lwd=1.5)
lines(x=st_coordinates(highway6$geometry)[,1], y=st_coordinates(highway6$geometry)[,2], col="red", lty=1, lwd=1.5)
lines(x=st_coordinates(highway7$geometry)[,1], y=st_coordinates(highway7$geometry)[,2], col="red", lty=1, lwd=1.5)
#lines(x=st_coordinates(desmaures$geometry)[,1], y=st_coordinates(desmaures$geometry)[,2], col="red", lty=2)


text(x=288.3291,y=46.97688, label="Lake Clair", cex=0.85)
text(x=288.6048,y=46.76888, label="Lake Saint-Augustin", cex=0.85)
text(x=288.6022,y=46.97064, label="Lake Saint-Charles", cex=0.85)
text(x=288.6401,y=46.95668, label="Lake Clément", adj=0, cex=0.85)

text(x=288.7550,y=46.80924, label="Quebec", col="black", font=2)
text(x=288.6691,y=46.98719, label="Stoneham-et-Tewkesbury", col="black", font=2, cex=0.8)
text(x=288.4867,y=46.75542, label="Saint-Augustin-", col="black", font=2, cex=0.8)
text(x=288.4867,y=46.74542, label="de-desmaures", col="black", font=2, cex=0.8)
text(x=288.3796,y=46.92, label="Sainte-Catherine-de-", col="black", font=2, cex=0.8)
text(x=288.3796,y=46.91, label="la-Jacques-Cartier", col="black", font=2, cex=0.8)

text(x=288.5400,y=46.71476, label="St-Lawrence River", srt=15)

rect(xleft=288.74, ybottom=46.707, xright=288.86, ytop=46.714, col="black", border="black")
rect(xleft=288.78, ybottom=46.707, xright=288.82, ytop=46.714, col="white", border="black")
text(x=c(288.74, 288.78, 288.82, 288.86), y=rep(46.72,4), labels=c(0,5,10,"15 km"), adj=0.2,
     cex=0.9)

legend(x=288.74, y=46.705, legend="Highways", col="red", pch=NA, lty=1, 
       horiz=TRUE, lwd=2, bty="n", x.intersp=0.8)

north.arrow(x=288.85,y=46.78,h=0.015, lab="N", lab.pos="below")

## north america 
pathLogo <- "map/North_America.jpg"
pchImage(x=288.39, y=46.83, file=pathLogo, cex.x=3.2, cex.y=3.2, atcenter=TRUE)
rect(288.2896, 46.77632, 288.5006, 46.88439, lwd=1, lty=1)
points(x=288.4469, y=46.82969, pch=-0x2605, col="red")

dev.off()