require(raster)

d_cembroides = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_cembroides_01.img")
d_discolor = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_discolor_01.img")

png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/problems_distribution/discolor_cembroides.png", width = 800, height = 400)
par(mfcol=c(1,2))
plot(crop(d_cembroides, c(-120,-90,10,40)), main="P. cembroides")
plot(crop(d_discolor, c(-120,-90,10,40)), main="P. discolor")
dev.off()

