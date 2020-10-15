#Code for calculating buffer around distribution and compare with euforgen data. The rationales of this is that we create the buffers to cover possible presences that Critchfield maps does not cover. We have used species with data in Euforgen as a test. If the buffer cover almost all distribution area according to euforgen, the buffer works. Of course, we only use speces for which there is data in euforgen. 


###definimos el directorio de trabajo
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus")

##########################################################
#####Calculating of widht of buffer#######################
##########################################################

#Librerias
require(raster)

#Load only a variable for having the background gray
bio1 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio1.asc")

#load atlas`s distribucion range of all species with data on euforgen
distr_halepensis = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_halepensis_01.img")
distr_brutia = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_brutia_01.img")
distr_cembra = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_cembra_01.img")
distr_nigra = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_nigra_01.img")
distr_pinaster = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_pinaster_01.img")
distr_pinea = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_pinea_01.img")
distr_sylvestris = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_sylvestris_01.img")
distr_peuce = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_peuce_01.img")
distr_heldreichii = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_heldreichii_01.img")

###make the crop
#halepensis
plot(distr_halepensis) 
distr_crop_halepensis =  crop(distr_halepensis,c(-13.18, 40, 30, 47))
plot(distr_crop_halepensis) 
#brutia
plot(distr_brutia) 
distr_crop_brutia =  crop(distr_brutia,c(17, 55, 25, 50))
plot(distr_crop_brutia) 
#cembra
plot(distr_cembra) 
distr_crop_cembra =  crop(distr_cembra,c(-5, 50, 30, 55))
plot(distr_crop_cembra)
#nigra
plot(distr_nigra) 
distr_crop_nigra =  crop(distr_nigra,c(-13.18, 45, 30, 50))
plot(distr_crop_nigra)
#pinaster
plot(distr_pinaster) 
distr_crop_pinaster =  crop(distr_pinaster,c(-13.18, 40, 30, 47))
plot(distr_crop_pinaster)
#pinea
plot(distr_pinea) 
distr_crop_pinea =  crop(distr_pinea,c(-13.18, 47, 30, 47))
plot(distr_crop_pinea)
#sylvestris
plot(distr_sylvestris) 
distr_crop_sylvestris =  crop(distr_sylvestris,c(-13.18, 140, 30, 75))
plot(distr_crop_sylvestris)
#peuce
plot(distr_peuce) 
distr_crop_peuce =  crop(distr_peuce,c(10, 30, 30, 50))
plot(distr_crop_peuce)
#heldreichii
plot(distr_heldreichii) 
distr_crop_heldreichii =  crop(distr_heldreichii,c(10, 30, 30, 50))
plot(distr_crop_heldreichii)


####convert rasters in polygons
#halepensis
polygon_halepensis_atlas = rasterToPolygons(distr_crop_halepensis, fun=function(x){x==1}, n=16) #esta funcion de raster te transforma un raster completo o una parte del mismo en un poliogno. En nuestro caso solo queremos las celdas con valor=1, es decir, presencias. Por eso ponemos x==1. 
new.CRS = CRS("+init=epsg:4326")
polygon_halepensis_atlas  = spTransform(polygon_halepensis_atlas, new.CRS) #cambiamos el CRS del poliogno para que coincida con el de euforgen. He comprobado que tanto los mapas de critfield rasterizados y los shapefile de euforgen caen donde deben y así es, además los he comparado con lo correspondientes mapas de origen. Lo único llamativo es la perdida de resolución brutal en los de critfield cuando se rasterizó.
#comprobamos is ha ido bien
plot(distr_crop_halepensis)
plot(polygon_halepensis_atlas, add=TRUE) #el poligono coincide con el area verde del raster. 
#brutia
polygon_brutia_atlas = rasterToPolygons(distr_crop_brutia, fun=function(x){x==1}, n=16)
polygon_brutia_atlas  = spTransform(polygon_brutia_atlas, new.CRS)
plot(distr_crop_brutia)
plot(polygon_brutia_atlas, add=TRUE) 
#cembra
polygon_cembra_atlas = rasterToPolygons(distr_crop_cembra, fun=function(x){x==1}, n=16)
polygon_cembra_atlas  = spTransform(polygon_cembra_atlas, new.CRS)
plot(distr_crop_cembra)
plot(polygon_cembra_atlas, add=TRUE) 
#nigra
polygon_nigra_atlas = rasterToPolygons(distr_crop_nigra, fun=function(x){x==1}, n=16)
polygon_nigra_atlas  = spTransform(polygon_nigra_atlas, new.CRS)
plot(distr_crop_nigra)
plot(polygon_nigra_atlas, add=TRUE) 
#pinaster
polygon_pinaster_atlas = rasterToPolygons(distr_crop_pinaster, fun=function(x){x==1}, n=16)
polygon_pinaster_atlas  = spTransform(polygon_pinaster_atlas, new.CRS)
plot(distr_crop_pinaster)
plot(polygon_pinaster_atlas, add=TRUE)
#pinea
polygon_pinea_atlas = rasterToPolygons(distr_crop_pinea, fun=function(x){x==1}, n=16)
polygon_pinea_atlas  = spTransform(polygon_pinea_atlas, new.CRS)
plot(distr_crop_pinea)
plot(polygon_pinea_atlas, add=TRUE)
#sylvestris
polygon_sylvestris_atlas = rasterToPolygons(distr_crop_sylvestris, fun=function(x){x==1}, n=16)
polygon_sylvestris_atlas  = spTransform(polygon_sylvestris_atlas, new.CRS)
plot(distr_crop_sylvestris)
plot(polygon_sylvestris_atlas, add=TRUE)
#peuce
polygon_peuce_atlas = rasterToPolygons(distr_crop_peuce, fun=function(x){x==1}, n=16)
polygon_peuce_atlas  = spTransform(polygon_peuce_atlas, new.CRS)
plot(distr_crop_peuce)
plot(polygon_peuce_atlas, add=TRUE)
#heldreichii
polygon_heldreichii_atlas = rasterToPolygons(distr_crop_heldreichii, fun=function(x){x==1}, n=16)
polygon_heldreichii_atlas  = spTransform(polygon_heldreichii_atlas, new.CRS)
plot(distr_crop_heldreichii)
plot(polygon_heldreichii_atlas, add=TRUE)


####Load the polygon of natural distribution from euforgem
#halepensis
polygon_halepensis_euforgen = shapefile("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/shapefiles/shapefiles/Pinus_halepensis/Pinus halepensis.shp")
polygon_halepensis_euforgen #es un shapefile lo que nos da euforgen
#Change to the CRS of my variables
new.CRS = CRS("+init=epsg:4326")
polygon_halepensis_euforgen  = spTransform(polygon_halepensis_euforgen, new.CRS) #cambiamos el CRS 
#brutia
polygon_brutia_euforgen = shapefile("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/shapefiles/shapefiles/Pinus_brutia/Pinus brutia.shp")
polygon_brutia_euforgen  = spTransform(polygon_brutia_euforgen, new.CRS) 
#cembra
polygon_cembra_euforgen = shapefile("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/shapefiles/shapefiles/Pinus_cembra/Pinus cembra.shp")
polygon_cembra_euforgen  = spTransform(polygon_cembra_euforgen, new.CRS) 
#nigra
polygon_nigra_euforgen = shapefile("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/shapefiles/shapefiles/Pinus_nigra/Pinus nigra.shp")
polygon_nigra_euforgen  = spTransform(polygon_nigra_euforgen, new.CRS) 
#pinaster
polygon_pinaster_euforgen = shapefile("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/shapefiles/shapefiles/Pinus_pinaster/Pinus pinaster.shp")
polygon_pinaster_euforgen  = spTransform(polygon_pinaster_euforgen, new.CRS) 
#pinea
polygon_pinea_euforgen = shapefile("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/shapefiles/shapefiles/Pinus_pinea/Pinus pinea.shp")
polygon_pinea_euforgen  = spTransform(polygon_pinea_euforgen, new.CRS) 
#sylvestris
polygon_sylvestris_euforgen = shapefile("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/shapefiles/shapefiles/Pinus_sylvestris/Pinus Sylvestris.shp")
polygon_sylvestris_euforgen  = spTransform(polygon_sylvestris_euforgen, new.CRS) 
#peuce
polygon_peuce_euforgen = shapefile("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/shapefiles/shapefiles/Pinus_peuce/Pinus_peuce_EUFORGEN.shp")
polygon_peuce_euforgen  = spTransform(polygon_peuce_euforgen, new.CRS) 
#heldreichii
polygon_heldreichii_euforgen = shapefile("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/shapefiles/shapefiles/Pinus_heldreichii/Pinus_leucodermis_EUFORGEN.shp")
polygon_heldreichii_euforgen  = spTransform(polygon_heldreichii_euforgen, new.CRS) 

####Calculate overlap between polygons of euforgen and atlas`s map
#halepensis
overlap_halepensis <- intersect(polygon_halepensis_atlas, polygon_halepensis_euforgen)
plot(overlap_halepensis)
##other way with the same results in relation to area
require(rgeos)
overlap2 = gIntersection(polygon_halepensis_atlas, polygon_halepensis_euforgen)
plot(overlap2)
#brutia
overlap_brutia <- intersect(polygon_brutia_atlas, polygon_brutia_euforgen)
plot(overlap_brutia)
#cembra
overlap_cembra <- intersect(polygon_cembra_atlas, polygon_cembra_euforgen)
plot(overlap_cembra)
#nigra
overlap_nigra <- intersect(polygon_nigra_atlas, polygon_nigra_euforgen)
plot(overlap_nigra)
#pinaster
overlap_pinaster <- intersect(polygon_pinaster_atlas, polygon_pinaster_euforgen)
plot(overlap_pinaster)
#pinea
overlap_pinea <- intersect(polygon_pinea_atlas, polygon_pinea_euforgen)
plot(overlap_pinea)
#sylvestris
overlap_sylvestris <- intersect(polygon_sylvestris_atlas, polygon_sylvestris_euforgen)
plot(overlap_sylvestris)
#peuce
overlap_peuce <- intersect(polygon_peuce_atlas, polygon_peuce_euforgen)
plot(overlap_peuce)
#heldreichii
overlap_heldreichii <- intersect(polygon_heldreichii_atlas, polygon_heldreichii_euforgen)
plot(overlap_heldreichii)


#plot with all polygons of atlas
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/atlas.png", width=1200, height=1200, pointsize=30)
par(mfcol=c(3,3))
plot(distr_crop_halepensis)
plot(distr_crop_brutia)
plot(distr_crop_cembra)
plot(distr_crop_nigra)
plot(distr_crop_pinaster)
plot(distr_crop_pinea)
plot(distr_crop_sylvestris)
plot(distr_crop_peuce)
plot(distr_crop_heldreichii)
dev.off()

#plot with all polygons of euforgen
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/euforgen.png", width=1200, height=1200, pointsize=30)
par(mfcol=c(3,3))
plot(polygon_halepensis_euforgen)
plot(polygon_brutia_euforgen)
plot(polygon_cembra_euforgen)
plot(polygon_nigra_euforgen)
plot(polygon_pinaster_euforgen)
plot(polygon_pinea_euforgen)
plot(polygon_sylvestris_euforgen)
plot(polygon_peuce_euforgen)
plot(polygon_heldreichii_euforgen)
dev.off()

#plot atlas`s polygons sobre euforgen`s polygons
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/atlas_sobre_euforgen.png", width=1200, height=1200, pointsize=30)
par(mfcol=c(3,3))
plot(polygon_halepensis_euforgen)
plot(polygon_halepensis_atlas, add=T)
plot(polygon_brutia_euforgen)
plot(polygon_brutia_atlas, add=T)
plot(polygon_cembra_euforgen)
plot(polygon_cembra_atlas, add=T)
plot(polygon_nigra_euforgen)
plot(polygon_nigra_atlas, add=T)
plot(polygon_pinaster_euforgen)
plot(polygon_pinaster_atlas, add=T)
plot(polygon_pinea_euforgen)
plot(polygon_pinea_atlas, add=T)
plot(polygon_sylvestris_euforgen)
plot(polygon_sylvestris_atlas, add=T)
plot(polygon_peuce_euforgen)
plot(polygon_peuce_atlas, add=T)
plot(polygon_heldreichii_euforgen)
plot(polygon_heldreichii_atlas, add=T)
dev.off()

#plot euforgen`s polygons sobre atlas`s polygons 
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/euforgen_sobre_atlas.png", width=1200, height=1200, pointsize=30)
par(mfcol=c(3,3))
plot(polygon_halepensis_atlas)
plot(polygon_halepensis_euforgen, add=T)
plot(polygon_brutia_atlas)
plot(polygon_brutia_euforgen, add=T)
plot(polygon_cembra_atlas)
plot(polygon_cembra_euforgen, add=T)
plot(polygon_nigra_atlas)
plot(polygon_nigra_euforgen, add=T)
plot(polygon_pinaster_atlas)
plot(polygon_pinaster_euforgen, add=T)
plot(polygon_pinea_atlas)
plot(polygon_pinea_euforgen, add=T)
plot(polygon_sylvestris_atlas)
plot(polygon_sylvestris_euforgen, add=T)
plot(polygon_peuce_atlas)
plot(polygon_peuce_euforgen, add=T)
plot(polygon_heldreichii_atlas)
plot(polygon_heldreichii_euforgen, add=T)
dev.off()

###plot species by species
#halepensis
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_atlas_euforgen/overlap_halepensis.png", width=1200, height=1200, pointsize=30)
plot(crop(bio1, c(-13.18, 40, 30, 47)), col="gray80", main="P. halepensis")
plot(polygon_halepensis_euforgen, col="red", add=T)
plot(polygon_halepensis_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.8)
dev.off()
#brutia
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_atlas_euforgen/overlap_brutia.png", width=1200, height=1200, pointsize=30)
plot(crop(bio1, c(17, 55, 25, 50)), col="gray80", main="P. brutia")
plot(polygon_brutia_euforgen, col="red", add=T)
plot(polygon_brutia_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.8)
dev.off()
#cembra
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_atlas_euforgen/overlap_cembra.png", width=1200, height=1200, pointsize=30)
plot(crop(bio1, c(-5, 50, 30, 55)), col="gray80", main="P. cembra")
plot(polygon_cembra_euforgen, col="red", add=T)
plot(polygon_cembra_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.8)
dev.off()
#nigra
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_atlas_euforgen/overlap_nigra.png", width=1200, height=1200, pointsize=30)
plot(crop(bio1,c(-13.18, 45, 30, 50)), col="gray80", main="P. nigra")
plot(polygon_nigra_euforgen, col="red", add=T)
plot(polygon_nigra_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.8)
dev.off()
#pinaster
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_atlas_euforgen/overlap_pinaster.png", width=1200, height=1200, pointsize=30)
plot(crop(bio1,c(-13.18, 40, 30, 47)), col="gray80", main="P. pinaster")
plot(polygon_pinaster_euforgen, col="red", add=T)
plot(polygon_pinaster_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.8)
dev.off()
#pinea
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_atlas_euforgen/overlap_pinea.png", width=1200, height=1200, pointsize=30)
plot(crop(bio1, c(-13.18, 47, 30, 47)), col="gray80", main="P. pinea")
plot(polygon_pinea_euforgen, col="red", add=T)
plot(polygon_pinea_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.8)
dev.off()
#sylvestris
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_atlas_euforgen/overlap_sylvestris.png", width=1200, height=1200, pointsize=30)
plot(crop(bio1, c(-13.18, 140, 30, 75)), col="gray80", main="P. sylvestris")
plot(polygon_sylvestris_euforgen, col="red", add=T)
plot(polygon_sylvestris_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.8)
dev.off()
#peuce
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_atlas_euforgen/overlap_peuce.png", width=1200, height=1200, pointsize=30)
plot(crop(bio1, c(10, 30, 30, 50)), col="gray80", main="P. peuce")
plot(polygon_peuce_euforgen, col="red", add=T)
plot(polygon_peuce_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.8)
dev.off()
#heldreichii
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_atlas_euforgen/overlap_heldreichii.png", width=1200, height=1200, pointsize=30)
plot(crop(bio1, c(10, 30, 30, 50)), col="gray80", main="P. heldreichii")
plot(polygon_heldreichii_euforgen, col="red", add=T)
plot(polygon_heldreichii_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.8)
dev.off()

####Calcualte the areas of each polygon and percentage of overlap
#There is a problem with the projection and gArea, but it does not matter because the decisión of the number of cells around the critfield distribution was made by Nick in basis on the figures with and without buffer. With two cells of buffer most of euforgen areas are covered. 
#halepensis
percent_overlap_halepensis_atlas = (gArea(overlap_halepensis)/gArea(polygon_halepensis_atlas))*100 #percentage of atlas area (critfield) overlapped with euforgen
percent_overlap_halepensis_euforgen = (gArea(overlap_halepensis)/gArea(polygon_halepensis_euforgen))*100 #percentage of atlas euforgen overlapped with area (critfield)
percent_overlap_halepensis_atlas
percent_overlap_halepensis_euforgen

#brutia
percent_overlap_brutia_atlas = (gArea(overlap_brutia)/gArea(polygon_brutia_atlas))*100 #percentage of atlas area (critfield) overlapped with euforgen
percent_overlap_brutia_euforgen = (gArea(overlap_brutia)/gArea(polygon_brutia_euforgen))*100 #percentage of atlas euforgen overlapped with area (critfield)
percent_overlap_brutia_atlas
percent_overlap_brutia_euforgen

#cembra
percent_overlap_cembra_atlas = (gArea(overlap_cembra)/gArea(polygon_cembra_atlas))*100 #percentage of atlas area (critfield) overlapped with euforgen
percent_overlap_cembra_euforgen = (gArea(overlap_cembra)/gArea(polygon_cembra_euforgen))*100 #percentage of atlas euforgen overlapped with area (critfield)
percent_overlap_cembra_atlas
percent_overlap_cembra_euforgen

#nigra
percent_overlap_nigra_atlas = (gArea(overlap_nigra)/gArea(polygon_nigra_atlas))*100 #percentage of atlas area (critfield) overlapped with euforgen
percent_overlap_nigra_euforgen = (gArea(overlap_nigra)/gArea(polygon_nigra_euforgen))*100 #percentage of atlas euforgen overlapped with area (critfield)
percent_overlap_nigra_atlas
percent_overlap_nigra_euforgen

#pinaster
percent_overlap_pinaster_atlas = (gArea(overlap_pinaster)/gArea(polygon_pinaster_atlas))*100 #percentage of atlas area (critfield) overlapped with euforgen
percent_overlap_pinaster_euforgen = (gArea(overlap_pinaster)/gArea(polygon_pinaster_euforgen))*100 #percentage of atlas euforgen overlapped with area (critfield)
percent_overlap_pinaster_atlas
percent_overlap_pinaster_euforgen

#pinea
percent_overlap_pinea_atlas = (gArea(overlap_pinea)/gArea(polygon_pinea_atlas))*100 #percentage of atlas area (critfield) overlapped with euforgen
percent_overlap_pinea_euforgen = (gArea(overlap_pinea)/gArea(polygon_pinea_euforgen))*100 #percentage of atlas euforgen overlapped with area (critfield)
percent_overlap_pinea_atlas
percent_overlap_pinea_euforgen

#sylvestris
percent_overlap_sylvestris_atlas = (gArea(overlap_sylvestris)/gArea(polygon_sylvestris_atlas))*100 #percentage of atlas area (critfield) overlapped with euforgen
percent_overlap_sylvestris_euforgen = (gArea(overlap_sylvestris)/gArea(polygon_sylvestris_euforgen))*100 #percentage of atlas euforgen overlapped with area (critfield)
percent_overlap_sylvestris_atlas
percent_overlap_sylvestris_euforgen

#peuce
percent_overlap_peuce_atlas = (gArea(overlap_peuce)/gArea(polygon_peuce_atlas))*100 #percentage of atlas area (critfield) overlapped with euforgen
percent_overlap_peuce_euforgen = (gArea(overlap_peuce)/gArea(polygon_peuce_euforgen))*100 #percentage of atlas euforgen overlapped with area (critfield)
percent_overlap_peuce_atlas
percent_overlap_peuce_euforgen

#heldreichii
percent_overlap_heldreichii_atlas = (gArea(overlap_heldreichii)/gArea(polygon_heldreichii_atlas))*100 #percentage of atlas area (critfield) overlapped with euforgen
percent_overlap_heldreichii_euforgen = (gArea(overlap_heldreichii)/gArea(polygon_heldreichii_euforgen))*100 #percentage of atlas euforgen overlapped with area (critfield)
percent_overlap_heldreichii_atlas
percent_overlap_heldreichii_euforgen


percent_overlaps = data.frame(c("halepensis", "brutia", "cembra", "nigra", "pinaster", "pinea", "sylvestris", "peuce", "heldreichii"), c(percent_overlap_halepensis_atlas, percent_overlap_brutia_atlas ,percent_overlap_cembra_atlas , percent_overlap_nigra_atlas ,percent_overlap_pinaster_atlas,percent_overlap_pinea_atlas,percent_overlap_sylvestris_atlas, percent_overlap_peuce_atlas, percent_overlap_heldreichii_atlas), c(percent_overlap_halepensis_euforgen, percent_overlap_brutia_euforgen ,percent_overlap_cembra_euforgen , percent_overlap_nigra_euforgen ,percent_overlap_pinaster_euforgen,percent_overlap_pinea_euforgen,percent_overlap_sylvestris_euforgen, percent_overlap_peuce_euforgen, percent_overlap_heldreichii_euforgen))
colnames(percent_overlaps)[1] <- "species"
colnames(percent_overlaps)[2] <- "% overlap critfield"
colnames(percent_overlaps)[3] <- "% overlap euforgen"
percent_overlaps


####Conclusions
##Nick saw the overlap of two polygons for the seven species and he decided to use an buffer of two cells. We will lose little areas like for example north africa of pinaster and halepensis, but it woths because in this way we are sure that we only use point included in the natural distribution of species. 

################################################################
#####Creating the buffer for species with euforgen data#########
################################################################

#load rgeos
require(rgeos)

#create the buffers for each speceis
polygon_halepensis_atlas_buffer = gBuffer(polygon_halepensis_atlas, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)
polygon_brutia_atlas_buffer = gBuffer(polygon_brutia_atlas, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)
polygon_cembra_atlas_buffer = gBuffer(polygon_cembra_atlas, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)
polygon_nigra_atlas_buffer = gBuffer(polygon_nigra_atlas, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)
polygon_pinaster_atlas_buffer = gBuffer(polygon_pinaster_atlas, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)
polygon_pinea_atlas_buffer = gBuffer(polygon_pinea_atlas, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)
polygon_sylvestris_atlas_buffer = gBuffer(polygon_sylvestris_atlas, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)
polygon_peuce_atlas_buffer = gBuffer(polygon_peuce_atlas, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)
polygon_heldreichii_atlas_buffer = gBuffer(polygon_heldreichii_atlas, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

#plot buffer on euforgen data for each specie
#halepensis
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_buffer_euforgen/overlap_buffer_euforgen_halepensis.png", width=1500, height=600, pointsize=30)
par(mfcol=c(1,2))
plot(crop(bio1, c(-13.18, 40, 30, 47)), col="gray80", main="P. halepensis")
plot(polygon_halepensis_euforgen, col="red", add=T)
plot(polygon_halepensis_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.5)
plot(crop(bio1, c(-13.18, 40, 30, 47)), col="gray80", main="P. halepensis")
plot(polygon_halepensis_atlas_buffer, add=T)
plot(polygon_halepensis_euforgen, col="red", add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield + buffer"), fill=c("red", "white"), cex=0.5)
dev.off()
#brutia
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_buffer_euforgen/overlap_buffer_euforgen_brutia.png", width=1500, height=600, pointsize=30)
par(mfcol=c(1,2))
plot(crop(bio1, c(17, 55, 25, 50)), col="gray80", main="P. brutia")
plot(polygon_brutia_euforgen, col="red", add=T)
plot(polygon_brutia_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.5)
plot(crop(bio1, c(17, 55, 25, 50)), col="gray80", main="P. brutia")
plot(polygon_brutia_atlas_buffer, add=T)
plot(polygon_brutia_euforgen, col="red", add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield + buffer"), fill=c("red", "white"), cex=0.5)
dev.off()
#cembra
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_buffer_euforgen/overlap_buffer_euforgen_cembra.png", width=1500, height=600, pointsize=30)
par(mfcol=c(1,2))
plot(crop(bio1, c(-5, 50, 30, 55)), col="gray80", main="P. cembra")
plot(polygon_cembra_euforgen, col="red", add=T)
plot(polygon_cembra_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.5)
plot(crop(bio1, c(-5, 50, 30, 55)), col="gray80", main="P. cembra")
plot(polygon_cembra_atlas_buffer, add=T)
plot(polygon_cembra_euforgen, col="red", add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield + buffer"), fill=c("red", "white"), cex=0.5)
dev.off()
#nigra
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_buffer_euforgen/overlap_buffer_euforgen_nigra.png", width=1500, height=600, pointsize=30)
par(mfcol=c(1,2))
plot(crop(bio1, c(-13.18, 45, 30, 50)), col="gray80", main="P. nigra")
plot(polygon_nigra_euforgen, col="red", add=T)
plot(polygon_nigra_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.5)
plot(crop(bio1, c(-13.18, 45, 30, 50)), col="gray80", main="P. nigra")
plot(polygon_nigra_atlas_buffer, add=T)
plot(polygon_nigra_euforgen, col="red", add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield + buffer"), fill=c("red", "white"), cex=0.5)
dev.off()
#pinaster
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_buffer_euforgen/overlap_buffer_euforgen_pinaster.png", width=1500, height=600, pointsize=30)
par(mfcol=c(1,2))
plot(crop(bio1, c(-13.18, 40, 30, 47)), col="gray80", main="P. pinaster")
plot(polygon_pinaster_euforgen, col="red", add=T)
plot(polygon_pinaster_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.5)
plot(crop(bio1, c(-13.18, 40, 30, 47)), col="gray80", main="P. pinaster")
plot(polygon_pinaster_atlas_buffer, add=T)
plot(polygon_pinaster_euforgen, col="red", add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield + buffer"), fill=c("red", "white"), cex=0.5)
dev.off()
#pinea
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_buffer_euforgen/overlap_buffer_euforgen_pinea.png", width=1500, height=600, pointsize=30)
par(mfcol=c(1,2))
plot(crop(bio1, c(-13.18, 47, 30, 47)), col="gray80", main="P. pinea")
plot(polygon_pinea_euforgen, col="red", add=T)
plot(polygon_pinea_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.5)
plot(crop(bio1, c(-13.18, 47, 30, 47)), col="gray80", main="P. pinea")
plot(polygon_pinea_atlas_buffer, add=T)
plot(polygon_pinea_euforgen, col="red", add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield + buffer"), fill=c("red", "white"), cex=0.5)
dev.off()
#sylvestris
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_buffer_euforgen/overlap_buffer_euforgen_sylvestris.png", width=1500, height=600, pointsize=30)
par(mfcol=c(1,2))
plot(crop(bio1, c(-13.18, 140, 30, 75)), col="gray80", main="P. sylvestris")
plot(polygon_sylvestris_euforgen, col="red", add=T)
plot(polygon_sylvestris_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.5)
plot(crop(bio1, c(-13.18, 140, 30, 75)), col="gray80", main="P. sylvestris")
plot(polygon_sylvestris_atlas_buffer, add=T)
plot(polygon_sylvestris_euforgen, col="red", add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield + buffer"), fill=c("red", "white"), cex=0.5)
dev.off()
#peuce
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_buffer_euforgen/overlap_buffer_euforgen_peuce.png", width=1500, height=600, pointsize=30)
par(mfcol=c(1,2))
plot(crop(bio1, c(10, 30, 30, 50)), col="gray80", main="P. peuce")
plot(polygon_peuce_euforgen, col="red", add=T)
plot(polygon_peuce_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.5)
plot(crop(bio1, c(10, 30, 30, 50)), col="gray80", main="P. peuce")
plot(polygon_peuce_atlas_buffer, add=T)
plot(polygon_peuce_euforgen, col="red", add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield + buffer"), fill=c("red", "white"), cex=0.5)
dev.off()
#heldreichii
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/overlap_buffer_euforgen/overlap_buffer_euforgen_heldreichii.png", width=1500, height=600, pointsize=30)
par(mfcol=c(1,2))
plot(crop(bio1, c(10, 30, 30, 50)), col="gray80", main="P. heldreichii")
plot(polygon_heldreichii_euforgen, col="red", add=T)
plot(polygon_heldreichii_atlas, add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield"), fill=c("red", "white"), cex=0.5)
plot(crop(bio1, c(10, 30, 30, 50)), col="gray80", main="P. heldreichii")
plot(polygon_heldreichii_atlas_buffer, add=T)
plot(polygon_heldreichii_euforgen, col="red", add=T)
legend(x="topright", legend=c("EUFORGEN", "Critfield + buffer"), fill=c("red", "white"), cex=0.5)
dev.off()

##figure paper
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/buffer/figure_paper.pdf")
#set graphic parameters
par(mfcol=c(5,2), mar=c(0,2,0,4)+0.1, oma = c(0,0.5,2,0), mai=c(0.4,0.5,0.1,0.1))

#brutia
plot(crop(bio1, c(20, 47, 32.5, 47.5)), col="gray80", legend=FALSE)
mtext(text=expression(bolditalic("P. brutia")), line=0.3)
plot(polygon_brutia_atlas_buffer, add=T)
plot(polygon_brutia_euforgen, col="black", add=T)
legend(x="topleft", legend=c("EUFORGEN", "Extended Critchfield"), fill=c("black", "white"), bg="white", cex=0.8)

#cembra
plot(crop(bio1, c(4.5, 30, 42.5, 51)), col="gray80", legend=FALSE)
mtext(text=expression(bolditalic("P. cembra")), line=0.3)
plot(polygon_cembra_atlas_buffer, add=T)
plot(polygon_cembra_euforgen, col="black", add=T)

#halepensis
plot(crop(bio1, c(-13.18, 40, 30, 47)), col="gray80", legend=FALSE)
mtext(text=expression(bolditalic("P. halepensis")), line=0.3)
plot(polygon_halepensis_atlas_buffer, add=T)
plot(polygon_halepensis_euforgen, col="black", add=T)

#heldreichii
plot(crop(bio1, c(12, 27.5, 37.5, 45.5)), col="gray80", legend=FALSE)
mtext(text=expression(bolditalic("P. heldreichii")), line=0.3)
plot(polygon_heldreichii_atlas_buffer, add=T)
plot(polygon_heldreichii_euforgen, col="black", add=T)

#nigra
plot(crop(bio1, c(-13.18, 45, 34, 50)), col="gray80", legend=FALSE)
mtext(text=expression(bolditalic("P. nigra")), line=0.3)
plot(polygon_nigra_atlas_buffer, add=T)
plot(polygon_nigra_euforgen, col="black", add=T)

#peuce
plot(crop(bio1, c(16, 27, 39, 45)), col="gray80", legend=FALSE)
mtext(text=expression(bolditalic("P. peuce")), line=0.3)
plot(polygon_peuce_atlas_buffer, add=T)
plot(polygon_peuce_euforgen, col="black", add=T)

#pinaster
plot(crop(bio1, c(-13.18, 14, 30, 49)), col="gray80", legend=FALSE)
mtext(text=expression(bolditalic("P. pinaster")), line=0.3)
plot(polygon_pinaster_atlas_buffer, add=T)
plot(polygon_pinaster_euforgen, col="black", add=T)

#pinea
plot(crop(bio1, c(-13.18, 43, 32.5, 47)), col="gray80", legend=FALSE)
mtext(text=expression(bolditalic("P. pinea")), line=0.3)
plot(polygon_pinea_atlas_buffer, add=T)
plot(polygon_pinea_euforgen, col="black", add=T)

#sylvestris
plot(crop(bio1, c(-13.18, 140, 30, 75)), col="gray80", legend=FALSE)
mtext(text=expression(bolditalic("P. sylvestris")), line=0.3)
plot(polygon_sylvestris_atlas_buffer, add=T)
plot(polygon_sylvestris_euforgen, col="black", add=T)
dev.off()

#save work.space
save.image("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/rdata/buffer_distribution_area.RData")
require(raster)