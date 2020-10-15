#Code for drop the area in the west russia for pinus sibirica. See note book for more information about the rationale of this

#Cargamos area distribucion
require(raster)
sibi_distribution = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/datos_brutos/Maps/p_sibirica_01.img") #load the distributon from the initial path, without changes

#select the area to drop 
plot(sibi_distribution)
#e = drawExtent() #directly on the map 
e = extent(29.97157, 33.84258, 67.72879, 69.48786) #these values have been obtained directly from the map using drawExtent(). With them and extent() function we can create a extent object, whici will be our interest area. 

#see in the map 
plot(sibi_distribution)
plot(e, add=T)

#select the cells of the raster, INSIDE our interest area with value=1
drop.idx <- which(sibi_distribution[cellsFromExtent(sibi_distribution, e)] %in% 1) #cellsFromExtent return the number of cells of a raster inside a extent object. Select these cells from our raster and from them, look for cells with value = 1, it is to say, cells with distribution. 

#give NA value to the cells with value=1 inside our interest area
sibi_distribution[cellsFromExtent(sibi_distribution, e)[drop.idx]] <- 0 #from the cells inside our interest area "cellsFromExtent(sibi_distribution, e)", select the cells with value=1. These cells will have 0, it is to say NO presence.  

#test 
plot(sibi_distribution)
plot(e, add=TRUE) #no distribution inside our interest area


#write the raster
writeRaster(sibi_distribution, "/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_sibirica_01.img", overwrite=TRUE)