####code for changinf pumila buffer
#P. pumila is the only species that has its buffer cutted at the end the map, specific, the mapper for calculating ranges (plotting) and the PA buffer

#required packages
require(raster) #for work with rasters
require(rgeos) #for creating the buffer and the centroids of the cells without gbif points

#load distribution map of pumila
distribution_pumila = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_pumila_01.img")

#create a polygon from distribution raster
polygon_distribution_pumila = rasterToPolygons(distribution_pumila, fun=function(x){x==1}, dissolve=TRUE)

#Create distribution buffer around
polygon_distribution_buffer_pumila = gBuffer(polygon_distribution_pumila, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

##Remove see areas of the buffer
#We want that the polygon that serves as model for the buffers will be exact the extended distribution of pine species, i.e. the extended distribution with the 1 width buffer but removing see areas. In that way when we apply a buffer of 20 width, is exactly around the extended distribution of the species (not counting see areas).
#first rasterize the buffer
###create a raster from buffer+bianca area converting the buffer+bianca polygon in a raster (rasterize) 
raster_buffer_distribution = raster() #create a empty raster
extent(raster_buffer_distribution) = extent(distribution_pumila) #give it a extention similar to distribution area of pinus halepensis polygon
res(raster_buffer_distribution) = res(distribution_pumila) #give it a resolution similar to he biancas maps because the polygon tha we will rasterize was obtained from this maps (I can`t extract resolution from polygon_buffer). In this way, we maintenain the same resolution. 
raster_buffer_distribution  = rasterize(polygon_distribution_buffer_pumila,raster_buffer_distribution) #convert the polygon of the halepensis buffer in a raster
distribution_pumila #take a raster of pumila distribution (without buffer) for obtaining coast limits
raster_buffer_distribution = distribution_pumila*raster_buffer_distribution #multiply both raster for creating a raster with NAs as much as these rasters, it is to say, we want a raster with NA in pixel where some of these raster have a NA. The zones without cells in distribution are places of sea only, if we would have used a environmental variable, terrestrial cells without data of the variable also would be computed as sea cells. 
raster_buffer_distribution[!is.na(raster_buffer_distribution)] <- 1 #give value of 1 to all cells without NA, in this way include cell of distribution raster with 1 and cell outside the distribution raster but inside the raster buffer (with 0).  

#create the new polygon halepensis distribution without sea areas
polygon_distribution_buffer_pumila = rasterToPolygons(raster_buffer_distribution, fun=function(x){x==1}, dissolve=TRUE)

#plot
plot(crop(raster_buffer_distribution, polygon_distribution_buffer_pumila))
plot(polygon_distribution_buffer_pumila, add=T)#as you can see, the raster and the polygon are almost the same, except in the joints (small difference) and of course the sea areas. The buffers for variables selection are also calculated with rasters of the PA buffers.

#compare with the initial distribution
plot(crop(distribution_pumila, polygon_distribution_buffer_pumila))
plot(polygon_distribution_buffer_pumila, add=T)

#save raster of distribution buffer
writeRaster(raster_buffer_distribution, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/pumila_distribution_buffer.asc", format="ascii", overwrite=TRUE)

#save kml file to see in google earth
plotKML::kml(polygon_distribution_pumila, file.name = paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/kml_files/pumila_polygon_distribution.kml", sep=""))
plotKML::kml(polygon_distribution_buffer_pumila, file.name = paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/kml_files/pumila_polygon_distribution_buffer.kml", sep=""))



##########################
####### PA buffer ########
##########################

#Create PA buffer around extended distribution
polygon_PA_buffer_pumila = gBuffer(polygon_distribution_buffer_pumila, byid=FALSE, id=NULL, width=22.5, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

#plot the distribution and the buffer
plot(distribution_pumila)
plot(polygon_PA_buffer_pumila, add=T)

#crop the PA_buffer using the distribution map, in that way we remove the areas outside the map
polygon_PA_buffer_pumila_cropped = crop(polygon_PA_buffer_pumila, distribution_pumila)
plot(polygon_PA_buffer_pumila, lty=1)
plot(polygon_PA_buffer_pumila_cropped, add=T, lty=2)
legend("bottomright", legend=c("complete PA buffer", "reduced PA buffer"), lty=c(1,2))

#calculate the difference between both buffers
extreme_PA_buffer = gDifference(polygon_PA_buffer_pumila, polygon_PA_buffer_pumila_cropped)

#plot polyogns
plot(polygon_PA_buffer_pumila)
plot(extreme_PA_buffer, add=T, col="red")
legend("bottomright", legend=c("complete PA buffer", "extreme PA buffer"), fill=c("white", "red"))

#plot the both polygons on the distribution map
plot(distribution_pumila)
plot(polygon_PA_buffer_pumila, add=T, lwd=4)
plot(extreme_PA_buffer, add=T, border="red")#the right extreme is not plotted, we have to changed the coordinates
legend("bottomright", legend=c("complete PA buffer", "extreme PA buffer"), fill=c("black", "red"))

#change the coordinates of the right extreme polygon, only the second polygon of that part, which should be in the left side of the map
extreme_PA_buffer@polygons[[1]]@Polygons[[2]]@coords[,1]#these are the long (x) coordinates of the second polygon
#we substract from each coordinate 180, the result will be substracted from 180 and the changed the sign
extreme_PA_buffer@polygons[[1]]@Polygons[[2]]@coords[,1] <- -(180 - (extreme_PA_buffer@polygons[[1]]@Polygons[[2]]@coords[,1] - 180))#Example to understand: The maximum x coordinate is 180, so all the points in the second polygon from 180 to 199 should be in the other extreme of the map. Imagine a point at 183 longitude degrees. The difference from 180, that it is 3 (183-180), should be on the other side, i.e. should be rested to 180 -> 177. Therefore: - (180 - (183-180)) -> -177. With 180 would be: - (180 - (180-180)) = -180, which is ok because we are exactly in the end of one size of the map (180) and beginning of the other side (-180).

#plot
plot(distribution_pumila)
plot(polygon_PA_buffer_pumila, add=T, lwd=4)
plot(extreme_PA_buffer, add=T, border="red")
legend("bottomright", legend=c("PA buffer inside the map", "the rest of PA buffer"), fill=c("black", "red"))

#bind both polygons
final_PA_buffer = union(extreme_PA_buffer, polygon_PA_buffer_pumila)

#remove the areas outside the map. Up for both polygons, left only for the polygon_PA_buffer_pumila. The left part of the  polygon_PA_buffer_pumila has been put in the other side, whilst the up part is not interesting because there is no earth in that
final_PA_buffer = crop(final_PA_buffer, distribution_pumila) #NO pongo la perte del buffer que sale por arriba porque: 1) No HAY tierra en esa zona, por tanto no se van a poner PAs; 2) Ese área sobre el globo (google earth) es super pequeña, un círculo sobre el polo.
    #Si quiseras añadir esa parte del polígono:
    #extreme_PA_buffer@polygons[[1]]@Polygons[[1]]@coords[,1:2]
    #extreme_PA_buffer@polygons[[1]]@Polygons[[1]]@coords[,2] <- (90 - (abs(extreme_PA_buffer@polygons[[1]]@Polygons[[1]]@coords[,2] - 90)))
    #extreme_PA_buffer@polygons[[1]]@Polygons[[1]]@coords[,1] <- -(extreme_PA_buffer@polygons[[1]]@Polygons[[1]]@coords[,1])

#plot
plot(distribution_pumila)
plot(polygon_PA_buffer_pumila, add=T, lwd=4)
plot(final_PA_buffer, add=T, border="red")
legend("bottomright", legend=c("initial PA buffer", "final_PA_buffer"), fill=c("black", "red"))

#save the buffers into kml files for opening in goole earth
plotKML::kml(polygon_PA_buffer_pumila, file.name = paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/kml_files/pumila_polygon_PA_buffer.kml", sep=""))
plotKML::kml(final_PA_buffer, file.name = paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/kml_files/pumila_final_PA_buffer.kml", sep=""))#The cutted parts match!


##Remove see areas of the buffer and save as a raster
#We want that the polygon that serves as model for the buffers will be exact the extended distribution of pine species, i.e. the extended distribution with the 1 width buffer but removing see areas. In that way when we apply a buffer of 20 width, is exactly around the extended distribution of the species (not counting see areas).
#first rasterize the buffer
###create a raster from buffer+bianca area converting the buffer+bianca polygon in a raster (rasterize) 
raster_buffer_PA = raster() #create a empty raster
extent(raster_buffer_PA) = extent(final_PA_buffer) #give it a extention similar to distribution area of pinus halepensis polygon
res(raster_buffer_PA) = res(distribution_pumila) #give it a resolution similar to he biancas maps because the polygon tha we will rasterize was obtained from this maps (I can`t extract resolution from polygon_buffer). In this way, we maintenain the same resolution. 
raster_buffer_PA  = rasterize(final_PA_buffer,raster_buffer_PA) #convert the polygon of the halepensis buffer in a raster
distribution_pumila #take a raster of pumila distribution (without buffer) for obtaining coast limits
raster_buffer_PA = distribution_pumila*raster_buffer_PA #multiply both raster for creating a raster with NAs as much as these rasters, it is to say, we want a raster with NA in pixel where some of these raster have a NA. The zones without cells in distribution are places of sea only, if we would have used a environmental variable, terrestrial cells without data of the variable also would be computed as sea cells. 
raster_buffer_PA[!is.na(raster_buffer_PA)] <- 1 #give value of 1 to all cells without NA, in this way include cell of distribution raster with 1 and cell outside the distribution raster but inside the raster buffer (with 0).  

#create the new polygon halepensis distribution without sea areas
final_PA_buffer = rasterToPolygons(raster_buffer_PA, fun=function(x){x==1}, dissolve=TRUE)

#plot
plot(crop(raster_buffer_PA, final_PA_buffer))
plot(final_PA_buffer, add=T)#as you can see, the raster and the polygon are almost the same, except in the joints (small difference) and of course the sea areas. The buffers for variables selection are also calculated with rasters of the PA buffers.

#compare with the initial distribution
plot(crop(distribution_pumila, final_PA_buffer))
plot(final_PA_buffer, add=T)

#save raster of distribution buffer
writeRaster(raster_buffer_PA, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences/pumila_PA_buffer.asc", format="ascii", overwrite=TRUE)

##save the final polygon as a shape file
#for that we need to covert the SpatialPolygons to a SpatialPolygonsDataFrame object. We have to create an empty data.frame with same number row than ID values we have in the polygon (1 in our case, we don't have IDs). See "https://stat.ethz.ch/pipermail/r-sig-geo/2006-August/001252.html"
IDs <- sapply(slot(final_PA_buffer, "polygons"), function(x) slot(x, "ID"))
df <- data.frame(rep(0, length(IDs)), row.names=IDs)
final_PA_buffer_poly_df <- SpatialPolygonsDataFrame(final_PA_buffer, df)#df: object of class ‘data.frame’; the number of rows in ‘data’ should equal the number of Polygons-class objects in ‘Sr’. We only have one polygon
final_PA_buffer@polygons[[1]]
class(final_PA_buffer_poly_df)
summary(final_PA_buffer_poly_df)
#write it
require(rgdal)
writeOGR(final_PA_buffer_poly_df, layer = 'PA_buffer_pumila', '/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/final_pumila_buffers/', driver="ESRI Shapefile", overwrite_layer=TRUE)
#check
final_PA_buffer_poly_df_check = readOGR("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/final_pumila_buffers/PA_buffer_pumila.shp")
plot(distribution_pumila)
plot(polygon_PA_buffer_pumila, add=T, lwd=4)
plot(final_PA_buffer_poly_df_check, add=T, border="red")
plotKML::kml(final_PA_buffer_poly_df_check, file.name = paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/kml_files/pumila_final_PA_buffer_poly_df_check.kml", sep=""))#perfect



##################################
####### Range calc buffer ########
##################################

#Create PA buffer around extended distribution
polygon_range_calc_pumila = gBuffer(polygon_distribution_buffer_pumila, byid=FALSE, id=NULL, width=12.5, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

#plot the distribution and the buffer
plot(distribution_pumila)
plot(polygon_range_calc_pumila, add=T)

#crop the range_calc using the distribution map, in that way we remove the areas outside the map
polygon_range_calc_pumila_cropped = crop(polygon_range_calc_pumila, distribution_pumila)
plot(polygon_range_calc_pumila, lty=1)
plot(polygon_range_calc_pumila_cropped, add=T, lty=2)
legend("bottomright", legend=c("complete calc buffer", "reduced calc buffer"), lty=c(1,2))

#calculate the difference between both buffers
extreme_range_calc = gDifference(polygon_range_calc_pumila, polygon_range_calc_pumila_cropped)

#plot polyogns
plot(polygon_range_calc_pumila)
plot(extreme_range_calc, add=T, col="red")
legend("bottomright", legend=c("complete calc buffer", "extreme calc buffer"), fill=c("white", "red"))

#plot the both polygons on the distribution map
plot(distribution_pumila)
plot(polygon_range_calc_pumila, add=T, lwd=4)
plot(extreme_range_calc, add=T, border="red")#the right extreme is not plotted, we have to changed the coordinates
legend("bottomright", legend=c("complete calc buffer", "extreme calc buffer"), fill=c("black", "red"))

#change the coordinates of the right extreme polygon. In this case we only have one polgyon (the buffer is not big enough to get to the up end of the map).
extreme_range_calc@polygons[[1]]@Polygons[[1]]@coords[,1]#these are the long (x) coordinates of the second polygon
#we substract from each coordinate 180, the result will be substracted from 180 and the changed the sign
extreme_range_calc@polygons[[1]]@Polygons[[1]]@coords[,1] <- -(180 - (extreme_range_calc@polygons[[1]]@Polygons[[1]]@coords[,1] - 180))#Example to understand: The maximum x coordinate is 180, so all the points in the second polygon from 180 to 199 should be in the other extreme of the map. Imagine a point at 183 longitude degrees. The difference from 180, that it is 3 (183-180), should be on the other side, i.e. should be rested to 180 -> 177. Therefore: - (180 - (183-180)) -> -177. With 180 would be: - (180 - (180-180)) = -180, which is ok because we are exactly in the end of one size of the map (180) and beginning of the other side (-180).

#plot
plot(distribution_pumila)
plot(extreme_range_calc, add=T, border="red")
plot(polygon_range_calc_pumila, add=T, lwd=4)
legend("bottomright", legend=c("calc buffer inside the map", "the rest of calc buffer"), fill=c("black", "red"))

#bind both polygons
final_polygon_range_calc = union(extreme_range_calc, polygon_range_calc_pumila)

#remove the areas outside the map. Up for both polygons, left only for the polygon_range_calc_pumila The left part of the  polygon_range_calc_pumila has been put in the other side, whilst the up part is not interesting because there is no earth in that
final_polygon_range_calc = crop(final_polygon_range_calc, distribution_pumila)

#plot
plot(distribution_pumila)
plot(polygon_range_calc_pumila, add=T, lwd=4)
plot(final_polygon_range_calc, add=T, border="red")
legend("bottomright", legend=c("initial PA buffer", "final_PA_buffer"), fill=c("black", "red"))

#save the buffers into kml files for opening in goole earth
plotKML::kml(polygon_range_calc_pumila, file.name = paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/kml_files/pumila_polygon_range_calc.kml", sep=""))
plotKML::kml(final_polygon_range_calc, file.name = paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/kml_files/pumila_final_polygon_range_calc.kml", sep="")) #The cutted parts match!

##Remove see areas of the buffer and save as a raster
#We want that the polygon that serves as model for the buffers will be exact the extended distribution of pine species, i.e. the extended distribution with the 1 width buffer but removing see areas. In that way when we apply a buffer of 20 width, is exactly around the extended distribution of the species (not counting see areas).
#first rasterize the buffer
###create a raster from buffer+bianca area converting the buffer+bianca polygon in a raster (rasterize) 
raster_buffer_range_calc = raster() #create a empty raster
extent(raster_buffer_range_calc) = extent(final_polygon_range_calc) #give it a extention similar to distribution area of pinus halepensis polygon
res(raster_buffer_range_calc) = res(distribution_pumila) #give it a resolution similar to he biancas maps because the polygon tha we will rasterize was obtained from this maps (I can`t extract resolution from polygon_buffer). In this way, we maintenain the same resolution. 
raster_buffer_range_calc  = rasterize(final_polygon_range_calc,raster_buffer_range_calc) #convert the polygon of the halepensis buffer in a raster
distribution_pumila #take a raster of pumila distribution (without buffer) for obtaining coast limits
raster_buffer_range_calc = distribution_pumila*raster_buffer_range_calc #multiply both raster for creating a raster with NAs as much as these rasters, it is to say, we want a raster with NA in pixel where some of these raster have a NA. The zones without cells in distribution are places of sea only, if we would have used a environmental variable, terrestrial cells without data of the variable also would be computed as sea cells. 
raster_buffer_range_calc[!is.na(raster_buffer_range_calc)] <- 1 #give value of 1 to all cells without NA, in this way include cell of distribution raster with 1 and cell outside the distribution raster but inside the raster buffer (with 0).  

#create the new polygon halepensis distribution without sea areas
final_polygon_range_calc = rasterToPolygons(raster_buffer_range_calc, fun=function(x){x==1}, dissolve=TRUE)

#plot
plot(crop(raster_buffer_range_calc, final_polygon_range_calc))
plot(final_polygon_range_calc, add=T)#as you can see, the raster and the polygon are almost the same, except in the joints (small difference) and of course the sea areas. The buffers for variables selection are also calculated with rasters of the PA buffers.

#compare with the initial distribution
plot(crop(distribution_pumila, final_polygon_range_calc))
plot(final_polygon_range_calc, add=T)

#save raster of distribution buffer
writeRaster(raster_buffer_range_calc, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/buffers_calc_ranges/pumila_buffer_range_calc.asc", format="ascii", overwrite=TRUE)

#save the final polygon as a shape file
#for that we need to covert the SpatialPolygons to a SpatialPolygonsDataFrame object. We have to create an empty data.frame with same number row than ID values we have in the polygon (1 in our case, we don't have IDs). See "https://stat.ethz.ch/pipermail/r-sig-geo/2006-August/001252.html"
IDs <- sapply(slot(final_polygon_range_calc, "polygons"), function(x) slot(x, "ID"))
df <- data.frame(rep(0, length(IDs)), row.names=IDs)
final_range_calc_poly_df <- SpatialPolygonsDataFrame(final_polygon_range_calc, df)#df: object of class ‘data.frame’; the number of rows in ‘data’ should equal the number of Polygons-class objects in ‘Sr’. We only have one polygon
final_polygon_range_calc@polygons[[1]]
class(final_range_calc_poly_df)
summary(final_range_calc_poly_df)
#write it
require(rgdal)
writeOGR(final_range_calc_poly_df, layer = 'range_calc_pumila', '/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/final_pumila_buffers/', driver="ESRI Shapefile", overwrite_layer=TRUE)
#check
final_range_calc_poly_df_check = readOGR("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/final_pumila_buffers/range_calc_pumila.shp")
plot(distribution_pumila)
plot(polygon_range_calc_pumila, add=T, lwd=4)
plot(final_range_calc_poly_df_check, add=T, border="red")
plotKML::kml(final_range_calc_poly_df_check, file.name = paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/kml_files/pumila_final_range_calc_poly_df_check.kml", sep=""))#perfect