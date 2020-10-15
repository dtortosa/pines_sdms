#set wroking directory
setwd("/Users/dsalazar/nicho_pinus/")

#require packages
require(raster)
require(rgeos)

#load bio1 for using it as a background
bio1 = raster("/Users/dsalazar/nicho_pinus/data/climate/finals/bio1.asc")
albicaulis_distribution = raster("/Users/dsalazar/nicho_pinus/data/MAPS/p_albicaulis_01.img") #load ablicaulis buffer to get resolution of distribution maps
bio1 = resample(bio1, albicaulis_distribution, method="bilinear") #reduce resolution. There is no problem with the resampling because cells with zero are actually zero, not like on phylo rasters in which NA cells were set as zero. This would artificially decrease the suitability of cells surrounding from false-zero cells. This is not the case. 
bio1[which(getValues(bio1) >= min(getValues(bio1), na.rm = TRUE))] <- 0 #set all continent areas as 0

#### problematic species
selected_epi = "pumila"

##### ocurrence buffer
#load [i] buffer as raster
ocurrences_buffer = raster(paste("/Users/dsalazar/nicho_pinus/data/buffers/", selected_epi, "_distribution_buffer", ".asc", sep=""))

#convert NAs into 0 to avoid problems in the sum
ocurrences_buffer[which(is.na(getValues(ocurrences_buffer)))] <- 0

#### migration buffer
#convert [i] buffer in a polygon 
ocurrences_buffer_polygon = rasterToPolygons(ocurrences_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon

#set species name as name of this raster
names(ocurrences_buffer) <- selected_epi

#create a migration buffer
polygon_migration_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=10, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) #byid=FALSE indicate that the function should be applied to the entire object or to sub-geometries

#plot migration buffer on bio1
pdf("/Users/dsalazar/nicho_pinus/data/problem_migration_pumila/migration_buffer_on_bio1.pdf")
plot(bio1)
plot(polygon_migration_buffer,add=T)
dev.off() #we lost the right extreme 

#select the area of migration buffer outside of the map
buffer_outside_map = crop(polygon_migration_buffer, c(xmax(bio1),xmax(polygon_migration_buffer),ymin(polygon_migration_buffer),ymax(polygon_migration_buffer)))
buffer_outside_map

#cehck that the extraction was correct
pdf("/Users/dsalazar/nicho_pinus/data/problem_migration_pumila/peak_on_migration_buffer.pdf")
plot(rbind(polygon_migration_buffer, buffer_outside_map, makeUniqueIDs = TRUE))
dev.off() #rbind bind two polygons: https://gis.stackexchange.com/questions/155328/merging-multiple-spatialpolygondataframes-into-1-spdf-in-r

#extract coordinates of buffer_outside_map
x_coords_buffer_outside_map = buffer_outside_map@polygons[[1]]@Polygons[[1]]@coords[,1] 
#convert these coordinate to other extreme of the map. If for example a point is at x=185, and the maximum of the map es 180, this 5 is the problem, and shoulb be in the other side of the map. So we sum that 5 to the minimun x of the map (-180) to get the real position, -175
new_coords=NULL
for(i in 1:length(x_coords_buffer_outside_map)){

    #select the [i] coordinate
    selected_coord = x_coords_buffer_outside_map[i]

    #calculate how much that point is otuisde of the map
    differ_x_bio1 = selected_coord-xmax(bio1)

    #add that difference to the minumu value of x for bio1. That is the newe coord
    new_coords = append(new_coords, xmin(bio1) + differ_x_bio1)
}


#add the new coord to the polygon
buffer_outside_map@polygons[[1]]@Polygons[[1]]@coords[,1] <- new_coords

#bind both polygons
final_polygon_migration_buffer = rbind(polygon_migration_buffer, buffer_outside_map, makeUniqueIDs = TRUE) #rbind bind two polygons: https://gis.stackexchange.com/questions/155328/merging-multiple-spatialpolygondataframes-into-1-spdf-in-r


#plot migration buffer on bio1
pdf("/Users/dsalazar/nicho_pinus/data/problem_migration_pumila/migration_buffer_on_bio1_corrected.pdf")
plot(bio1)
plot(final_polygon_migration_buffer,add=T)
dev.off() #we lost the right extreme 

#create a raster from that polygon
final_migration_buffer = bio1
final_migration_buffer[final_polygon_migration_buffer]<-1

#plot it
pdf("/Users/dsalazar/nicho_pinus/data/problem_migration_pumila/plot_final_migration_raster.pdf")
plot(final_migration_buffer)
plot(final_polygon_migration_buffer,add=T)
dev.off() #we lost the right extreme 

#save.it
writeRaster(final_migration_buffer, "/Users/dsalazar/nicho_pinus/data/problem_migration_pumila/final_migration_buffer_pumila.asc")

