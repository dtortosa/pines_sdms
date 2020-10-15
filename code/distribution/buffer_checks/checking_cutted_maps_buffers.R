####Code for plotting and checking that the buffer of occurrences, global figures and pseudoabsences do not reach the end of the maps. 

#required packages
require(raster) #for work with rasters
require(rgeos) #for creating the buffer and the centroids of the cells without gbif points

#load list of species
list_species = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/species.txt", sep="\t", header=T)

#extract epithet from species list
epithet_species_list = NULL
for(i in 1:nrow(list_species)){

    #selected species
    selected_species = as.vector(list_species[i,])

    #extract epithet
    epithet_species_list = append(epithet_species_list, strsplit(selected_species, split=" ")[[1]][2])
}
summary(is.na(epithet_species_list)) #all false
#check
require(tidyverse)
paste("Pinus", epithet_species_list, sep=" ") == str_trim(as.vector(list_species[,1])) #the seperated epithet more Pinus are equal to list_species? We use str_trim from tidyverse to remove end spaces in each element of list_species

#remove tecunumanii, jaliscana y discolor. These species are not used for Niche paper. The two first because we should downloadad gbif data now, so we would mix gbif data form 2016 and 2019. The third one was impossible to differentiate from P. cembroides
epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("tecunumanii", "jaliscana", "discolor"))]
#check
c("tecunumanii", "jaliscana", "discolor") %in% epithet_species_list



###############################################################################
####### CHECK THAT VARIATION BECAUSE OF SHAPE OF EARTH DOES NOT AFFECT ########
###############################################################################
#Earth is an ellipsoid. This can produce alterations in the calculation of buffers (see these links for further information "https://stackoverflow.com/questions/25411251/buffer-geospatial-points-in-r-with-gbuffer" "https://www.esri.com/news/arcuser/0111/geodesic.html"). These variations can produce that the same buffer would be altered when projected in different part of the planet. This is specially relevant when the buffer enconpases ver large areas in different UTM zones See "https://seethedatablog.wordpress.com/2017/08/03/euclidean-vs-geodesic-buffering-in-r/".

#I'm going to check that species with big distributions at high and low latitudes.
species_check_ellipsoid = c("contorta", "pseudostrobus", "merkusii", "sylvestris")

#for each species
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/check_buffer_ellipsoid.pdf")
for(i in 1:length(species_check_ellipsoid)){

    #select the [i] species
    selected_species = species_check_ellipsoid[i]

    #load the distribution maps of [i] species
    distribution = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_", selected_species, "_01.img", sep=""))    

    #add the csr data in the case of sylvestris and merkussi. It seems that Bianca did not add thes einfo when modified the map. The previous maps (previous with errors) have this information. In addition, all Bianca maps tend o have similar datum, I have checed species form distance latitudes and longitudes
    if(selected_species %in% c("merkusii", "sylvestris")){
 
        #add the datum to the new sylvestris map
        proj4string(distribution) <- proj4string(raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_", selected_species, "_01_error.img", sep="")))
    }

    #create a polygon from distribution raster
    polygon_distribution = rasterToPolygons(distribution, fun=function(x){x==1}, dissolve=TRUE)

    #Create a buffer around the current distribution
    polygon_distribution_buffer = gBuffer(polygon_distribution, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

    #Create another buffer around the extended distribution. This will be for calculating range loss and range
    polygon_range_calcs_buffer = gBuffer(polygon_distribution_buffer, byid=FALSE, id=NULL, width=12.5, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

    #Create another buffer around the extended distribution. This will be for calculating PAs
    polygon_PA_buffer = gBuffer(polygon_distribution_buffer, byid=FALSE, id=NULL, width=22.5, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

    #save the buffer into kml files for opening in goole earth
    plotKML::kml(polygon_distribution, file.name = paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/kml_files/", selected_species, "_polygon_distribution.kml", sep=""))
    plotKML::kml(polygon_distribution_buffer, file.name = paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/kml_files/", selected_species, "_polygon_distribution_buffer.kml", sep=""))
    plotKML::kml(polygon_range_calcs_buffer, file.name = paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/kml_files/", selected_species, "_polygon_range_calcs_buffer.kml", sep=""))
    plotKML::kml(polygon_PA_buffer, file.name = paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/kml_files/", selected_species, "_polygon_PA_buffer.kml", sep=""))

    #plot
    plot(crop(distribution, polygon_PA_buffer), main=selected_species)
    plot(polygon_distribution_buffer, add=T)
    plot(polygon_range_calcs_buffer, add=T)
    plot(polygon_PA_buffer, add=T)
}
dev.off()

#El error que da gbuffer se refiere a que le estamos dando poligonos proyectados sobre latitud longitud, cuando espera una projecton sobre el globo terraqueo. Eso puede dar lugar a que se hagan buffers más pequeños de lo esperado en los polos ("https://stackoverflow.com/questions/25411251/buffer-geospatial-points-in-r-with-gbuffer")

#LAS UNIDADES DE LOS BUFFERS SON GRADOS DECIMALES!!!

#En general las lineas pasan exactamente por donde se indnca en el mapa plano, el porblema es que como la tierra es esferica, a latitudes elevadas el buffer es más pequeño. Se reocge lo que se dice en el mapa plano, pro los buffers son menos anchos a latitudes elevadas. Mira "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/contorta_buffers.pdf"
    #hay tres buffers: i) El de distribución, el más pequeño (1º; ~100 km en el ecuador) que extiende la distribución original de Critchfield; ii) El mediano (12.5º; ~1250 km en el ecuador), dentro del cual se calcularán los cambios de idoneidad entre las condiciones actuales y futuras; iii) El más grande (22.5º; ~2250 km en el ecuador), dentro del cual se seleccionan pseudoausencias. 

#Pasa exactamente lo mismo con las celdas de bianca: Por ejemplo, 5 celdas en la parte sur de la distribución de Contorta (Colorado) suman 200 km, mientras que las mismas 5 celdas en Alaska suman 127 km. Esto me hace reafirmarme en la idea de que los suizos son conscientes de esto pero les da un poco igual dada la escala a la que estamos trabajando. 5 celdas puestas una al lado de la otra suman 200 y no 250 km, porque no son exactamente 50x50, son 50 de largo pero ~43 ancho (en el ecuador).

#Hemos expandido un poco los dos buffers grandes (10 y 20 a 12.5 y 22.5) para reducir la influencia del estrechamiento al norte y meter así más areas relevantes, pero sin colarnos y meter areas poco poco releveantes. De esta forma para contorata entra casi toda alaska en el buffer de PAs y sin embargo sylvestris solo toca de refilón Groenlandia. Puede caer alguna PAs, pero la probabilidad es baja, a sbemos que el peso más grande de los modelos está en las presencias (ten en cuenta que antes se metía TODA groenlandia..). Parece un buen compromiso
    #Puede que para algunas especies mu pequeñas como maestrensis pueda faltar un pooc, en ese caso casi no se llelga al desierto de chihuahan, pero esto es un problema insitriseco de la falta de distribución de estas especies, y que no puede afectar al resto. Ese el mejor compromiso

#Puedes mirar estos buffers en google earth abriendo los archivos kml

##################################################
####### ATTEMPT TO CREATE A MANUAL BUFFER ########
##################################################
#I have tried to create a buffer add cels around the bianca's maps (with the function boundary). However the results were not very good, because as I said above, the cells of bianca's map also modify across latitude. 

#for each species
if(FALSE){
for(i in 1:length(epithet_species_list)){

    #select the [i] species
    selected_species = epithet_species_list[i]

    #load the distribution maps of [i] species
    distribution = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_", selected_species, "_01.img", sep=""))    

    #create a polygon from distribution raster
    polygon_distribution = rasterToPolygons(distribution, fun=function(x){x==1}, dissolve=TRUE)

    #we calculate several buffers using as reference the previous map calculated. In that way the size of the map increases
    boundary_rasters = list()
    boundary_polygons = list()    
    for(j in 1:5){

        #if this is the firs map 
        if(j == 1){
            #the reference is the distribution                
            reference_raster = distribution

            #remove areas with NA (see areas)
            reference_raster[which(is.na(getValues(reference_raster)))] <- 0

        } else {
            #if not the reference map is the maps previously calculated (j-1 in the list)
            reference_raster = boundary_rasters[[j-1]]
        }

        #calculate the boundary in the reference maps (cells between 0 and 1)
        boundary_distribution = boundaries(reference_raster, classes=TRUE, directions=8) 

        #all the areas as 1 in the reference maps should be included (we are increasing the size)
        boundary_distribution[which(getValues(reference_raster) == 1)] <- 1 

        #convert to polygon
        boundary_polygon = rasterToPolygons(boundary_distribution, fun=function(x){x==1}, dissolve=TRUE)

        #plot 
        plot(crop(boundary_distribution, boundary_polygon))
        plot(boundary_polygon, add=T)
        plot(polygon_distribution, add=T)

        #save the raster and polygon of the buffer
        boundary_rasters[[j]] =  boundary_distribution
        boundary_polygons[[j]] = boundary_polygon
    }
}
} #It does not work. As you can see in "/Users/diegosalazar/Desktop/attemp_manual_buffer.png", the buffer of contorar is reduced at higher latitudes. This is expected given the variation across latitude of Bianca's cells.
    
    #Además, en los casso con una sola celda que no es distributon pero que queda rodeada, no se termina incluyendo nunca y se queda el hueco ahí...

    #FAIL

##################################################################
####### CHECK THAT BUFFER DO NOT REACH THE END OF THE MAP ########
##################################################################


#plot the different buffers: Info about buffer size inside the loop
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/buffer_checks/plot_buffers.pdf")
for(i in 1:length(epithet_species_list)){

    #select the [i] species
    selected_species = epithet_species_list[i]

    #load the distribution maps of [i] species
    distribution = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_", selected_species, "_01.img", sep=""))    

    #create a polygon from distribution raster
    polygon_distribution = rasterToPolygons(distribution, fun=function(x){x==1}, dissolve=TRUE)
    
    #Create a buffer around the current distribution
    polygon_distribution_buffer = gBuffer(polygon_distribution, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

    #Create another buffer around the extended distribution. This will be for calculating range loss and range
    polygon_range_calcs_buffer = gBuffer(polygon_distribution_buffer, byid=FALSE, id=NULL, width=12.5, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

    #Create another buffer around the extended distribution. This will be for calculating PAs
    polygon_PA_buffer = gBuffer(polygon_distribution_buffer, byid=FALSE, id=NULL, width=22.5, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

    #crop the distribution raster with the polygon PA buffer (the largest)
    distribution_crop =  crop(distribution, polygon_PA_buffer) 

    #plot all buffers
    par(mfcol=c(2,1), mai = c(0.5, 0.5, 0.4, 0.1))

    #first using the whole map
    plot(distribution, main=selected_species)
    plot(polygon_distribution_buffer, add=T, lty=1)
    plot(polygon_range_calcs_buffer, add=T, lty=2)
    plot(polygon_PA_buffer, add=T, lty=3)
    #legend
    legend("topright", legend=c("Distribution buffer (extended distribution; 1º)", "Range calculation buffer (12.5º)", "Pseudoabsences buffer (22.5º)"), lty=c(1,2,3), cex=0.5)

    #second, using the croped map
    plot(distribution_crop)
    plot(polygon_distribution_buffer, add=T, lty=1)
    plot(polygon_range_calcs_buffer, add=T, lty=2)
    plot(polygon_PA_buffer, add=T, lty=3)
}
dev.off()


##########################################################################
####### ESTIMATE THE DISTANCE IN KMS OF THE BUFFER AT THE EQUATOR ########
##########################################################################

#load distribution map of halepensis
distribution_halepensis = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_halepensis_01.img")

#create a polygon from distribution raster
polygon_distribution_halepensis = rasterToPolygons(distribution_halepensis, fun=function(x){x==1}, dissolve=TRUE)

#Create distribution buffer around
polygon_distribution_buffer_halepensis = gBuffer(polygon_distribution_halepensis, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

#Create PA buffer around extended distribution
polygon_PA_buffer_halepensis = gBuffer(polygon_distribution_buffer_halepensis, byid=FALSE, id=NULL, width=20, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

#plot the distribution buffer of halepensis and a smaller one of 0.5 (50 km)
plot(crop(distribution_halepensis, polygon_distribution_buffer_halepensis))
plot(gBuffer(polygon_distribution_halepensis, byid=FALSE, id=NULL, width=0.5, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0), add=T)
plot(gBuffer(polygon_distribution_halepensis, byid=FALSE, id=NULL, width=1, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0), add=T)
    #we can see that an increment of with in 0.5 is like a cell of 50 km, hence increments of width in 1 are equal to 100 km. The dimensions are 50x43 km, bercause cells are not perfect squared.

#plot the distribution of halepensis and add buffer of 100 km each time to 2000 km. The fist one is the distribution buffer (extended distribution) without any add. From there, we add the buffers
plot(crop(distribution_halepensis, polygon_PA_buffer_halepensis))
for(j in 0:20){

    plot(gBuffer(polygon_distribution_buffer_halepensis, byid=FALSE, id=NULL, width=j, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0), add=T)
} #width of 20 is equal to 20 buffers of width=1 (i.e. 100), hence 20*100=2000 km. A width of 10 would be equal to 1000 km. This is in EQUATOR!! BECAUSE OF THE SHAPE OF THE EARTH SEE SECTION "CHECK THAT VARIATION BECAUSE OF SHAPE OF EARTH DOES NOT AFFECT"