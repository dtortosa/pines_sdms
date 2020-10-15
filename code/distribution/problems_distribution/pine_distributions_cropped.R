#loop that cleans ocurrences and save them for modelling

###definimos el directorio de trabajo
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus")

#Librerias
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

#plot cropped distribution of all species
for (i in 1:length(epithet_species_list)){

    #select the [i] species
    selected_species = epithet_species_list[i]

    #Cargamos area distribucion
    distribution = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(selected_species, "01.img", sep="_") ,sep="_")) #select the path of the distribution file of the corresponding species
    
    #create a polygon of the dsitribution
    polygon = rasterToPolygons(distribution, fun=function(x){x==1}) #convertimos en poligono the raster using the cell with values=1

    #create a buffer
    polygon_buffer = gBuffer(polygon, byid=FALSE, width=15) #aumentamos el area del poligono en 2 celdas (width=1), según nos indico el test de las especies con datos de euforgen. byid determining if the function should be applied across subgeometries (TRUE) or the entire object (FALSE).  

    #crop the distribution using the polygon (more detail)
    distribution_cropped = crop(distribution, polygon_buffer)

    #plot
    png(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/distribution_pines_cropped/distribution_pines_cropped_", selected_species, ".png", sep=""))
    plot(distribution_cropped)
    dev.off()
}    

#bind all png file into one single pdf using imagemagick (convert command)
system("cd /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/datos/raw_ocurrences/distribution_pines_cropped; convert *.png full.pdf")
