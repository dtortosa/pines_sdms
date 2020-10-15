#Code for make comprobations in relaton distribution maps created by Bianca two years ago.

###definimos el directorio de trabajo
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus")

#required packages
require(raster)
require(rgeos)

#### Notes about maps archives #####
##ALL the chage will made in my folder of the Bianca`s data (/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS), we will not change the folder with data of Bianca without changes. 
##Bianca's thesis (https://www.wsl.ch/staff/niklaus.zimmermann/research/Thesis_Saladin_2013.pdf)
##critfield and little maps: "https://archive.org/details/geographicdistri991crit"

#p_tabuliformis have two files, p_tabuliformis.img and p_tabuliformis.img TRUE . gri. The second has a litle bit more area and it was called as TRUE, thus we will select this and delete the other. 
par(mfcol=c(1,2))
plot(crop(raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/datos_brutos/Maps/p_tabulaeformis_01.img"), c(60,130,10,55)), main="p_tabuliformis.img")
plot(crop(raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/datos_brutos/Maps/p_tabulaeformis_01.img TRUE.gri"), c(60,130,10,55)), main="p_tabulaeformis_01.img TRUE.gri")

#actualization!!! 26/04/18: There was a problem with this species, the area included in the second page of this species was in other proijection , thus areas in this page would be difficult to locate. Bianca has dropped in the new map,  I have changd its maps, so "/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_tabuliformis_01.img" is the correct now, and "p_tabuliformis_01_error.img" is the ancient with the error. 

#p_sylvestris have two files: p_sylvestris.img and p_sylvestris_michale.img. The difference between them is that the second include area of P. sibirica as P.sylvesitrs, on the contrary of the first file. Because we are considering sibirica as a different species of P. sylvestris, we will use the first file. 
    #In addition, we have "p_sylvestris_01_sin_siberia.img", a raster of sylvestris distribution without the area of siberia, this raster is wrong, the correct is "p_sylvestris_01.img", and it has been used in all analyses. 
par(mfcol=c(1,2))
plot(crop(raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_sylvestris_01_sin_siberia.img"), c(-13.18, 140, 30, 75)), main="Sylvestris without Siberia")
plot(crop(raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_sylvestris_01.img"), c(-13.18, 140, 30, 75)), main="Sylvestris with Siberia")
dev.off()

#p_devoniana is called as p_devoniana.michoacana_01.img. La verdad es que no lo encuentro en los mapas de critfield, but the distribution is similar to the Bianca's thesis and wikipedia maps ("https://en.wikipedia.org/wiki/Pinus_devoniana")
plot(raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_devoniana_01.img"), c(-120,-80, 0,50))

#p_kesiya is called as p_kesiya.insularis_01.img. But the distribution is similar to the Bianca's thesis. 
plot(raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_kesiya_01.img"))

#p_muricata is called as p_murricata_01.img. It is an error. But the distribution is similar to the Bianca's thesis. 
plot(raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_muricata_01.img"))

#p_tabuliformis is called as p_tabulaeformis_01.img. It is an error. But the distribution is similar to the Bianca's thesis. 
plot(raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_tabuliformis_01.img"))

#p_wallichiana is called as p_wallichiana.griffithii_01.img. But the distribution is similar to the Bianca's thesis. 
plot(raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_wallichiana_01.img"))


#############LOAD  REVISION DATA
#revison done species per species showed in the critfield maps: "https://archive.org/stream/geographicdistri991crit"

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

#load revision data
revision = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/distribution/revised_pines.txt", sep=",", header=TRUE)
str(revision) #epithet; included or not in critchfield maps; there is problem?; notes

#all species from revision included in epithet_species_list
length(which(revision$epithet_species %in% epithet_species_list)) == nrow(revision) #discolor left

#all species analysed
length(which(revision$epithet_species %in% epithet_species_list)) == length(epithet_species_list)

#searh notes for a given species
revision$note[which(revision$epithet_species=="discolor")] #e.g. for discolor

#species do not included in critfield and species included but with problems
species_not_included = revision$epithet_species[which(revision$critfield == "NO")]

#plot cropped distribution of all species
for (i in 1:length(species_not_included)){

    #select the [i] species
    selected_species = species_not_included[i]

    #Cargamos area distribucion
    distribution = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(selected_species, "01.img", sep="_") ,sep="_")) #select the path of the distribution file of the corresponding species
    
    #create a polygon of the dsitribution
    polygon = rasterToPolygons(distribution, fun=function(x){x==1}) #convertimos en poligono the raster using the cell with values=1

    #create a buffer
    polygon_buffer = gBuffer(polygon, byid=FALSE, width=15) #aumentamos el area del poligono for seeing well the distribution 

    #crop the distribution using the polygon (more detail)
    distribution_cropped = crop(distribution, polygon_buffer)

    #plot
    png(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/distribution_pines_cropped/problematic_species/distribution_pines_cropped_", selected_species, ".png", sep=""))
    plot(distribution_cropped, main=paste("Pinus", selected_species, sep=" "))
    dev.off()
}    

#bind all png file into one single pdf using imagemagick (convert command)
system("cd /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/datos/raw_ocurrences/distribution_pines_cropped/problematic_species; convert *.png full.pdf")

#final problematic species (after manual revision)
final_problematic_species = revision$epithet_species[which(revision$problem == "YES")] #In my opinion, the only important case is discolor, which distribution is identicol to cembroides, it is a pseudoreplic and we have to remove it

#bind all png file of FINAL problematic (made at hand the selection) into one single pdf using imagemagick (convert command)
system("
    cd /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/datos/raw_ocurrences/final_problematic_species;
    rm full.pdf; 
    convert * full.pdf")


#####change of BIANCA
species_solved = c("amamiana", "aristata", "bungeana", "dalatensis", "densiflora", "douglasiana", "fenzeliana", "hwangshanensis", "kesiya", "latteri", "luchuensis", "maximartinezii", "merkusii", "monophylla", "patula", "peuce", "squamata", "sylvestris", "maximinoi", "densata", "kwangtungensis")
#plot cropped distribution of all species
for (i in 1:length(species_solved)){

    #select the [i] species
    selected_species = species_solved[i]

    #Cargamos area distribucion buena y mala
    old_distribution = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(selected_species, "01_error.img", sep="_") ,sep="_")) #select the path of the distribution file of the corresponding species
    new_distribution = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(selected_species, "01.img", sep="_") ,sep="_")) #select the path of the distribution file of the corresponding species    

    #if the [i] species is not sylvestris create a buffer and crop (with sylvestris is very slow and in general changes can be seen without cropping)
    if(!selected_species == "sylvestris"){
        #create a polygon of the dsitribution
        polygon = rasterToPolygons(new_distribution, fun=function(x){x==1}) #convertimos en poligono the raster using the cell with values=1

        #create a buffer
        polygon_buffer = gBuffer(polygon, byid=FALSE, width=15) #aumentamos el area del poligono for seeing well the distribution 

        #crop the distribution using the polygon (more detail)
        old_distribution_cropped = crop(old_distribution, polygon_buffer)
        new_distribution_cropped = crop(new_distribution, polygon_buffer)
    } else {
        #if sylvestris, do nothing
        old_distribution_cropped = old_distribution
        new_distribution_cropped = new_distribution
    }

    #plot
    png(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/distribution_pines_cropped/problematic_species/comparison_solutions_bianca/comparison_old_new_", selected_species, ".png", sep=""))
    par(mfrow=c(2,1))
    plot(old_distribution_cropped, main=paste("Pinus", selected_species, "old", sep=" "))
    plot(new_distribution_cropped, main=paste("Pinus", selected_species, "new", sep=" "))    
    dev.off()
} 

####new speceis from BIANCA
new_species = c("tecunumanii", "jaliscana")
#plot cropped distribution of all species
for (i in 1:length(new_species)){

    #select the [i] species
    selected_species = new_species[i]

    #Cargamos area distribucion buena nueva
    new_distribution = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(selected_species, "01.img", sep="_") ,sep="_")) #select the path of the distribution file of the corresponding species    

    #create a polygon of the dsitribution
    polygon = rasterToPolygons(new_distribution, fun=function(x){x==1}) #convertimos en poligono the raster using the cell with values=1

    #create a buffer
    polygon_buffer = gBuffer(polygon, byid=FALSE, width=15) #aumentamos el area del poligono for seeing well the distribution 

    #crop the distribution using the polygon (more detail)
    new_distribution_cropped = crop(new_distribution, polygon_buffer)

    #plot
    png(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/distribution_pines_cropped/problematic_species/new_species/new_species", selected_species, ".png", sep=""))
    plot(new_distribution_cropped, main=paste("Pinus", selected_species, "new", sep=" "))    
    dev.off()
}