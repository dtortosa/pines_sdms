#code for plotting distribution maps of Bianca for all species
require(raster)

#list species
list_species = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/species.txt", sep="\t", header=T)
str(list_species)
summary(list_species)

#extract epithet from species list
epithet_species_list = NULL
for(i in 1:nrow(list_species)){

    #selected species
    selected_species = as.vector(list_species[i,])

    #extract epithet
    epithet_species_list = append(epithet_species_list, strsplit(selected_species, split=" ")[[1]][2])
}
summary(is.na(epithet_species_list)) #all false

#####################################################################
########Plot distribution of all species ############################
#####################################################################

pdf("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/ocurrences/distribution_pines.pdf", width=12, height=12)
for (i in length(epithet_species_list)){

    #select the [i] species
    selected_species = epithet_species_list[i]

    #load the raster of the [i] species
    distribution = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(selected_species, "01.img", sep="_") ,sep="_"))
    
    #plot it
    plot(distribution, main=paste("Pinus", selected_species, sep=" ")) 
}
dev.off()