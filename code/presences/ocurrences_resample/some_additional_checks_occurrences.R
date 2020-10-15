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
#check
require(tidyverse)
paste("Pinus", epithet_species_list, sep=" ") == str_trim(as.vector(list_species[,1])) #the seperated epithet more Pinus are equal to list_species? We use str_trim from tidyverse to remove end spaces in each element of list_species


#remove tecunumanii, jaliscana y discolor. These species are not used for Niche paper. The two first because we should downloadad gbif data now, so we would mix gbif data form 2016 and 2019. The third one was impossible to differentiate from P. cembroides
epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("tecunumanii", "jaliscana", "discolor"))]
#check
c("tecunumanii", "jaliscana", "discolor") %in% epithet_species_list


#make the loop for calculating the number of occurrences per species and checking that all species have lon-lat indicated in that name
test_lon_lat = NULL
ocurrences_per_species = data.frame(species=NA, number_ocurrences=NA, n_unknown_basisOfRecord=NA)
for (i in 1:length(epithet_species_list)){

    #select the [i] species
    species = epithet_species_list[i]

    #load presencias raw
    presencia_raw<-read.csv(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus", paste(species, "csv", sep="."), sep="_"), header=T, fill=TRUE, check.names=TRUE, stringsAsFactors=FALSE) 

    #check that lon-lat are indicated with these names
    test_lon_lat = append(test_lon_lat, length(which(colnames(presencia_raw) %in% c("lon", "lat"))) == 2)

    #calculate the number of occurrences with basis of record = "unknown"
    n_unknown_basisOfRecord = length(which(presencia_raw$basisOfRecord == "UNKNOWN"))

    #load the presences of the [i] species
    presences = read.table(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences", paste(species, "final.presences.csv", sep="_"), sep="/"), sep=",", header=T)

    #extract the number of rows (i.e. the number of occurrences)
    number_ocurrences = nrow(presences)

    #save it
    ocurrences_per_species = rbind.data.frame(ocurrences_per_species, cbind.data.frame(species, number_ocurrences, n_unknown_basisOfRecord))
}

#remove the first row with NA from ocurrences_per_species
ocurrences_per_species = ocurrences_per_species[-1,]
#check that you have 112 nrows
nrow(ocurrences_per_species) == 112
#take a look
ocurrences_per_species

#calculate the ratio between presenceis with unknown basis of record and the total number of presences
ocurrences_per_species$unknown_presence_ratio = ocurrences_per_species$n_unknown_basisOfRecord / ocurrences_per_species$number_ocurrences

#save as csv
write.table(ocurrences_per_species, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/ocurrences_per_species.csv", sep=",", col.names = TRUE, row.names = FALSE)

#check that the test of long-lat is ok
summary(test_lon_lat)