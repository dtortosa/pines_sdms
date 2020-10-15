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
ratio_PA_occurrence_df = data.frame(selected_species=NA, n_row_PAs=NA, n_row_occurr=NA, ratio_PA_occurrence=NA)
test_PA_weight = NULL
for (i in 1:length(epithet_species_list)){

    #select the [i] species
    selected_species = epithet_species_list[i]

    #load complete presences
    complete_presences = read.table(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences", paste(selected_species, "complete.presences.csv", sep="_"),sep="/"), sep=",", header=T)

    #extract number of rows with PAs and occurrences
    n_row_PAs = nrow(complete_presences[which(complete_presences$presence == 0),])
    n_row_occurr = nrow(complete_presences[which(complete_presences$presence == 1),])

    #calculate the ratio
    ratio_PA_occurrence = n_row_PAs/n_row_occurr

    #save it
    ratio_PA_occurrence_df = rbind.data.frame(ratio_PA_occurrence_df, cbind.data.frame(selected_species, n_row_PAs, n_row_occurr, ratio_PA_occurrence))

    #calculate PA weight with the final data
    calculate_PA_weight = sum(complete_presences[which(complete_presences$presence == 1),]$precision_weight) / n_row_PAs

    #extract the PA weight calculate when the PAs were calculated
    current_PA_weight = unique(complete_presences[which(complete_presences$presence == 0),]$precision_weight)

    #check that both are equal and save
    test_PA_weight = append(test_PA_weight, round(calculate_PA_weight, 4) == round(current_PA_weight, 4))
}

#remove the first row with NA from ratio_PA_occurrence_df
ratio_PA_occurrence_df = ratio_PA_occurrence_df[-1,]
#check that you have 112 nrows
nrow(ratio_PA_occurrence_df) == 112
#take a look
ratio_PA_occurrence_df

#see speceis whose PA/occurrence ratio is not 10
ratio_PA_occurrence_df[ratio_PA_occurrence_df$ratio_PA_occurrence != 10,]#No problem, all speices with less than 120 occurrences, so it possible that n_occrruence*10 does not reach to 30 PAs per strata in these species

#check the test
summary(test_PA_weight)#ALL true. The PA weight is ok taking into account all the dataset, but this does not take into account the partition of the data!! This will be done in model_loop_v4 (or higher) in Rafa pro.