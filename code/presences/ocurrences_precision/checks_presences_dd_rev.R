#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
    #https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
    #https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file




#########################################################################
########## SCRIPTS FOR CALCULATING PLOTS OF REGULATORY DENSITY ##########
#########################################################################

#In this script, we make a check about the comment 21 of the first revision in Diversity and Distribution.




###################################################
###### CHANGES RESPECT TO PREVIOUS VERSIONS #######
###################################################

#Version 1




##########################################
########## REMOVE PREVIOUS WORKSPACE #####
##########################################
remove(list=ls(all=TRUE))



##########################################
########## SET THE WORKING DIRECTORY #####
##########################################

#main folder nicho_pinus
setwd("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus")




#################################
###### REQUIRED PACKAGES ########
#################################

require(tidyverse) #for a check




#################################
###### EXTRACT PINE NAMES #######
#################################

#load list of species
list_species = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/code/presences/species.txt", sep="\t", header=T)

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
if(FALSE){
    require(tidyverse)
    paste("Pinus", epithet_species_list, sep=" ") == str_trim(as.vector(list_species[,1])) #the seperated epithet more Pinus are equal to list_species? We use str_trim from tidyverse to remove end spaces in each element of list_species
}#it is in false because loading tidyverse lead to load several packages that have a function name "extract", and this gives problems with the extract function of raster. If you want to check run these lines manually

#remove tecunumanii, jaliscana y discolor. These species are not used for Niche paper. The two first because we should downloadad gbif data now, so we would mix gbif data form 2016 and 2019. The third one was impossible to differentiate from P. cembroides
epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("tecunumanii", "jaliscana", "discolor"))]
#check
!TRUE %in% c("tecunumanii", "jaliscana", "discolor") %in% epithet_species_list




###################################################
###### CALCULATE HIGH PRECISION OCCURRENCES #######
###################################################

#make a loop for calculating the number of high precision points according to both precision variables in each species

#open an empty data.frame
number_high_precision_occurrences = data.frame(selected_species=NA, high_precision_presences_coord_uncertainty=NA, high_precision_presences_coord_precision=NA, check_1=NA, check_2=NA)

#open the loop
for(i in 1:length(epithet_species_list)){

    #select the [i] species
    selected_species = epithet_species_list[i]

    #load the occurrences
    species_presences <- read.csv(paste("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_", selected_species, ".csv", sep=""), header=TRUE, fill=TRUE, check.names=TRUE, stringsAsFactors=FALSE) 
        #stringsAsFactors=FALSE sirve para procesar mejor la columna coordinatePrecision, que lleva mezclados números y caracteres
        #we need the raw data to get coordinate precision/coordinate uncertainty. In the final data.frame obtained from loop_cleaning_occurrences, there are only 4 columns: longitude, latitude, precision_weight and presence.

    #if there is a column called coordinateUncertaintyInMeters
    if(length(which(colnames(species_presences) == "coordinateUncertaintyInMeters")) > 0){

        #Esta variable indica la distancia horizontal en metros desde la latitud decimal dada y la longitud decimal dada describring the smallest circle containing the whole of the location. Se deja empty si la incertidumbre es desconocida, no puede ser estimada o no es aplicable (porque no hay coordenadas). Zero no es valor valido para este termino. Debe ser mayor y diferente de cero además de menor de 5000000 (5000 km). 
        #points with precision lower than 4 km (4000 m) should receive a weight of 1, whilst point with coarser precision will receive a weight of 0.5.
            #info obtained from loop_cleaning_occurrences_v2.R
    
        #select those presences with coordinateUncertaintyInMeters lower than 4000
        subset_high_precision_presences_coord_uncertainty = species_presences[which(species_presences$coordinateUncertaintyInMeters < 4000),]

        #check that the subset worked well
        check_1 = length(which(subset_high_precision_presences_coord_uncertainty$coordinateUncertaintyInMeters > 4000)) == 0

        #extract the number of high precision points according to this variable
        high_precision_presences_coord_uncertainty = nrow(subset_high_precision_presences_coord_uncertainty)
    } else { #if not

        #no value for this variable exist for the [i] species
        high_precision_presences_coord_uncertainty = NA
    }

    #if there is a column called coordinatePrecision
    if(length(which(colnames(species_presences) == "coordinatePrecision")) > 0){

        #Hasta abril de 2016 este termino era el radio de error que tiene cada coordenada. La unidad es en metros, por tanto un valor de 25 indica un radio de 25x25 metros, y eso es lo que tenemos nosotros. 
        #Desde abril de 2016 ha cambiado. Una representación lineal de la precisión de las coordenadas dada en latitud y longitud decimal. Debe estar entre 0 y 1. Nosotros NO tenemos esto, por que la peña sigue subiendo los datos para esta variable como metros. 
        #Mis variables tenian una resolucion de 4x4 km  y aqui tenemos valores por debajo (1, 10, 25 = 1.85 km) y algunos por encima (0 =  18 km), los calculos están explicados en niche_questions.txt. Points with 1,10,25 will receive a weight of 1, whilst points with 0 will receive a weight of 0.5.             
            #info obtained from loop_cleaning_occurrences_v2.R
    
        #select those presences with coordinatePrecision >= 1 AND <= 25
        subset_high_coord_precision = species_presences[which(species_presences$coordinatePrecision>=1 & species_presences$coordinatePrecision<=25),]

        #check that the subset worked well
        check_2 = length(which(subset_high_coord_precision$coordinatePrecision < 1 | subset_high_coord_precision$coordinatePrecision > 25)) == 0

        #extract the number of high precision points according to this variable
        high_precision_presences_coord_precision = nrow(subset_high_coord_precision)
    } else {

        #no value for this variable exist for the [i] species
        high_precision_presences_coord_precision = NA
    }

    #save everything
    number_high_precision_occurrences = rbind.data.frame(number_high_precision_occurrences, cbind.data.frame(selected_species, high_precision_presences_coord_uncertainty, high_precision_presences_coord_precision, check_1, check_2))
}

#remove first row with NAs
number_high_precision_occurrences = number_high_precision_occurrences[-which(rowSums(is.na(number_high_precision_occurrences)) == ncol(number_high_precision_occurrences)),]

#see the results
number_high_precision_occurrences