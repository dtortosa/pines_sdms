#checks comment 21

###definimos el directorio de trabajo

#Librerias

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
c("tecunumanii", "jaliscana", "discolor") %in% epithet_species_list


#make the loop for calculating the number of high precision points according to both precision variables in each species
number_high_precision_occurrences = data.frame(selected_species=NA, high_precision_presences_coord_uncertainty=NA, high_precision_presences_coord_precision=NA, check_1=NA, check_2=NA)
for (i in 1:length(epithet_species_list)){

    #select the [i] species
    selected_species = epithet_species_list[i]

    #load the occurrences
    species_presences <- read.csv(paste("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_", selected_species, ".csv", sep=""), header=TRUE, fill=TRUE, check.names=TRUE, stringsAsFactors=FALSE) 
        #stringsAsFactors=FALSE sirve para procesar mejor la columna coordinatePrecision, que lleva mezclados números y caracteres
        #we need the raw data to get coordinate precision/coordinate uncertainty. In the final data.frame obtained from loop_cleaning_occurrences, there are only 4 columns: longitude, latitude, precision_weight and presence.

    #if there is any column about coordinate precision
    if(length(which(colnames(species_presences) == "coordinateUncertaintyInMeters")) > 0){

        #Esta variable indica la distancia horizontal en metros desde la latitud decimal dada y la longitud decimal ada describrined the smallest circle containing the whole of the location. Se deja empty si la incertidumbre es desconocida, no puede ser estimada o no es aplicable (porque no hay coordenadas). Zero no es valor valido para este termino. Debe ser mayor y diferente de cero además de menor de 5000000 (5000 km). 
        #points with precision lower than 4 km (4000 m) will receive a weight of 1, whilts point with coarser precision will receive a weight of 0.5.
            #info obtained from loop_cleaning_occurrences_v2.R
    

        subset_high_precision_presences_coord_uncertainty = species_presences[which(species_presences$coordinateUncertaintyInMeters < 4000),]

        check_1 = length(which(subset_high_precision_presences_coord_uncertainty$coordinateUncertaintyInMeters > 4000)) == 0

        #
        high_precision_presences_coord_uncertainty = nrow(subset_high_precision_presences_coord_uncertainty)
    } else {

        #
        high_precision_presences_coord_uncertainty = NA
    }

    #if there is any column about coordinate precision
    if(length(which(colnames(species_presences) == "coordinatePrecision")) > 0){

        #Hasta abril de 2016 este termino era el radio de error que tiene cada coordenada. La unidad es en metros, por tanto un valor de 25 indica un radio de 25x25 metros, y eso es lo que tenemos nosotros. 
        #Desde abril de 2016 ha cambiado. Una representación lineal de la precisión de las coordenadas dada en latitud y longitud decimal. Debe estar entre 0 y 1. Nosotros NO tenemos esto, por que la peña sigue subiendo los datos para esta variable como metros. 
        #Mis variables tenian una resolucion de 4x4 km  y aqui tenemos valores por debajo (1,10, 25 = 1.85 km) y algunos por encima (0 =  18 km), los calculos están explicados en niche_questions.txt. Points with 1,10,25 will receive a weight of 1, whilst points with 0 will receive a weight of 0.5.             
            #info obtained from loop_cleaning_occurrences_v2.R
    
        subset_high_coord_precision = species_presences[which(species_presences$coordinatePrecision>=1 & species_presences$coordinatePrecision<=25),]

        check_2 = length(which(subset_high_coord_precision$coordinatePrecision < 1 | subset_high_coord_precision$coordinatePrecision > 25)) == 0

        #
        high_precision_presences_coord_precision = nrow(subset_high_coord_precision)
    } else {

        #
        high_precision_presences_coord_precision = NA
    }

    number_high_precision_occurrences = rbind.data.frame(number_high_precision_occurrences, cbind.data.frame(selected_species, high_precision_presences_coord_uncertainty, high_precision_presences_coord_precision, check_1, check_2))
}

#REVISA ESTE CODIGO Y PREPARA LA RESPUESTA CON niche_questions.txt.

#se pierden casos de high precision al solo coger coordinate Precision. revisa niche_questions.txt., di en la respuesta que cuando lo hicimos en 2016, coordinate uncertaining se acaba de aplicar, y vimos en varias especies problemas. No encontramos una forma clara de separar ocurrencias adecuadas de no adecuadas usando esa variable, como defines de forma precisa que ocurrencias vienen de un raster/atlas? esto reduciría la reproducibilidad, así que directamente nos quedamos con la variable original.
