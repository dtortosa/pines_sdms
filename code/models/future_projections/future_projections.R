#SEWAL. Code for modelling and project into the future. It is prepared for run in sewall.  

###################################
#ESTABLECE EL DIRECTORIO DE TRABAJO
###################################
#DIRECTORIO DE TRABAJO
setwd("/Users/dsalazar/nicho_pinus/")


#function for make future projections
projections = function(species){

    #required libraries
    #library(raster)
    #library(dismo)
    #library(randomForest) #for predict rf models
    #library(gam) #for predict gam model
    #library(rgdal) #for problem in ensamble predictions

    #begin the species
    print(paste("begin",species))

    #load cluster number for each species
    group_species = read.csv("/Users/dsalazar/nicho_pinus/data/climate/complete_2_10_g_2.csv", header=TRUE)[,c("species", "groups")] 

    #select the group of the corresponding species
    variables_cluster = group_species[which(group_species$species == species),]$groups 

    #Load data with presences and values of environmental variables
    data = read.csv(paste("/Users/dsalazar/nicho_pinus/data/data_prepare_modelling", paste(species, "csv", sep="."), sep="/"), header=TRUE)

    #select the number cluster of this species 
    rasters_list = list.files("/Users/dsalazar/nicho_pinus/data/climate/finals", pattern=".asc", full.names=TRUE) #list the corresponding group of variables
    rasters_names = list.files("/Users/dsalazar/nicho_pinus/data/climate/finals", pattern=".asc", full.names=FALSE) #list names of raster
    names_variables = NULL #loop for separate raster names from extension ".asc"
    for (i in rasters_names){
        names_variables = append(names_variables, strsplit(i, split=".asc")[[1]])
    }
    variables_stack = stack(rasters_list) #stack them 
    names(variables_stack) = names_variables #give the names to layers of the stack

    #load names of selected variables  
    load("/Users/dsalazar/nicho_pinus/data/climate/list_selected_variables.rda") 

    #take selected variables from the list of variables
    selected_variables = ultimate_variables[[variables_cluster]]

    #calculate the number of ocurrences for each species 
    number_ocurrences = read.csv("/Users/dsalazar/nicho_pinus/data/ocurrences/ocurrences_per_species.csv", header=TRUE)
    n_ocurrence = number_ocurrences[number_ocurrences$species==species,]$number_ocurrences

    if (n_ocurrence<length(selected_variables)*10){ #if there is not 10 ocurrences for each variable 
        
        #load the names of selected variables for los number ocurrences species
        load("/Users/dsalazar/nicho_pinus/data/climate/finals/final_variables_low_number_ocurrence_species.rda")

        #select the selecte variables for low number ocurrences species
        variables_stack = variables_stack[[final_variables_low_number_ocurrence_species_new[[species]]]]    
    
    } else {

        variables_stack = variables_stack[[selected_variables]]
    
    }       
 
    #extract names of the current variables for select soil variable in the loop
    names_variables_stack = names(variables_stack)

    #load raster of PA buffer
    raster_PA_buffer = raster(paste("/Users/dsalazar/nicho_pinus/data/pa_buffers", paste(species, "PA_buffer.asc", sep="_"), sep="/"))

    #crop variables using distribution buffer
    variables_stack = crop(variables_stack, raster_PA_buffer)

    #create a vector with all bio variables
    bio = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

    #select from selected variables only bio variables, which will be replaced by future variables
    if (n_ocurrence<length(selected_variables)*10){
        future_selected = final_variables_low_number_ocurrence_species_new[[species]][final_variables_low_number_ocurrence_species_new[[species]] %in% bio]        
    } else {
        future_selected = selected_variables[selected_variables %in% bio]
    }

    #list of continuos projections and binary projections for all scenarios and climatic models
    climatic_scenarios = c("bc26", "bc45", "bc60", "bc85", "cc26", "cc45", "cc60", "cc85", "gs26", "gs45", "gs60", "gs85", "he26", "he45", "he60", "he85", "ip26", "ip45", "ip60", "ip85", "mg26", "mg45", "mg60", "mg85", "mr26", "mr45", "mr60", "mr85")
    
    #load models 
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/models", paste(species, "glm_model.rda", sep="_"), sep="/"))
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/models", paste(species, "gam_model.rda", sep="_"), sep="/"))
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/models", paste(species, "rf_model.rda", sep="_"), sep="/"))

    #load thresholds
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/threshold", paste(species, "glm_threshold.rda", sep="_"), sep="/"))
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/threshold", paste(species, "gam_threshold.rda", sep="_"), sep="/"))

    #loop for run 12 glms, 12 gams and 12 rf for each climatic_model*IPCC_scenario
    for (i in climatic_scenarios){

        #load future climatic variables 
        list_future_climate = list.files("/Users/dsalazar/nicho_pinus/data/climate/future", pattern=i, full.names=TRUE) #list all raster of the [i] modelo*scenario 
        list_names_climatic_variables = list.files("/Users/dsalazar/nicho_pinus/data/climate/future", pattern=i) #list names 
        names_climatic_variables = NULL #drop .asc extension
        for (k in list_names_climatic_variables){
            names_climatic_variables = append(names_climatic_variables, strsplit(k, split=".asc")[[1]])
        }

        #stack all of them 
        future_climate = stack(list_future_climate)
        names(future_climate) = names_climatic_variables #give names to the layers

        #crop the future variable with the PA buffer (only predict in the interesting area)
        future_climate = crop(future_climate, raster_PA_buffer)

        #resample the future climatic variables 
        future_climate = resample(future_climate, variables_stack[[1]], method="bilinear")
        
        #stack future climatic variables and the current soil variables
        if(nlayers(variables_stack)>1){
            final_future_climate = stack(
                future_climate[[paste(future_selected, i, sep="_")]], #select the selected variables using i (model*scenario) and "bioXX"
                variables_stack[[names_variables_stack[names_variables_stack %in% c("ph", "cec", "carbon", "depth", "sand", "silt", "clay")]]]) #select from the stack the only the SOIL selected variables, becuase of this we use a vector with all soil variables names
        } else {
            final_future_climate = stack(future_climate[[paste(future_selected, i, sep="_")]])
        }   

        #loop for change the names of variables and match with variable names in the models   
        require(stringr) #require for nchar (calculate the number of characters) and str_split_fixed

        final_names = NULL 
        for (k in names(final_future_climate)){ #for each name of the variables in the final stack

            name_splitted = str_split_fixed(k, "_", 2)
           
            if (nchar(name_splitted[2])>0){ #if the number has more than 0 characters
               
                final_names =  append(final_names, name_splitted[1]) #select the first part of the split (bio name)

            } else { #if not and thus is a soil variable 

                 final_names = append(final_names, k) #save exactly the same name 
            }
        }

        #change the names of the variables in the stack 
        names(final_future_climate) = final_names

        #list of continuos projections for each model 
        glm_projections = list()
        gam_projections = list()
        rf_projections = list()

        #list of binary projections for each model 
        glm_projection_bin = list()
        gam_projection_bin = list()
        rf_projection_bin = list()

        for (k in 1:12){

            glm_projections[[k]] = predict(final_future_climate, glm_resample[[k]], type="response")
            gam_projections[[k]] = predict(final_future_climate, gam_resample[[k]], type="response")
            rf_projections[[k]] = predict(final_future_climate, rf_resample[[k]], type="response")

            #convert each projection in binary using the threshold calculated with fit_eval_models. 
            #glm
            glm_projection_bin[[k]] = glm_projections[[k]] #copy the raster 
            glm_projection_bin[[k]][glm_projection_bin[[k]]>=(glm_threshold[[k]][2,2]),] <- 1 #give 1 to the pixels with a predicted value higher or equal than the threshold (TSS)
            glm_projection_bin[[k]][glm_projection_bin[[k]]<(glm_threshold[[k]][2,2]),] <- 0 #give 0 to the pixels with a predicted value lower than the threshold (TSS)

            #gam
            gam_projection_bin[[k]] = gam_projections[[k]] #copy the raster 
            gam_projection_bin[[k]][gam_projection_bin[[k]]>=(gam_threshold[[k]][2,2]),] <- 1 #give 1 to the pixels with a predicted value higher or equal than the threshold (TSS)
            gam_projection_bin[[k]][gam_projection_bin[[k]]<(gam_threshold[[k]][2,2]),] <- 0 #give 0 to the pixels with a predicted value lower than the threshold (TSS)

            #rf
            rf_projection_bin[[k]] = rf_projections[[k]] #It does not change  because it is already a binary projection 
        }

        #save continous projections of glm and gam
        writeRaster(stack(glm_projections), filename=paste(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/continuous_projections", paste("continuous", paste("projection_glm", i, sep="_"), sep="_"), sep="/"), paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
        writeRaster(stack(gam_projections), filename=paste(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/continuous_projections", paste("continuous", paste("projection_gam", i, sep="_"), sep="_"), sep="/"), paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)        

        #save binary projections of gam and glm
        writeRaster(stack(glm_projection_bin), filename=paste(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_projections", paste("binary", paste(paste("projection", "glm", sep="_"), i, sep="_"), sep="_"), sep="/"), paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
        writeRaster(stack(gam_projection_bin), filename=paste(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_projections", paste("binary", paste(paste("projection", "gam", sep="_"), i, sep="_"), sep="_"), sep="/"), paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)

        #save binary projectiosn of rf
        writeRaster(stack(rf_projection_bin), filename=paste(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_projections", paste("binary", paste(paste("projection", "rf", sep="_"), i, sep="_"), sep="_"), sep="/"), paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
        
        #rm projections in the final scenario to release memory space before ensambling
        if(i == tail(climatic_scenarios, 1)){
            rm(
                glm_projections,
                gam_projections,
                glm_projection_bin,
                gam_projection_bin,
                rf_projection_bin)
        }          
    }

    ##reate the ensamble
    #load all projections
    list_projections = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_projections", pattern=glob2rx(paste("*", species, "*", "tif", sep="")), full.names=TRUE)
    all_projections = stack(list_projections)

    #calculate the percentage of models for which a pixel is suitable
    ensamble_projections_bin = calc(all_projections, function(x) (sum(x)*100)/nlayers(all_projections))

    #save ensamble
    writeRaster(ensamble_projections_bin, filename=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin", paste("ensamble_projections_bin", paste(species, "tif", sep="."), sep="_"), sep="/"), overwrite=TRUE)

    #########################################
    ######ZIP AND DELETE MODELS #############
    ######################################### 
    #list continuous projections created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files    
        models = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/models", pattern=glob2rx(paste("*", species, "*", "rda", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        models = models[!grepl("pseudostrobus", models)]
    } else {
        models = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/models", pattern=glob2rx(paste("*", species, "*", "rda", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)        
    } 

    #zip all continuous projections for lack of space in disk
    zip(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/models/models", paste(species, "zip", sep="."), sep="_"), models, flags="-j") #j indicate that you don't want all the directory structure

    #delete continuous projections
    file.remove(models) #binary will be used in the ensamble, because of this we don't delete them now    

    #########################################
    ######ZIP AND DELETE THRESHOLDS #########
    ######################################### 
    #list continuous projections created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files        
        threshold = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/threshold", pattern=glob2rx(paste("*", species, "*", "rda", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        threshold = threshold[!grepl("pseudostrobus", threshold)]
    } else {
        threshold = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/threshold", pattern=glob2rx(paste("*", species, "*", "rda", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
    }

    #zip all continuous projections for lack of space in disk
    zip(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/threshold/threshold", paste(species, "zip", sep="."), sep="_"), threshold, flags="-j") #j indicate that you don't want all the directory structure

    #delete continuous projections
    file.remove(threshold)  

    ########################################################
    ######ZIP AND DELETE BIN AND CONTINUOUS PROJECTIONS ####
    ######################################################## 
    #list continuous projections created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files       
        continuous_projections = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/continuous_projections", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        continuous_projections = continuous_projections[!grepl("pseudostrobus", continuous_projections)]
    } else {
        continuous_projections = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/continuous_projections", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)        
    }
    #zip all continuous projections for lack of space in disk
    zip(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/continuous_projections/continuous_projections", paste(species, "zip", sep="."), sep="_"), continuous_projections, flags="-j") #j indicate that you don't want all the directory structure

    #delete continuous projections
    file.remove(continuous_projections) #binary will be used in the ensamble, because of this we don't delete them now    

    #list binary projections created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files           
        binary_projections = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_projections", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        binary_projections = binary_projections[!grepl("pseudostrobus", binary_projections)]
    } else {
        binary_projections = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_projections", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)        
    }
    #zip all binary projections for lack of space in disk
    zip(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_projections/binary_projections", paste(species, "zip", sep="."), sep="_"), binary_projections, flags="-j") #j indicate that you don't want all the directory structure

    #delete binary projections
    file.remove(binary_projections)

    #name of species
    print(paste(species, "ended"))  
}


########Paralelize the process######
require(foreach)
require(doParallel) #for parallel

###load species names: First species with low and high precision points, then the rest of species ordered in basis on the difference between correct and incorrect PA weight (more ratio, more difference first).

#list species
list_species = read.table("/Users/dsalazar/nicho_pinus/data/list_species.txt", sep="\t", header=T)
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


#load data about precision occurrence
results_high_vs_low_occurrences = read.table("/Users/dsalazar/nicho_pinus/data/results_high_vs_low_occurrences.txt", sep="\t", header=TRUE)

#select those species with both types od precision points
pines_high_low = results_high_vs_low_occurrences[which(results_high_vs_low_occurrences$high_low_precision == TRUE),]$selected_species

#reorder being sylvestris (biggest dataset) the first, then species with more effect of phylocorrection. 
pines_high_low = pines_high_low[c(7,2,4,5,6,3,1)]

#load data about differences between calculations of PA weight
results_weight_PA = read.table("/Users/dsalazar/nicho_pinus/data/results_weight_PA.txt", sep="\t", header=TRUE)

#reorder the df in basis on the differences between both PA weight calculations: Species with more differences first
results_weight_PA_order = results_weight_PA[order(results_weight_PA$correct_incorrect_ratio, decreasing=TRUE),]

#remove those species with both types of occurrences (hig an low), which has been included in pines_high_low
results_weight_PA_order = results_weight_PA_order[-which(results_weight_PA_order$selected_species %in% pines_high_low),]
pines_high_low %in% results_weight_PA_order$selected_species

#save the rest of species
rest_pines_order_by_weight_PA = results_weight_PA_order$selected_species

#bind all
species_list = c(as.vector(pines_high_low), as.vector(rest_pines_order_by_weight_PA))
length(species_list) == 113
summary(epithet_species_list %in% species_list)

#extract species for which we have binary projections
species_with_projections_path = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin", pattern=".tif")

#if else por si no hay ninguna especie con ensamble
if(length(species_with_projections_path) > 0){

    #extract species names
    species_with_projections = NULL
    for(i in 1:length(species_with_projections_path)){
        species_with_projections = append(species_with_projections, strsplit(strsplit(species_with_projections_path[i], split="_")[[1]][4], split=".tif")[[1]])
    }

    #analyse only species withput final ensamble
    species_list = species_list[which(!species_list %in% species_with_projections)]
    length(species_list) + length(species_with_projections) == 113
}

# set up cluster
clust <- makeCluster(2) 
registerDoParallel(clust)

###########################
#########PROJECTIONS#######
###########################
#project to the future for non_problematic species
foreach(i = species_list, .packages=c("raster", "dismo", "randomForest", "gam", "rgdal")) %dopar% { 
    projections(species = i)
} 


#stop the cluster 
stopCluster(clust)