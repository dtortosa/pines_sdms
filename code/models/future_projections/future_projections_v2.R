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
    #library(rgeos)

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

    #load distribution map
    distribution = raster(paste("/Users/dsalazar/nicho_pinus/data/MAPS/p_", species, "_01.img", sep=""))
    
    #read raster distribution buffer
    raster_buffer = raster(paste("/Users/dsalazar/nicho_pinus/data/buffers", paste(species, "distribution_buffer.asc", sep="_"), sep="/"))

    #create a polygon from that raster
    polygon_buffer = rasterToPolygons(raster_buffer, fun=function(x){x==1}, dissolve=TRUE)
    
    ###Create the pseudoabsence buffer for crop the environmental variables. The difference between this buffer and the distribution_buffer will be the area with PAs.
    if(species == "pumila"){#if the species is pumila we load the PA buffer previously created with buffer in both sides of the map
        #load the polygon with polygons in both sides of the map previously created
        raster_range_calc_buffer = raster("/Users/dsalazar/nicho_pinus/data/buffers_range_calc/pumila_buffer_range_calc.asc")

        #create the new polygon halepensis distribution without sea areas
        polygon_range_calc_buffer = rasterToPolygons(raster_range_calc_buffer, fun=function(x){x==1}, dissolve=TRUE)        
    } else {#if not, we create the buffer here
        #create a polygon buffer around the distribution buffer
        polygon_range_calc_buffer = gBuffer(polygon_buffer, byid=FALSE, id=NULL, width=12.5, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

        #crop the distribution raster with the polygon PA buffer
        distribution_crop =  crop(distribution, polygon_range_calc_buffer) 

        #create a raster from the PA buffer polygon 
        raster_range_calc_buffer = raster() 
        extent(raster_range_calc_buffer) = extent(distribution_crop) 
        res(raster_range_calc_buffer) = res(distribution) 
        raster_range_calc_buffer  = rasterize(polygon_range_calc_buffer, raster_range_calc_buffer) 
        #plot(crop(bio6, raster_range_calc_buffer), col="gray80")
        #plot(raster_range_calc_buffer, add=T)
        #plot(polygon_range_calc_buffer, add=T)
        #plot(polygon_buffer, add=T)

        #drop the sea areas
        raster_range_calc_buffer = distribution_crop*raster_range_calc_buffer 
        raster_range_calc_buffer[!is.na(raster_range_calc_buffer)] <- 1 
        #plot(crop(bio6, raster_range_calc_buffer), col="gray80")
        #plot(raster_range_calc_buffer, add=T, col="yellow")
        #plot(raster_buffer, add=T, col="blue")
        #legend(10, 57, legend=c("Pseudo-Absences buffer", "Presences buffer"), fill=c("yellow", "blue"))
        #We mantenain some water bodies inside the continents because we can't multiply our soil variables with our PA buffer, they have different resolution. Because of this, we will use part of vaule of environmental variables over water bodies for variables selection, I think it is not very important.
           
        #create the new polygon halepensis distribution without sea areas
        polygon_range_calc_buffer = rasterToPolygons(raster_range_calc_buffer, fun=function(x){x==1}, dissolve=TRUE)

        #write the PA buffer without water bodies.  
        writeRaster(raster_range_calc_buffer, paste("/Users/dsalazar/nicho_pinus/data/buffers_range_calc", paste(species, "range_calc_buffer.asc", sep="_"), sep="/"), format="ascii", overwrite=TRUE)              
    }

    #crop variables using calc range buffer to reduce the extension of the map (reduce running time)
    variables_stack = crop(variables_stack, polygon_range_calc_buffer)
    
    #crop variables using range calc buffer to reduce all areas outside the buffer (we use crop before because mask not reduce the extension of the map)
    variables_stack = mask(variables_stack, polygon_range_calc_buffer)#no problem with sea areas inside the buffer because these will be removed when range calculations be done

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

        #crop variables using range calc buffer to reduce all areas outside the buffer (we use crop before because mask not reduce the extension of the map)
        future_climate = crop(future_climate, polygon_range_calc_buffer)

        #crop variables using range calc buffer to reduce all areas outside the buffer (we use crop before because mask not reduce the extension of the map)
        future_climate = mask(future_climate, polygon_range_calc_buffer)#no problem with sea areas inside the buffer because these will be removed when range calculations be done

        #resample the future climatic variables 
        future_climate = resample(future_climate, variables_stack[[1]], method="bilinear")
            #IMPORTANT: THE LAST RUN OF THE MODELS WAS MADE RUNNING THIS LINE, BUT I DO NOT UNDERSTAND WHY, FUTURE CLIMATE IS ALREADY AT THE SAME RESOLUTION THAN CURRENT CLIMATE. THE RESULTS IS THE SAME MAP, BUT ADITIONAL CELLS AROUND COAST AREAS, Y EXTRMELY SIMILAR, BUT GIVE PROBLEMS IN RANGE LOSS CALCULATIONS OF SPECIES WITHOUT SOIL VARIABLES IN THEIR MODELS. FOR THE FIGURES IS NOT AN ISSUE BECAUSE WE ADD A RASTER WITH WATER BODIES, THIS RASTER WAS OBTAINED USING CURREN SUITABILITY DATA, FOR WHICH WE DO NOT RESAMPLE CLIMATE, SO THIS PROBLEM DOES NOT EXIST.
            #IF YOU HAVE TO RE-RUN, CHECK IF YOU CAN AVOID THIS LINE.

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
#check
if(FALSE){
    require(tidyverse)
    paste("Pinus", epithet_species_list, sep=" ") == str_trim(as.vector(list_species[,1])) #the seperated epithet more Pinus are equal to list_species? We use str_trim from tidyverse to remove end spaces in each element of list_species
}#it is in false because loading tidyverse lead to load several packages that have a function name "extract", and this gives problems with the extract function of raster. If you want to check run these lines manually

#remove tecunumanii, jaliscana y discolor. These species are not used for Niche paper. The two first because we should downloadad gbif data now, so we would mix gbif data form 2016 and 2019. The third one was impossible to differentiate from P. cembroides
epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("tecunumanii", "jaliscana", "discolor"))]
#check
c("tecunumanii", "jaliscana", "discolor") %in% epithet_species_list


# set up cluster
clust <- makeCluster(2) 
registerDoParallel(clust)

###########################
#########PROJECTIONS#######
###########################
#project to the future for non_problematic species
foreach(i = epithet_species_list, .packages=c("raster", "dismo", "randomForest", "gam", "rgdal", "rgeos")) %dopar% { 
    projections(species = i)
} 


if(FALSE){
####do the remaining species. We separate sibirica form sylvestris, and add the the closest species those smaller (low RAM)
group_1 = c("sibirica", "torreyana", "tropicalis", "tabuliformis", "taeda", "teocote", "virginiana",  "yecorensis", "yunnanensis")
group_2 = c("sylvestris", "taiwanensis", "wallichiana", "washoensis", "thunbergii", "squamata", "strobiformis", "strobus")

#project to the future for non_problematic species
foreach(i = group_1, .packages=c("raster", "dismo", "randomForest", "gam", "rgdal", "rgeos")) %dopar% { 
    projections(species = i)
} 

#project to the future for non_problematic species
foreach(i = group_2, .packages=c("raster", "dismo", "randomForest", "gam", "rgdal", "rgeos")) %dopar% { 
    projections(species = i)
} 
}

#stop the cluster 
stopCluster(clust)