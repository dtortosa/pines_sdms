#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
    #https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
    #https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



#########################################################
####### COMPARISON OF SEVERAL ENSEMBLE THRESHOLDS ####### #########################################################

#This script compares the results of considering different ensemble thresholds.



###################################################
##### DIFFERENCES RESPECT TO PREVIOUS VERSION #####
###################################################

#Respect to version 1:



########################
##### BEGIN SCRIPT #####
########################

#set wroking directory
#setwd("/home/dftortosa/diego_docs/science/phd/nicho_pinus")

#make a folders
system("mkdir -p ./results/global_figures/final_global_figures/threshold_comparisons")
system("mkdir -p ./results/global_figures/final_global_figures/threshold_comparisons/range_change_loss")
system("mkdir -p ./results/global_figures/final_global_figures/threshold_comparisons/pine_richness_change")
system("mkdir -p ./results/global_figures/final_global_figures/threshold_comparisons/suitability_stacks")
system("mkdir -p ./results/global_figures/final_global_figures/threshold_comparisons/raster_range_calc")

#pre-defined functions
plot_sin=function(input){
    jpeg("./singularity_plot.jpeg", height=2000, width=2000, res=300)
    plot(input)
    dev.off()
}

#require packages
require(raster)
require(sf)

#load species names
list_species = read.table("code/presences/species.txt", sep="\t", header=TRUE)

#extract epithet from species list
epithet_species_list = NULL
for(i in 1:nrow(list_species)){

    #selected species
    selected_species = as.vector(list_species[i,])

    #extract epithet
    epithet_species_list = append(epithet_species_list, strsplit(selected_species, split=" ")[[1]][2])
}
#check there is no NA
summary(!is.na(epithet_species_list))
#check
if(FALSE){
    require(tidyverse)
    paste("Pinus", epithet_species_list, sep=" ") == str_trim(as.vector(list_species[,1])) #the seperated epithet more Pinus are equal to list_species? We use str_trim from tidyverse to remove end spaces in each element of list_species
}#it is in false because loading tidyverse lead to load several packages that have a function name "extract", and this gives problems with the extract function of raster. If you want to check run these lines manually

#remove tecunumanii, jaliscana y discolor. These species are not used for Niche paper. The two first because we should downloadad gbif data now, so we would mix gbif data form 2016 and 2019. The third one was impossible to differentiate from P. cembroides
epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("tecunumanii", "jaliscana", "discolor"))]
#check these species are not present
!c("tecunumanii", "jaliscana", "discolor") %in% epithet_species_list

#load environment variables for using them as a background
clay = raster("datos/finals/clay.asc")
bio1 = raster("datos/finals/bio1.asc")
environment_var = clay*bio1
environment_var[which(getValues(environment_var) >= min(getValues(environment_var), na.rm = TRUE))] <- 0 #set all continent areas as 0. These are areas with data, that is not NA. All continent areas will be zero, while the rest would be zero.

#load buffer albicaulis to get a reduced resolution version of environment_var to mask the distribution buffers used for the sum of distribution
buffer_albicaulis = raster(paste("results/ocurrences/albicaulis_distribution_buffer", ".asc", sep=""))

#resample environment_var
environment_var_low_res = resample(environment_var, buffer_albicaulis, method="bilinear")

#It's key that you remove all areas outside the range_calc_buffer and the water bodies for ALL rasters, because these areas would enter into the calculations. Because of this, I have carefully masked and cropped all the predictions (current, future)




###########################
##### DEFINE FUNCTION #####
###########################

#set the thresholds and species to be tested
thresholds_to_test = seq(0,100,1)
species_to_test = epithet_species_list[1:length(epithet_species_list)]

#write the function
#species="albicaulis"
master_processor=function(species){

    #load distribution buffer
    ocurrences_buffer = raster(paste("results/ocurrences/", species, "_distribution_buffer", ".asc", sep=""))

    #drop sea areas inside the ocurrences_buffer
    ocurrences_buffer = mask(ocurrences_buffer, environment_var_low_res, inverse=FALSE)     

    #convert NAs into 0 to avoid problems in the sum
    ocurrences_buffer[which(is.na(getValues(ocurrences_buffer)))] <- 0

    #load the polygon used for calculations of changes of suitability (calc_ranges)
    if(!species=="pumila"){
        raster_range_calc = raster(paste("results/global_figures/buffers_calc_ranges/", species, "_range_calc_buffer.asc", sep=""))
    } else {
        raster_range_calc = raster(paste("results/global_figures/buffers_calc_ranges/", species, "_buffer_range_calc.asc", sep=""))        
    }             
    polygon_range_calc = rasterToPolygons(raster_range_calc, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon

    #crop and mask clay (this raster will be used to remove sea areas from polygon_range_calc)
    environment_var_cropped = crop(environment_var, polygon_range_calc)
    environment_var_cropped = mask(environment_var_cropped, polygon_range_calc)
    environment_var_low_res_cropped = crop(environment_var_low_res, polygon_range_calc)
    environment_var_low_res_cropped = mask(environment_var_low_res_cropped, polygon_range_calc)

    #we want to save raster_range_calc for all species in a stack to have the area terrestrial area considered in the range calculations so we can crop the predictions across the globe
    #crop current suitability to reduce map size
    raster_range_calc = crop(raster_range_calc, polygon_range_calc)

    #mask current suitability to remove all areas outside the buffer calc range
    raster_range_calc = mask(raster_range_calc, polygon_range_calc)

    #mask with clay to remove water bodies
    raster_range_calc = mask(raster_range_calc, environment_var_low_res_cropped)

    #extend the raster
    raster_range_calc = extend(raster_range_calc, environment_var_low_res)

    #convert NAs into 0 to avoid problems in the sum
    raster_range_calc[which(is.na(getValues(raster_range_calc)))] <- 0

    #save raster_range_calc. This will be used to create a polygon with the whole distribution of all pines globally
    writeRaster(raster_range_calc, paste("./results/global_figures/final_global_figures/threshold_comparisons/raster_range_calc/raster_range_calc_", species, sep=""), options="COMPRESS=LZW", overwrite=TRUE)

    #load predicted suitability
    current_suit = raster(paste("results/ensamble_predictions_bin/ensamble_predictions_bin_", species, ".tif", sep=""))

    #crop current suitability to reduce map size
    current_suit = crop(current_suit, polygon_range_calc)

    #mask current suitability to remove all areas outside the buffer calc range
    current_suit = mask(current_suit, polygon_range_calc)

    #mask with clay to remove water bodies
    current_suit = mask(current_suit, environment_var_cropped)

    #load projected suitability
    projected_suit = raster(paste("results/ensamble_projections_bin/ensamble_projections_bin_", species, ".tif", sep=""))

    #crop current suitability to reduce map size
    projected_suit = crop(projected_suit, polygon_range_calc)
    
    #mask current suitability to remove all areas outside the buffer calc range
    projected_suit = mask(projected_suit, polygon_range_calc)

    #mask with clay to remove water bodies
    projected_suit = mask(projected_suit, environment_var_cropped)
        #in the case of pumila, masking with the polygon_range_calc buffer leaves the sea of japan with zero instead of NA. This is caused when removing sea areas from the raster of that buffer, the two extremes of the Japan's sea almost touch and the area inside is included. This is not a problem because after that, we mask with the environmnetal varaible (bio1 and clay), so sea areas are removed. In species with several polygons of distribution is not a problem because: 1) The calc_range_buffer is ver big, so in almost all cases all polygons are included within it. If sea areas inside of them they will be removedd with environment_var. 

    #open data frame to save metrics of suitability change
    suitability_changes = data.frame(species=NA, selected_threshold=NA, current_suitable_area=NA, future_suitable_area_inside_current_range=NA, future_suitable_area_elsewhere=NA, range_change=NA, range_loss=NA)

    #open stacks for saving binary raster with current and future suitability
    selected_current_suit_stack = stack()
    selected_projected_suit_stack = stack()

    #for each threshold
    #selected_threshold=50
    for(selected_threshold in thresholds_to_test){

        ##obtain maps with zero and ones from suitability maps (1 means suitable) using the selected threshold
        #current suitability
        selected_current_suit = current_suit
        selected_current_suit[which(getValues(selected_current_suit) < selected_threshold)] <- 0 #set as zero those areas with suitability lower than 75
        selected_current_suit[which(getValues(selected_current_suit) > 0)] <- 1 #set as 1 all areas with suitability higher than zero (i.e. all with suit equal or higher than k)

        #future suitability
        selected_projected_suit = projected_suit
        selected_projected_suit[which(getValues(selected_projected_suit) < selected_threshold)] <- 0 #set as zero those areas with suitability lower than 75
        selected_projected_suit[which(getValues(selected_projected_suit) > 0)] <- 1 #set as 1 all areas with suitability higher than zero (i.e. all with suit equal or higher than k)

        #extract suitability under future conditions from areas that are currently suitables
        projected_suit_inside_range = selected_projected_suit
        projected_suit_inside_range[which(!getValues(selected_current_suit)==1)] <- 0

        #extract size of area suitable
        current_suitable_area = length(which(getValues(selected_current_suit)==1))
        future_suitable_area_inside_current_range = length(which(getValues(projected_suit_inside_range)==1))
        future_suitable_area_elsewhere = length(which(getValues(selected_projected_suit)==1))

        #we should have the same cells with 1 in projected_suit_inside_range and selected_current_suit
        if(sum(!which(getValues(projected_suit_inside_range)==1) %in% which(getValues(selected_current_suit)==1)) > 0){
            stop("ERROR! PROBLEM! WE HAVE A PROBLEM WHEN CALCULATING THE FUTURE SUITABLE AREAS IN CURRENT SUITABLE REGIONS")
        }

        #check that suitability outside current range and under future conditons is at least equal to the suitability inside areas that are currently suitable
        if((future_suitable_area_elsewhere < future_suitable_area_inside_current_range) | (future_suitable_area_inside_current_range > current_suitable_area)){
            stop("ERROR! FALSE! PROBLEM WITH CALCULATION SUITABLE ARE IN FUTURE")
        }

        #calculate range loss as (current suitable area - nº cells of that area that remain suitable ) / current suitable area, then multiplied by 100
        range_loss = ((current_suitable_area - future_suitable_area_inside_current_range) / current_suitable_area) * 100

        #calculate range change as (nº cells of that area that are suitable across the whole calc_range_buffer - current suitable areas) / current suitable area, then multiplied by 100. Here we consider future suitability of both areas that are suitable or unsuitable currently
        range_change = ((future_suitable_area_elsewhere - current_suitable_area) / current_suitable_area ) * 100

        #save metrics of suitability changes
        suitability_changes = rbind.data.frame(suitability_changes, cbind.data.frame(species, selected_threshold, current_suitable_area, future_suitable_area_inside_current_range, future_suitable_area_elsewhere, range_change, range_loss))

        #save suitability maps
        #first, update layer names
        names(selected_current_suit)=paste("threshold_", selected_threshold, "_", species, sep="")
        names(selected_projected_suit)=paste("threshold_", selected_threshold, "_", species, sep="")
        #then save
        selected_current_suit_stack = stack(selected_current_suit_stack, selected_current_suit)
        selected_projected_suit_stack = stack(selected_projected_suit_stack, selected_projected_suit)
    }

    #check we have the correct layer names
    if(
        (!identical(names(selected_current_suit_stack), paste("threshold_", thresholds_to_test, "_", species, sep=""))) | 
        (!identical(names(selected_projected_suit_stack), paste("threshold_", thresholds_to_test, "_", species, sep="")))){
        stop(paste("ERROR! FALSE! WE HAVE A PROBLEM CALCULATING THE PREDICTIONS ACROSS THRESHOLDS, WE DO NOT HAVE ALL THE CORRECT LAYER NAMES FOR ", species, sep="")) 
    }

    #check we have a layer per each threshold in both stacks
    if((nlayers(selected_current_suit_stack) != length(thresholds_to_test)) | (nlayers(selected_projected_suit_stack) != length(thresholds_to_test))){
        stop(paste("ERROR! FALSE! WE HAVE A PROBLEM CALCULATING THE PREDICTIONS ACROSS THRESHOLDS, WE DO NOT HAVE ALL THRESHOLDS FOR ", species, sep=""))
    }

    #save the stacks
    writeRaster(selected_current_suit_stack, paste("./results/global_figures/final_global_figures/threshold_comparisons/suitability_stacks/stack_current_suit_", species, sep=""), options="COMPRESS=LZW", overwrite=TRUE)
    writeRaster(selected_projected_suit_stack, paste("./results/global_figures/final_global_figures/threshold_comparisons/suitability_stacks/stack_future_suit_", species, sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        #https://stackoverflow.com/questions/42041695/writeraster-output-file-size
        #CHECK COMPRESSION

    #remove first row without NAs
    suitability_changes = suitability_changes[-which(rowSums(is.na(suitability_changes)) == ncol(suitability_changes)),]

    #check that the suitable area is equal or lower always as the threshold increases. Note that we can have higher decreases of current suitable area than future suitable area as the threshold increases, leading to less range loss with a higher threshold, like threshold 66 vs 67 in albicaulis. This is ok. We calculate range loss as (current-future)/current. If the current area decreases more than the future, this means that the denominator is smaller while the numerator is larger, so the total is larger. This is an expected behaviour because we have proportionally higher future suitability respect to the current suitability, as current has decreased more.
    #the important thing here is that always the current suitability is lower than in the previous threshold, and the same for the future suitability. This is what we are going to check here.
    #suit_var="current_suitable_area"
    for(suit_var in c("current_suitable_area", "future_suitable_area_inside_current_range", "future_suitable_area_elsewhere")){

        #select the column corresponding with the selected variable 
        selected_variable=suitability_changes[,which(colnames(suitability_changes)==suit_var)]
        
        #sort the variable from higher to lower
        selected_variable_sorted=sort(selected_variable, decreasing=TRUE)

        #check whether the variable is the same after sorting
        if(identical(selected_variable_sorted, selected_variable) == FALSE){
            stop("ERROR! FALSE! WE HAVE A PROBLEM WITH THE CALCULATION OF THE SIZE OF SUITABLE AREA")
        }
    }

    #return only the table with range change and loss, the rest of results are saved as stacks
    return(suitability_changes)
}

#run it for just one species
#master_processor(species="albicaulis")




#######################
##### PARALLELIZE #####
#######################
require(foreach)
require(doParallel) #for parallel

#set up cluster
clust <- makeCluster(20, outfile="")
    #only 20 to avoid memory explosion
    #You can usually figure out why the worker died by using the makeCluster "outfile" option so that the error message generated by the worker isn't thrown away. I usually recommend using outfile=""
        #https://stackoverflow.com/a/24352032
registerDoParallel(clust)

#run the function in parallel
threshold_results_df = foreach(i=species_to_test, .packages=c("raster", "sf"), .combine="rbind.data.frame") %dopar% { 
    master_processor(species=i)
}
    #.combine: 
        #function that is used to process the tasks results as they generated.  This can be specified as either a function or a non-empty character string naming the function. Specifying 'c' is useful for concatenating the results into a vector, for example.  The values 'cbind' and 'rbind' can combine vectors into a matrix.
        #we use rbind.data.frame to combine each vector as a row in a data.frame

#stop the cluster 
stopCluster(clust)

#see the results
print(head(threshold_results_df))
print(summary(threshold_results_df))

#check we do not have NANs
if(!identical(na.omit(threshold_results_df), threshold_results_df)){
    stop("ERROR! FALSE! WE HAVE NANS IN THE DATA FRAME WITH RANGE LOSS AND CHANGE ACROSS THRESHOLDS AND SPECIES")
}

#check we have all species and the correct number of rows
check_1=nrow(threshold_results_df) == length(species_to_test)*101
    #we have calculated range loss/change per each of the 112 species across 101 thresholds
check_2=identical(unique(threshold_results_df$species), species_to_test)
if(check_1 & check_2){
    print("GOOD TO GO! We have the correct number of rows")
} else {
    stop("ERROR! FALSE! WE HAVE A PROBLEM WITH THE CALCULATIONS ACROSS THREHOLDS, WE DO NOT HAVE ALL SPECIES OR THE EXPECTED NUMBER OF ROWS")
}

#save the table
write.table(threshold_results_df, "./results/global_figures/final_global_figures/threshold_comparisons/range_change_loss/range_change_loss_thresholds.tsv", sep="\t", row.names=FALSE, col.names=TRUE)
#threshold_results_df=read.table("./results/global_figures/final_global_figures/threshold_comparisons/range_change_loss_thresholds.tsv", sep="\t", header=TRUE)





########################################################
##### PLOT RANGE CHANGE AND LOSS ACROSS THRESHOLDS #####
########################################################

#calculate median of range change/loss per threshold
library(dplyr)
mean_range_loss <- group_by(threshold_results_df, selected_threshold) %>% 
    summarise(range_loss=median(range_loss), na.rm=TRUE)
mean_range_change <- group_by(threshold_results_df, selected_threshold) %>% 
    summarise(range_change=median(range_change), na.rm=TRUE)

#open the plot
jpeg("./results/global_figures/final_global_figures/threshold_comparisons/range_change_loss/range_change_loss_thresholds.jpeg", height=2000, width=2000, res=300)
par(mfcol=c(2,1), mar=c(4.1, 4, 1, 4))


#plot range loss against the threshold value for each species
plot(x=threshold_results_df$selected_threshold, y=threshold_results_df$range_loss, ylab="", xlab="", type="n", lwd=0.1, cex=0.5)
#species="albicaulis"
for(species in unique(threshold_results_df$species)){
    lines(
        x=threshold_results_df[threshold_results_df$species==species,]$selected_threshold, 
        y=threshold_results_df[threshold_results_df$species==species,]$range_loss, 
        lwd=0.2, cex=0.5)
}
title(ylab="Range loss (%)", line=2.5, cex.lab=1.2)
lines(x=mean_range_loss$selected_threshold, y=mean_range_loss$range_loss, lwd=3, col="red")

#plot range change
plot(x=threshold_results_df$selected_threshold, y=threshold_results_df$range_change, ylab="", xlab="", type="n", lwd=0.1, cex=0.5)
#species="albicaulis"
for(species in unique(threshold_results_df$species)){
    lines(
        x=threshold_results_df[threshold_results_df$species==species,]$selected_threshold, 
        y=threshold_results_df[threshold_results_df$species==species,]$range_change, 
        lwd=0.2, cex=0.5)
}
title(xlab="Threshold (%)", ylab="Range change (%)", line=2.5, cex.lab=1.2)
lines(x=mean_range_change$selected_threshold, y=mean_range_change$range_change, lwd=3, col="red")
dev.off()
detach("package:dplyr", unload=TRUE)
    #some function names are overlapped with raster and other packages




#######################################################################
##### CALCULATE THE DIFFERENCE IN PINE RICHNESS IN EACH THRESHOLD #####
#######################################################################


##obtain a polygon with the combined area used to male calculations in all pine species
#load the polygons used for calculations of changes of suitability (calc_ranges) in all species
raster_range_calc_stack = stack(list.files("./results/global_figures/final_global_figures/threshold_comparisons/raster_range_calc", full.names=TRUE, pattern=".grd"))

#check we have all species
if(nlayers(raster_range_calc_stack) != length(species_to_test)){
    stop("ERROR! FALSE! WE DO NOT HAVE ALL SPECIES IN raster_range_calc_stack")
}

#sum the calculation area of all species
raster_range_calc_stack_sum = calc(raster_range_calc_stack, function(x) (sum(x)))

#convert raster_range_calc_stack_sum to zero-one
raster_range_calc_stack_sum[which(getValues(raster_range_calc_stack_sum) > 0)] <- 1
    #now every cell being considered in the calculation of at least 1 species will be 1

#convert to polygon
polygon_range_calc_stack_sum = rasterToPolygons(raster_range_calc_stack_sum, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon
    #use that polygon for masking predictions rasters and remove areas outside the buffers calc ranges


##define a function to calculate the difference in pine richness between current and future predictions for a given threshold
#n_layer=50
#we use the number of the layer instead of the name of the threshold. Threshold 0 is layer number 1. Threshold 100 is layer number 101. We can extract layer 1 from a stack, but not number 0.
stack_pred_threshold = function(n_layer){
    
    #starting
    print(paste("STARTING THRESHOLD ", n_layer-1, sep=""))

    #open stacks to save the pine richness of each species for the same layer (threshold) across the different threshold calculations
    current_suit_stack=stack()
    future_suit_stack=stack()

    #species=species_to_test[1]
    for(species in species_to_test){

        #load the specific layer within the current and future predictions for the selected species, which includes the different thresholds
        #we use the number of layer to extract it, getting a stack of 1 layer. Then we have to extract it to get the raster layer instead of a stack
        current_suit_species_threhold = stack(paste("./results/global_figures/final_global_figures/threshold_comparisons/suitability_stacks/stack_current_suit_", species, ".grd", sep=""), bands=n_layer)[[1]]
        future_suit_species_threhold = stack(paste("./results/global_figures/final_global_figures/threshold_comparisons/suitability_stacks/stack_future_suit_", species, ".grd", sep=""), bands=n_layer)[[1]]

        #check we have selected the correct layer name
        if(
            (names(current_suit_species_threhold)!=paste("threshold_", n_layer-1, "_", species, sep="")) |
            (names(future_suit_species_threhold)!=paste("threshold_", n_layer-1, "_", species, sep=""))){
            stop(paste("ERROR! FALSE! WE HAVE A PROBLEM WITH THE SELECTION OF THE SPECIFIC LAYER FOR THRESHOLD ", n_layer-1, " in ", species, sep=""))

        }

        #extend the extent of the predictions to the whole globe
        current_suit_species_threhold = extend(current_suit_species_threhold, environment_var)
        future_suit_species_threhold = extend(future_suit_species_threhold, environment_var)

        #set NAs as zero to avoid propagation of NAs in the sum
        current_suit_species_threhold[which(is.na(getValues(current_suit_species_threhold)))] <- 0      
        future_suit_species_threhold[which(is.na(getValues(future_suit_species_threhold)))] <- 0  

        #remove the areas outside the global distribution of pines
        current_suit_species_threhold = mask(current_suit_species_threhold, polygon_range_calc_stack_sum)
        future_suit_species_threhold = mask(future_suit_species_threhold, polygon_range_calc_stack_sum)

        #save the raster
        current_suit_stack=stack(current_suit_stack, current_suit_species_threhold)
        future_suit_stack=stack(future_suit_stack, future_suit_species_threhold)
    }

    #check we have a layer per each threshold in both stacks
    if((nlayers(current_suit_stack) != length(species_to_test)) | (nlayers(future_suit_stack) != length(species_to_test))){
        stop(paste("ERROR! FALSE! WE HAVE A PROBLEM CALCULATING PINE RICHNESS ACROSS THRESHOLDS, WE DO NOT HAVE ALL SPECIES FOR THRESHOLD ", n_layer-1, sep=""))
    }

    #check we have selected the correct layers
    if( 
        (!identical(names(current_suit_stack), paste("threshold_", n_layer-1, "_", species_to_test, sep=""))) |
        (!identical(names(future_suit_stack), paste("threshold_", n_layer-1, "_", species_to_test, sep="")))){
        stop(paste("ERROR! FALSE! WE HAVE NOT SELECTED THE CORRECT LAYERS FOR THRESHOLD ", n_layer-1, sep=""))
    }

    #sum the presence of all pines into one single raster for current and future conditions, respectively
    current_suit_stack_sum = sum(current_suit_stack)
    future_suit_stack_sum = sum(future_suit_stack)
    
    #calculate the difference in pine richness between current and future conditions for the selected threshold
    diff_suit = future_suit_stack_sum -current_suit_stack_sum 

    #add name to the layer using the threshold number
    names(diff_suit) = paste("threshold_", n_layer-1, sep="")

    #ending
    print(paste("ENDING THRESHOLD ", n_layer-1, sep=""))
    return(diff_suit)
}
#stack_pred_threshold(1)

##parallelize the function
#load packages
require(foreach)
require(doParallel)

#set up cluster
clust <- makeCluster(15, outfile="")
    #5 less cores to avoid memory explosion
    #You can usually figure out why the worker died by using the makeCluster "outfile" option so that the error message generated by the worker isn't thrown away. I usually recommend using outfile=""
        #https://stackoverflow.com/a/24352032
registerDoParallel(clust)

#run the function in parallel
diff_suit_stack = foreach(i=thresholds_to_test+1, .packages=c("raster", "sf"), .combine="stack") %dopar% { 
    stack_pred_threshold(n_layer=i)
}
    #thresholds_to_test+1 is the number of layers. Threshold 0 is layer 1, threshold 1 is layer 2, and so on...
    #the output is a raster, so use stack to save the different raster into a stack

#stop the cluster 
stopCluster(clust)

#reorder the layers to ensure we have thresholds in increasing order
diff_suit_stack=diff_suit_stack[[paste("threshold_", thresholds_to_test, sep="")]]
    #if a threshold is missing, you will get a warning and a missing layer

#check
if(
    (class(diff_suit_stack)[1] != "RasterStack") | 
    (nlayers(diff_suit_stack) != length(thresholds_to_test)) |
    (!identical(names(diff_suit_stack), paste("threshold_", thresholds_to_test, sep="")))){
    stop("ERROR! FALSE! WE HAVE A PROBLEM CALCULATING THE DIFFERENCE IN PINE RICHNESS ACROSS THRESHOLDS")
}


##calculate the percentage of thresholds for which a cell have less or more pine species
#copy the stack with the difference
diff_suit_stack_positive=diff_suit_stack
diff_suit_stack_negative=diff_suit_stack

#in each layer, convert to zero those cases that have negative or positive change in richness for the positive and negative rasters, respectively
#layer=1
for(layer in 1:nlayers(diff_suit_stack)){
    diff_suit_stack_positive[[layer]][which(getValues(diff_suit_stack_positive[[layer]])<0)] <- 0
    diff_suit_stack_positive[[layer]][which(getValues(diff_suit_stack_positive[[layer]])>0)] <- 1
    diff_suit_stack_negative[[layer]][which(getValues(diff_suit_stack_negative[[layer]])>0)] <- 0
    diff_suit_stack_negative[[layer]][which(getValues(diff_suit_stack_negative[[layer]])<0)] <- 1
}
    #IMPORTANT: this approach loses information about cases where pine richness increases in more than 1, but I think this is ok. If you have an area where richness increases by 2 up to 75% of uncertainty, and then other area increases richness by 1 from 1 to 90% of uncertainty. We have more certainty that richness is going to increase in the second scenario even if the number of species is reduced.

#calculate the proportion of thresholds that have positive and negative change in pine richness, respectively
positive_pine_richness_across_thresholds = (sum(diff_suit_stack_positive)/nlayers(diff_suit_stack_positive))*100
negative_pine_richness_across_thresholds = (sum(diff_suit_stack_negative)/nlayers(diff_suit_stack_negative))*100

#save them
writeRaster(positive_pine_richness_across_thresholds, "./results/global_figures/final_global_figures/threshold_comparisons/pine_richness_change/positive_pine_richness_across_thresholds", options="COMPRESS=LZW", overwrite=TRUE)
writeRaster(negative_pine_richness_across_thresholds, "./results/global_figures/final_global_figures/threshold_comparisons/pine_richness_change/negative_pine_richness_across_thresholds", options="COMPRESS=LZW", overwrite=TRUE)


##make two plots with percentage of thresholds gaining or losing species 
#set extent of the raster to be plotted
plot_extent = c(-180,180,-10,90) #if you change the extent of the plot, some sea border could change (only difference of 1 cell). For example, when you plot the whole globe, save the plot as pdf and then zoom to the Canary Islands, the shape of the island is not very accurate. If you crop the maps to P. canariensis buffer, then the shape of the islands in the pdf is more similar to the reality.

##define the color palette
require(RColorBrewer)
#We selected from Colorbrewer a single hue pallete with green. As we are going to plot only one variable in each plot.
green_palette <-brewer.pal(9,"Greens")
    #Names taken from "http://colorbrewer2.org/#type=sequential&scheme=Greens&n=9"
    #All works for anomalous trychromacy and dychromacy ("http://www.color-blindness.com/coblis-color-blindness-simulator/")
#these palletes are used in colorRampPalette to create a function that can create a great number of colors 
colfunc_green <- colorRampPalette(green_palette)


##open the plot
jpeg("./results/global_figures/final_global_figures/threshold_comparisons/pine_richness_change/change_pine_richness.jpeg", height=2000, width=2000, res=300)
par(mfcol=c(2,1), mai=c(0,0.4,0,0.5), oma=c(0,0,2,1))

#upper plot
plot(crop(environment_var, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, axes=FALSE, box=FALSE, main="") #higher values in argument start of gray colors lead to brighter gray
plot(positive_pine_richness_across_thresholds, add=TRUE, col=colfunc_green(161), breaks=seq(0,100,1), axes=FALSE, box=FALSE, axis.args=list(at=seq(0,100,20), cex.axis=1.3), legend=TRUE, legend.shrink=0.8, legend.args=list(text=expression(bold(paste('% Thresholds + richness'))), side=4, font=2, line=3.4, cex=1.2)) 
    #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"
    #breaks indicate the number of partitions between colors, whilst "at" indicate the numbers in the legend. The number of colors have to be EQUAL to the number o breaks. Breaks should encompass the RANGE of values of the raster

#lower plot
plot(crop(environment_var, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, axes=FALSE, box=FALSE, main="") #higher values in argument start of gray colors lead to brighter gray
plot(negative_pine_richness_across_thresholds, add=TRUE, col=colfunc_green(161), breaks=seq(0,100,1), axes=FALSE, box=FALSE, axis.args=list(at=seq(0,100,20), cex.axis=1.3), legend=TRUE, legend.shrink=0.8, legend.args=list(text=expression(bold(paste('% Thresholds - richness'))), side=4, font=2, line=3.4, cex=1.2))
dev.off()




#######################################################################
##### FINISH THE SCRIPT #####
#######################################################################

#finish the script
print("## FINISH ##")
