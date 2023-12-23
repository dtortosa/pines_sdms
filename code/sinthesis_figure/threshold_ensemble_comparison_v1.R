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
setwd("/home/dftortosa/diego_docs/science/phd/nicho_pinus")

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


##HERE OPEN STACKS TO SAVE RASTER TO COMBINE PLOTS
#open stacks for saving binary raster with current and future suitability
if(FALSE){ #Right now we are not interested in saving the rasters
    current_suit_stack = stack()
    projected_suit_inside_range_stack = stack()
    projected_suit_stack = stack()
    sum_distributions = stack()
    raster_range_calc_stack = stack()
}

#It's key that you remove all areas outside the range_calc_buffer and the water bodies for ALL rasters, because these areas would enter into the calculations. Because of this, I have carefully masked and cropped all the predictions (current, future)



#################################################################
##### CALCULATE RANGE LOSS AND CHANGE AND STACK PREDICTIONS #####
#################################################################

#species="albicaulis"
master_processor=function(species){

    #load distribution buffer
    ocurrences_buffer = raster(paste("results/ocurrences/", species, "_distribution_buffer", ".asc", sep=""))

    #drop sea areas inside the ocurrences_buffer
    ocurrences_buffer = mask(ocurrences_buffer, environment_var_low_res, inverse=FALSE)     

    #convert NAs into 0 to avoid problems in the sum
    ocurrences_buffer[which(is.na(getValues(ocurrences_buffer)))] <- 0

    #save it
    if(FALSE){ #Right now, we are not interesting in saving. 
        sum_distributions = stack(sum_distributions, ocurrences_buffer)
    }

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

    #we want to save raster_range_calc for all species in a stack to have the area terrestrial area considered in the range calculations
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

    #save raster_range_calc
    if(FALSE){ #Right now, we are not interested in saving rasters now
        raster_range_calc_stack = stack(raster_range_calc_stack, raster_range_calc)
    }

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
    suitability_changes = data.frame(species=NA, selected_threshold=NA, range_change=NA, range_loss=NA)

    #for each threshold
    #selected_threshold=50
    for(selected_threshold in seq(0,100,1)){

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

        #check that suitability outside current range and under future conditons is at least equal to the suitability inside areas that are currently suitable
        if(future_suitable_area_elsewhere < future_suitable_area_inside_current_range){
            stop("ERROR! FALSE! PROBLEM WITH CALCULATION SUITABLE ARE IN FUTURE")
        }

        #calculate range loss as (current suitable area - nº cells of that area that remain suitable ) / current suitable area, then multiplied by 100
        range_loss = ((current_suitable_area - future_suitable_area_inside_current_range) / current_suitable_area ) * 100

        #calculate range change as (nº cells of that area that are suitable across the whole calc_range_buffer - current suitable areas) / current suitable area, then multiplied by 100. Here we consider future suitability of both areas that are suitable or unsuitable currently
        range_change = ((future_suitable_area_elsewhere - current_suitable_area) / current_suitable_area ) * 100

        #save metrics of suitability changes
        suitability_changes = rbind.data.frame(suitability_changes, cbind.data.frame(species, selected_threshold, current_suitable_area, future_suitable_area_inside_current_range, future_suitable_area_elsewhere, range_change, range_loss))
    }

    suitability_changes=suitability_changes[-1,]

    return(suitability_changes)
}

#master_processor(species="albicaulis")

require(foreach)
require(doParallel) #for parallel

threshold_results_df = foreach(i=epithet_species_list[1:20], .packages=c("raster", "sf"), .combine="rbind.data.frame") %dopar% { 
    master_processor(species=i)
}

##I do not see difference cores working!!!

print(threshold_results_df)

system("mkdir -p ./results/global_figures/final_global_figures/threshold_comparisons")

write.table(threshold_results_df, "./results/global_figures/final_global_figures/threshold_comparisons/eso.tsv", sep="\t", row.names=FALSE, col.names=TRUE)




if(FALSE){


    threshold_results_df=read.table("./results/global_figures/final_global_figures/threshold_comparisons/eso.tsv", sep="\t", header=TRUE)

    #CHECK THAT THE SUITABLE AREA CURRENT IS LOWER IN EACH STEP, AND THE SAME FOR FUTURE suitaiblity inside and outisde
    #you can have higher decreases of current suitable area than future suitable area as the threshold increase, leading to less range loss with a higher threshold, like threshold 66 vs 67 in albicaulis. This is ok.
    #the important thing is that always the current suitability is lower than in the previous threshold, and the same for the future suitability 

    library(dplyr)
    mean_data <- group_by(threshold_results_df, selected_threshold) %>% summarise(range_loss = mean(range_loss, na.rm = TRUE))

    pdf("./results/global_figures/final_global_figures/threshold_comparisons/eso.pdf")

    plot(x=threshold_results_df$selected_threshold, y=threshold_results_df$range_loss, cex=0.5)
    lines(x=mean_data$selected_threshold, y=mean_data$range_loss, lwd=2, col="red")    dev.off()





    ##do some operations with the rasters so we can save all of them in stacks
    #extend the extent of the predictions to the whole globe
    current_suit = extend(current_suit, environment_var)
    projected_suit_inside_range = extend(projected_suit_inside_range, environment_var)
    phylo_ensamble_inside_range_intersect_projected_suit_inside_range = extend(phylo_ensamble_inside_range_intersect_projected_suit_inside_range, environment_var)
    projected_suit = extend(projected_suit, environment_var)
    phylo_ensamble_intersect_projected_suit = extend(phylo_ensamble_intersect_projected_suit, environment_var)

    #set NAs as zero to avoid propagation of NAs in the sum
    current_suit[which(is.na(getValues(current_suit)))] <- 0      
    projected_suit_inside_range[which(is.na(getValues(projected_suit_inside_range)))] <- 0      
    phylo_ensamble_inside_range_intersect_projected_suit_inside_range[which(is.na(getValues(phylo_ensamble_inside_range_intersect_projected_suit_inside_range)))] <- 0      
    projected_suit[which(is.na(getValues(projected_suit)))] <- 0      
    phylo_ensamble_intersect_projected_suit[which(is.na(getValues(phylo_ensamble_intersect_projected_suit)))] <- 0      

    #save rasters
    if(FALSE){ #Right now we are not interested in saving the rasters
        current_suit_stack = stack(current_suit_stack, current_suit)
        projected_suit_inside_range_stack = stack(projected_suit_inside_range_stack, projected_suit_inside_range)
        phylo_ensamble_inside_range_intersect_projected_suit_inside_range_stack = stack(phylo_ensamble_inside_range_intersect_projected_suit_inside_range_stack, phylo_ensamble_inside_range_intersect_projected_suit_inside_range)
        projected_suit_stack = stack(projected_suit_stack, projected_suit)
        phylo_ensamble_intersect_projected_suit_stack = stack(phylo_ensamble_intersect_projected_suit_stack, phylo_ensamble_intersect_projected_suit)
    }


    #remove first row without NAs
    suitability_changes = suitability_changes[-which(rowSums(is.na(suitability_changes)) == ncol(suitability_changes)),]

    #check all species are included in the table
    nrow(suitability_changes) == length(epithet_species_list)

    #save the table
    write.table(suitability_changes, "results/global_figures/initial_global_figures/suitability_changes_phylo_non_scaled_v1.csv", sep=",", row.names=FALSE, col.names=TRUE)
        #suitability_changes = read.table("results/global_figures/initial_global_figures/suitability_changes_phylo_non_scaled_v1.csv", sep=",", header=TRUE)

    #check all species are included in the stacks
    if(FALSE){
        nlayers(current_suit_stack) == 112
        nlayers(projected_suit_inside_range_stack) == 112
        nlayers(phylo_ensamble_inside_range_intersect_projected_suit_inside_range_stack) == 112
        nlayers(projected_suit_stack) == 112
        nlayers(phylo_ensamble_intersect_projected_suit_stack) == 112
        nlayers(raster_range_calc_stack)
        nlayers(sum_distributions)
    }

    #sum predictions under current and future conditions
    if(FALSE){ #Right now we are not interested in saving the rasters
        current_suit_stack_sum = calc(current_suit_stack, function(x) (sum(x)))
        projected_suit_stack_sum = calc(projected_suit_stack, function(x) (sum(x)))
        projected_suit_phylo_stack_sum = calc(phylo_ensamble_intersect_projected_suit_stack, function(x) (sum(x)))
        sum_distributions_sum = calc(sum_distributions, function(x) (sum(x)))
        raster_range_calc_stack_sum = calc(raster_range_calc_stack, function(x) (sum(x)))
    }

    #save the rasters
    if(FALSE){ #Right now we are not interested in saving the rasters
        #writeRaster(current_suit_stack_sum, "results/global_figures/initial_global_figures/current_suit_stack_sum.asc", overwrite=TRUE)
        #writeRaster(projected_suit_stack_sum, "results/global_figures/initial_global_figures/projected_suit_stack_sum.asc", overwrite=TRUE)
        writeRaster(projected_suit_phylo_stack_sum, "results/global_figures/initial_global_figures/projected_suit_phylo_nonscaled_stack_sum.asc", overwrite=TRUE)
        #writeRaster(sum_distributions_sum, "results/global_figures/initial_global_figures/sum_distributions_sum.asc", overwrite=TRUE)
        #writeRaster(raster_range_calc_stack_sum, "results/global_figures/initial_global_figures/raster_range_calc_stack_sum.asc", overwrite=TRUE)
    }


    ##Calculate the differences between phylo - no phylo
    differ_percent = data.frame(selected_species=NA, range_loss_no_phylo=NA, range_loss_phylo=NA, differ_range_loss=NA, range_change_no_phylo=NA, range_change_phylo=NA, differ_range_change=NA)
    for(i in 1:length(epithet_species_list)){

        #selected species
        selected_species = epithet_species_list[i]

        #select the [i] row
        selected_row = suitability_changes[which(suitability_changes$species==selected_species),]

        #extract percentage
        range_loss_no_phylo = selected_row$range_loss_no_phylo
        range_loss_phylo = selected_row$range_loss_phylo   
        range_change_no_phylo = selected_row$range_change_no_phylo
        range_change_phylo = selected_row$range_change_phylo  

        #calculate absolute difference
        differ_range_loss = abs(range_loss_phylo-range_loss_no_phylo)
        differ_range_change = abs(range_change_phylo-range_change_no_phylo)

        #save it
        differ_percent = rbind.data.frame( differ_percent, cbind.data.frame(selected_species, range_loss_no_phylo, range_loss_phylo, differ_range_loss, range_change_no_phylo, range_change_phylo, differ_range_change))
    }
    differ_percent = differ_percent[-which(rowSums(is.na(differ_percent)) == ncol(differ_percent)),]

    #check all species are included in the table
    nrow(differ_percent) == length(epithet_species_list)

    #calculate median and the interquartile range
    median_values = cbind.data.frame("global median", rbind.data.frame(apply(differ_percent[,-1], 2, median)))
    first_quartile_values = cbind.data.frame("global first quartile", rbind.data.frame(apply(differ_percent[,-1], 2, quantile, 0.25)))
    third_quartile_values = cbind.data.frame("global third quartile", rbind.data.frame(apply(differ_percent[,-1], 2, quantile, 0.75)))
    interquartile_range = cbind.data.frame("global interquartile range", third_quartile_values[-1]-first_quartile_values[-1]) #we have to add [-1] to avoid selecting the first column with the name
    #check
    interquartile_range[,-1] == third_quartile_values[-1] - first_quartile_values[-1]
    #add column names
    names(median_values) <- c("selected_species", colnames(differ_percent[,-1]))
    names(first_quartile_values) <- c("selected_species", colnames(differ_percent[,-1]))
    names(third_quartile_values) <- c("selected_species", colnames(differ_percent[,-1]))
    names(interquartile_range) <- c("selected_species", colnames(differ_percent[,-1]))
        #It makes sense to use the mean and the standard deviation as measures of center and spread only for distributions that are reasonably symmetric with a central peak. When outliers are present, the mean and standard deviation are not a good choice (unless you want that these outliers influence the summary statistic). An outlier can decrease/increase the mean so that the mean is too low or too high to be representative of whole sample The outlier can also decrease/increase the standard deviation, which gives the impression of a wide variability in the variable This makes sense because the standard deviation measures the average deviation of the data from the mean. So a point that has a large deviation from the mean will increase the average of the deviations. This is the case for several of the variables presented in this table (see annotated plot line). Therefore, we are going to use the median and the interquartile range.
            #par(mfcol=c(3,2)); for(i in 2:ncol(differ_percent)){plot(differ_percent[,i])}
            #Link very useful and simple: https://courses.lumenlearning.com/wmopen-concepts-statistics/chapter/standard-deviation-4-of-4/
        # I guess it does not make sense to use median +- IQR because the distance between the median and the second quartile could be different from the distance to the third quartile. 
            #For example, consider the following sequence: c(1,2,3,4,5,6,7,8,9,9,9,10). The second quartile would be 3.750, the median would be 6.500, and the third quartile would be 9.000. 9.000-6.500=2.5, while 6.500-3.750=2.75. Therefore, the median is closer to the second quartile.

    #bind to the data.frame
    differ_percent = rbind.data.frame(median_values, first_quartile_values, third_quartile_values, interquartile_range, differ_percent)
    str(differ_percent)

    #save it
    write.table(differ_percent, "results/global_figures/final_global_figures/differ_phylo_inside_nonscaled_v1.csv", sep=",", col.names = TRUE, row.names = FALSE)
        #differ_percent = read.table("results/global_figures/final_global_figures/differ_phylo_inside_nonscaled_v1.csv", header=TRUE, sep=",")



    ##############################################
    ##### COMPARE PHYLO SCALE AND NON-SCALED #####
    ##############################################

    #copy differ_percent as phylo_non_scaled to avoid confusion 
    phylo_nonscaled = differ_percent
    #remove rows of summary metrics (median, IQR...)
    phylo_nonscaled = phylo_nonscaled[which(!phylo_nonscaled$selected_species %in% c("global median", "global first quartile", "global third quartile", "global interquartile range")),]

    #load results phylo scaled
    phylo_scaled = read.table("results/global_figures/final_global_figures/differ_phylo_inside_v3.csv", header=TRUE, sep=",")
    #remove rows of summary metrics (median, IQR...)
    phylo_scaled = phylo_scaled[which(!phylo_scaled$selected_species %in% c("global median", "global first quartile", "global third quartile", "global interquartile range")),]

    #merge both dataframes
    phylo_scaled_nonscaled = merge(x=phylo_nonscaled, y=phylo_scaled, by="selected_species", suffixes = c("_p_nonscaled", "_p_scaled"), all.x = TRUE, all.y = TRUE)
        #suffixes: a character vector of length 2 specifying the suffixes to be used for making unique the names of columns in the result which are not used for merging (appearing in ‘by’ etc).
        #all.x: logical; if ‘TRUE’, then extra rows will be added to the output, one for each row in ‘x’ that has no matching row in ‘y’.  These rows will have ‘NA’s in those columns that are usually filled with values from ‘y’.  The default is ‘FALSE’, so that only rows with data from both ‘x’ and ‘y’ are included in the output.
        #all.y: logical; analogous to ‘all.x’.


    ##plot range loss phylo scaled vs non-scaled
    #open the pdf
    pdf("results/global_figures/final_global_figures/phylo_scaled_vs_non_scaled.pdf", height = 6, width = 12)
    par(mfrow=c(1,2),  mar=c(6.5, 4, 2, 2) +0.1)

    #make the plot
    plot(phylo_scaled_nonscaled$range_loss_phylo_p_nonscaled, phylo_scaled_nonscaled$range_loss_phylo_p_scaled, type="p", xlab="Range loss - Phylo non-scaled", ylab="Range loss - Phylo scaled", cex.lab=1.5)

    #make the correlation
    tests_rl = cor.test(phylo_scaled_nonscaled$range_loss_phylo_p_nonscaled, phylo_scaled_nonscaled$range_loss_phylo_p_scaled, method="spearman")

    #extract and plot the results of the correlation
    if(tests_rl$p.value < 2.2e-16){
        tests_rl_p = bquote(italic(p.value) < .(format(2.2e-16)))
    }else{
        tests_rl_p = bquote(italic(p.value) == .(format(tests_rl$p.value, digits = 3)))
    }
    text(x=20, y=60, labels = tests_rl_p, cex=1.3)
    tests_rl_s = bquote(italic(S) == .(format(tests_rl$statistic, digits = 3)))
    text(x=20, y=55, labels = tests_rl_s, cex=1.3)
    tests_rl_rho = bquote(italic(rho) == .(format(tests_rl$estimate, digits = 3)))
    text(x=20, y=50, labels = tests_rl_rho, cex=1.3)


    ##plot range change phylo scaled vs. non-scaled
    #make the plot
    plot(phylo_scaled_nonscaled$range_change_phylo_p_nonscaled, phylo_scaled_nonscaled$range_change_phylo_p_scaled, type="p", xlab="Range change - Phylo non-scaled", ylab="Range change - Phylo scaled", cex.lab=1.5)

    #make the correlation
    tests_rg = cor.test(phylo_scaled_nonscaled$range_change_phylo_p_nonscaled, phylo_scaled_nonscaled$range_change_phylo_p_scaled, method="spearman")

    #extract and plot the results of the correlation
    if(tests_rg$p.value < 2.2e-16){
        tests_rg_p = bquote(italic(p.value) < .(format(2.2e-16)))
    }else{
        tests_rg_p = bquote(italic(p.value) == .(format(tests_rg$p.value, digits = 3)))
    }
    text(x=-35, y=47, labels = tests_rg_p, cex=1.3)
    tests_rg_s = bquote(italic(S) == .(format(tests_rg$statistic, digits = 3)))
    text(x=-35, y=37, labels = tests_rg_s, cex=1.3)
    tests_rg_rho = bquote(italic(rho) == .(format(tests_rg$estimate, digits = 3)))
    text(x=-35, y=27, labels = tests_rg_rho, cex=1.3)

    #add the title plot
    #mtext("Online supplementary figure ...", side=1, font=2, cex=2, adj=0.015, padj=1.5, outer=TRUE, line=-3)

    #close the pdf
    dev.off()


    ## see the data reordered based on the difference between phylo scaled and non-scaled
    #range loss
    phylo_scaled_nonscaled[order(abs(phylo_scaled_nonscaled$range_loss_phylo_p_nonscaled - phylo_scaled_nonscaled$range_loss_phylo_p_scaled), decreasing = TRUE), c("selected_species", "range_loss_phylo_p_nonscaled", "range_loss_phylo_p_scaled")]
    #range change
    phylo_scaled_nonscaled[order(abs(phylo_scaled_nonscaled$range_change_phylo_p_nonscaled - phylo_scaled_nonscaled$range_change_phylo_p_scaled), decreasing = TRUE), c("selected_species", "range_change_phylo_p_nonscaled", "range_change_phylo_p_scaled")]

    #Only a few species show some differences between scaled and non-scaled. Given this and the great correlation between the two approaches across the whole genus, it seems that the scaling is not influencing results at the scale of the whole genus.
}