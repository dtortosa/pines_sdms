#set wroking directory
setwd("/Users/dsalazar/nicho_pinus")

#require packages
require(raster)
require(rgeos)

#load species names
list_species = read.table("/Users/dsalazar/nicho_pinus/data/list_species.txt", sep="\t", header=TRUE)

#extract epithet from species list
epithet_species_list = NULL
for(i in 1:nrow(list_species)){

    #selected species
    selected_species = as.vector(list_species[i,])

    #extract epithet
    epithet_species_list = append(epithet_species_list, strsplit(selected_species, split=" ")[[1]][2])
}
summary(is.na(epithet_species_list)) #all false

#drop discolor (problem tazonomy, no diferetiaced from cembriodes)
if("discolor" %in% epithet_species_list){
    epithet_species_list = epithet_species_list[-which(epithet_species_list == "discolor")]
}

#check it
!"discolor" %in% epithet_species_list
length(epithet_species_list) == 112

#load clay for using it as a background
clay = raster("/Users/dsalazar/nicho_pinus/data/climate/finals/clay.asc") #better clay than bioclim variables because some water bodies inside continent have climate data and then we would obtain the percentage of pines with suitability there, whcih would be zero, but this is not a real data.
albicaulis_distribution = raster("/Users/dsalazar/nicho_pinus/data/MAPS/p_albicaulis_01.img") #load ablicaulis buffer to get resolution of distribution maps
clay = resample(clay, albicaulis_distribution, method="bilinear") #reduce resolution. There is no problem with the resampling because cells with zero are actually zero, not like on phylo rasters in which NA cells were set as zero. This would artificially decrease the suitability of cells surrounding from false-zero cells. This is not the case. 
clay[which(getValues(clay) >= min(getValues(clay), na.rm = TRUE))] <- 0 #set all continent areas as 0

#load all distributions into a stack and also create migration buffers
stack_distribution = stack() #stack of distributions
stack_migration = stack() #stack of migration buffer WITHOUT distribution area inside
stack_migration_distrib = stack() #stack of migration buffer WITH distribution area inside
for(i in 1:length(epithet_species_list)){
    
    #select the [i] epithet
    selected_epi = epithet_species_list[i]

    #print species name
    print(selected_epi)

    ##### ocurrence buffer
    #load [i] buffer as raster
    ocurrences_buffer = raster(paste("/Users/dsalazar/nicho_pinus/data/buffers/", selected_epi, "_distribution_buffer", ".asc", sep=""))

    #drop sea areas inside the ocurrences_buffer
    ocurrences_buffer = mask(ocurrences_buffer, clay, inverse=FALSE)     

    #convert NAs into 0 to avoid problems in the sum
    ocurrences_buffer[which(is.na(getValues(ocurrences_buffer)))] <- 0

    #### migration buffer
    #convert [i] buffer in a polygon 
    ocurrences_buffer_polygon = rasterToPolygons(ocurrences_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon

    #set species name as name of this raster
    names(ocurrences_buffer) <- selected_epi

    #create a migration buffer
    polygon_migration_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=10, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) #byid=FALSE indicate that the function should be applied to the entire object or to sub-geometries

    #copy ocurrences_buffer as migration_buffer
    migration_buffer = ocurrences_buffer

    #set as 1 those areas included in the polygon migration buffer
    migration_buffer[polygon_migration_buffer] <- 1

    #drop sea areas inside the migration buffer
    migration_buffer = mask(migration_buffer, clay, inverse=FALSE) 

    #convert NAs into 0 to avoid problems in the sum
    migration_buffer[which(is.na(getValues(migration_buffer)))] <- 0

    #set species name as name of this raster
    names(migration_buffer) <- selected_epi

    #save it before remove areas inside of the buffer included in the distribution
    stack_migration_distrib = stack(stack_migration_distrib, migration_buffer)

    #drop areas included in the distribution
    migration_buffer = mask(migration_buffer, ocurrences_buffer_polygon, inverse=TRUE)

    #convert NAs into 0 to avoid problems in the sum, because the new masked areas inside the distribution are indicated with NA, so we have to change NA by 0
    migration_buffer[which(is.na(getValues(migration_buffer)))] <- 0

    #save into a raster
    stack_distribution = stack(stack_distribution, ocurrences_buffer)
    stack_migration = stack(stack_migration, migration_buffer)
}
nlayers(stack_distribution) == 112 #without idscolor
nlayers(stack_migration) == 112 #without idscolor
nlayers(stack_migration_distrib) == 112 #without idscolor

#sum distributions to get the number of pines in each cell
sum_distributions = calc(stack_distribution, function(x) (sum(x)))
sum_distributions[which(getValues(sum_distributions)==0)] <- NA #drop areas without pines

#sum the same but for migration buffer
sum_buffers = calc(stack_migration, function(x) (sum(x)))
sum_buffers[which(getValues(sum_buffers)==0)] <- NA #drop areas without pines

#sum migration_distribution buffers
sum_buffers_distrib = calc(stack_migration_distrib, function(x) (sum(x)))
sum_buffers_distrib[which(getValues(sum_buffers_distrib)==0)] <- NA #drop areas without pines

#save both stacks
writeRaster(stack_distribution, filename="/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/stack_distribution.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(stack_migration, filename="/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/stack_migration.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(stack_migration_distrib, filename="/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/stack_migration_distrib.tif", options="INTERLEAVE=BAND", overwrite=TRUE)

#save all rasters 
writeRaster(sum_distributions, filename="/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/sum_distributions.asc", overwrite=TRUE) #sum of distributions
writeRaster(sum_buffers, filename="/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/sum_buffers.asc", overwrite=TRUE) #sum of migration buffers without including distribution areas
writeRaster(sum_buffers_distrib, filename="/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/sum_buffers_distrib.asc", overwrite=TRUE) #sum of migration buffers including distribution areas 

#read them
sum_distributions=raster("/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/sum_distributions.asc")
sum_buffers=raster("/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/sum_buffers.asc")
sum_buffers_distrib=raster("/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/sum_buffers_distrib.asc")

#in case we want a fast test with Mediterranean pines
mediterranean_epithet = epithet_species_list[which(epithet_species_list %in% c("halepensis", "pinaster", "pinea", "nigra", "sylvestris", "mugo"))]

#calculate sum*100 of cells with loss of suitability
synthetic_plots = function(type, phylo_cor="no"){

    #loop for calculating suitability loss across species
    if(!type=="histogram"){
        suit_change_stack = stack()
    } else {
        percent_loss = data.frame(species=NA, percent_loss=NA)
    }        
    for(i in 1:length(epithet_species_list)){
    
        #select the [i] epithet
        selected_epi = epithet_species_list[i]

        #print the name of [i] species
        print(selected_epi)

        #load [i] buffer as raster
        ocurrences_buffer = raster(paste("/Users/dsalazar/nicho_pinus/data/buffers/", selected_epi, "_distribution_buffer", ".asc", sep=""))

        #convert [i] buffer in a polygon 
        ocurrences_buffer_polygon = rasterToPolygons(ocurrences_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon

        #load the [i] current suitability
        current_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_predictions_bin/ensamble_predictions_bin_", selected_epi, ".tif", sep=""))

        #reduce resolution
        current_suit = resample(current_suit, sum_distributions, method="bilinear") #there is no problem with the resampling in this case because the cases with zero are real data, are areas with no suitability. The problem occurred in phylo-suitability rasters because all NAs were converted into zero, so phylo-suitable areas close to false-zeros were artificially reduced its suitability. 

        #mask the current suitability with the polyong distribution
        current_suit_masked = mask(current_suit, ocurrences_buffer_polygon)

        ##According to the plot
        #% of reduction of suitability per species
        if(type=="histogram"){

            #load the [i] future suitability
            future_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep="")) #there is no problem with the resampling in this case because the cases ith zero are real data, are areas with no suitability. The problem occured in phylo-suitability rasters because all NAs were converted into zero, so phylo-suitable areas close to false-zeros were artificially reduced its suitbaility. 
        
            #reduce resolution
            future_suit = resample(future_suit, sum_distributions, method="bilinear")

            #mask the future suitability with the polyong distribution
            future_suit_masked = mask(future_suit, ocurrences_buffer_polygon)

            #number of cells inside the buffer without NA
            total_cells = length(which(!is.na(getValues(current_suit_masked))))
            #total_cells == length(which(!is.na(getValues(future_suit_masked)))) #both current and future raster have the same number of cells

            #cells suitable currently
            suitable_cells_current = length(which(getValues(current_suit_masked) >= 75))

            #if we want to include phylocorrected data
            if(phylo_cor=="no"){

                #cells suitable future
                suitable_cells_future = length(which(getValues(future_suit_masked) >= 75))
            } else {

                #set as NA all areas with lower suitability than 75
                future_suit_masked[which(getValues(future_suit_masked)<75)] <- NA
                #set the rest as 1
                future_suit_masked[which(!is.na(getValues(future_suit_masked)))] <- 1

                #load again future_suit because we want it without resampling. This raster will be used to select those areas of the phylo raster that are considered as intermediate, and set NAs all outside of that range. This is important to avoid FALSE zeros in the resampling process, these zeros would artificially reduce the phylo suitaiblity
                future_suit_to_subset = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))

                #load phylo niche
                phylo_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/phylo_ensamble/with_proportions/", selected_epi, "_phylo_ensamble_with_proportions.asc", sep="")) 

                #create a raster with only 1 in areas with suitability uncertainity in future_suit_to_subset, the rest NA
                phylo_mask = phylo_suit
                phylo_mask[which(getValues(future_suit_to_subset) > 25 & getValues(future_suit_to_subset) < 75)] <- 1 #select all uncertainty areas
                phylo_suit[which(getValues(phylo_mask)==0)] <- NA #the rest of areas will be NA in phylo_suit. In this way, we remove the false zeros (all phylo_suit is cover by zero except in those areas where phylo-suit is higher than zero), an those zeros area not used in the resampling. Think that if you have close two cells one with 0.3 and other of zero, the bigger cells under lower resolution cannot be 0.3 neither 0, it should be 0.15.... this is ok if the zero is real, i.e. an areas with uncertainity according to SDMs that it is not included in the phylo range, but this is a artificially reduction of phylo suitability if that zero is false, an area with higher than 75 or lower than 25 regular suitability that was set as zero for using it as background. 

                #reduce resolution
                phylo_suit = resample(phylo_suit, sum_distributions, method="bilinear")#there is a reduction of suitability, considering zeros is very very low, and here we take all positive values as 1, so it doesn't matter.
                    #mean values of suitability before and after of resampling change a litle bit

                #mask the future suitability with the polyong distribution
                phylo_suit_masked = mask(phylo_suit, ocurrences_buffer_polygon)

                #select only areas inside of the phylo range with probability higher than 0.25. We are sure that they are not very bad, a higher threshold is not possible given the low suitability
                phylo_suit_masked[which(getValues(phylo_suit_masked) < 0.25 )] <- NA
                phylo_suit_masked[which(!is.na(getValues(phylo_suit_masked)))] <- 1

                #sum suitable cells included in the phylo niche but not in the regular niche, and sum them to the suitable cells of the regular niche
                sum_suitable_phylo_no_phylo = length(which(!is.na(getValues(phylo_suit_masked)) & is.na(getValues(future_suit_masked)))) + length(which(getValues(future_suit_masked) > 0))

                #add to the suitable areas, the suitable areas according to phylo niche
                future_suit_masked[which(getValues(phylo_suit_masked)==1)] <- 1

                #cells suitable future
                suitable_cells_future = length(which(getValues(future_suit_masked) > 0))

                #check if cells included in the phylo niche but not in the regular niche, summed to suitable areas of the regular niche gives the same number than suitable cells in the final future_suit_masked.
                print("Test is ok:")
                print(suitable_cells_future == sum_suitable_phylo_no_phylo)
                print("")
            } 

            #percentage of reduction (current-future/total)*100
            percent_species = cbind.data.frame(selected_epi, ((suitable_cells_current-suitable_cells_future)/total_cells)*100)
            colnames(percent_species) <- c("species", "percent_loss") #For each species, suitability loss was calculated as the difference between suitable area under current and future conditions, and then divided by the total area of species range (it is expressed in percentage). The rationale of this is that we want to know the differences between current and future, e.g. from 80 to 40. This entails a 50% of reduction, but we also want to correct by the size of the range, this reduction of 40 how much is it from the total range? If the total range area is 100, then a reduction of 40 entails 40%, not 50. In that way we can compare across species with different range sizes.

            #save it
            percent_loss = rbind.data.frame(percent_loss, percent_species)
        }
        #sum of suitability across pine distributions
        if(type=="sum_suitability"){

            if(phylo_cor=="no"){
 
                #load the [i] future suitability
                future_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))
        
                #reduce resolution
                future_suit = resample(future_suit, sum_distributions, method="bilinear") #there is no problem with the resampling in this case because the cases ith zero are real data, are areas with no suitability. The problem occured in phylo-suitability rasters because all NAs were converted into zero, so phylo-suitable areas close to false-zeros were artificially reduced its suitbaility. 

                #mask the future suitability with the polyong distribution
                future_suit_masked = mask(future_suit, ocurrences_buffer_polygon)

                #set NAs as zero to avoid propagation of NAs in the sum
                future_suit_masked[which(is.na(getValues(future_suit_masked)))] <- 0 

                #save into an stack with the speceis name the suitability
                suit_change_stack = stack(suit_change_stack, future_suit_masked)
                names(suit_change_stack[[i]]) <- selected_epi   
            } else {
                return(print("ERROR: This plot can not be phylocorrected automatically")) #I want the phylosuitability in only one raster for inside and outside of the distribution
            }             
        } 
        #sum of suitability across pine distributions
        if(type=="sum_suitability_outside"){

            if(phylo_cor=="no"){
 
                #load the [i] future suitability
                future_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))
        
                #reduce resolution
                future_suit = resample(future_suit, sum_distributions, method="bilinear")#there is no problem with the resampling in this case because the cases ith zero are real data, are areas with no suitability. The problem occured in phylo-suitability rasters because all NAs were converted into zero, so phylo-suitable areas close to false-zeros were artificially reduced its suitbaility. 

                #create a polygon buffer around the distribution buffer
                polygon_migration_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=10, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) 

                #mask the future suitability with the polyong distribution
                future_suit_masked = mask(future_suit, polygon_migration_buffer)
                
                #drop areas inside distribution
                future_suit_masked = mask(future_suit_masked, ocurrences_buffer_polygon, inverse=TRUE)                

                #set NAs as zero to avoid propagation of NAs in the sum
                future_suit_masked[which(is.na(getValues(future_suit_masked)))] <- 0 

                #save into an stack with the speceis name the suitability
                suit_change_stack = stack(suit_change_stack, future_suit_masked)
                names(suit_change_stack[[i]]) <- selected_epi   
            } else {
                return(print("ERROR: This plot can not be phylocorrected automatically")) #I want the phylosuitability in only one raster for inside and outside of the distribution
            }             
        }                   
        #sum of phylo suitability inside and outside pine distribution
        if(type=="sum_suitability_phylo_inout"){
            if(phylo_cor=="yes"){

                #load again future_suit because we want it without resampling. This raster will be used to select those areas of the phylo raster that are considered as intermediate, and set NAs all outside of that range. This is important to avoid FALSE zeros in the resampling process, these zeros would artificially reduce the phylo suitaiblity
                future_suit_to_subset = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))

                #load phylo niche
                phylo_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/phylo_ensamble/with_proportions/", selected_epi, "_phylo_ensamble_with_proportions.asc", sep="")) 

                #create a raster with only 1 in areas with suitability uncertainity in future_suit_to_subset, the rest NA
                phylo_mask = phylo_suit
                phylo_mask[which(getValues(future_suit_to_subset) > 25 & getValues(future_suit_to_subset) < 75)] <- 1 #select all uncertainity areas
                phylo_suit[which(getValues(phylo_mask)==0)] <- NA #the rest of areas will be NA in phylo_suit. In this way, we remove the false zeros (all phylo_suit is cover by zero except in those areas where phylo-suit is higher than zero), an those zeros area not used in the resampling. Think that if you have close two cells one with 0.3 and other of zero, the bigger cells under lower resolution cannot be 0.3 neither 0, it should be 0.15.... this is ok if the zero is real, i.e. an areas with uncertainity according to SDMs that it is not included in the phylo range, but this is a artificially reduction of phylo suitability if that zero is false, an area with higher than 75 or lower than 25 regular suitability that was set as zero for using it as background. 

                #reduce resolution
                phylo_suit = resample(phylo_suit, sum_distributions, method="bilinear")

                #create a polygon buffer around the distribution buffer
                polygon_migration_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=10, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) 

                #mask the future suitability with the polyong distribution
                phylo_suit_masked = mask(phylo_suit, polygon_migration_buffer)

                #set NAs as zero to avoid that NAs se expandan and get huecos
                phylo_suit_masked[which(is.na(getValues(phylo_suit_masked)))] <- 0
                
                #save into an stack with the speceis name the phylo suitability
                suit_change_stack = stack(suit_change_stack, phylo_suit_masked)
                names(suit_change_stack[[i]]) <- selected_epi                  
            } else { 
                return(print("ERROR: This plot cannob be calulated witout phylocorrection")) #this raster is only phylocorrection, in-out, it is not make sense to get a non-phylocorrected version
            }
        }                
        #suitable areas in the future inside distribution (categoric)
        if(type=="total_more_75"){
 
            #load the [i] future suitability
            future_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))
        
            #reduce resolution
            future_suit = resample(future_suit, sum_distributions, method="bilinear") #there is no problem with the resampling in this case because the cases ith zero are real data, are areas with no suitability. The problem occured in phylo-suitability rasters because all NAs were converted into zero, so phylo-suitable areas close to false-zeros were artificially reduced its suitbaility. 

            #mask the future suitability with the polyong distribution
            future_suit_masked = mask(future_suit, ocurrences_buffer_polygon)

            #set as NA cells with suit lower than 75 and as 1 cells with higher or equal than 75. 1 for unsuitable
            future_suit_masked[which(getValues(future_suit_masked) < 75)] <- NA
            future_suit_masked[which(!is.na(getValues(future_suit_masked)))] <- 1

            #if we don't want phylo niche
            if(phylo_cor=="yes"){
 
                #load again future_suit because we want it without resampling. This raster will be used to select those areas of the phylo raster that are considered as intermediate, and set NAs all outside of that range. This is important to avoid FALSE zeros in the resampling process, these zeros would artificially reduce the phylo suitaiblity
                future_suit_to_subset = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))

                #load phylo niche
                phylo_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/phylo_ensamble/with_proportions/", selected_epi, "_phylo_ensamble_with_proportions.asc", sep="")) 

                #create a raster with only 1 in areas with suitability uncertainity in future_suit_to_subset, the rest NA
                phylo_mask = phylo_suit
                phylo_mask[which(getValues(future_suit_to_subset) > 25 & getValues(future_suit_to_subset) < 75)] <- 1 #select all uncertainity areas
                phylo_suit[which(getValues(phylo_mask)==0)] <- NA #the rest of areas will be NA in phylo_suit. In this way, we remove the false zeros (all phylo_suit is cover by zero except in those areas where phylo-suit is higher than zero), an those zeros area not used in the resampling. Think that if you have close two cells one with 0.3 and other of zero, the bigger cells under lower resolution cannot be 0.3 neither 0, it should be 0.15.... this is ok if the zero is real, i.e. an areas with uncertainity according to SDMs that it is not included in the phylo range, but this is a artificially reduction of phylo suitability if that zero is false, an area with higher than 75 or lower than 25 regular suitability that was set as zero for using it as background. 

                #reduce resolution
                phylo_suit = resample(phylo_suit, sum_distributions, method="bilinear")

                #mask the future suitability with the polyong distribution
                phylo_suit_masked = mask(phylo_suit, ocurrences_buffer_polygon)

                #select only areas inside of the phylo range
                phylo_suit_masked[which(getValues(phylo_suit_masked)==0)] <- NA
                phylo_suit_masked[which(!is.na(getValues(phylo_suit_masked)))] <- 1

                #sub suitable cells included in the phyloniche but not in the regular niche, and sum them to the suitable cells of the regular niche
                sum_suitable_phylo_no_phylo = length(which(!is.na(getValues(phylo_suit_masked)) & is.na(getValues(future_suit_masked)))) + length(which(getValues(future_suit_masked) > 0))

                #add to the suitable areas, the suitable areas according to phylo niche
                future_suit_masked[which(getValues(phylo_suit_masked)==1)] <- 1

                #cells suitable future
                suitable_cells_future = length(which(getValues(future_suit_masked) > 0))

                #check if cells included in the phylo niche but not in the regular niche, summed to suitable areas of the regular niche gives the same number than suitable cells in the final future_suit_masked.
                print("Test is ok:")
                print(suitable_cells_future == sum_suitable_phylo_no_phylo)
                print("")
            }

            #if there are cells with low of suitability in future
            if(length(which(!is.na(getValues(future_suit_masked)))) > 0){

                #create a polygon from suit_change
                polygon_suit_change = rasterToPolygons(future_suit_masked, fun=function(x){x==1}, n=16, dissolve = TRUE) #convertimos en poligono the raster using the cell with values=1

                #copy clay
                clay_to_mask = clay

                #set as 1 those areas included in the polygon of loss
                clay_to_mask[polygon_suit_change] <- 1
            } else { #if not, i.e, non suitable areas in the futuere are non suitable yet today

                #all remains as zero
                clay_to_mask = clay
            } 
            #save into an stack with the speceis name
            suit_change_stack = stack(suit_change_stack, clay_to_mask)
            names(suit_change_stack[[i]]) <- selected_epi 
        }
        #loss of suitability below 75% to less than 25 inside distribution (categoric)
        if(type=="more_75_less_25"){
            #set as NA cells with suit lower than 75 and as 1 cells with higher or equal than 75. 1 for suitable.
            current_suit_masked[which(getValues(current_suit_masked) < 75)] <- NA 
            current_suit_masked[which(!is.na(getValues(current_suit_masked)))] <- 1 

            #load the [i] future suitability
            future_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))
        
            #reduce resolution
            future_suit = resample(future_suit, sum_distributions, method="bilinear")#there is no problem with the resampling in this case because the cases ith zero are real data, are areas with no suitability. The problem occured in phylo-suitability rasters because all NAs were converted into zero, so phylo-suitable areas close to false-zeros were artificially reduced its suitbaility. 

            #mask the future suitability with the polyong distribution
            future_suit_masked = mask(future_suit, ocurrences_buffer_polygon)

            #set as NA cells with suit lower or equal than 75 and as 1 cells with higher than 75. 1 for unsuitble
            future_suit_masked[which(getValues(future_suit_masked) > 25)] <- NA
            future_suit_masked[which(!is.na(getValues(future_suit_masked)))] <- 1

            #if we don't want phylo niche
            if(phylo_cor=="yes"){

                #load again future_suit because we want it without resampling. This raster will be used to select those areas of the phylo raster that are considered as intermediate, and set NAs all outside of that range. This is important to avoid FALSE zeros in the resampling process, these zeros would artificially reduce the phylo suitaiblity
                future_suit_to_subset = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))

                #load phylo niche
                phylo_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/phylo_ensamble/with_proportions/", selected_epi, "_phylo_ensamble_with_proportions.asc", sep="")) 

                #create a raster with only 1 in areas with suitability uncertainity in future_suit_to_subset, the rest NA
                phylo_mask = phylo_suit
                phylo_mask[which(getValues(future_suit_to_subset) > 25 & getValues(future_suit_to_subset) < 75)] <- 1 #select all uncertainity areas
                phylo_suit[which(getValues(phylo_mask)==0)] <- NA #the rest of areas will be NA in phylo_suit. In this way, we remove the false zeros (all phylo_suit is cover by zero except in those areas where phylo-suit is higher than zero), an those zeros area not used in the resampling. Think that if you have close two cells one with 0.3 and other of zero, the bigger cells under lower resolution cannot be 0.3 neither 0, it should be 0.15.... this is ok if the zero is real, i.e. an areas with uncertainity according to SDMs that it is not included in the phylo range, but this is a artificially reduction of phylo suitability if that zero is false, an area with higher than 75 or lower than 25 regular suitability that was set as zero for using it as background. 

                #reduce resolution
                phylo_suit = resample(phylo_suit, sum_distributions, method="bilinear")

                #mask the future suitability with the polyong distribution
                phylo_suit_masked = mask(phylo_suit, ocurrences_buffer_polygon)

                #select only areas inside of the phylo range
                phylo_suit_masked[which(getValues(phylo_suit_masked)==0)] <- NA
                phylo_suit_masked[which(!is.na(getValues(phylo_suit_masked)))] <- 1

                #select unsuitable cells included in the phylo niche and considered as unsuitable in the regular niche, and rest them to the unsuitable cells of the regular niche
                sum_suitable_phylo_no_phylo = abs(length(which(!is.na(getValues(phylo_suit_masked)) & !is.na(getValues(future_suit_masked)))) - length(which(getValues(future_suit_masked) > 0))) #the difference will be the the cells non-suitable for both regular and phylo niche.

                #drop to the suitable areas, the suitable areas according to phylo niche
                future_suit_masked[which(getValues(phylo_suit_masked)==1)] <- NA

                #cells suitable future
                suitable_cells_future = length(which(getValues(future_suit_masked) > 0))

                #check if cells included in the phylo niche but considered as unsuitable by the regular niche, rested to unsuitable areas of the regular niche gives the same number than unsuitable cells in the final future_suit_masked.
                print("Test is ok:")
                print(suitable_cells_future == sum_suitable_phylo_no_phylo)
                print("")
            }

            #if there are cells with low of suitability in future
            if(length(which(!is.na(getValues(future_suit_masked)))) > 0){ #the conditional is only applied to future suit masked because we are sure that current climate inside distribution is suitable for species (occurrences inside distribution lead to that), so is unlikely to get zero cell with any data under current climate.

                #multiply to get areas suitable currently but not suitable in the future (map not included in any of these maps will bi removed; NA propagation)
                suit_change = current_suit_masked*future_suit_masked

                #if the resulting crossing of both raster has numbers, i.e. areas suitable currently than are predicted to be below 25 (not all is na)
                if(length(which(!is.na(getValues(suit_change)))) > 0){

                    #create a polygon from suit_change
                    polygon_suit_change = rasterToPolygons(suit_change, fun=function(x){x==1}, n=16, dissolve = TRUE) #convertimos en poligono the raster using the cell with values=1

                    #copy clay
                    clay_to_mask = clay

                    #set as 1 those areas included in the polygon of loss
                    clay_to_mask[polygon_suit_change] <- 1
                } else { #if not, i.e, non suitable areas in the futuere are non suitable yet today

                    #all remains as zero
                    clay_to_mask = clay
                }       
            } else { #if not

                #all remains as zero
                clay_to_mask = clay
            }
            #save into an stack with the speceis name
            suit_change_stack = stack(suit_change_stack, clay_to_mask)
            names(suit_change_stack[[i]]) <- selected_epi
        }
        #loss of suitability below 75%
        if(type=="less_75"){
            #set as NA cells with suit lower than 75 and as 1 cells with higher or equal than 75. 1 for suitable.
            current_suit_masked[which(getValues(current_suit_masked) < 75)] <- NA 
            current_suit_masked[which(!is.na(getValues(current_suit_masked)))] <- 1 

            #load the [i] future suitability
            future_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))
        
            #reduce resolution
            future_suit = resample(future_suit, sum_distributions, method="bilinear")#there is no problem with the resampling in this case because the cases ith zero are real data, are areas with no suitability. The problem occured in phylo-suitability rasters because all NAs were converted into zero, so phylo-suitable areas close to false-zeros were artificially reduced its suitbaility. 

            #mask the future suitability with the polyong distribution
            future_suit_masked = mask(future_suit, ocurrences_buffer_polygon)

            #set as NA cells with suit lower or equal than 75 and as 1 cells with higher than 75. 1 for unsuitble
            future_suit_masked[which(getValues(future_suit_masked) > 75)] <- NA
            future_suit_masked[which(!is.na(getValues(future_suit_masked)))] <- 1 

            #if we don't want phylo niche
            if(phylo_cor=="yes"){

                #load again future_suit because we want it without resampling. This raster will be used to select those areas of the phylo raster that are considered as intermediate, and set NAs all outside of that range. This is important to avoid FALSE zeros in the resampling process, these zeros would artificially reduce the phylo suitaiblity
                future_suit_to_subset = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))

                #load phylo niche
                phylo_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/phylo_ensamble/with_proportions/", selected_epi, "_phylo_ensamble_with_proportions.asc", sep="")) 

                #create a raster with only 1 in areas with suitability uncertainity in future_suit_to_subset, the rest NA
                phylo_mask = phylo_suit
                phylo_mask[which(getValues(future_suit_to_subset) > 25 & getValues(future_suit_to_subset) < 75)] <- 1 #select all uncertainity areas
                phylo_suit[which(getValues(phylo_mask)==0)] <- NA #the rest of areas will be NA in phylo_suit. In this way, we remove the false zeros (all phylo_suit is cover by zero except in those areas where phylo-suit is higher than zero), an those zeros area not used in the resampling. Think that if you have close two cells one with 0.3 and other of zero, the bigger cells under lower resolution cannot be 0.3 neither 0, it should be 0.15.... this is ok if the zero is real, i.e. an areas with uncertainity according to SDMs that it is not included in the phylo range, but this is a artificially reduction of phylo suitability if that zero is false, an area with higher than 75 or lower than 25 regular suitability that was set as zero for using it as background. 

                #reduce resolution
                phylo_suit = resample(phylo_suit, sum_distributions, method="bilinear")

                #mask the future suitability with the polyong distribution
                phylo_suit_masked = mask(phylo_suit, ocurrences_buffer_polygon)

                #select only areas inside of the phylo range
                phylo_suit_masked[which(getValues(phylo_suit_masked)==0)] <- NA
                phylo_suit_masked[which(!is.na(getValues(phylo_suit_masked)))] <- 1

                #sum suitable cells included in the phyloniche but considered as unsuitable in the regular niche, and rest them to the unsuitable cells of the regular niche
                sum_suitable_phylo_no_phylo = abs(length(which(!is.na(getValues(phylo_suit_masked)) & !is.na(getValues(future_suit_masked)))) - length(which(getValues(future_suit_masked) > 0)))#the difference will be the the cells non-suitable for both regular and phylo niche.

                #drop to the suitable areas, the suitable areas according to phylo niche
                future_suit_masked[which(getValues(phylo_suit_masked)==1)] <- NA

                #cells suitable future
                suitable_cells_future = length(which(getValues(future_suit_masked) > 0))

                #check if cells included in the phylo niche but considered as unsuitable by the regular niche, rested to unsuitable areas of the regular niche gives the same number than unsuitable cells in the final future_suit_masked.
                print("Test is ok:")
                print(suitable_cells_future == sum_suitable_phylo_no_phylo)
                print("")
            }

            #if there are cells with low of suitability in future
            if(length(which(!is.na(getValues(future_suit_masked)))) > 0){#the conditional is only applied to future suit masked because we are sure that current climate inside distribution is suitable for species (occurrences inside distribution lead to that), so is unlikely to get zero cell with any data under current climate.

                #multiply to get areas suitables currently but not suitable in the future (map not included in any of these maps will bi removed; NA propagation)
                suit_change = current_suit_masked*future_suit_masked 

                #if the resulting crossing of both raster has numbers, i.e. areas suitable currently than are predicted to be lost (not all is na)
                if(length(which(!is.na(getValues(suit_change)))) > 0){

                    #create a polygon from suit_change
                    polygon_suit_change = rasterToPolygons(suit_change, fun=function(x){x==1}, n=16, dissolve = TRUE) #convertimos en poligono the raster using the cell with values=1

                    #copy clay
                    clay_to_mask = clay

                    #set as 1 those areas included in the polygon of loss
                    clay_to_mask[polygon_suit_change] <- 1
                } else { #if not, i.e, non suitable areas in the futuere are non suitable yet today

                    #all remains as zero
                    clay_to_mask = clay
                }       
            } else { #if not

                #all remains as zero
                clay_to_mask = clay
            }
            #save into an stack with the speceis name
            suit_change_stack = stack(suit_change_stack, clay_to_mask)
            names(suit_change_stack[[i]]) <- selected_epi
        }
        #increase of suitability over 75%
        if(type=="more_75"){
            #set as NA cells with suit higher or equal than 75 and as 1 cells with lower than 75. 1 for unsuitable.
            current_suit_masked[which(getValues(current_suit_masked) >= 75)] <- NA
            current_suit_masked[which(!is.na(getValues(current_suit_masked)))] <- 1

            #load the [i] future suitability
            future_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))
        
            #reduce resolution
            future_suit = resample(future_suit, sum_distributions, method="bilinear")#there is no problem with the resampling in this case because the cases ith zero are real data, are areas with no suitability. The problem occured in phylo-suitability rasters because all NAs were converted into zero, so phylo-suitable areas close to false-zeros were artificially reduced its suitbaility. 

            #mask the future suitability with the polyong distribution
            future_suit_masked = mask(future_suit, ocurrences_buffer_polygon)

            #set as NA cells with suit lower than 75 and as 1 cells with higher than 75. 1 for unsuitable
            future_suit_masked[which(getValues(future_suit_masked) < 75)] <- NA 
            future_suit_masked[which(!is.na(getValues(future_suit_masked)))] <- 1

            #if we don't want phylo niche
            if(phylo_cor=="yes"){
 
                #load again future_suit because we want it without resampling. This raster will be used to select those areas of the phylo raster that are considered as intermediate, and set NAs all outside of that range. This is important to avoid FALSE zeros in the resampling process, these zeros would artificially reduce the phylo suitaiblity
                future_suit_to_subset = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))

                #load phylo niche
                phylo_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/phylo_ensamble/with_proportions/", selected_epi, "_phylo_ensamble_with_proportions.asc", sep="")) 

                #create a raster with only 1 in areas with suitability uncertainity in future_suit_to_subset, the rest NA
                phylo_mask = phylo_suit
                phylo_mask[which(getValues(future_suit_to_subset) > 25 & getValues(future_suit_to_subset) < 75)] <- 1 #select all uncertainity areas
                phylo_suit[which(getValues(phylo_mask)==0)] <- NA #the rest of areas will be NA in phylo_suit. In this way, we remove the false zeros (all phylo_suit is cover by zero except in those areas where phylo-suit is higher than zero), an those zeros area not used in the resampling. Think that if you have close two cells one with 0.3 and other of zero, the bigger cells under lower resolution cannot be 0.3 neither 0, it should be 0.15.... this is ok if the zero is real, i.e. an areas with uncertainity according to SDMs that it is not included in the phylo range, but this is a artificially reduction of phylo suitability if that zero is false, an area with higher than 75 or lower than 25 regular suitability that was set as zero for using it as background. 

                #reduce resolution
                phylo_suit = resample(phylo_suit, sum_distributions, method="bilinear")

                #mask the future suitability with the polyong distribution
                phylo_suit_masked = mask(phylo_suit, ocurrences_buffer_polygon)

                #select only areas inside of the phylo range
                phylo_suit_masked[which(getValues(phylo_suit_masked)==0)] <- NA
                phylo_suit_masked[which(!is.na(getValues(phylo_suit_masked)))] <- 1

                #sub suitable cells included in the phyloniche but not in the regular niche, and sum them to the suitable cells of the regular niche
                sum_suitable_phylo_no_phylo = length(which(!is.na(getValues(phylo_suit_masked)) & is.na(getValues(future_suit_masked)))) + length(which(getValues(future_suit_masked) > 0))

                #add to the suitable areas, the suitable areas according to phylo niche
                future_suit_masked[which(getValues(phylo_suit_masked)==1)] <- 1

                #cells suitable future
                suitable_cells_future = length(which(getValues(future_suit_masked) > 0))

                #check if cells included in the phylo niche but not in the regular niche, summed to suitable areas of the regular niche gives the same number than suitable cells in the final future_suit_masked.
                print("Test is ok:")
                print(suitable_cells_future == sum_suitable_phylo_no_phylo)
                print("")
            }

            #if there are cells with high of suitability in future
            if(length(which(!is.na(getValues(future_suit_masked)) & !is.na(getValues(current_suit_masked)))) > 0){ #in this case we are removing suitable areas from current maps, so maybe it is possible that the full range is suitable today. Because of this we include the conditional for current. 

                #multiply to get areas suitable currently but not suitable in the future
                suit_change = current_suit_masked*future_suit_masked

                #if the resulting crossing of both raster has numbers, i.e. areas non-suitable currently than are predicted to be suitable in future (not all is na)
                if(length(which(!is.na(getValues(suit_change)))) > 0){

                    #create a polygon from suit_change
                    polygon_suit_change = rasterToPolygons(suit_change, fun=function(x){x==1}, n=16, dissolve = TRUE) #convertimos en poligono the raster using the cell with values=1

                    #copy clay
                    clay_to_mask = clay

                    #set as 1 those areas included in the polygon of loss
                    clay_to_mask[polygon_suit_change] <- 1
                } else { #if not, i.e, non suitable areas in the futuere are non suitable yet today

                    #all remains as zero
                    clay_to_mask = clay
                } 
            } else { #if not

                #all remains as zero
                clay_to_mask = clay
            }
            #save into an stack with the speceis name
            suit_change_stack = stack(suit_change_stack, clay_to_mask)
            names(suit_change_stack[[i]]) <- selected_epi
        }
        #loss of suitability below 75% between 75-25.
        if(type=="more_75_intermediate"){
            #set as NA cells with suit lower than 75 and as 1 cells with higher or equal than 75. 1 for suitable.
            current_suit_masked[which(getValues(current_suit_masked) < 75)] <- NA 
            current_suit_masked[which(getValues(current_suit_masked) > 75)] <- 1 

            #load the [i] future suitability
            future_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))
        
            #reduce resolution
            future_suit = resample(future_suit, sum_distributions, method="bilinear")#there is no problem with the resampling in this case because the cases ith zero are real data, are areas with no suitability. The problem occured in phylo-suitability rasters because all NAs were converted into zero, so phylo-suitable areas close to false-zeros were artificially reduced its suitbaility.

            #mask the future suitability with the polyong distribution
            future_suit_masked = mask(future_suit, ocurrences_buffer_polygon)

            #set as NA cells with suit higher than 75 or lower than 25
            future_suit_masked[which(getValues(future_suit_masked) >= 75 | getValues(future_suit_masked) <= 25)] <- NA
            #set as 1 the rest of cells (intermediates)
            future_suit_masked[which(!is.na(getValues(future_suit_masked)))] <- 1

            #if we want phylo niche
            if(phylo_cor=="yes"){

                #load again future_suit because we want it without resampling. This raster will be used to select those areas of the phylo raster that are considered as intermediate, and set NAs all outside of that range. This is important to avoid FALSE zeros in the resampling process, these zeros would artificially reduce the phylo suitaiblity
                future_suit_to_subset = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))

                #load phylo niche
                phylo_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/phylo_ensamble/with_proportions/", selected_epi, "_phylo_ensamble_with_proportions.asc", sep="")) 

                #create a raster with only 1 in areas with suitability uncertainity in future_suit_to_subset, the rest NA
                phylo_mask = phylo_suit
                phylo_mask[which(getValues(future_suit_to_subset) > 25 & getValues(future_suit_to_subset) < 75)] <- 1 #select all uncertainity areas
                phylo_suit[which(getValues(phylo_mask)==0)] <- NA #the rest of areas will be NA in phylo_suit. In this way, we remove the false zeros (all phylo_suit is cover by zero except in those areas where phylo-suit is higher than zero), an those zeros area not used in the resampling. Think that if you have close two cells one with 0.3 and other of zero, the bigger cells under lower resolution cannot be 0.3 neither 0, it should be 0.15.... this is ok if the zero is real, i.e. an areas with uncertainity according to SDMs that it is not included in the phylo range, but this is a artificially reduction of phylo suitability if that zero is false, an area with higher than 75 or lower than 25 regular suitability that was set as zero for using it as background. 

                #reduce resolution
                phylo_suit = resample(phylo_suit, sum_distributions, method="bilinear")

                #mask the future suitability with the polyong distribution
                phylo_suit_masked = mask(phylo_suit, ocurrences_buffer_polygon)

                #select only areas inside of the phylo range
                phylo_suit_masked[which(getValues(phylo_suit_masked)==0)] <- NA
                phylo_suit_masked[which(!is.na(getValues(phylo_suit_masked)))] <- 1

                #sum suitable cells included in the phyloniche but considered as unsuitable in the regular niche, and rest them to the unsuitable cells of the regular niche
                sum_suitable_phylo_no_phylo = abs(length(which(!is.na(getValues(phylo_suit_masked)) & !is.na(getValues(future_suit_masked)))) - length(which(getValues(future_suit_masked) > 0)))

                #drop to the suitable areas, the suitable areas according to phylo niche
                future_suit_masked[which(getValues(phylo_suit_masked)==1)] <- NA

                #cells suitable future
                suitable_cells_future = length(which(getValues(future_suit_masked) > 0))

                #check if cells included in the phylo niche but considered as unsuitable by the regular niche, rested to unsuitable areas of the regular niche gives the same number than unsuitable cells in the final future_suit_masked.
                print("Test is ok:")
                print(suitable_cells_future == sum_suitable_phylo_no_phylo)
                print("")
            }

            #if there are cells with intermediate of suitability in future
            if(length(which(!is.na(getValues(future_suit_masked)))) > 0){

                #multiply to get areas suitables currently but not suitable in the future (intermediate)
                suit_change = current_suit_masked*future_suit_masked

                #if the resulting crossing of both raster has numbers, i.e. areas suitable currently than are predicted to be lost (not all is na)
                if(length(which(!is.na(getValues(suit_change)))) > 0){

                    #create a polygon from suit_change
                    polygon_suit_change = rasterToPolygons(suit_change, fun=function(x){x==1}, n=16, dissolve = TRUE) #convertimos en poligono the raster using the cell with values=1

                    #copy clay
                    clay_to_mask = clay

                    #set as 1 those areas included in the polygon of loss
                    clay_to_mask[polygon_suit_change] <- 1
                } else { #if not, i.e, non suitable areas in the futuere are non suitable yet today

                    #all remains as zero
                    clay_to_mask = clay
                }
            } else { #if not

                #all remains as zero
                clay_to_mask = clay
            }
            #save into an stack with the speceis name
            suit_change_stack = stack(suit_change_stack, clay_to_mask)
            names(suit_change_stack[[i]]) <- selected_epi
        }
        #increase of suitability from 25% to intermediate.
        if(type=="less_25_intermediate"){
            #set as NA cells with suit higher or equal than 25 and as 1 cells with lower than 75. 1 for non-suitable.
            current_suit_masked[which(getValues(current_suit_masked) > 25)] <- NA 
            current_suit_masked[which(!is.na(getValues(current_suit_masked)))] <- 1 

            #load the [i] future suitability
            future_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))
        
            #reduce resolution
            future_suit = resample(future_suit, sum_distributions, method="bilinear")#there is no problem with the resampling in this case because the cases ith zero are real data, are areas with no suitability. The problem occured in phylo-suitability rasters because all NAs were converted into zero, so phylo-suitable areas close to false-zeros were artificially reduced its suitbaility.

            #mask the future suitability with the polyong distribution
            future_suit_masked = mask(future_suit, ocurrences_buffer_polygon)

            #set as NA cells with suit higher than 75 or lower than 25
            future_suit_masked[which(getValues(future_suit_masked) >= 75 | getValues(future_suit_masked) <= 25)] <- NA
            #set as 1 the rest of cells (intermediate)
            future_suit_masked[which(!is.na(getValues(future_suit_masked)))] <- 1

            #if we don't want phylo niche
            if(phylo_cor=="yes"){

                #load again future_suit because we want it without resampling. This raster will be used to select those areas of the phylo raster that are considered as intermediate, and set NAs all outside of that range. This is important to avoid FALSE zeros in the resampling process, these zeros would artificially reduce the phylo suitaiblity
                future_suit_to_subset = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))

                #load phylo niche
                phylo_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/phylo_ensamble/with_proportions/", selected_epi, "_phylo_ensamble_with_proportions.asc", sep="")) 

                #create a raster with only 1 in areas with suitability uncertainity in future_suit_to_subset, the rest NA
                phylo_mask = phylo_suit
                phylo_mask[which(getValues(future_suit_to_subset) > 25 & getValues(future_suit_to_subset) < 75)] <- 1 #select all uncertainity areas
                phylo_suit[which(getValues(phylo_mask)==0)] <- NA #the rest of areas will be NA in phylo_suit. In this way, we remove the false zeros (all phylo_suit is cover by zero except in those areas where phylo-suit is higher than zero), an those zeros area not used in the resampling. Think that if you have close two cells one with 0.3 and other of zero, the bigger cells under lower resolution cannot be 0.3 neither 0, it should be 0.15.... this is ok if the zero is real, i.e. an areas with uncertainity according to SDMs that it is not included in the phylo range, but this is a artificially reduction of phylo suitability if that zero is false, an area with higher than 75 or lower than 25 regular suitability that was set as zero for using it as background. 

                #reduce resolution
                phylo_suit = resample(phylo_suit, sum_distributions, method="bilinear")

                #mask the future suitability with the polyong distribution
                phylo_suit_masked = mask(phylo_suit, ocurrences_buffer_polygon)

                #select only areas inside of the phylo range
                phylo_suit_masked[which(getValues(phylo_suit_masked)==0)] <- NA
                phylo_suit_masked[which(!is.na(getValues(phylo_suit_masked)))] <- 1

                #sum suitable cells included in the phyloniche but considered as unsuitable in the regular niche, and rest them to the unsuitable cells of the regular niche
                sum_suitable_phylo_no_phylo = abs(length(which(!is.na(getValues(phylo_suit_masked)) & !is.na(getValues(future_suit_masked)))) - length(which(getValues(future_suit_masked) > 0)))

                #drop to the suitable areas, the suitable areas according to phylo niche
                future_suit_masked[which(getValues(phylo_suit_masked)==1)] <- NA

                #cells suitable future
                suitable_cells_future = length(which(getValues(future_suit_masked) > 0))

                #check if cells included in the phylo niche but considered as unsuitable by the regular niche, rested to unsuitable areas of the regular niche gives the same number than unsuitable cells in the final future_suit_masked.
                print("Test is ok:")
                print(suitable_cells_future == sum_suitable_phylo_no_phylo)
                print("")
            }

            #if there are cells with intermediate of suitability in future
            if(length(which(!is.na(getValues(future_suit_masked)) & !is.na(getValues(current_suit_masked)))) > 0){#in this case we are removing suitable areas from current maps, so maybe it is possible that the full range is suitable today. Because of this we include the conditional for current.

                #multiply to get areas suitables currently but not suitable in the future
                suit_change = current_suit_masked*future_suit_masked

                #if the resulting crossing of both raster has numbers, i.e. areas suitable currently than are predicted to be lost (not all is na)
                if(length(which(!is.na(getValues(suit_change)))) > 0){

                    #create a polygon from suit_change
                    polygon_suit_change = rasterToPolygons(suit_change, fun=function(x){x==1}, n=16, dissolve = TRUE) #convertimos en poligono the raster using the cell with values=1

                    #copy clay
                    clay_to_mask = clay

                    #set as 1 those areas included in the polygon of loss
                    clay_to_mask[polygon_suit_change] <- 1
                } else { #if not, i.e, non suitable areas in the futuere are non suitable yet today

                    #all remains as zero
                    clay_to_mask = clay
                }
            } else { #if not

                #all remains as zero
                clay_to_mask = clay
            }

            #save into an stack with the speceis name
            suit_change_stack = stack(suit_change_stack, clay_to_mask)
            names(suit_change_stack[[i]]) <- selected_epi
        }
        #possibility of migration outside distribution
        if(type=="migration"){

            #load the [i] future suitability
            future_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))
        
            #reduce resolution
            future_suit = resample(future_suit, sum_distributions, method="bilinear")#there is no problem with the resampling in this case because the cases ith zero are real data, are areas with no suitability. The problem occured in phylo-suitability rasters because all NAs were converted into zero, so phylo-suitable areas close to false-zeros were artificially reduced its suitbaility.

            #create a polygon buffer around the distribution buffer
            polygon_migration_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=10, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) 

            #mask the future suitability with the polyong distribution
            future_suit_masked = mask(future_suit, polygon_migration_buffer)

            #drop areas inside distribution
            future_suit_masked = mask(future_suit_masked, ocurrences_buffer_polygon, inverse=TRUE)

            #set as NA cells with suit lower than 75 and as 1 cells with higher or equal than 75. 1 for suitable
            future_suit_masked[which(getValues(future_suit_masked) < 75)] <- NA
            future_suit_masked[which(!is.na(getValues(future_suit_masked)))] <- 1 

            #if we don't want phylo niche
            if(phylo_cor=="yes"){
 
                #load again future_suit because we want it without resampling. This raster will be used to select those areas of the phylo raster that are considered as intermediate, and set NAs all outside of that range. This is important to avoid FALSE zeros in the resampling process, these zeros would artificially reduce the phylo suitaiblity
                future_suit_to_subset = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))

                #load phylo niche
                phylo_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/phylo_ensamble/with_proportions/", selected_epi, "_phylo_ensamble_with_proportions.asc", sep="")) 

                #create a raster with only 1 in areas with suitability uncertainity in future_suit_to_subset, the rest NA
                phylo_mask = phylo_suit
                phylo_mask[which(getValues(future_suit_to_subset) > 25 & getValues(future_suit_to_subset) < 75)] <- 1 #select all uncertainity areas
                phylo_suit[which(getValues(phylo_mask)==0)] <- NA #the rest of areas will be NA in phylo_suit. In this way, we remove the false zeros (all phylo_suit is cover by zero except in those areas where phylo-suit is higher than zero), an those zeros area not used in the resampling. Think that if you have close two cells one with 0.3 and other of zero, the bigger cells under lower resolution cannot be 0.3 neither 0, it should be 0.15.... this is ok if the zero is real, i.e. an areas with uncertainity according to SDMs that it is not included in the phylo range, but this is a artificially reduction of phylo suitability if that zero is false, an area with higher than 75 or lower than 25 regular suitability that was set as zero for using it as background. 

                #reduce resolution
                phylo_suit = resample(phylo_suit, sum_distributions, method="bilinear")

                #mask the future suitability with the polyong distribution
                phylo_suit_masked = mask(phylo_suit, polygon_migration_buffer)

                #drop areas inside distribution
                phylo_suit_masked = mask(phylo_suit_masked, ocurrences_buffer_polygon, inverse=TRUE)                

                #select only areas inside of the phylo range
                phylo_suit_masked[which(getValues(phylo_suit_masked)==0)] <- NA
                phylo_suit_masked[which(!is.na(getValues(phylo_suit_masked)))] <- 1

                #sub suitable cells included in the phyloniche but not in the regular niche, and sum them to the suitable cells of the regular niche
                sum_suitable_phylo_no_phylo = length(which(!is.na(getValues(phylo_suit_masked)) & is.na(getValues(future_suit_masked)))) + length(which(getValues(future_suit_masked) > 0))

                #add to the suitable areas, the suitable areas according to phylo niche
                future_suit_masked[which(getValues(phylo_suit_masked)==1)] <- 1

                #cells suitable future
                suitable_cells_future = length(which(getValues(future_suit_masked) > 0))

                #check if cells included in the phylo niche but not in the regular niche, summed to suitable areas of the regular niche gives the same number than suitable cells in the final future_suit_masked.
                print("Test is ok:")
                print(suitable_cells_future == sum_suitable_phylo_no_phylo)
                print("")
            }

            #if there are cells with low of suitability in future
            if(length(which(!is.na(getValues(future_suit_masked)))) > 0){

                #create a polygon from future_suit_masked
                polygon_future_suit = rasterToPolygons(future_suit_masked, fun=function(x){x==1}, n=16, dissolve = TRUE) #convertimos en poligono the raster using the cell with values=1

                #copy clay
                clay_to_mask = clay

                #set as 1 those areas included in the polygon of loss
                clay_to_mask[polygon_future_suit] <- 1
            } else { #if not

                #all remains as zero
                clay_to_mask = clay
            }
            #save into an stack with the speceis name
            suit_change_stack = stack(suit_change_stack, clay_to_mask)
            names(suit_change_stack[[i]]) <- selected_epi
        }
        #possibility of migration considering also suitability inside distribution
        if(type=="total_more_75_distri_migra"){
 
            #load the [i] future suitability
            future_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))
        
            #reduce resolution
            future_suit = resample(future_suit, sum_distributions, method="bilinear")#there is no problem with the resampling in this case because the cases ith zero are real data, are areas with no suitability. The problem occured in phylo-suitability rasters because all NAs were converted into zero, so phylo-suitable areas close to false-zeros were artificially reduced its suitbaility.

            #create a polygon buffer around the distribution buffer
            polygon_migration_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=10, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

            #mask the future suitability with the polyong distribution
            future_suit_masked = mask(future_suit, polygon_migration_buffer)

            #set as NA cells with suit lower than 75 and as 1 cells with higher or equal than 75. 1 for suitable
            future_suit_masked[which(getValues(future_suit_masked) < 75)] <- NA
            future_suit_masked[which(!is.na(getValues(future_suit_masked)))] <- 1

            #if we don't want phylo niche
            if(phylo_cor=="yes"){
 
                #load again future_suit because we want it without resampling. This raster will be used to select those areas of the phylo raster that are considered as intermediate, and set NAs all outside of that range. This is important to avoid FALSE zeros in the resampling process, these zeros would artificially reduce the phylo suitaiblity
                future_suit_to_subset = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", selected_epi, ".tif", sep=""))

                #load phylo niche
                phylo_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/phylo_ensamble/with_proportions/", selected_epi, "_phylo_ensamble_with_proportions.asc", sep="")) 

                #create a raster with only 1 in areas with suitability uncertainity in future_suit_to_subset, the rest NA
                phylo_mask = phylo_suit
                phylo_mask[which(getValues(future_suit_to_subset) > 25 & getValues(future_suit_to_subset) < 75)] <- 1 #select all uncertainity areas
                phylo_suit[which(getValues(phylo_mask)==0)] <- NA #the rest of areas will be NA in phylo_suit. In this way, we remove the false zeros (all phylo_suit is cover by zero except in those areas where phylo-suit is higher than zero), an those zeros area not used in the resampling. Think that if you have close two cells one with 0.3 and other of zero, the bigger cells under lower resolution cannot be 0.3 neither 0, it should be 0.15.... this is ok if the zero is real, i.e. an areas with uncertainity according to SDMs that it is not included in the phylo range, but this is a artificially reduction of phylo suitability if that zero is false, an area with higher than 75 or lower than 25 regular suitability that was set as zero for using it as background. 

                #reduce resolution
                phylo_suit = resample(phylo_suit, sum_distributions, method="bilinear")

                #mask the future suitability with the polyong distribution
                phylo_suit_masked = mask(phylo_suit, polygon_migration_buffer)

                #select only areas inside of the phylo range
                phylo_suit_masked[which(getValues(phylo_suit_masked)==0)] <- NA
                phylo_suit_masked[which(!is.na(getValues(phylo_suit_masked)))] <- 1

                #sub suitable cells included in the phyloniche but not in the regular niche, and sum them to the suitable cells of the regular niche
                sum_suitable_phylo_no_phylo = length(which(!is.na(getValues(phylo_suit_masked)) & is.na(getValues(future_suit_masked)))) + length(which(getValues(future_suit_masked) > 0))

                #add to the suitable areas, the suitable areas according to phylo niche
                future_suit_masked[which(getValues(phylo_suit_masked)==1)] <- 1

                #cells suitable future
                suitable_cells_future = length(which(getValues(future_suit_masked) > 0))

                #check if cells included in the phylo niche but not in the regular niche, summed to suitable areas of the regular niche gives the same number than suitable cells in the final future_suit_masked.
                print("Test is ok:")
                print(suitable_cells_future == sum_suitable_phylo_no_phylo)
                print("")
            }

            #if there are cells with low of suitability in future
            if(length(which(!is.na(getValues(future_suit_masked)))) > 0){

                #create a polygon from suit_change
                polygon_suit_change = rasterToPolygons(future_suit_masked, fun=function(x){x==1}, n=16, dissolve = TRUE) #convertimos en poligono the raster using the cell with values=1

                #copy clay
                clay_to_mask = clay

                #set as 1 those areas included in the polygon of loss
                clay_to_mask[polygon_suit_change] <- 1
            } else { #if not, i.e, non suitable areas in the futuere are non suitable yet today

                #all remains as zero
                clay_to_mask = clay
            } 
            #save into an stack with the speceis name
            suit_change_stack = stack(suit_change_stack, clay_to_mask)
            names(suit_change_stack[[i]]) <- selected_epi 
        }        
    }

    #save results
    if(!type=="histogram"){
        #nlayers of the stack
        nlayers(suit_change_stack) ==112 #without discolor

        #sum the suitability loss across species
        suit_change_stack_sum = calc(x=suit_change_stack, fun=function(x) (sum(x)))

        #if the type plot is not migration
        if(!type %in% c("migration", "total_more_75_distri_migra")){
            
            if(!type %in% c("sum_suitability_outside", "sum_suitability_phylo_inout", "sum_suitability")){
                #set relative to the number of pines in each site (percentage) using the sum of distribution buffers
                suit_change_stack_sum_proport = (suit_change_stack_sum*100)/sum_distributions

                #drop areas without pines from suit_change_stack_sum 
                suit_change_stack_sum[which(is.na(getValues(sum_distributions)))] <- NA
                suit_change_stack_sum_proport[which(is.na(getValues(sum_distributions)))] <- NA                
            } else {
                if(type=="sum_suitability"){
                    #set relative to the number of pines in each site (percentage) using the sum of distribution buffers
                    suit_change_stack_sum_proport = (suit_change_stack_sum*100)/(sum_distributions*100) #if for example we have 6 species in an area, the maxmium suitaiblity is 6*100=600. Therefore, if 600 is 100, 300 will be x; (300*100)/600

                    #drop areas without pines from suit_change_stack_sum 
                    suit_change_stack_sum[which(is.na(getValues(sum_distributions)))] <- NA #we have zero over all the globe, so we have to drop areas without pines
                    suit_change_stack_sum_proport[which(is.na(getValues(sum_distributions)))] <- NA                    
                } 
                if(type=="sum_suitability_phylo_inout"){
                    #set relative to the number of pines in each site (percentage) using the sum of migration buffers (including distribution)
                    suit_change_stack_sum_proport = (suit_change_stack_sum*100)/(sum_buffers_distrib) #if for example we have 6 species in an area, the maxmium phylo suitaiblity is 6. Therefore, if 6 is 100, 3 will be x; (3*100)/6

                    #drop areas without pines from suit_change_stack_sum 
                    suit_change_stack_sum[which(is.na(getValues(sum_buffers_distrib)))] <- NA #we have zero over all the globe, so we have to drop areas without pines
                    suit_change_stack_sum_proport[which(is.na(getValues(sum_buffers_distrib)))] <- NA                                          
                }
                if(type=="sum_suitability_outside"){
                    #set relative to the number of pines in each site (percentage) using the sum of migration buffers (without distribution)
                    suit_change_stack_sum_proport = (suit_change_stack_sum*100)/(sum_buffers*100) #if for example we have 6 species in an area, the maxmium suitaiblity is 6*100=600. Therefore, if 600 is 100, 300 will be x; (300*100)/600

                    #drop areas without pines from suit_change_stack_sum 
                    suit_change_stack_sum[which(is.na(getValues(sum_buffers)))] <- NA #we have zero over all the globe, so we have to drop areas without pines
                    suit_change_stack_sum_proport[which(is.na(getValues(sum_buffers)))] <- NA   
                } 
            }    
        } else { #if not

            if(type %in% c("migration")){
                #set relative to the number of pines in each site (percentage) using the sum of migration buffers (not including distribution)
                suit_change_stack_sum_proport = (suit_change_stack_sum*100)/sum_buffers

                #drop areas without pines from suit_change_stack_sum 
                suit_change_stack_sum[which(is.na(getValues(sum_buffers)))] <- NA
                suit_change_stack_sum_proport[which(is.na(getValues(sum_buffers)))] <- NA                
            } else{ #and thus it is total_more_75_distri_migra (migration including distribution area)
            
                #set relative to the number of pines in each site (percentage) using the sum of migration buffers (including distributionarea)
                suit_change_stack_sum_proport = (suit_change_stack_sum*100)/sum_buffers_distrib

                #drop areas without pines from suit_change_stack_sum 
                suit_change_stack_sum[which(is.na(getValues(sum_buffers_distrib)))] <- NA
                suit_change_stack_sum_proport[which(is.na(getValues(sum_buffers_distrib)))] <- NA                
            }
        }

        #plot absolute and relative results
        pdf(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/globalPlot_", type, "_phylo_cor_", phylo_cor,  ".pdf", sep=""), width=14, height=6)
        par(mfcol=c(1,2), mai=c(0.5,0.5,0.5,1))
        plot(crop(clay, c(-180,180,-10,90)), col="gray", legend=FALSE, main="Absolute number of species")
        plot(crop(suit_change_stack_sum, c(-180,180,0,90)), add=TRUE, col=colorRampPalette(c("yellow", "red"))(40))
        plot(crop(clay, c(-180,180,-10,90)), col="gray", legend=FALSE, main="Percentage of species")
        plot(crop(suit_change_stack_sum_proport, c(-180,180,0,90)), add=TRUE, col=colorRampPalette(c("yellow", "red"))(40))
        dev.off()

        #save rasters
        writeRaster(suit_change_stack, filename=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/suit_change_stack_", type, "_phylo_cor_", phylo_cor,  ".tif", sep=""), options="INTERLEAVE=BAND", overwrite=TRUE)
        writeRaster(suit_change_stack_sum, filename=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/suit_change_stack_abs_", type, "_phylo_cor_", phylo_cor,  ".asc", sep=""), overwrite=TRUE)
        writeRaster(suit_change_stack_sum_proport, filename=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/suit_change_stack_sum_proport_", type, "_phylo_cor_", phylo_cor,  ".asc", sep=""), overwrite=TRUE)  
    } else {

        #removes the first row with NAs
        percent_loss = percent_loss[-1,]
        nrow(percent_loss) == 112 #without discolor

        #save table
        write.table(percent_loss, paste("/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/percentage_data_phylo_cor_", phylo_cor, ".csv", sep=""), sep=",", row.names = FALSE, col.names=TRUE)

        #plot
        pdf(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/globalPlot_", type, "_phylo_cor_", phylo_cor, ".pdf", sep=""), width=8, height=8)
        hist(percent_loss$percent_loss, xlab="Suitability loss (%)", main="", cex.lab=1.4)
        dev.off()
    }
}

########Paralelize the process######
require(foreach)
require(doParallel) #for parallel

#create a vector with species names
type_plots = c(
    "histogram",
    "sum_suitability",
    "sum_suitability_outside",
    "sum_suitability_phylo_inout",
    "total_more_75",
    "more_75_less_25",
    "less_75",
    "more_75",
    "more_75_intermediate",
    "less_25_intermediate",
    "migration",
    "total_more_75_distri_migra")

# set up cluster
clust <- makeCluster(11) 
registerDoParallel(clust)

#plot all type of plots with and without phylogenetic correction
foreach(i=type_plots, .packages=c("raster", "rgeos")) %dopar% { 

    if(!i %in% c("sum_suitability", "sum_suitability_outside", "sum_suitability_phylo_inout")){ #except in the case of these three plots
        synthetic_plots(type=i, phylo_cor="no")    
        synthetic_plots(type=i, phylo_cor="yes")
    } else {
        if(i %in% c("sum_suitability", "sum_suitability_outside")){ #the first two only with no phylo, and the second with yes phylo
            synthetic_plots(type=i, phylo_cor="no")    
        }else{
            synthetic_plots(type=i, phylo_cor="yes")            
        }
    } 
}

#stop the cluster 
stopCluster(clust)