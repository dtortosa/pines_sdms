###########################################################################
######### Extraction and preparation of climatic and species data #########
###########################################################################

#load packages
library(raster)
library(rgeos)

#function to perform phylo analyses
phylo_analyses = function(species){

    #set function to check if cell values is between ancestral and current value
    is.between <- function(cell_value, ancestral, current) {
        
        #empty vector to save results
        result = NULL
    
        #loop for extractinb results
        for(v in 1:length(cell_value)){ #for each value of cell_value vector
    
            #select the [v] cell_value
            selected_value = cell_value[v]
    
            #test if the [v] cell value is between ancestral an current values
            test = (selected_value - ancestral)  *  (current - selected_value) #code taken from "https://stat.ethz.ch/pipermail/r-help/2008-August/170749.html". The order is irrelevant, ancestral can be higher or lower than current value. Idem for the sign of numbers, it works with only negative, only positive and negative-positive numbers.  
 
            #if test is not zero 
            if(!test == 0){
    
                #test if test is lower or higher than zero to know is the [v] cell value is between current and ancestral values. then save
                result = append(result, test > 0)
    
            } else { #if not, then [v] cell value is equal to the current or ancestral value, but we only want TRUE if the value is equal to the current value. 
    
                #If the [v] cell value is equal to the current value
                if(current == selected_value){
    
                    #result is TRUE
                    result = append(result, TRUE)
                } else { #if not
    
                    #if the [v] cell value is equal to ancestral value
                    if(ancestral == selected_value){
    
                        #result is FALSE
                        result = append(result, FALSE)
    
                    } else {
    
                        #result is NA, problem
                        result = append(result, NA)
    
                    }
                }
            }
        }
    
        #return results
        return(result)
    } #Is very important to add ancestral first, and second current value, becuase TRUE will be returned if the cell_values is equal to "current" (second argument), but FALSE if it is equal to ancestral (first argument)

    #set function to covert the suitibalityi phylo correct to a proportion from 0 to 1
    phylo_proportion = function(x, ancestral_value, current_value){
    
        #calculate the maximum distance to the ancestal value (i.e the current value)
        range_length = abs(current_value - ancestral_value)
    
        #calculate between the cell value and the ancestral value
        distance_to_ancestral = abs(x - ancestral_value)
    
        #if range_length is the 1, distance_to_ancestral will be x; so x = (distance_to_ancestral*1)/range_length 
        distance_to_ancestral/range_length
    } #Like in the latter function, the order is key. The proportion will have 1 as value is close to the second argument (current value).

    #load bio4 currently
    bio4 = raster("/Users/dsalazar/nicho_pinus/data/climatic_data_phylo/bio4.asc")
    res(bio4)
    
    #load bio17 currently
    bio17 = raster("/Users/dsalazar/nicho_pinus/data/climatic_data_phylo/bio17.asc")
    res(bio17)
    
    #list of continuos projections and binary projections for all scenarios and climatic models
    climatic_scenarios = c("bc26", "bc45", "bc60", "bc85", "cc26", "cc45", "cc60", "cc85", "gs26", "gs45", "gs60", "gs85", "he26", "he45", "he60", "he85", "ip26", "ip45", "ip60", "ip85", "mg26", "mg45", "mg60", "mg85", "mr26", "mr45", "mr60", "mr85")

    #load projections of variables for all scenarios
    stack_bio17 = stack(paste("/Users/dsalazar/nicho_pinus/data/climate_proj_phylo/bio17_", climatic_scenarios, ".asc", sep=""))
    stack_bio4 = stack(paste("/Users/dsalazar/nicho_pinus/data/climate_proj_phylo/bio4_", climatic_scenarios, ".asc", sep=""))
    
    #check that the order of scenarios is correct
    names(stack_bio4) == paste("bio4", climatic_scenarios, sep="_") #all TRUE
    names(stack_bio17) == paste("bio17", climatic_scenarios, sep="_") #all TRUE
    
    ## load ancestral reconstruction
    final_anc_ou_bio4 = read.table("/Users/dsalazar/nicho_pinus/data/final_recons/final_anc_ou_bio4.csv", sep=",", header=TRUE)
    final_anc_ou_bio17 = read.table( "/Users/dsalazar/nicho_pinus/data/final_recons/final_anc_ou_bio17.csv", sep=",", header=TRUE)
    final_anc_bm_bio4 = read.table("/Users/dsalazar/nicho_pinus/data/final_recons/final_anc_bm_bio4.csv", sep=",", header=TRUE)
    final_anc_bm_bio17 = read.table( "/Users/dsalazar/nicho_pinus/data/final_recons/final_anc_bm_bio17.csv", sep=",", header=TRUE)
    
    #list of evolution models: ONLY BM
    #list_models_bio4 = c("final_anc_ou_bio4", "final_anc_bm_bio4")
    #list_models_bio17 = c("final_anc_ou_bio17", "final_anc_bm_bio17")
    list_models_bio4 = c("final_anc_bm_bio4")
    list_models_bio17 = c("final_anc_bm_bio17")

    #load calc range buffer to crop climatic variables. This buffer is the one used for migration space
    if(!species=="pumila"){
        calc_range_buffer = raster(paste("/Users/dsalazar/nicho_pinus/data/buffers_range_calc/", species, "_range_calc_buffer.asc", sep=""))
    } else {
        calc_range_buffer = raster(paste("/Users/dsalazar/nicho_pinus/data/buffers_range_calc/", species, "_buffer_range_calc.asc", sep=""))        
    }     
    calc_range_polygon_buffer = rasterToPolygons(calc_range_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon

    #crop climatic variables with calc range buffer
    stack_bio17_cropped = crop(stack_bio17, calc_range_polygon_buffer)
    stack_bio17_cropped = mask(stack_bio17_cropped, calc_range_polygon_buffer)    
    stack_bio4_cropped = crop(stack_bio4, calc_range_polygon_buffer)
    stack_bio4_cropped = mask(stack_bio4_cropped, calc_range_polygon_buffer)

    #load species distribution (withput buffer)
    distri_raster = raster(paste("/Users/dsalazar/nicho_pinus/data/MAPS/p", paste(species, "01.img", sep="_"), sep="_"))
    
    #load ensamble of binary suitability
    ensamble_suitability = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", species, ".tif", sep=""))
        
    #loop for comparing projected climate and phylo range
    phylo_rasters_bio17 = stack()
    phylo_rasters_bio17_proportion = stack()
    for(s in 1:length(climatic_scenarios)){

        #select the climatic scenario
        selected_scenario = climatic_scenarios[s]

        #create a empty raster with the same extent and resolution of the raster layer of the [s] IPCC scenario
        raster_subsetted = raster(extent(stack_bio17_cropped[[s]]), resolution = res(stack_bio17_cropped[[s]]))

        #add to the empty raster those cells of the [s] IPCC raster in which the habitat suitability is higher than 25 and lower than 75
        raster_subsetted[which(getValues(ensamble_suitability) > 25 & getValues(ensamble_suitability) < 75)] <- stack_bio17_cropped[[s]][which(getValues(ensamble_suitability) > 25 & getValues(ensamble_suitability) < 75)] #suitability map has the same resolution and extent than IPCC maps because the suitability map was obtained from these raster of ICPP scenarios. 

        #for each evolution model 
        for(m in 1:length(list_models_bio17)){

            #select the [m] evolution model 
            selected_model = list_models_bio17[m]

            #extract data of [m] model
            model = get(selected_model)

            #select the row of the corresponding species
            model = model[which(model$species == paste("Pinus_", species, sep="")),]

            #extract all cell values from the raster with climatic data of the [s] scenario only in those areas with uncertainty (raster_subsetted)
            cell_values = getValues(raster_subsetted)

            #extract ID of those cells without NA
            cells_withot_NA = which(!is.na(cell_values))

            #extract, from all cells withput NA, those whose value is inside the phylogenetic range (including the current value but not including the ancestral). For that we used is.between function, created by me. 
            cell_inside_phylo_range = which(is.between(cell_value = na.omit(getValues(raster_subsetted)), ancestral = model$ace, current = model$current_value))

            #from ID of cells without NA, select the ID of those whose vale is inside of the phylo range
            final_cells = cells_withot_NA[cell_inside_phylo_range]

            #create a empty raster with the same extent and resolution than the [s] IPCC raster
            final_raster = raster(extent(stack_bio17_cropped[[s]]), resolution = res(stack_bio17_cropped[[s]]))
            final_raster_proportion = raster(extent(stack_bio17_cropped[[s]]), resolution = res(stack_bio17_cropped[[s]]))

            #fill the raster with zeros
            final_raster[] <- 0
            final_raster_proportion[] <- 0

            #add to these final cells a value of suitability without and with proportion
            final_raster[final_cells] <- 1
            final_raster_proportion[final_cells] <- phylo_proportion(x=raster_subsetted[final_cells], ancestral_value=model$ace, current_value=model$current_value)

            #add the name of the raster
            names(final_raster) <- paste(selected_scenario, "_", strsplit(selected_model, split="_")[[1]][3], "_", strsplit(selected_model, split="_")[[1]][4], sep="")
            names(final_raster_proportion) <- paste(selected_scenario, "_", strsplit(selected_model, split="_")[[1]][3], "_", strsplit(selected_model, split="_")[[1]][4], sep="")        

            #save the raster into a stack
            phylo_rasters_bio17 = stack(phylo_rasters_bio17, final_raster)
            phylo_rasters_bio17_proportion = stack(phylo_rasters_bio17_proportion, final_raster_proportion)
        }    
    }

    #check that all scenarios has been included
    nlayers(phylo_rasters_bio17) == length(climatic_scenarios)
    nlayers(phylo_rasters_bio17_proportion) == length(climatic_scenarios)

    #sum suitability across IPCC scenarios for bio17
    sum_phylo_bio17 = calc(phylo_rasters_bio17, function(x) (sum(x)))

    #loop for comparing projected climate and phylo range
    phylo_rasters_bio4 = stack()
    phylo_rasters_bio4_proportion = stack()
    for(s in 1:length(climatic_scenarios)){

        #select the climatic scenario
        selected_scenario = climatic_scenarios[s]

        #create a empty raster with the same extent and resolution of the raster layer of the [s] IPCC scenario
        raster_subsetted = raster(extent(stack_bio4_cropped[[s]]), resolution = res(stack_bio4_cropped[[s]]))

        #add to the empty raster those cells of the [s] IPCC raster in which the habitat suitability is higher than 25 and lower than 75
        raster_subsetted[which(getValues(ensamble_suitability) > 25 & getValues(ensamble_suitability) < 75)] <- stack_bio4_cropped[[s]][which(getValues(ensamble_suitability) > 25 & getValues(ensamble_suitability) < 75)] #suitability map has the same resolution and extent than IPCC maps because the suitability map was obtained from these raster of ICPP scenarios. 

        #for each evolution model 
        for(m in 1:length(list_models_bio4)){

            #select the [m] evolution model 
            selected_model = list_models_bio4[m]

            #extract data of [m] model
            model = get(selected_model)

            #select the row of the corresponding species
            model = model[which(model$species == paste("Pinus_", species, sep="")),]

            #extract all cell values from the raster with climatic data of the [s] scenario only in those areas with uncertainty (raster_subsetted)
            cell_values = getValues(raster_subsetted)

            #extract ID of those cells without NA
            cells_withot_NA = which(!is.na(cell_values))

            #extract, from all cells withput NA, those whose value is inside the phylogenetic range (including the current value but not including the ancestral). For that we used is.between function, created by me. 
            cell_inside_phylo_range = which(is.between(cell_value = na.omit(getValues(raster_subsetted)), ancestral = model$ace, current = model$current_value))

            #from ID of cells without NA, select the ID of those whose vale is inside of the phylo range
            final_cells = cells_withot_NA[cell_inside_phylo_range]

            #create a empty raster with the same extent and resolution than the [s] IPCC raster
            final_raster = raster(extent(stack_bio4_cropped[[s]]), resolution = res(stack_bio4_cropped[[s]]))
            final_raster_proportion = raster(extent(stack_bio4_cropped[[s]]), resolution = res(stack_bio4_cropped[[s]]))

            #fill the raster with zeros
            final_raster[] <- 0
            final_raster_proportion[] <- 0

            #add to these final cells a value of suitability without and with proportion
            final_raster[final_cells] <- 1
            final_raster_proportion[final_cells] <- phylo_proportion(x=raster_subsetted[final_cells], ancestral_value=model$ace, current_value=model$current_value)

            #add the name of the raster
            names(final_raster) <- paste(selected_scenario, "_", strsplit(selected_model, split="_")[[1]][3], "_", strsplit(selected_model, split="_")[[1]][4], sep="")
            names(final_raster_proportion) <- paste(selected_scenario, "_", strsplit(selected_model, split="_")[[1]][3], "_", strsplit(selected_model, split="_")[[1]][4], sep="")        

            #save the raster into a stack
            phylo_rasters_bio4 = stack(phylo_rasters_bio4, final_raster)
            phylo_rasters_bio4_proportion = stack(phylo_rasters_bio4_proportion, final_raster_proportion)
        }    
    }    

    #check that all scenarios has been included
    nlayers(phylo_rasters_bio4) == length(climatic_scenarios)
    nlayers(phylo_rasters_bio4_proportion) == length(climatic_scenarios)
    
    #sum suitability across IPCC scenarios for bio4
    sum_phylo_bio4 = calc(phylo_rasters_bio4, function(x) (sum(x)))
    
    #calculate intersection between sum both suitaiblity maps (bio17, bio4)
    intersection_ensamble_phylo = sum_phylo_bio17 * sum_phylo_bio4 #this will be used to exclude areas not shared between variables from the final ensamble. 
    
    #bind bio4 and bio17 rasters without proportions
    phylo_rasters = stack(phylo_rasters_bio4, phylo_rasters_bio17)
    nlayers(phylo_rasters) == length(climatic_scenarios)*2
    
    #bind bio4 and bio17 rasters with proportions
    phylo_rasters_proportion = stack(phylo_rasters_bio4_proportion, phylo_rasters_bio17_proportion)
    nlayers(phylo_rasters_proportion) == length(climatic_scenarios)*2
    
    #calculate the proportion of cells that fall inside the differents phylo ranges across IPCC scenarios and global circulation models without and with proportions
    ensamble_phylo = calc(phylo_rasters, function(x) sum(x)/nlayers(phylo_rasters)) 
    ensamble_phylo_proportion = calc(phylo_rasters_proportion, function(x) sum(x)/nlayers(phylo_rasters_proportion))
    
    #check that raster without proportions has a equal or higher suitability than the raster with proportions
    median(getValues(ensamble_phylo)[which(!getValues(ensamble_phylo)==0)]) >= median(getValues(ensamble_phylo_proportion)[which(!getValues(ensamble_phylo_proportion)==0)])

    #select only those areas suitable for at least one scenario for each variables in the suitability without proportions
    #bio17 and bio4 suitaiblity don't overlap in areas with zero in intersection_ensamble_phylo
    ensamble_phylo[intersection_ensamble_phylo == 0] <- 0
    ensamble_phylo_proportion[intersection_ensamble_phylo == 0] <- 0
    
    #save second ensables
    writeRaster(ensamble_phylo, paste("/Users/dsalazar/nicho_pinus/results/phylo_ensamble/without_proportions/", species, "_phylo_ensamble_without_proportions.asc", sep=""), overwrite=TRUE)
    writeRaster(ensamble_phylo_proportion, paste("/Users/dsalazar/nicho_pinus/results/phylo_ensamble/with_proportions/", species, "_phylo_ensamble_with_proportions.asc", sep=""), overwrite=TRUE)
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

#transform to a vector
epithet_species_list = as.vector(epithet_species_list)

# set up cluster
clust <- makeCluster(5)
registerDoParallel(clust)

#paralelize
foreach(species = epithet_species_list, .packages=c("raster", "rgeos")) %dopar% {
    phylo_analyses(species = species)
} 

#stop the cluster 
stopCluster(clust)