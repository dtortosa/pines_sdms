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
#check
if(FALSE){
    require(tidyverse)
    paste("Pinus", epithet_species_list, sep=" ") == str_trim(as.vector(list_species[,1])) #the seperated epithet more Pinus are equal to list_species? We use str_trim from tidyverse to remove end spaces in each element of list_species
}#it is in false because loading tidyverse lead to load several packages that have a function name "extract", and this gives problems with the extract function of raster. If you want to check run these lines manually

#remove tecunumanii, jaliscana y discolor. These species are not used for Niche paper. The two first because we should downloadad gbif data now, so we would mix gbif data form 2016 and 2019. The third one was impossible to differentiate from P. cembroides
epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("tecunumanii", "jaliscana", "discolor"))]
#check
c("tecunumanii", "jaliscana", "discolor") %in% epithet_species_list

#load environment variables for using them as a background
clay = raster("/Users/dsalazar/nicho_pinus/data/climate/finals/clay.asc") 
bio1 = raster("/Users/dsalazar/nicho_pinus/data/climate/finals/bio1.asc") 
environment_var = clay*bio1 
environment_var[which(getValues(environment_var) >= min(getValues(environment_var), na.rm = TRUE))] <- 0 #set all continent areas as 0

#load buffer albicaulis to get a reduced resolution version of environment_var to mask the distribution buffers used for the sum of distribution
buffer_albicaulis = raster(paste("/Users/dsalazar/nicho_pinus/data/buffers/albicaulis_distribution_buffer", ".asc", sep=""))

#resample environment_var
environment_var_low_res = resample(environment_var, buffer_albicaulis, method="bilinear")

#open stacks for saving binary raster with current and future suitability
current_suit_stack = stack()
projected_suit_inside_range_stack = stack()
phylo_ensamble_inside_range_intersect_projected_suit_inside_range_stack = stack()
projected_suit_stack = stack()
phylo_ensamble_intersect_projected_suit_stack = stack()
sum_distributions = stack()
raster_range_calc_stack = stack()

#It's key that you remove all areas outside the range_calc_buffer and the water bodies for ALL rasters, because these areas would enter into the calculations. Becasue of this, I have carefully masked and cropped all the predicions (current, future and phylogenetic)

#open data frame to save metrics of suitability change
suitability_changes = data.frame(species=NA, range_loss_no_phylo=NA, range_loss_phylo=NA, range_change_no_phylo=NA, range_change_phylo=NA)

#for each species
for(i in 1:length(epithet_species_list)){
    
    #select the [i] epithet
    species = epithet_species_list[i]

    #print species name
    print(species)

    #load distribution buffer
    ocurrences_buffer = raster(paste("/Users/dsalazar/nicho_pinus/data/buffers/", species, "_distribution_buffer", ".asc", sep=""))

    #drop sea areas inside the ocurrences_buffer
    ocurrences_buffer = mask(ocurrences_buffer, environment_var_low_res, inverse=FALSE)     

    #convert NAs into 0 to avoid problems in the sum
    ocurrences_buffer[which(is.na(getValues(ocurrences_buffer)))] <- 0

    #save it
    sum_distributions = stack(sum_distributions, ocurrences_buffer)

    #load the polygon used for calculations of changes of substitutability (calc_ranges)
    if(!species=="pumila"){
        raster_range_calc = raster(paste("/Users/dsalazar/nicho_pinus/data/buffers_range_calc/", species, "_range_calc_buffer.asc", sep=""))
    } else {
        raster_range_calc = raster(paste("/Users/dsalazar/nicho_pinus/data/buffers_range_calc/", species, "_buffer_range_calc.asc", sep=""))        
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
    raster_range_calc_stack = stack(raster_range_calc_stack, raster_range_calc)

    #load predicted suitability
    current_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_predictions_bin/ensamble_predictions_bin_", species, ".tif", sep=""))    

    #crop current suitability to reduce map size
    current_suit = crop(current_suit, polygon_range_calc)
    
    #mask current suitability to remove all areas outside the buffer calc range
    current_suit = mask(current_suit, polygon_range_calc)

    #mask with clay to remove water bodies
    current_suit = mask(current_suit, environment_var_cropped)

    #load projected suitability
    projected_suit = raster(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_projections_bin/ensamble_projections_bin_", species, ".tif", sep=""))

    #crop current suitability to reduce map size
    projected_suit = crop(projected_suit, polygon_range_calc)
    
    #mask current suitability to remove all areas outside the buffer calc range
    projected_suit = mask(projected_suit, polygon_range_calc)

    #mask with clay to remove water bodies
    projected_suit = mask(projected_suit, environment_var_cropped)

    #load phylo ensamble with proportion
    phylo_ensamble = raster(paste("/Users/dsalazar/nicho_pinus/results/phylo_ensamble/with_proportions/", species, "_phylo_ensamble_with_proportions.asc", sep=""))

    #crop phylo ensamble to reduce map size
    phylo_ensamble = crop(phylo_ensamble, polygon_range_calc)
    
    #mask phylo ensamble to remove all areas outside the buffer calc range
    phylo_ensamble = mask(phylo_ensamble, polygon_range_calc)

    #mask with clay to remove water bodies
    phylo_ensamble = mask(phylo_ensamble, environment_var_cropped) #in the case of pumila, masking with the polygon_range_calc buffer leaves the sea of japan with zero instead of NA. This is caused when removing sea areas from the raster of that buffer, the two extremes of the Japan's sea almost touch and the area inside is included. This is not a problem because after that, we mask with the environmnetal varaible (bio1 and clay), so sea areas are removed. In species with several polygons of disitribution is not a problem because: 1) The calc_range_buffer is ver big, so in almost all cases all polygons are included within it. If sea areas inside of them they will be removedd with environment_var. 

    #obtain maps with zero and ones from suitability maps (1 means suitable). Threshold for pre and projections of suitabiity is equal or higher 75. Threshold for phylosuitability higher than 0.25 (we remove areas very close to the ancestors or very variable across scenairos, being in the phylo range of both climatic variables (only areas inside both were considered, which is strict) is itself an indicative of suitability.)
    #current suitability
    current_suit[which(getValues(current_suit) < 75)] <- 0 #set as zero those areas with suitability lower than 75
    current_suit[which(getValues(current_suit) > 0)] <- 1 #set as 1 all areas with suitability higher than zero (i.e. all with suit equal or higher than 75)

    #future suitability
    projected_suit[which(getValues(projected_suit) < 75)] <- 0 #set as zero those areas with suitability lower than 75
    projected_suit[which(getValues(projected_suit) > 0)] <- 1 #set as 1 all areas with suitability higher than zero (i.e. all with suit equal or higher than 75)

    #current suitability
    phylo_ensamble[which(getValues(phylo_ensamble) < 0.1)] <- 0 #set as zero those areas with phylo suitability lower than 0.25
    phylo_ensamble[which(getValues(phylo_ensamble) > 0)] <- 1 #set as 1 all areas with suitability higher than zero (i.e. all with suit equal or higher than 0.1)

    #extract suitability under future conditions from areas that are currently suitables
    projected_suit_inside_range = projected_suit
    projected_suit_inside_range[which(!getValues(current_suit)==1)] <- 0

    #extract phylo suitability under future conditions from areas that are currently suitables
    phylo_ensamble_inside_range = phylo_ensamble
    phylo_ensamble_inside_range[which(!getValues(current_suit)==1)] <- 0

    #calculate the intersection between phylo_ensamble and projected_suit to get all suitable areas under future conditions according to both, SDMs and phylogenetic correction
    #empty raster for suitability elsewhere
    phylo_ensamble_intersect_projected_suit = raster()
    extent(phylo_ensamble_intersect_projected_suit) = extent(projected_suit)
    res(phylo_ensamble_intersect_projected_suit) = res(projected_suit)
    #empty raster for suitability inside current suitable areas    
    phylo_ensamble_inside_range_intersect_projected_suit_inside_range = raster()
    extent(phylo_ensamble_inside_range_intersect_projected_suit_inside_range) = extent(projected_suit_inside_range)
    res(phylo_ensamble_inside_range_intersect_projected_suit_inside_range) = res(projected_suit_inside_range)

    #set as 1 those cells with 1 in projected_suit and the phylo_ensamble raster
    phylo_ensamble_intersect_projected_suit[c(which(getValues(projected_suit) == 1), which(getValues(phylo_ensamble) == 1))] <- 1
    phylo_ensamble_inside_range_intersect_projected_suit_inside_range[c(which(getValues(projected_suit_inside_range) == 1), which(getValues(phylo_ensamble_inside_range) == 1))] <- 1

    #crop phylo ensamble to reduce map size
    phylo_ensamble_intersect_projected_suit = crop(phylo_ensamble_intersect_projected_suit, polygon_range_calc)
    phylo_ensamble_inside_range_intersect_projected_suit_inside_range = crop(phylo_ensamble_inside_range_intersect_projected_suit_inside_range, polygon_range_calc)

    #mask phylo ensamble to remove all areas outside the buffer calc range
    phylo_ensamble_intersect_projected_suit = mask(phylo_ensamble_intersect_projected_suit, polygon_range_calc)
    phylo_ensamble_inside_range_intersect_projected_suit_inside_range = mask(phylo_ensamble_inside_range_intersect_projected_suit_inside_range, polygon_range_calc)

    #mask with environment_var to remove water bodies
    phylo_ensamble_intersect_projected_suit = mask(phylo_ensamble_intersect_projected_suit, environment_var_cropped)    
    phylo_ensamble_inside_range_intersect_projected_suit_inside_range = mask(phylo_ensamble_inside_range_intersect_projected_suit_inside_range, environment_var_cropped)

    #check that all cells with 1 with and without phylo are included in the phinal raster
    #inside-outside current range
    print(summary(which(getValues(projected_suit) == 1) %in% which(getValues(phylo_ensamble_intersect_projected_suit) == 1)))
    print(summary(which(getValues(phylo_ensamble) == 1) %in% which(getValues(phylo_ensamble_intersect_projected_suit) == 1)))
    print(length(which(getValues(projected_suit) == 1) %in% which(getValues(phylo_ensamble_intersect_projected_suit) == 1)) + length(which(getValues(phylo_ensamble) == 1) %in% which(getValues(phylo_ensamble_intersect_projected_suit) == 1)) == length(na.omit(getValues(phylo_ensamble_intersect_projected_suit))))
    #inside current range
    print(summary(which(getValues(projected_suit_inside_range) == 1) %in% which(getValues(phylo_ensamble_inside_range_intersect_projected_suit_inside_range) == 1)))
    print(summary(which(getValues(phylo_ensamble_inside_range) == 1) %in% which(getValues(phylo_ensamble_inside_range_intersect_projected_suit_inside_range) == 1)))
    print(length(which(getValues(projected_suit_inside_range) == 1) %in% which(getValues(phylo_ensamble_inside_range_intersect_projected_suit_inside_range) == 1)) + length(which(getValues(phylo_ensamble_inside_range) == 1) %in% which(getValues(phylo_ensamble_inside_range_intersect_projected_suit_inside_range) == 1)) == length(na.omit(getValues(phylo_ensamble_inside_range_intersect_projected_suit_inside_range))))

    #extract size of area suitable
    current_suitable_area = length(which(getValues(current_suit)==1))
    future_suitable_area_no_phylo_inside_current_range = length(which(getValues(projected_suit_inside_range)==1))
    future_suitable_area_phylo_inside_current_range = length(which(getValues(phylo_ensamble_inside_range_intersect_projected_suit_inside_range)==1))    
    future_suitable_area_no_phylo_elsewhere = length(which(getValues(projected_suit)==1))
    future_suitable_area_phylo_elsewhere = length(which(getValues(phylo_ensamble_intersect_projected_suit)==1))

    #calculate range loss as (current suitable area - nº cells of that area that remain suitable ) / current suitable area, then multiplied by 100
    range_loss_no_phylo = ((current_suitable_area - future_suitable_area_no_phylo_inside_current_range) / current_suitable_area ) * 100
    range_loss_phylo = ((current_suitable_area - future_suitable_area_phylo_inside_current_range) / current_suitable_area ) * 100    

    #check that phylosuitaiblity increase suitable areas
    print(future_suitable_area_phylo_elsewhere >= future_suitable_area_no_phylo_elsewhere)
    print(future_suitable_area_phylo_inside_current_range >= future_suitable_area_no_phylo_inside_current_range)

    #check that suitability outside current range and under future conditons is at least equal to the suitability inside areas that are currently suitable
    print(future_suitable_area_no_phylo_elsewhere >= future_suitable_area_no_phylo_inside_current_range)
    print(future_suitable_area_phylo_elsewhere >= future_suitable_area_phylo_inside_current_range)

    #calculate range change as (nº cells of that area that are suitable across the whole calc_range_buffer - current suitable areas) / current suitable area, then multiplied by 100. Here we consider future suitability of both areas that are suitable or unsuitable currently
    range_change_no_phylo = ((future_suitable_area_no_phylo_elsewhere - current_suitable_area) / current_suitable_area ) * 100
    range_change_phylo = ((future_suitable_area_phylo_elsewhere - current_suitable_area) / current_suitable_area ) * 100   

    #save metrics of suitability changes
    suitability_changes = rbind.data.frame(suitability_changes, cbind.data.frame(species, range_loss_no_phylo, range_loss_phylo, range_change_no_phylo,range_change_phylo))

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
    current_suit_stack = stack(current_suit_stack, current_suit)
    projected_suit_inside_range_stack = stack(projected_suit_inside_range_stack, projected_suit_inside_range)
    phylo_ensamble_inside_range_intersect_projected_suit_inside_range_stack = stack(phylo_ensamble_inside_range_intersect_projected_suit_inside_range_stack, phylo_ensamble_inside_range_intersect_projected_suit_inside_range)
    projected_suit_stack = stack(projected_suit_stack, projected_suit)
    phylo_ensamble_intersect_projected_suit_stack = stack(phylo_ensamble_intersect_projected_suit_stack, phylo_ensamble_intersect_projected_suit)
}

#remove first row without NAs
suitability_changes = suitability_changes[-1,]

#check all species are included in the table
nrow(suitability_changes) == length(epithet_species_list)

#check all species are included in the stacks
nlayers(current_suit_stack) == 112
nlayers(projected_suit_inside_range_stack) == 112
nlayers(phylo_ensamble_inside_range_intersect_projected_suit_inside_range_stack) == 112
nlayers(projected_suit_stack) == 112
nlayers(phylo_ensamble_intersect_projected_suit_stack) == 112
nlayers(raster_range_calc_stack)
nlayers(sum_distributions)

#sum predictions under current and future conditions
current_suit_stack_sum = calc(current_suit_stack, function(x) (sum(x)))
projected_suit_stack_sum = calc(projected_suit_stack, function(x) (sum(x)))
projected_suit_phylo_stack_sum = calc(phylo_ensamble_intersect_projected_suit_stack, function(x) (sum(x)))
sum_distributions_sum = calc(sum_distributions, function(x) (sum(x)))
raster_range_calc_stack_sum = calc(raster_range_calc_stack, function(x) (sum(x)))

#save the rasters
writeRaster(current_suit_stack_sum, "/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/current_suit_stack_sum.asc", overwrite=TRUE)
writeRaster(projected_suit_stack_sum, "/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/projected_suit_stack_sum.asc", overwrite=TRUE)
writeRaster(projected_suit_phylo_stack_sum, "/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/projected_suit_phylo_stack_sum.asc", overwrite=TRUE)
writeRaster(sum_distributions_sum, "/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/sum_distributions_sum.asc", overwrite=TRUE)
writeRaster(raster_range_calc_stack_sum, "/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/raster_range_calc_stack_sum.asc", overwrite=TRUE)

#save the table
write.table(suitability_changes, "/Users/dsalazar/nicho_pinus/results/final_analyses/synthesis_figure/suitability_changes.csv", sep=",", row.names=FALSE, col.names=TRUE)