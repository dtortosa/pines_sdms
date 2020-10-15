##code for checking our hypothess about the phylo correction
require(raster)

##############################
####### Mediterranean ########
##############################

#extension europe
extension = c(-10,40,30,70)

#background for ployying
background = crop(bio4, extension)
background[which(!is.na(getValues(background)))] <- 0

#species
species = c("brutia", "cembra", "halepensis", "heldreichii", "mugo", "nigra", "peuce", "pinaster", "pinea", "sylvestris")

#loop for plotting
for(s in 1:length(species)){

    #selected_species
    selected_species = paste("Pinus_", species[s], sep="")
    
    #extract phylo range
    range_bio_4 = as.numeric(final_anc_bm_bio4[which(final_anc_bm_bio4$species==selected_species), which(colnames(final_anc_bm_bio4) %in% c("current_value", "ace"))])    
    range_bio_17 = as.numeric(final_anc_bm_bio17[which(final_anc_bm_bio17$species==selected_species), which(colnames(final_anc_bm_bio17) %in% c("current_value", "ace"))])

    #extract extremes of the range    
    max_bio4 = range_bio_4[which(range_bio_4 == max(range_bio_4))]
    min_bio4 = range_bio_4[which(range_bio_4 == min(range_bio_4))]
    max_bio17 = range_bio_17[which(range_bio_17 == max(range_bio_17))]
    min_bio17 = range_bio_17[which(range_bio_17 == min(range_bio_17))]
    
    
    #set climatic scenarios    
    circ_models = c("bc26", "cc26", "gs26", "he26", "ip26", "mg26", "mr26", "bc45", "cc45", "gs45", "he45", "ip45", "mg45", "mr45", "bc60", "cc60", "gs60", "he60", "ip60", "mg60", "mr60", "bc85", "cc85", "gs85", "he85", "ip85", "mg85", "mr85")
    
    #for each scenario plot areas inside the phylo range in bio4    
    pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/explanation_evol_history/mediterranean/plots_bio4_", selected_species, ".pdf", sep=""))
    for(i in 1:length(circ_models)){
        
        #selected scenario
        circ_models_selected = circ_models[i]
        
        #extract the corresponding raster cropped
        bio4_crop = crop(raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/climate_proj_phylo/bio4_", circ_models_selected, ".asc", sep="")),extension)
        
        #remove areas outside the range
        bio4_crop[which(getValues(bio4_crop)<min_bio4 | getValues(bio4_crop)>max_bio4)]<- NA

        #plot
        plot(crop(background, extension), col="gray", legend=FALSE, main=paste(selected_species, circ_models_selected, sep="_"))
        plot(bio4_crop, add=TRUE)
    }
    dev.off()

    #for each scenario plot areas inside the phylo range in bio17
    pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/explanation_evol_history/mediterranean/plots_bio17_", selected_species, ".pdf", sep=""))    
    for(i in 1:length(circ_models)){
        
        #selected scenario
        circ_models_selected = circ_models[i]
        
        #extract the corresponding raster cropped
        bio17_crop = crop(raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/climate_proj_phylo/bio17_", circ_models_selected, ".asc", sep="")),extension)

        #remove areas outside the range
        bio17_crop[which(getValues(bio17_crop)<min_bio17 | getValues(bio17_crop)>max_bio17)]<- NA

        #plot
        plot(crop(background, extension), col="gray", legend=FALSE, main=paste(selected_species, circ_models_selected, sep="_"))
        plot(bio17_crop, add=TRUE)
    }
    dev.off()    
}   


#######################
####### Hudson ########
#######################

#extension europe
extension = c(-100,-50,30,70)

#background for ployying
background = crop(bio4, extension)
background[which(!is.na(getValues(background)))] <- 0

#species
species = c("banksiana", "resinosa", "strobus")

#loop for plotting
for(s in 1:length(species)){

    #selected_species
    selected_species = paste("Pinus_", species[s], sep="")
    
    #extract phylo range
    range_bio_4 = as.numeric(final_anc_bm_bio4[which(final_anc_bm_bio4$species==selected_species), which(colnames(final_anc_bm_bio4) %in% c("current_value", "ace"))])    
    range_bio_17 = as.numeric(final_anc_bm_bio17[which(final_anc_bm_bio17$species==selected_species), which(colnames(final_anc_bm_bio17) %in% c("current_value", "ace"))])

    #extract extremes of the range    
    max_bio4 = range_bio_4[which(range_bio_4 == max(range_bio_4))]
    min_bio4 = range_bio_4[which(range_bio_4 == min(range_bio_4))]
    max_bio17 = range_bio_17[which(range_bio_17 == max(range_bio_17))]
    min_bio17 = range_bio_17[which(range_bio_17 == min(range_bio_17))]
    
    
    #set climatic scenarios    
    circ_models = c("bc26", "cc26", "gs26", "he26", "ip26", "mg26", "mr26", "bc45", "cc45", "gs45", "he45", "ip45", "mg45", "mr45", "bc60", "cc60", "gs60", "he60", "ip60", "mg60", "mr60", "bc85", "cc85", "gs85", "he85", "ip85", "mg85", "mr85")
    
    #for each scenario plot areas inside the phylo range in bio4    
    pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/explanation_evol_history/hudson/plots_bio4_", selected_species, ".pdf", sep=""))
    for(i in 1:length(circ_models)){
        
        #selected scenario
        circ_models_selected = circ_models[i]
        
        #extract the corresponding raster cropped
        bio4_crop = crop(raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/climate_proj_phylo/bio4_", circ_models_selected, ".asc", sep="")),extension)
        
        #remove areas outside the range
        bio4_crop[which(getValues(bio4_crop)<min_bio4 | getValues(bio4_crop)>max_bio4)]<- NA

        #plot
        plot(crop(background, extension), col="gray", legend=FALSE, main=paste(selected_species, circ_models_selected, sep="_"))
        plot(bio4_crop, add=TRUE)
    }
    dev.off()

    #for each scenario plot areas inside the phylo range in bio17
    pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/explanation_evol_history/hudson/plots_bio17_", selected_species, ".pdf", sep=""))    
    for(i in 1:length(circ_models)){
        
        #selected scenario
        circ_models_selected = circ_models[i]
        
        #extract the corresponding raster cropped
        bio17_crop = crop(raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/climate_proj_phylo/bio17_", circ_models_selected, ".asc", sep="")),extension)

        #remove areas outside the range
        bio17_crop[which(getValues(bio17_crop)<min_bio17 | getValues(bio17_crop)>max_bio17)]<- NA

        #plot
        plot(crop(background, extension), col="gray", legend=FALSE, main=paste(selected_species, circ_models_selected, sep="_"))
        plot(bio17_crop, add=TRUE)
    }
    dev.off()    
}
