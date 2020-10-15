######## COMPARISON of models with and without soil variables. This is the only difference, we only remove soil variable and did not add any climate variable. 

#required packages
require(raster)

## select species for which we have current ensamble withput soil
# Not all current suit without soil analyses are done
# path of ensambles
species_current_suit_n_soil_path = list.files("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_without_soil/ensamble_predictions_bin_without_soil")
# names of species
species_current_suit_n_soil = NULL
for(i in 1:length(species_current_suit_n_soil_path)){
    species_current_suit_n_soil = append(species_current_suit_n_soil, strsplit(strsplit(species_current_suit_n_soil_path[i], split="_")[[1]][4], split=".tif")[[1]])
}

## plot current suitability with and without soil
# open pdf
pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/soil_comprobations/plots/current_soil.pdf", sep=""), width=12, height = 6)
par(oma=c(0,0,2.7,2))

#for each of the species with current suit without soil
for(i in 1:length(species_current_suit_n_soil)){

    #select the [i] species
    species = species_current_suit_n_soil[i]

    #load current suitability without soil variables 
    current_suit_without_soil = raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_without_soil/ensamble_predictions_bin_without_soil/ensamble_predictions_bin_", species, ".tif", sep=""))

    #load current suitability with soil variables 
    current_suit_with_soil = raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_final/ensamble_predictions_bin/ensamble_predictions_bin_", species, ".tif", sep=""))

    #two pannels 
    par(mfcol=c(1,2))
    plot(current_suit_without_soil, main="Without soil")
    plot(current_suit_with_soil, main="With soil")

    #### main title
    title("") #this is used because mtext only does not work (see http://stackoverflow.com/questions/12895783/r-language-mtext-not-working-with-image-plot-array)
    mtext(bquote(italic('Pinus') ~italic(.(species))), outer = TRUE, cex = 2.5, font=2, padj=0.5) #bquote is used to convert to italic the specific epithet (see http://stackoverflow.com/questions/27266398/using-italics-in-the-title-on-an-object-from-a-dataframe)

}
dev.off()

## select species for which we have future ensamble withput soil
# Not all future suit without soil analyses are done
# path of ensambles
species_future_suit_n_soil_path = list.files("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_without_soil/ensamble_projections_bin_without_soil")
# names of species
species_future_suit_n_soil = NULL
for(i in 1:length(species_future_suit_n_soil_path)){
    species_future_suit_n_soil = append(species_future_suit_n_soil, strsplit(strsplit(species_future_suit_n_soil_path[i], split="_")[[1]][4], split=".tif")[[1]])
}

## plot future suitability with and without soil
# open pdf
pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/soil_comprobations/plots/future_soil.pdf", sep=""), width=12, height = 6)
par(oma=c(0,0,2.7,2))

#for each of the species with future suit without soil
for(i in 1:length(species_future_suit_n_soil)){

    #select the [i] species
    species = species_future_suit_n_soil[i]

    #load future suitability without soil variables 
    future_suit_without_soil = raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_without_soil/ensamble_projections_bin_without_soil/ensamble_projections_bin_", species, ".tif", sep=""))

    #load future suitability with soil variables 
    future_suit_with_soil = raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_final/ensamble_projections_bin/ensamble_projections_bin_", species, ".tif", sep=""))

    #two pannels 
    par(mfcol=c(1,2))
    plot(future_suit_without_soil, main="Without soil")
    plot(future_suit_with_soil, main="With soil")

    #### main title
    title("") #this is used because mtext only does not work (see http://stackoverflow.com/questions/12895783/r-language-mtext-not-working-with-image-plot-array)
    mtext(bquote(italic('Pinus') ~italic(.(species))), outer = TRUE, cex = 2.5, font=2, padj=0.5) #bquote is used to convert to italic the specific epithet (see http://stackoverflow.com/questions/27266398/using-italics-in-the-title-on-an-object-from-a-dataframe)

}
dev.off()