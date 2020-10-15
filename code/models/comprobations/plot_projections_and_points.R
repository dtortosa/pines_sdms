##Code for plot all projections with gbif and critfield distribution points in sewall

#function for plot current distribution and current/future predictions of the models 
plot_sdm = function(species){
    
    #libraries
    #library(raster)

    #load the PA_buffer to crop the naturla distribution
    PA_buffer = raster(paste("/home/dsalazar/data/pa_buffers", paste(species, "PA_buffer.asc", sep="_"), sep="/"))

    #load natural distribution
    nat_distribution =raster(paste("/home/dsalazar/data/maps/p", paste(species, "01.img", sep="_") ,sep="_")) 
    nat_distribution = crop(nat_distribution, PA_buffer) #crop with PA_buffer
    
    #load predicted habitat suitability currently
    current_prediction = raster(paste("/home/dsalazar/modelos/ensamble_predictions_bin/ensamble_predictions_bin", paste(species, "tif", sep=".") ,sep="_"))
    
    #load predicted habitat suitablity in 2070
    future_prediction = raster(paste("/home/dsalazar/modelos/ensamble_projections_bin/ensamble_projections_bin", paste(species, "tif", sep=".") ,sep="_"))

    #load ocurrences
    ocurrence_data = read.csv(paste("/home/dsalazar/data/ocurrences/", species, "_complete.presences.csv", sep=""), header=TRUE)

    #subset high precision presences
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){
        high_precision_ocurrences = ocurrence_data[ocurrence_data$presence==1 & ocurrence_data$precision_weight==1,]
    } else {
        high_precision_ocurrences = data.frame()
    }

    #subset low precision presences
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){
        low_precision_ocurrences = ocurrence_data[ocurrence_data$presence==1 & ocurrence_data$precision_weight==0.5,]
    } else {
        low_precision_ocurrences = data.frame()
    }

    #plot all
    pdf(paste("/home/dsalazar/modelos/final_plots_sdm/plot_sdm_", species, ".pdf", sep=""), width=12, height = 6)
    par(oma=c(0,0,2.7,2))
    par(mfrow=c(1,2))
    #plot current predicted habitat suitability
    plot(current_prediction, main="Predicted habitat suitability currently")
    #plot low precision points
    if(nrow(low_precision_ocurrences)>0){
        points(low_precision_ocurrences$longitude, low_precision_ocurrences$latitude, cex=0.3, col="blue", lwd=0.5, type="p", pch=".")
    }
    #plot high precision points
    if(nrow(high_precision_ocurrences)>0){
        points(high_precision_ocurrences$longitude, high_precision_ocurrences$latitude, cex=0.3, col="red", lwd=0.1, type="p", pch=".")
    } 
    #plot the legend    
    legend("topright", legend=c("High precision points", "Low precision points"), fill=c("red", "blue"), cex=0.7)
    #plot the future predictions of habitat suitability
    plot(future_prediction, main="Predicted habitat suitability in 2070")
    #plot the main title of the plot
    title("") #this is used because mtext only does not work (see http://stackoverflow.com/questions/12895783/r-language-mtext-not-working-with-image-plot-array)
    mtext(bquote(italic('Pinus') ~italic(.(species))), outer = TRUE, cex = 2) #bquote is used to convert to italic the specific epithet (see http://stackoverflow.com/questions/27266398/using-italics-in-the-title-on-an-object-from-a-dataframe)
    dev.off()
}

########Paralelize the process######
require(foreach)
require(doParallel) #for parallel

#create a vector with species names
list_species = read.table("/home/dsalazar/data/species.txt", sep="", header=T)
species = as.vector(list_species$specific_epithet)

#load problematic species
list_problematic_species = read.csv("/home/dsalazar/data/problematic_species.csv", header=TRUE)
problematic_species = as.vector(list_problematic_species$specific_epithet)

#select non_problematic species from the pool of species
non_problematic_species = species[!species %in% problematic_species]

# set up cluster
clust <- makeCluster(6) 
registerDoParallel(clust)

###########################
#########PLOTS#############
###########################
#project to the future for non_problematic species
foreach(i = non_problematic_species, .packages="raster") %dopar% { 
    plot_sdm(species = i)
}

#project to the future for problematic species
foreach(i = problematic_species, .packages="raster") %dopar% { 
    plot_sdm(species = i)
}  

#stop the cluster 
stopCluster(clust)