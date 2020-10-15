#required packages
require(raster)
require(rgeos)

#function
compare_ancient_new_results = function(species, proportions, legend=TRUE){

    #load species distribution (withput buffer)
    distri_raster = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(species, "01.img", sep="_"), sep="_"))

    #create a polygon from distributon: DISTRIBUTION BUFFER USED INSTEAD OF THIS (see above)
    #distri_polygon = rasterToPolygons(distri_raster, fun=function(x){x==1}, n=16, dissolve=TRUE) #esta funcion de raster te transforma un raster completo o una parte del mismo en un poliogno. En nuestro caso solo queremos las celdas con valor=1, es decir, presencias. Por eso ponemos x==1. 
        #dissolve = TRUE for dissolve limit inside the polygon, only external

    #create a polygon from distributon + buffer
    ocurrences_buffer_path = unzip("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences.zip",  list=FALSE, exdir = "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/phlyo/buffers", files = paste("ocurrences/", species, "_distribution_buffer.asc", sep=""), junkpaths = TRUE) #unzip species distribution with buffer
    ocurrences_buffer_raster = raster(ocurrences_buffer_path) #load it
    ocurrences_buffer_polygon = rasterToPolygons(ocurrences_buffer_raster, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon

    #create a polygon buffer around the distribution buffer to reduce the area in plots
    polygon_plot_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=15, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) 
    
    #create a raster from the plot buffer polygon 
    raster_plot_buffer = raster() 
    extent(raster_plot_buffer) = extent(distri_raster) 
    res(raster_plot_buffer) = res(distri_raster) 
    raster_plot_buffer  = rasterize(polygon_plot_buffer, raster_plot_buffer) 
    
    #drop the sea areas
    raster_plot_buffer = distri_raster*raster_plot_buffer 
    raster_plot_buffer[!is.na(raster_plot_buffer)] <- 1 
    
    #write the plot buffer without water bodies.  
    writeRaster(raster_plot_buffer, paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/plot_buffers", paste("final_figures", species, "plot_buffer.asc", sep="_"), sep="/"), format="ascii", overwrite=TRUE)
    
    #load ensamble of binary suitability
    projected_suit = raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensamble_projections_bin/ensamble_projections_bin_", species, ".tif", sep=""))

    #load ensamble of binary suitability ancient
    projected_suit_ancient = raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensamble_projections_bin_ancient/ensamble_projections_bin_", species, ".tif", sep=""))

    #crop projected suitability
    projected_suit = crop(projected_suit, polygon_plot_buffer)

    #crop projected suitability ancient
    projected_suit_ancient = crop(projected_suit_ancient, polygon_plot_buffer)

    #create a raster with acuatic bodies
    aquatic_bodies =  raster(extent(projected_suit), resolution=res(projected_suit))
    aquatic_bodies[which(is.na(getValues(projected_suit)))] <- 1

    #load ocurrences
    ocurrence_data = read.csv(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ocurrences/", species, "_complete.presences.csv", sep=""), header=TRUE)
    
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
    
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){    
        #create a spatial point data frame
        coors_high_precision_points = high_precision_ocurrences[,which(colnames(high_precision_ocurrences) %in% c("longitude", "latitude"))]
        coordinates(coors_high_precision_points) = c("longitude", "latitude")
        
        #extrac id cells from that spatial point data frame
        id.cell <- extract(ocurrences_buffer_raster, coors_high_precision_points, cellnumbers=TRUE)[,1]
        
        #create a polygon with those id cell, which will include all high precision points
        raster_high_precision_points = raster(extent(ocurrences_buffer_raster), resolution = res(ocurrences_buffer_raster)) #empty raster with the same extent and res than ocurrences_buffer_raster
        raster_high_precision_points[id.cell] <- ocurrences_buffer_raster[id.cell] #fill the raster with the values from the cells ocurrences_buffer_raster with presence of the species
        polygon_high_precision_points = rasterToPolygons(raster_high_precision_points, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to polygon
    }

    #load phylo ensamble with and without proportions
    second_ensamble_phylo = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/phylo_ensamble/without_proportions/", species, "_phylo_ensamble_without_proportions.asc", sep=""))
    second_ensamble_phylo_proportion = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/phylo_ensamble/with_proportions/", species, "_phylo_ensamble_with_proportions.asc", sep=""))

    #crop with polygon_plot_buffer
    second_ensamble_phylo = crop(second_ensamble_phylo, polygon_plot_buffer)
    second_ensamble_phylo_proportion = crop(second_ensamble_phylo_proportion, polygon_plot_buffer)

    #creamos un raster con todas las áreas que caen en algún ranog filo para algún escenario para el raster sin proporciones. PARA plotear
    final_ensamble_phylo = raster(extent(second_ensamble_phylo), resolution = res(second_ensamble_phylo))
    final_ensamble_phylo[]<-0
    final_ensamble_phylo[which(getValues(second_ensamble_phylo) > 0)] <- 1
    
    #creamos un raster con todas las áreas que caen en algún ranog filo para algún escenario para el raster sin proporciones. PARA plotear. 
    final_ensamble_phylo_proportion = raster(extent(second_ensamble_phylo_proportion), resolution = res(second_ensamble_phylo_proportion))
    final_ensamble_phylo_proportion[]<-0
    final_ensamble_phylo_proportion[which(getValues(second_ensamble_phylo_proportion) > 0)] <- 1

    #load phylo ensamble with and without proportions ancient
    second_ensamble_phylo_ancient = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures_ancient/phylo_ensamble/without_proportions/", species, "_phylo_ensamble_without_proportions.asc", sep=""))
    second_ensamble_phylo_proportion_ancient = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures_ancient/phylo_ensamble/with_proportions/", species, "_phylo_ensamble_with_proportions.asc", sep=""))

    #crop with polygon_plot_buffer
    second_ensamble_phylo_ancient = crop(second_ensamble_phylo_ancient, polygon_plot_buffer)
    second_ensamble_phylo_proportion_ancient = crop(second_ensamble_phylo_proportion_ancient, polygon_plot_buffer)

    #creamos un raster con todas las áreas que caen en algún ranog filo para algún escenario para el raster sin proporciones. PARA plotear
    final_ensamble_phylo_ancient = raster(extent(second_ensamble_phylo_ancient), resolution = res(second_ensamble_phylo))
    final_ensamble_phylo_ancient[]<-0
    final_ensamble_phylo_ancient[which(getValues(second_ensamble_phylo_ancient) > 0)] <- 1
    
    #creamos un raster con todas las áreas que caen en algún ranog filo para algún escenario para el raster sin proporciones. PARA plotear. 
    final_ensamble_phylo_proportion_ancient = raster(extent(second_ensamble_phylo_proportion_ancient), resolution = res(second_ensamble_phylo_proportion_ancient))
    final_ensamble_phylo_proportion_ancient[]<-0
    final_ensamble_phylo_proportion_ancient[which(getValues(second_ensamble_phylo_proportion_ancient) > 0)] <- 1

    ###PLOTS
    par(oma=c(0,0,2.7,2))
    par(mfcol=c(1,2))

    #plot new results
    plot(projected_suit, main="New results", axis.args=list(cex.axis=1.5))
    #mtext(text=bquote(italic('Pinus') ~italic(.(species))), side=3, line=2, outer=FALSE, cex=1.4, font = 2)
    #add acuatic bodies
    plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
    #add lines around high precision points
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){
        plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=2)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){    
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=2)
    }
    #add info phlo
    if(!proportions == TRUE){ #if we don't want proportions
        plot(final_ensamble_phylo, col="#00ffe9", alpha=second_ensamble_phylo, add=TRUE, legend=FALSE)
            #col="#0000FF": color azul oscuro; ALTERNATIVA el negro ("#0000FF"). 
            #alpha=second_ensamble_phylo: transparencia dependiente de la incertidumbre (número casos en los que la celda ha caído dentro del rango filo)
            #add=TRUE: Add to the previous plot
            #legend=FALSE: Not add legend, in that way the legen of the first plot will not be hidden
    } else { #if yes
        plot(final_ensamble_phylo_proportion, col="#00ffe9", alpha=second_ensamble_phylo_proportion, add=TRUE, legend=FALSE) 
        #col="#0000FF": color azul oscuro; ALTERNATIVA el negro ("#0000FF"). 
        #alpha=second_ensamble_phylo: transparencia dependiente de la incertidumbre (número casos en los que la celda ha caído dentro del rango filo)
        #add=TRUE: Add to the previous plot
        #legend=FALSE: Not add legend, in that way the legen of the first plot will not be hidden
    }  

    #plot the legends if   
    if(legend==TRUE){
        legend("top", legend=c("High precision points", "Low precision points"), lty=c(1,3), cex=0.9, bty="o", lwd=2, horiz=TRUE, bg="white")    
        legend("bottom", legend="Phylo-corrected suitability", fill="#00ffe9", bty="o", cex=1.1, bg="white")
    }    

    #plot ancient results
    plot(projected_suit_ancient, main="Ancient results", axis.args=list(cex.axis=1.5))
    #mtext(text=bquote(italic('Pinus') ~italic(.(species))), side=3, line=2, outer=FALSE, cex=1.4, font = 2)
    #add acuatic bodies
    plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
    #add lines around high precision points
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){
        plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=2)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){    
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=2)
    }
    #add info phlo
    if(!proportions == TRUE){ #if we don't want proportions
        plot(final_ensamble_phylo_ancient, col="#00ffe9", alpha=second_ensamble_phylo_ancient, add=TRUE, legend=FALSE)
            #col="#0000FF": color azul oscuro; ALTERNATIVA el negro ("#0000FF"). 
            #alpha=second_ensamble_phylo: transparencia dependiente de la incertidumbre (número casos en los que la celda ha caído dentro del rango filo)
            #add=TRUE: Add to the previous plot
            #legend=FALSE: Not add legend, in that way the legen of the first plot will not be hidden
    } else { #if yes
        plot(final_ensamble_phylo_proportion_ancient, col="#00ffe9", alpha=second_ensamble_phylo_proportion_ancient, add=TRUE, legend=FALSE) 
        #col="#0000FF": color azul oscuro; ALTERNATIVA el negro ("#0000FF"). 
        #alpha=second_ensamble_phylo: transparencia dependiente de la incertidumbre (número casos en los que la celda ha caído dentro del rango filo)
        #add=TRUE: Add to the previous plot
        #legend=FALSE: Not add legend, in that way the legen of the first plot will not be hidden
    } 

    #title
    title("") #this is used because mtext only does not work (see http://stackoverflow.com/questions/12895783/r-language-mtext-not-working-with-image-plot-array)
    mtext(bquote(italic('Pinus') ~italic(.(species))), outer = TRUE, cex = 2.5, font=2) #bquote is used to convert to italic the specific epithet (see http://stackoverflow.com/questions/27266398/using-italics-in-the-title-on-an-object-from-a-dataframe) 
}

#extract species for which there is phylo ensamble
path_species_phylo = list.files("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/phylo_ensamble/without_proportions", pattern=".asc")

#select their names
species_phylo = NULL
for(i in 1:length(path_species_phylo)){
    species_phylo = append(species_phylo, strsplit(path_species_phylo[i], split="_")[[1]][1])
}


#plot all of them
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/comparison_ancient_new_results/comparison_ancient_new_results.pdf", width=12, height = 6)
for(i in 1:length(species_phylo)){
    compare_ancient_new_results(species_phylo[i], proportions=FALSE)
}
dev.off()


