#code for making final figures of niche paper

#required packages
require(raster)
require(rgeos)

#########################
##### SET FUNCTIONS #####
#########################
plot_suples = function(species, size="medium"){

    #load species distribution (withput buffer)
    distri_raster = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(species, "01.img", sep="_"), sep="_"))

    #create a polygon from distributon: DISTRIBUTION BUFFER USED INSTEAD OF THIS (see above)
    #distri_polygon = rasterToPolygons(distri_raster, fun=function(x){x==1}, n=16, dissolve=TRUE) #esta funcion de raster te transforma un raster completo o una parte del mismo en un poliogno. En nuestro caso solo queremos las celdas con valor=1, es decir, presencias. Por eso ponemos x==1. 
        #dissolve = TRUE for dissolve limit inside the polygon, only external

    #create a polygon from distributon + buffer
    ocurrences_buffer_raster = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/", species, "_distribution_buffer.asc", sep="")) #load it
    ocurrences_buffer_polygon = rasterToPolygons(ocurrences_buffer_raster, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon

    #load the polygon used for calculations of changes of substitutability (calc_ranges)
    if(!species=="pumila"){
        raster_plot_buffer = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/buffers_calc_ranges/", species, "_range_calc_buffer.asc", sep=""))
    } else {
        raster_plot_buffer = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/buffers_calc_ranges/", species, "_buffer_range_calc.asc", sep=""))        
    }                
    polygon_plot_buffer = rasterToPolygons(raster_plot_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon

    #load current suitability
    current_suit = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ensamble_predictions_bin/ensamble_predictions_bin_", species, ".tif", sep=""))
    
    #crop current suitability to reduce map size
    current_suit = crop(current_suit, polygon_plot_buffer)
    
    #mask current suitability to remove all areas outside the buffer calc range
    current_suit = mask(current_suit, polygon_plot_buffer)

    #load ensamble of binary suitability
    projected_suit = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ensamble_projections_bin/ensamble_projections_bin_", species, ".tif", sep=""))

    #crop current suitability to reduce map size
    projected_suit = crop(projected_suit, polygon_plot_buffer)
    
    #mask current suitability to remove all areas outside the buffer calc range
    projected_suit = mask(projected_suit, polygon_plot_buffer)

    #load ocurrences
    ocurrence_data = read.csv(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences/", species, "_complete.presences.csv", sep=""), header=TRUE)

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
    
    #create a polygon from high precision points
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

    #creamos un raster con todas las áreas que caen en algún ranog filo para algún escenario para el raster sin proporciones. PARA plotear
    final_ensamble_phylo = raster(extent(second_ensamble_phylo), resolution = res(second_ensamble_phylo))
    final_ensamble_phylo[]<-0
    final_ensamble_phylo[which(getValues(second_ensamble_phylo) > 0)] <- 1
    
    #creamos un raster con todas las áreas que caen en algún ranog filo para algún escenario para el raster sin proporciones. PARA plotear. 
    final_ensamble_phylo_proportion = raster(extent(second_ensamble_phylo_proportion), resolution = res(second_ensamble_phylo_proportion))
    final_ensamble_phylo_proportion[]<-0
    final_ensamble_phylo_proportion[which(getValues(second_ensamble_phylo_proportion) > 0)] <- 1

    #convert the proportion from 1 to 0.5. If 1 is 0.5, 0.9 would be X; x=(0.5*0.9)/1=0.45. Then sum 0.5 to have the values form 0.5 to 1 and increase in that way visibility in the plot
    third_ensamble_phylo = ((0.5*second_ensamble_phylo)/1)+0.2
    third_ensamble_phylo_proportion = ((0.5*second_ensamble_phylo_proportion)/1)+0.2
    #remove those areas with phylo predicted suitability equal to zero, as they will be 0.5 after the sum
    summary(which(getValues(second_ensamble_phylo)==0) == which(getValues(third_ensamble_phylo)==0.2))#check that zero cell are now 0.1  
    third_ensamble_phylo[which(getValues(second_ensamble_phylo)==0)] <- 0
    summary(which(getValues(second_ensamble_phylo_proportion)==0) == which(getValues(third_ensamble_phylo_proportion)==0.2))#check that zero cell are now 0.1  
    third_ensamble_phylo_proportion[which(getValues(second_ensamble_phylo_proportion)==0)] <- 0
        #these new values are comparable across species and between methods (proportion and not proportion), as the same modifications have been performed in all cases

    #load BIO1 and clay for complete the backgroun around predictions
    clay = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/clay.asc")
    bio1 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio1.asc")
    #combine them
    env_var = bio1*clay
    env_var = crop(env_var, polygon_plot_buffer)#crop
    #set all values higher than the lowest value (i.e. with data) as 1
    env_var[which(getValues(env_var) >= min(getValues(env_var), na.rm=TRUE))] <- 1

    #create a raster with acuatic bodies
    aquatic_bodies =  raster(extent(current_suit), resolution=res(current_suit))
    aquatic_bodies[which(is.na(getValues(env_var)))] <- 1
    aquatic_bodies[which(!is.na(getValues(env_var)))] <- 0

    #color pallete for water and terrestrial areas
    colfunc_water_terrestrial <- colorRampPalette(c(gray.colors(1, start=0.2), "steelblue2"))

    ###plot final con todos los paneles sin usar proporciones
    #ploteamos seas áreas sonbre la idoneidad de hábitat, pero con un nivel de transparencia dependiente del número casos en los que la celda ha caído dentro del rango filo (incertidumbre; alpha=second_ensamble_phylo). Por tanto, aquellas zonas que han caído dentro del rango para pocos scenarios ó solo para una de las varuables se ven poco (podría verse con solo una variable dentro del rango sin entrease en el rango bajo muuuchos escenarios).
    pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple/without_proportions/", species, "_without_proportions.pdf", sep=""), width=12, height = 12)
    par(oma=c(0,0,2.7,2))
    par(mfcol=c(2,2))

    ###Pannel 1
    plot(aquatic_bodies, col=colfunc_water_terrestrial(2), legend=FALSE)
    plot(current_suit,  main="", add=T, axis.args=list(cex.axis=1.5))
    mtext(text="Predicted-current habitat suitability", side=3,line=2, outer=FALSE, cex=1.4, font = 2)
    #add lines around high precision points
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){
        plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=1)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){       
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=1)
    }
    #plot the legend    
    legend("top", legend=c("High precision points", "Low precision points"), lty=c(1,3), cex=1.1, bty="o", lwd=2, horiz=TRUE, bg="white")
    
    #### add evaluation metrics
    ## load data of evaluation metrics
    median_metrics = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/evaluations/medians_evaluations_v3.csv", sep=",", header=TRUE)

    ##calculate median and sd of AUC and OOB for random forest 
    #kappa
    rf_kappa_median = round(median_metrics[which(median_metrics$species==species),]$rf_kappa_median, 3)
    rf_kappa_sd = round(median_metrics[which(median_metrics$species==species),]$rf_kappa_sd, 3)
    #tss
    rf_tss_median = round(median_metrics[which(median_metrics$species==species),]$rf_tss_median, 3)
    rf_tss_sd = round(median_metrics[which(median_metrics$species==species),]$rf_tss_sd, 3)

    ##calculate median and sd of AUC, kappa and TSS for gam 
    #auc
    gam_auc_median = round(median_metrics[which(median_metrics$species==species),]$gam_auc_median, 3)
    gam_auc_sd = round(median_metrics[which(median_metrics$species==species),]$gam_auc_sd, 3)
    #kappa
    gam_kappa_median = round(median_metrics[which(median_metrics$species==species),]$gam_kappa_median, 3)
    gam_kappa_sd = round(median_metrics[which(median_metrics$species==species),]$gam_kappa_sd, 3)
    #tss
    gam_tss_median = round(median_metrics[which(median_metrics$species==species),]$gam_tss_median, 3)
    gam_tss_sd = round(median_metrics[which(median_metrics$species==species),]$gam_tss_sd, 3)

    ##calculate median and sd of AUC, kappa and TSS for glm 
    #auc
    glm_auc_median = round(median_metrics[which(median_metrics$species==species),]$glm_auc_median, 3)
    glm_auc_sd = round(median_metrics[which(median_metrics$species==species),]$glm_auc_sd, 3)
    #kappa
    glm_kappa_median = round(median_metrics[which(median_metrics$species==species),]$glm_kappa_median, 3)
    glm_kappa_sd = round(median_metrics[which(median_metrics$species==species),]$glm_kappa_sd, 3)
    #tss
    glm_tss_median = round(median_metrics[which(median_metrics$species==species),]$glm_tss_median, 3)
    glm_tss_sd = round(median_metrics[which(median_metrics$species==species),]$glm_tss_sd, 3)


    ##plot the results using legend, because let us to include a background
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('RF: Kappa') == .(rf_kappa_median) %+-% .(rf_kappa_sd) ~ bold('TSS') == .(rf_tss_median) %+-% .(rf_tss_sd)),cex=.7, text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=-2)    
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('GAM: AUC') == .(gam_auc_median) %+-% .(gam_auc_sd) ~ bold('Kappa') == .(gam_kappa_median) %+-% .(gam_kappa_sd) ~ bold('TSS') == .(gam_tss_median) %+-% .(gam_tss_sd)),cex=.7, text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=-1) 
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('GLM: AUC') == .(glm_auc_median) %+-% .(glm_auc_sd) ~ bold('Kappa') == .(glm_kappa_median) %+-% .(glm_kappa_sd) ~ bold('TSS') == .(glm_tss_median) %+-% .(glm_tss_sd)),cex=.7, text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=0) 

    ###Pannel 2
    plot(aquatic_bodies, col=colfunc_water_terrestrial(2), legend=FALSE)    
    plot(current_suit,  main="", axis.args=list(cex.axis=1.5), add=T)
    mtext(text="Predicted-current habitat suitability", side=3,line=2, outer=FALSE, cex=1.4, font = 2)
    #plot low precision points
    if(nrow(low_precision_ocurrences)>0){
        points(low_precision_ocurrences$longitude, low_precision_ocurrences$latitude, cex=0.5, col="white", bg=NA, lwd=0.6, type="p", pch=21)
    }
    #plot high precision points
    if(nrow(high_precision_ocurrences)>0){
        points(high_precision_ocurrences$longitude, high_precision_ocurrences$latitude, cex=0.5, col="black", bg=NA, lwd=0.2, type="p", pch=21)
    } 
    #plot the legend    
    legend("top", legend=c("High precision points", "Low precision points"), fill=c("black", "white"), cex=1.1, bty="o", horiz=TRUE, bg="white")

    ###Pannel 3
    plot(aquatic_bodies, col=colfunc_water_terrestrial(2), legend=FALSE)     
    plot(projected_suit, main="", , axis.args=list(cex.axis=1.5), add=T)
    mtext(text="Predicted habitat suitability in 2070", side=3,line=2, outer=FALSE, cex=1.4, font = 2)
    #add lines around high precision points
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){
        plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=1)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){   
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=1)
    }
    #plot the legend    
    #legend("top", legend=c("High precision pts.", "Low precision pts."), lty=c(1,3), cex=1.455, bty="o", lwd=2, horiz=TRUE, bg="white", x.intersp=0.2) #We reduce the interspace between characters to reduce the space occupied by the legend and then we can increase the size of the legend

    ###Pannel 4
    plot(aquatic_bodies, col=colfunc_water_terrestrial(2), legend=FALSE)     
    plot(projected_suit, main="", axis.args=list(cex.axis=1.5), add=T)
    mtext(text="Predicted habitat suitability in 2070", side=3, line=3, outer=FALSE, cex=1.4, font = 2)
    mtext(text="+", side=3,line=2.2, outer=FALSE, cex=1, font = 2)
    mtext(text="Phylo-corrected suitability", side=3,line=1, outer=FALSE, cex=1.4, font = 2)
    #add lines around high precision points
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){
       plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=1)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){   
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=1)
    }
    #add info phlo
    plot(final_ensamble_phylo, col="#009900", alpha=third_ensamble_phylo, add=TRUE, legend=FALSE) 
        #col="#0000FF": color azul oscuro; ALTERNATIVA el negro ("#0000FF"). 
        #alpha=second_ensamble_phylo: transparencia dependiente de la incertidumbre (número casos en los que la celda ha caído dentro del rango filo)
        #add=TRUE: Add to the previous plot
        #legend=FALSE: Not add legend, in that way the legen of the first plot will not be hidden
    #plot the legend of precision    
    #legend("top", legend=c("High precision pts.", "Low precision pts."), lty=c(1,3), cex=1.455, bty="o", lwd=2, horiz=TRUE, bg="white", x.intersp=0.2) #We reduce the interspace between characters to reduce the space occupied by the legend and then we can increase the size of the legend 
    #add legend of acnestral
    #legend("bottom", legend="Phylo-corrected suitability", fill="#00ffe9", bty="o", cex=1.1, bg="white")
    
    #### main title
    title("") #this is used because mtext only does not work (see http://stackoverflow.com/questions/12895783/r-language-mtext-not-working-with-image-plot-array)
    mtext(bquote(italic('Pinus') ~italic(.(species))), outer = TRUE, cex = 2.5, font=2) #bquote is used to convert to italic the specific epithet (see http://stackoverflow.com/questions/27266398/using-italics-in-the-title-on-an-object-from-a-dataframe)
    dev.off()

    ###plot final con todos los paneles con proporciones
    #ploteamos seas áreas sonbre la idoneidad de hábitat, pero con un nivel de transparencia dependiente del número casos en los que la celda ha caído dentro del rango filo (incertidumbre; alpha=second_ensamble_phylo). Por tanto, aquellas zonas que han caído dentro del rango para pocos scenarios ó solo para una de las varuables se ven poco (podría verse con solo una variable dentro del rango sin entrease en el rango bajo muuuchos escenarios).
    pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple/with_proportions/", species, "_with_proportions.pdf", sep=""), width=12, height = 12)
    par(oma=c(0,0,2.7,2))
    par(mfcol=c(2,2))

    ###Pannel 1
    plot(aquatic_bodies, col=colfunc_water_terrestrial(2), legend=FALSE)     
    plot(current_suit,  main="", axis.args=list(cex.axis=1.5), add=T)
    mtext(text="Predicted-current habitat suitability", side=3,line=2, outer=FALSE, cex=1.4, font = 2)
    #add lines around high precision points
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){
        plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=1)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){       
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=1)
    }
    #plot the legend    
    legend("top", legend=c("High precision points", "Low precision points"), lty=c(1,3), cex=1.1, bty="o", lwd=2, horiz=TRUE, bg="white")
    
    #### add evaluation metrics
    ## load data of evaluation metrics
    median_metrics = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/evaluations/medians_evaluations_v3.csv", sep=",", header=TRUE)

    ##calculate median and sd of AUC and OOB for random forest 
    #kappa
    rf_kappa_median = round(median_metrics[which(median_metrics$species==species),]$rf_kappa_median, 3)
    rf_kappa_sd = round(median_metrics[which(median_metrics$species==species),]$rf_kappa_sd, 3)
    #tss
    rf_tss_median = round(median_metrics[which(median_metrics$species==species),]$rf_tss_median, 3)
    rf_tss_sd = round(median_metrics[which(median_metrics$species==species),]$rf_tss_sd, 3)

    ##calculate median and sd of AUC, kappa and TSS for gam 
    #auc
    gam_auc_median = round(median_metrics[which(median_metrics$species==species),]$gam_auc_median, 3)
    gam_auc_sd = round(median_metrics[which(median_metrics$species==species),]$gam_auc_sd, 3)
    #kappa
    gam_kappa_median = round(median_metrics[which(median_metrics$species==species),]$gam_kappa_median, 3)
    gam_kappa_sd = round(median_metrics[which(median_metrics$species==species),]$gam_kappa_sd, 3)
    #tss
    gam_tss_median = round(median_metrics[which(median_metrics$species==species),]$gam_tss_median, 3)
    gam_tss_sd = round(median_metrics[which(median_metrics$species==species),]$gam_tss_sd, 3)

    ##calculate median and sd of AUC, kappa and TSS for glm 
    #auc
    glm_auc_median = round(median_metrics[which(median_metrics$species==species),]$glm_auc_median, 3)
    glm_auc_sd = round(median_metrics[which(median_metrics$species==species),]$glm_auc_sd, 3)
    #kappa
    glm_kappa_median = round(median_metrics[which(median_metrics$species==species),]$glm_kappa_median, 3)
    glm_kappa_sd = round(median_metrics[which(median_metrics$species==species),]$glm_kappa_sd, 3)
    #tss
    glm_tss_median = round(median_metrics[which(median_metrics$species==species),]$glm_tss_median, 3)
    glm_tss_sd = round(median_metrics[which(median_metrics$species==species),]$glm_tss_sd, 3)


    ##plot the results using legend, because let us to include a background
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('RF: Kappa') == .(rf_kappa_median) %+-% .(rf_kappa_sd) ~ bold('TSS') == .(rf_tss_median) %+-% .(rf_tss_sd)),cex=.7, text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=-2)   
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('GAM: AUC') == .(gam_auc_median) %+-% .(gam_auc_sd) ~ bold('Kappa') == .(gam_kappa_median) %+-% .(gam_kappa_sd) ~ bold('TSS') == .(gam_tss_median) %+-% .(gam_tss_sd)),cex=.7, text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=-1) 
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('GLM: AUC') == .(glm_auc_median) %+-% .(glm_auc_sd) ~ bold('Kappa') == .(glm_kappa_median) %+-% .(glm_kappa_sd) ~ bold('TSS') == .(glm_tss_median) %+-% .(glm_tss_sd)),cex=.7, text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=0) 

    ###Pannel 2
    plot(aquatic_bodies, col=colfunc_water_terrestrial(2), legend=FALSE)     
    plot(current_suit,  main="", axis.args=list(cex.axis=1.5), add=T)
    mtext(text="Predicted-current habitat suitability", side=3,line=2, outer=FALSE, cex=1.4, font = 2)
    #plot low precision points
    if(nrow(low_precision_ocurrences)>0){
        points(low_precision_ocurrences$longitude, low_precision_ocurrences$latitude, cex=0.5, col="white", bg=NA, lwd=0.6, type="p", pch=21)
    }
    #plot high precision points
    if(nrow(high_precision_ocurrences)>0){
        points(high_precision_ocurrences$longitude, high_precision_ocurrences$latitude, cex=0.5, col="black", bg=NA, lwd=0.2, type="p", pch=21)
    } 
    #plot the legend    
    legend("top", legend=c("High precision points", "Low precision points"), fill=c("black", "white"), cex=1.1, bty="o", horiz=TRUE, bg="white")

    ###Pannel 3
    plot(aquatic_bodies, col=colfunc_water_terrestrial(2), legend=FALSE)     
    plot(projected_suit, main="", , axis.args=list(cex.axis=1.5), add=T)
    mtext(text="Predicted habitat suitability in 2070", side=3,line=2, outer=FALSE, cex=1.4, font = 2)
    #add lines around high precision points
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){
        plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=1)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){   
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=1)
    }
    #plot the legend    
    #legend("top", legend=c("High precision pts.", "Low precision pts."), lty=c(1,3), cex=1.455, bty="o", lwd=2, horiz=TRUE, bg="white", x.intersp=0.2) #We reduce the interspace between characters to reduce the space occupied by the legend and then we can increase the size of the legend

    ###Pannel 4
    plot(aquatic_bodies, col=colfunc_water_terrestrial(2), legend=FALSE)     
    plot(projected_suit, main="", axis.args=list(cex.axis=1.5), add=T)
    mtext(text="Predicted habitat suitability in 2070", side=3, line=3, outer=FALSE, cex=1.4, font = 2)
    mtext(text="+", side=3,line=2.2, outer=FALSE, cex=1, font = 2)
    mtext(text="Phylo-corrected suitability", side=3,line=1, outer=FALSE, cex=1.4, font = 2)
    #add lines around high precision points
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){
       plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=1)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){   
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=1)
    }
    #add info phlo
    plot(final_ensamble_phylo_proportion, col="#009900", alpha=third_ensamble_phylo_proportion, add=TRUE, legend=FALSE) 
        #col="#0000FF": color azul oscuro; ALTERNATIVA el negro ("#0000FF"). 
        #alpha=second_ensamble_phylo_proportion: transparencia dependiente de la incertidumbre (número casos en los que la celda ha caído dentro del rango filo)
        #add=TRUE: Add to the previous plot
        #legend=FALSE: Not add legend, in that way the legen of the first plot will not be hidden
    #plot the legend of precision    
    #legend("top", legend=c("High precision pts.", "Low precision pts."), lty=c(1,3), cex=1.455, bty="o", lwd=2, horiz=TRUE, bg="white", x.intersp=0.2) #We reduce the interspace between characters to reduce the space occupied by the legend and then we can increase the size of the legend 
    #add legend of acnestral
    #legend("bottom", legend="Phylo-corrected suitability", fill="#00ffe9", bty="o", cex=1.1, bg="white")
    
    #### main title
    title("") #this is used because mtext only does not work (see http://stackoverflow.com/questions/12895783/r-language-mtext-not-working-with-image-plot-array)
    mtext(bquote(italic('Pinus') ~italic(.(species))), outer = TRUE, cex = 2.5, font=2) #bquote is used to convert to italic the specific epithet (see http://stackoverflow.com/questions/27266398/using-italics-in-the-title-on-an-object-from-a-dataframe)
    dev.off()
}

##################
##### SUPLES #####
##################

#list species
list_species = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/species.txt", sep="\t", header=T)
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

#for each path, extract the specie names and plot the suple
for(i in 1:length(epithet_species_list)){

    #select the [i] speceis
    selected_name = epithet_species_list[i]

    #plot the suple
    plot_suples(selected_name, size="medium")
}

#bind all pdf generated into one single pdf file with pdftk
system("rm /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple/full/full_with_proportions.pdf;

    cd /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple/with_proportions; 

    pdftk *.pdf cat output /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple/full/full_with_proportions.pdf") #first delete the prevous full file in full directory; next change the drectory to the that with separate plots; finally run pdftk to bind all pdf in that directory and save the result in full directory

system("rm /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple/full/full_without_proportions.pdf;

    cd /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple/without_proportions; 

    pdftk *.pdf cat output /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple/full/full_without_proportions.pdf") #first delete the prevous full file in full directory; next change the drectory to the that with separate plots; finally run pdftk to bind all pdf in that directory and save the result in full directory