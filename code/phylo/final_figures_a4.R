#code for making final figures of niche paper

#required packages
require(raster)
require(rgeos)

#########################
##### SET FUNCTIONS #####
#########################
plot_suples = function(species, size="big"){

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
    if(size=="small"){ #for small species narrower buffer
        polygon_plot_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=2.5, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)
    }
    if(size=="medium") {
        polygon_plot_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=10, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) 
    }
    if(size=="big") {
        polygon_plot_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=15, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) 
    }    

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

    #load current suitability
    current_suit = raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_final/ensamble_predictions_bin/ensamble_predictions_bin_", species, ".tif", sep=""))
    
    #crop current suitability
    current_suit = crop(current_suit, polygon_plot_buffer)
    
    #create a raster with acuatic bodies
    aquatic_bodies =  raster(extent(current_suit), resolution=res(current_suit))
    aquatic_bodies[which(is.na(getValues(current_suit)))] <- 1
    
    #load ensamble of binary suitability
    projected_suit = raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_final/ensamble_projections_bin/ensamble_projections_bin_", species, ".tif", sep=""))

    #crop projected suitability
    projected_suit = crop(projected_suit, polygon_plot_buffer)

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

    #set range of PDF pages
    range_pages = seq(252, 251+length(epithet_species_list), 1)
    length(range_pages)

    ###plot final con todos los paneles sin usar proporciones
    #ploteamos seas áreas sonbre la idoneidad de hábitat, pero con un nivel de transparencia dependiente del número casos en los que la celda ha caído dentro del rango filo (incertidumbre; alpha=second_ensamble_phylo). Por tanto, aquellas zonas que han caído dentro del rango para pocos scenarios ó solo para una de las varuables se ven poco (podría verse con solo una variable dentro del rango sin entrease en el rango bajo muuuchos escenarios).
    pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/without_proportions/", species, "_without_proportions.pdf", sep=""), width=8.27, height=11.69)
    par(oma=c(0,1.5,6,2), mai = c(0.45, 0.45, 0.45, 0.45))
    par(mfcol=c(3,2))
    
    ###Pannel 1
    plot(current_suit,  main="", axis.args=list(cex.axis=1.5))
    mtext(text="Predicted-current habitat suitability", side=3,line=1.7, outer=FALSE, cex=1, font = 2, adj=0.45)
    #add acuatic bodies
    plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
    #add lines around high precision pts.
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){
        plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=1)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){    
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=1)
    }
    #plot the legend    
    legend("top", legend=c("High precision pts.", "Low precision pts."), lty=c(1,3), cex=1, bty="o", lwd=2, horiz=TRUE, bg="white")
    
    #### add evaluation metrics
    ## load data of evaluation metrics
    median_metrics = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/evaluations/medians_evaluations.csv", sep=",", header=TRUE)

    ##calculate median and sd of AUC and OOB for random forest 
    #auc
    rf_auc_median = round(median_metrics[which(median_metrics$species==species),]$rf_auc_median, 3)
    rf_auc_sd = round(median_metrics[which(median_metrics$species==species),]$rf_auc_sd, 3)
    #oob
    rf_oob_median = round(median_metrics[which(median_metrics$species==species),]$rf_oob_median, 3)
    rf_oob_sd = round(median_metrics[which(median_metrics$species==species),]$rf_oob_sd, 3)

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
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('RF: AUC') == .(rf_auc_median) %+-% .(rf_auc_sd) ~ bold('OOB') == .(rf_oob_median) %+-% .(rf_oob_sd) ~ '%'),cex=.6,text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=-2)  
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('GAM: AUC') == .(gam_auc_median) %+-% .(gam_auc_sd) ~ bold('Kappa') == .(gam_kappa_median) %+-% .(gam_kappa_sd) ~ bold('TSS') == .(gam_tss_median) %+-% .(gam_tss_sd)),cex=.6, text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=-1) 
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('GLM: AUC') == .(glm_auc_median) %+-% .(glm_auc_sd) ~ bold('Kappa') == .(glm_kappa_median) %+-% .(glm_kappa_sd) ~ bold('TSS') == .(glm_tss_median) %+-% .(glm_tss_sd)),cex=.6, text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=0) 
 
    ###Pannel 2
    plot(current_suit,  main="", axis.args=list(cex.axis=1.5))
    mtext(text="Predicted-current habitat suitability", side=3,line=1.7, outer=FALSE, cex=1, font = 2, adj=0.3)
    #add acuatic bodies
    plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
    #plot low precision pts.
    if(nrow(low_precision_ocurrences)>0){
        points(low_precision_ocurrences$longitude, low_precision_ocurrences$latitude, cex=0.5, col="white", bg=NA, lwd=0.6, type="p", pch=21)
    }
    #plot high precision pts.
    if(nrow(high_precision_ocurrences)>0){
        points(high_precision_ocurrences$longitude, high_precision_ocurrences$latitude, cex=0.5, col="black", bg=NA, lwd=0.2, type="p", pch=21)
    } 
    #plot the legend    
    legend("top", legend=c("High precision pts.", "Low precision pts."), fill=c("black", "white"), cex=1.1, bty="o", horiz=TRUE, bg="white", x.intersp=0.2) #We reduce the interspace between characters to reduce the space occupied by the legend and then we can increase the size of the legend


    ###Pannel 3
    plot(projected_suit, main="", , axis.args=list(cex.axis=1.5))
    mtext(text="Predicted habitat suitability in 2070", side=3,line=1.7, outer=FALSE, cex=1, font = 2)
    #add acuatic bodies
    plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
    #add lines around high precision pts.
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){       
        plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=1)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){       
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=1)
    }
    #plot the legend    
    #legend("top", legend=c("High precision pts.", "Low precision pts."), lty=c(1,3), cex=1.1, bty="o", lwd=2, horiz=TRUE, bg="white")
    
    ###Pannel 4
    plot(projected_suit, main="", axis.args=list(cex.axis=1.5))
    mtext(text="Predicted habitat suitability in 2070", side=3, line=2.6, outer=FALSE, cex=1, font = 2, adj=0.3)
    mtext(text="+", side=3,line=1.7, outer=FALSE, cex=1, font = 2)
    mtext(text="Phylo-corrected suitability", side=3,line=0.8, outer=FALSE, cex=1, font = 2, adj=0.3)
    #add acuatic bodies
    plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
    #add lines around high precision pts.
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){
        plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=1)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){   
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=1)
    }
    #add info phlo
    plot(final_ensamble_phylo, col="#009900", alpha=second_ensamble_phylo, add=TRUE, legend=FALSE) 
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
    number_figure = toString(i)#extract number figure as string. In that way the number can be printed in bold
    mtext(bquote(bold("Figure") ~ bold(.(number_figure)) ~ bold("Appendix S11:") ~ "Projections of habitat suitability under current and future"), outer = TRUE, cex = 1.2, font=2, line=-56.5, adj=0.3, family="Times") #bquote is used to convert to italic the specific epithet (see http://stackoverflow.com/questions/27266398/using-italics-in-the-title-on-an-object-from-a-dataframe)
        #* let us to add strings without spaces between them
    mtext(bquote("conditions for" ~ italic('P.') ~italic(.(species))*"."), outer = TRUE, cex = 1.2, font=2, line=-58, adj=0.0725, family="Times")

    #number of PDFs
    mtext(side = 1, text = range_pages[i], outer = TRUE, adj=0.88, line=-4.5, cex=0.8)
    dev.off()

    ###plot final con todos los paneles con proporciones
    #ploteamos seas áreas sonbre la idoneidad de hábitat, pero con un nivel de transparencia dependiente del número casos en los que la celda ha caído dentro del rango filo (incertidumbre; alpha=second_ensamble_phylo). Por tanto, aquellas zonas que han caído dentro del rango para pocos scenarios ó solo para una de las varuables se ven poco (podría verse con solo una variable dentro del rango sin entrease en el rango bajo muuuchos escenarios).
    pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/with_proportions/", species, "_with_proportions.pdf", sep=""), width=8.27, height=11.69)
    par(oma=c(0,1.5,6,2), mai = c(0.45, 0.45, 0.45, 0.45))
    par(mfcol=c(3,2))
    
    ###Pannel 1
    plot(current_suit,  main="", axis.args=list(cex.axis=1.5))
    mtext(text="Predicted-current habitat suitability", side=3,line=1.7, outer=FALSE, cex=1, font = 2, adj=0.45)
    #add acuatic bodies
    plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
    #add lines around high precision pts.
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){
        plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=1)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){    
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=1)
    }
    #plot the legend    
    legend("top", legend=c("High precision pts.", "Low precision pts."), lty=c(1,3), cex=1, bty="o", lwd=2, horiz=TRUE, bg="white")
    
    #### add evaluation metrics
    ## load data of evaluation metrics
    median_metrics = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/evaluations/medians_evaluations.csv", sep=",", header=TRUE)

    ##calculate median and sd of AUC and OOB for random forest 
    #auc
    rf_auc_median = round(median_metrics[which(median_metrics$species==species),]$rf_auc_median, 3)
    rf_auc_sd = round(median_metrics[which(median_metrics$species==species),]$rf_auc_sd, 3)
    #oob
    rf_oob_median = round(median_metrics[which(median_metrics$species==species),]$rf_oob_median, 3)
    rf_oob_sd = round(median_metrics[which(median_metrics$species==species),]$rf_oob_sd, 3)

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
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('RF: AUC') == .(rf_auc_median) %+-% .(rf_auc_sd) ~ bold('OOB') == .(rf_oob_median) %+-% .(rf_oob_sd) ~ '%'),cex=.6,text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=-2)  
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('GAM: AUC') == .(gam_auc_median) %+-% .(gam_auc_sd) ~ bold('Kappa') == .(gam_kappa_median) %+-% .(gam_kappa_sd) ~ bold('TSS') == .(gam_tss_median) %+-% .(gam_tss_sd)),cex=.6, text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=-1) 
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('GLM: AUC') == .(glm_auc_median) %+-% .(glm_auc_sd) ~ bold('Kappa') == .(glm_kappa_median) %+-% .(glm_kappa_sd) ~ bold('TSS') == .(glm_tss_median) %+-% .(glm_tss_sd)),cex=.6, text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=0) 
 
    ###Pannel 2
    plot(current_suit,  main="", axis.args=list(cex.axis=1.5))
    mtext(text="Predicted-current habitat suitability", side=3,line=1.7, outer=FALSE, cex=1, font = 2, adj=0.3)
    #add acuatic bodies
    plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
    #plot low precision pts.
    if(nrow(low_precision_ocurrences)>0){
        points(low_precision_ocurrences$longitude, low_precision_ocurrences$latitude, cex=0.5, col="white", bg=NA, lwd=0.6, type="p", pch=21)
    }
    #plot high precision pts.
    if(nrow(high_precision_ocurrences)>0){
        points(high_precision_ocurrences$longitude, high_precision_ocurrences$latitude, cex=0.5, col="black", bg=NA, lwd=0.2, type="p", pch=21)
    } 
    #plot the legend    
    legend("top", legend=c("High precision pts.", "Low precision pts."), fill=c("black", "white"), cex=1.1, bty="o", horiz=TRUE, bg="white", x.intersp=0.2) #We reduce the interspace between characters to reduce the space occupied by the legend and then we can increase the size of the legend


    ###Pannel 3
    plot(projected_suit, main="", , axis.args=list(cex.axis=1.5))
    mtext(text="Predicted habitat suitability in 2070", side=3,line=1.7, outer=FALSE, cex=1, font = 2)
    #add acuatic bodies
    plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
    #add lines around high precision pts.
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){       
        plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=1)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){       
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=1)
    }
    #plot the legend    
    #legend("top", legend=c("High precision pts.", "Low precision pts."), lty=c(1,3), cex=1.1, bty="o", lwd=2, horiz=TRUE, bg="white")
    
    ###Pannel 4
    plot(projected_suit, main="", axis.args=list(cex.axis=1.5))
    mtext(text="Predicted habitat suitability in 2070", side=3, line=2.6, outer=FALSE, cex=1, font = 2, adj=0.3)
    mtext(text="+", side=3,line=1.7, outer=FALSE, cex=1, font = 2)
    mtext(text="Phylo-corrected suitability", side=3,line=0.8, outer=FALSE, cex=1, font = 2, adj=0.3)
    #add acuatic bodies
    plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
    #add lines around high precision pts.
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==1,])>0){
        plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=1)
    }
    #add polygon of cirtifield distribution + buffer
    if(nrow(ocurrence_data[ocurrence_data$precision_weight==0.5,])>0){   
        plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=1)
    }
    #add info phlo
    plot(final_ensamble_phylo_proportion, col="#009900", alpha=second_ensamble_phylo_proportion, add=TRUE, legend=FALSE) 
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
    number_figure = toString(i)#extract number figure as string. In that way the number can be printed in bold
    mtext(bquote(bold("Figure") ~ bold(.(number_figure)) ~ bold("Appendix S11:") ~ "Projections of habitat suitability under current and future"), outer = TRUE, cex = 1.2, font=2, line=-56.5, adj=0.3, family="Times") #bquote is used to convert to italic the specific epithet (see http://stackoverflow.com/questions/27266398/using-italics-in-the-title-on-an-object-from-a-dataframe)
        #* let us to add strings without spaces between them
    mtext(bquote("conditions for" ~ italic('P.') ~italic(.(species))*"."), outer = TRUE, cex = 1.2, font=2, line=-58, adj=0.0725, family="Times")

    #number of PDFs
    mtext(side = 1, text = range_pages[i], outer = TRUE, adj=0.88, line=-4.5, cex=0.8) 
    dev.off()   
}
plot_current_main = function(species, legend=TRUE, small=FALSE, size="big", title_species=TRUE){

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
    if(size=="small"){ #for small species narrower buffer
        polygon_plot_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=2.5, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)
    }
    if(size=="medium") {
        polygon_plot_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=10, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) 
    }
    if(size=="big") {
        polygon_plot_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=15, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) 
    }  

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

    #load current suitability
    current_suit = raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_final/ensamble_predictions_bin/ensamble_predictions_bin_", species, ".tif", sep=""))

    #crop current suitability
    current_suit = crop(current_suit, polygon_plot_buffer)
    
    #create a raster with acuatic bodies
    aquatic_bodies =  raster(extent(current_suit), resolution=res(current_suit))
    aquatic_bodies[which(is.na(getValues(current_suit)))] <- 1        

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

    #plot        
    plot(current_suit,  main="", cex.axis=1.5, axis.args=list(cex.axis=1.8)) #axis.args=list(cex.axis=1.8) is for the size of the colour legend
    
    #add title
    if(title_species==FALSE){

    } else {
        mtext(text=bquote(italic('Pinus') ~italic(.(species))), side=3,line=2, outer=FALSE, cex=2.8, font = 2)
    }
        
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
    #plot the legend if
    if(legend==TRUE){  
        legend("top", legend=c("High precision pts.", "Low precision pts."), lty=c(1,3), cex=1.455, bty="o", lwd=2, horiz=TRUE, bg="white", x.intersp=0.2) #We reduce the interspace between characters to reduce the space occupied by the legend and then we can increase the size of the legend
    }  

    #### add evaluation metrics
    ## load data of evaluation metrics
    median_metrics = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/evaluations/medians_evaluations.csv", sep=",", header=TRUE)

    ##calculate median and sd of AUC and OOB for random forest 
    #auc
    rf_auc_median = round(median_metrics[which(median_metrics$species==species),]$rf_auc_median, 3)
    rf_auc_sd = round(median_metrics[which(median_metrics$species==species),]$rf_auc_sd, 3)
    #oob
    rf_oob_median = round(median_metrics[which(median_metrics$species==species),]$rf_oob_median, 3)
    rf_oob_sd = round(median_metrics[which(median_metrics$species==species),]$rf_oob_sd, 3)

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
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('RF: AUC') == .(rf_auc_median) %+-% .(rf_auc_sd) ~ bold('OOB') == .(rf_oob_median) %+-% .(rf_oob_sd) ~ '%'),cex=0.95,text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=-2, x.intersp=0.5)  
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('GAM: AUC') == .(gam_auc_median) %+-% .(gam_auc_sd) ~ bold('Kappa') == .(gam_kappa_median) %+-% .(gam_kappa_sd) ~ bold('TSS') == .(gam_tss_median) %+-% .(gam_tss_sd)),cex=0.95, text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=-1, x.intersp=0.42) 
    legend(xmin(current_suit), ymin(current_suit),legend=bquote(bold('GLM: AUC') == .(glm_auc_median) %+-% .(glm_auc_sd) ~ bold('Kappa') == .(glm_kappa_median) %+-% .(glm_kappa_sd) ~ bold('TSS') == .(glm_tss_median) %+-% .(glm_tss_sd)),cex=0.95, text.col="black", box.col="black",bg="white", adj=c(0.035,0.5), xjust=0, yjust=0, x.intersp=0.5) 
}
plot_proj_main = function(species, proportions, legend=TRUE, small=FALSE, size="big", title_species=TRUE){

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
    if(size=="small"){ #for small species narrower buffer
        polygon_plot_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=2.5, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)
    }
    if(size=="medium") {
        polygon_plot_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=10, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) 
    }
    if(size=="big") {
        polygon_plot_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=15, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) 
    }  
    
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
    projected_suit = raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_final/ensamble_projections_bin/ensamble_projections_bin_", species, ".tif", sep=""))

    #crop projected suitability
    projected_suit = crop(projected_suit, polygon_plot_buffer)

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

    #plot
    plot(projected_suit, main="", cex.axis=1.5, axis.args=list(cex.axis=1.8))
    
    #add title
    if(title_species==FALSE){

    } else {
        mtext(text=bquote(italic('Pinus') ~italic(.(species))), side=3, line=2, outer=FALSE, cex=2.8, font = 2)
    }
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

    #plot the legends if   
    if(legend==TRUE){
        legend("top", legend=c("High precision pts.", "Low precision pts."), lty=c(1,3), cex=1.455, bty="o", lwd=2, horiz=TRUE, bg="white", x.intersp=0.2) #We reduce the interspace between characters to reduce the space occupied by the legend and then we can increase the size of the legend  
        #legend("bottom", legend="Phylo-corrected suitability", fill="#00ffe9", bty="o", cex=1.1, bg="white")
    } 

    #add info phlo
    if(!proportions == TRUE){ #if we don't want proportions
        plot(final_ensamble_phylo, col="#009900", alpha=second_ensamble_phylo, add=TRUE, legend=FALSE)
            #col="#0000FF": color azul oscuro; ALTERNATIVA el negro ("#0000FF"). 
            #alpha=second_ensamble_phylo: transparencia dependiente de la incertidumbre (número casos en los que la celda ha caído dentro del rango filo)
            #add=TRUE: Add to the previous plot
            #legend=FALSE: Not add legend, in that way the legen of the first plot will not be hidden
    } else { #if yes
        plot(final_ensamble_phylo_proportion, col="#009900", alpha=second_ensamble_phylo_proportion, add=TRUE, legend=FALSE) 
        #col="#0000FF": color azul oscuro; ALTERNATIVA el negro ("#0000FF"). 
        #alpha=second_ensamble_phylo: transparencia dependiente de la incertidumbre (número casos en los que la celda ha caído dentro del rango filo)
        #add=TRUE: Add to the previous plot
        #legend=FALSE: Not add legend, in that way the legen of the first plot will not be hidden
    }    
}

##################
##### SUPLES #####
##################

#species with medium and small distribution
medium_species = c("amamiana", "cubensis", "culminicola", "dalatensis", "fragilissima", "krempfii", "kwangtungensis", "maestrensis", "maximartinezii", "morrisonicola", "quadrifolia", "radiata", "rzedowskii", "squamata", "taiwanensis", "torreyana", "tropicalis", "washoensis")
small_species = c("canariensis", "luchuensis")

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

#drop discolor (problem tazonomy, no diferetiaced from cembriodes). Also remove tecunumi y jaslicana, for which we didn't create maps in the beginning
epithet_species_list=epithet_species_list[which(!epithet_species_list %in% c("discolor", "jaliscana", "tecunumanii"))]

#check it
!c("discolor", "jaliscana", "tecunumanii") %in% epithet_species_list
length(epithet_species_list) == 112

#for each path, extract the specie names and plot the suple
for(i in 1:length(epithet_species_list)){

    #select the [i] speceis
    selected_name = epithet_species_list[i]

    #plot the suple
    if(selected_name %in% small_species){
        plot_suples(selected_name, size="small")
    } else{
        if(selected_name %in% medium_species){
            plot_suples(selected_name, size="medium")
        } else {
            plot_suples(selected_name, size="big")        
        }
    } 
}

#bind all pdf generated into one single pdf file with pdftk
system("rm /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/full/full_with_proportions.pdf;

    cd /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/with_proportions; 

    pdftk *.pdf cat output /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/full/full_with_proportions.pdf") #first delete the prevous full file in full directory; next change the drectory to the that with separate plots; finally run pdftk to bind all pdf in that directory and save the result in full directory

system("rm /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/full/full_without_proportions.pdf;

    cd /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/without_proportions; 

    pdftk *.pdf cat output /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple_a4/full/full_without_proportions.pdf") #first delete the prevous full file in full directory; next change the drectory to the that with separate plots; finally run pdftk to bind all pdf in that directory and save the result in full directory



####PROPUESTA DE SPECIES PARA EL MAIN TEXT 

###############################
##### CURRENT SUITABILITY #####
###############################
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/main_text/current_suitability.pdf", width=12, height = 12)
par(oma=c(0,0,2.7,2))
par(mfrow=c(2,2))
plot_current_main(species="nigra", legend=TRUE)
plot_current_main(species="banksiana", size="medium", legend=FALSE)
plot_current_main(species="wallichiana", legend=FALSE)
plot_current_main(species="cembroides", legend=FALSE)
dev.off()


################################################
##### PROJECT SUITABILITY with proportions #####
################################################
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/main_text/projected_suitability_with_proportions.pdf", width=12, height = 12)
par(oma=c(0,0,2.7,2))
par(mfrow=c(2,2))
plot_proj_main(species="nigra", proportions=TRUE, legend=TRUE)
plot_proj_main(species="banksiana", size="medium", proportions=TRUE, legend=FALSE)
plot_proj_main(species="wallichiana", proportions=TRUE, legend=FALSE)
plot_proj_main(species="cembroides", proportions=TRUE, legend=FALSE)
dev.off()

###################################################
##### PROJECT SUITABILITY without proportions #####
####################################3##############
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/main_text/projected_suitability_without_proportions.pdf", width=12, height = 12)
par(oma=c(0,0,2.7,2))
par(mfrow=c(2,2))
plot_proj_main(species="nigra", proportions=FALSE, legend=TRUE)
plot_proj_main(species="banksiana", size="medium", proportions=FALSE, legend=FALSE)
plot_proj_main(species="wallichiana", proportions=FALSE, legend=FALSE)
plot_proj_main(species="cembroides", proportions=FALSE, legend=FALSE)
dev.off()


#############################
##### Figures for talks #####
#############################

#For lab kathleen with proportions
pdf("/Volumes/GoogleDrive/My Drive/science/talks_across_labs/talks_USA_2018/fotos/figure_sdms.pdf", width=14, height = 6)
par(mfcol=c(1,2))
par(oma=c(0,0,2.7,2))
plot_current_main(species="nigra", legend=TRUE, title_species=FALSE)
plot_proj_main(species="nigra", proportions=TRUE, legend=FALSE, title_species=FALSE)
dev.off()

