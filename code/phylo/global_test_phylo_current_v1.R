#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
    #https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
    #https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



#########################################################
####### COMPARE CURRENT EXSITU PREDICTIONS WITH AND WITHOUT PHYLO ####### #########################################################

#This script compares the results of predicting current suitability under current conditions outside the natural distribution of the pines.



###################################################
##### DIFFERENCES RESPECT TO PREVIOUS VERSION #####
###################################################

#Respect to version 1:



########################
##### BEGIN SCRIPT #####
########################

#set the seed for reproducibility
set.seed(56756)

#create some folders
system("mkdir -p ./results/global_test_phylo_current/exsitu_occurrences")
system("mkdir -p ./results/global_test_phylo_current/env_predictors")
system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_no_phylo/", sep=""))

#pre-defined functions
plot_sin=function(input){
    jpeg("./singularity_plot.jpeg", height=2000, width=2000, res=300)
    plot(input)
    dev.off()
}

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
    #Multiply both variables for obtaining a raster with all NAs. We will use them to ensure that all presences have environmental data, and specific we want to ensure that 0.5 fall in areas with value of environmental variables. We use a soil variable because there is bioclim data for some water bodies inside the continents, and we don't want to get presences in these areas. 

#reduce resolution of environmental variable
#load buffer albicaulis to get the target resolution
buffer_albicaulis = raster(paste("results/ocurrences/albicaulis_distribution_buffer", ".asc", sep=""))
#prepare fact value for aggregate
fact_value=res(buffer_albicaulis)[1]/res(environment_var)[1]
    #in raster::aggregate, fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
    #we calculate it dividing the target (coarser) resolution by the original (finer) resolution
#aggregate cells
environment_var_low_res=raster::aggregate(environment_var, fact=fact_value, fun=mean)
    #We use mean as the function to aggregate values. This is the same than using resample with "bilinear" method, as you can see in the code of "resample". 
        #https://github.com/cran/raster/blob/4d218a7565d3994682557b8ae4d5b52bc2f54241/R/resample.R#L42
    #See also dummy example with both functions:
        #r <- raster(nrow=4, ncol=8)
        #r2 <- raster(nrow=2, ncol=4)
        #r <- setValues(r, values = 1:32)
        #r_agg <- aggregate(r, fact=2, fun=mean)
            #4/2=2 and 8/2=4, so r2 is half r. indeed resolution of r is 45x45 and resolution of r2 is 90x90
            #therefore, aggregating 2 cells in r should give r2
        #r_resam <- resample(r, r2, method="bilinear")
            #now we use directly r2 as source resolution
        #getValues(r_resam) == getValues(r_agg)
        #https://gis.stackexchange.com/a/255155
#check we have the correct resolution
if(res(environment_var_low_res)[1]!=0.5 | res(environment_var_low_res)[2]!=0.5){
    stop("ERROR! FALSE! WE HAVE A PROBLEM REDUCING THE RESOLUTION OF THE ENVIRONMENTAL VARIABLE")
}

#load elevation at 10x10 resolution
elev = raster("./datos/topografia/elev_low_resolution.asc")
    #we will calculate the percentile of altitude inside each cell of 50x50 inside the buffer, so we don't need great resolution. Only 10x10 is enough (0.08333334, 0.08333334)
#check elevation has a similar resolution of moisture bioclimvariables up to the 7th decimal
if((round(res(elev)[1], 7) != round(res(environment_var)[1], 7)) | (round(res(elev)[2], 7) != round(res(environment_var)[2], 7))){
    stop("ERROR! FALSE! WE HAVE A PROBLEM WITH THE ELEVATION VARIABLE")
}

#extract occurrence data
#we are using data from Perret et al 2018
    #the dataset of Perret is carefully curated, obtaining data from herbarium specimens and accounts in the literature. They only consider as naturalized those records with clear evidence of self-sustainability and geographical coordinates. This solves a very important limitation compared to just use our GBIF data, as we do not include occurrences like from gardens... This is also a great independent validation, as the data does not come from GBIF, it is completely independent data. Different regions and different data sources.
    #Naturalized distributions show that climatic disequilibrium is structured by niche size in pines (Pinus L.)

#extract and load the naturalized occurrences from zip file
system(paste("unzip -o ./datos/phlyo/method_validation/doi_10_5061_dryad_1hr1n52__v20181213.zip pinus_occurrences_fordryad_exoticonly.csv -d ./results/global_test_phylo_current/exsitu_occurrences/", sep=""))
naturalized_occurrences=read.csv("./results/global_test_phylo_current/exsitu_occurrences/pinus_occurrences_fordryad_exoticonly.csv", header=TRUE)
summary(naturalized_occurrences)
#check
if(nrow(naturalized_occurrences)!=597){
    stop("ERROR! FALSE! WE DO NOT HAVE ALL NATURALIZED OCCURRENCES WE SHOULD HAVE ACCORDING TO THE PERRET'S MANUSCRIPT")
}
if(sum(unique(naturalized_occurrences$species) %in% epithet_species_list)!= length(unique(naturalized_occurrences$species))){
    stop("ERROR! FALSE! WE DO NOT HAVE THE SAME SPECIES NAMES IN PERRET DATA")
}

#NOTE: In general, I have used extract to obtain any information about the points: elevation, environmental data and number of cell. Extract only consider that a point falls inside a cell if its center is inside that cell. It is important consider this. 




###################################################
##### SELECT OCCURRENCES OUTSIDE DISTRIBUTION #####
###################################################

#packages
require(raster)
require(sf)

#species="radiata"
exsitu_occurrences=function(species){

    ##stop if pumila, just in case, because this species has part of its distribution in a corner and we should check it carefully
    if(species=="pumila"){
        stop("ERROR! FALSE! YOU ARE TRYING TO ANALYZE PUMILA, AND WE HAVE NOT CHECKED THIS SCRIPT FOR THAT")
    }

    ##open folders
    system(paste("mkdir -p ./results/global_test_phylo_current/exsitu_occurrences/", species, sep=""))

    ##load the distribution buffer
    #read the raster
    distribution_buffer = raster(paste("results/ocurrences/", species, "_distribution_buffer.asc", sep=""))

    #convert to polygon
    polygon_distribution_buffer = rasterToPolygons(distribution_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon


    ##load the buffer used to create PAs.
    #This is out area of known absences due to edaphoclimatic conditions, i.e., no migration, etc... so we will exclude it from these analyses.
    #Reasons
        #It is more fair for the non-phylo models. 
            #Of course, our SDMs are not going to predict well naturalized occurrences in the PA buffer. This buffer is considered as the area with absences truly caused by edaphoclimatic conditions because it is outside but close to the natural distribution.
            #if we include these regions in the evaluations, the SDMs are going to be worse, and this could inflate the performance of the phylogenetic correction.
            #yes, it would be better if no naturalized presences is inside the PA buffer because that means we have selected an area that is truly unsuitable for the species, but in general I have not seen many presences inside, so I think we are ok, our PA buffer is not too big.
        #we want an actual validation of the model predictions, 
            #we need to use an area not used in any way during modeling and occurrences obtained from other sources.
    raster_pa_buffer = raster(paste("./results/pseudo_absences/", species, "_PA_buffer.asc", sep=""))  
    
    #convert to polygon 
    polygon_raster_pa_buffer = rasterToPolygons(raster_pa_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon


    ##obtain a polygon for all the world outside of the PA buffer
    #we will use the environmental raster, as it already includes areas with environmental data (both climatic and edaphic) across the globe and exclude water bodies
    #do an inverse mask of the environmental variable using the PA buffer as mask. This will give us all areas outside of the PA buffer.
    raster_outside_pa_buffer=mask(environment_var_low_res, polygon_raster_pa_buffer, inverse=TRUE)

    #convert to 1 all cells with no NA, i.e., with edaphoclimatic data
    raster_outside_pa_buffer[which(!is.na(getValues(raster_outside_pa_buffer)))] = 1

    #create a polygon with all rows having edaphoclimatic data
    raster_outside_pa_buffer_polygon = rasterToPolygons(raster_outside_pa_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE)

    #create the same polygon but with cell borders inside just for checking occurrences resampling in the plots
    raster_outside_pa_buffer_polygon_non_dissolved = rasterToPolygons(raster_outside_pa_buffer, fun=function(x){x==1}, n=16, dissolve = FALSE)

    #select naturalized occurrences for the selected species
    occurrences=naturalized_occurrences[which(naturalized_occurrences$species==species),]
    if(unique(occurrences$species)!=species){
        stop("ERROR! FALSE! WE HAVE A PROBLEM SELECTING THE NATURALIZED OCCURRENCES OF THE SPECIES")
    }
    occurrences=occurrences[,c("longitude", "latitude")]
    
    #drop rows with NA for longitude and latitude
    if(length(which(is.na(occurrences$longitude)))>0 | length(which(is.na(occurrences$latitude)))>0){
        stop("ERROR! FALSE! WE HAVE A PROBLEM SELECTING THE NATURALIZED OCCURRENCES OF THE SPECIES, WE HAVE NA FOR COORDINATES")  
    }
    occurrences<-occurrences[!(is.na(occurrences$longitude)),] 
    occurrences<-occurrences[!(is.na(occurrences$latitude)),] 

    #drop points without environmental variables.
    environment_presences = extract(x=environment_var, y=occurrences[c("longitude", "latitude")]) 
        #extract coge los puntos d epresencia, coge las varialbes, y para cada punto de presencai coge los valoers de las variables. Los puntos que quedan fuera del area de las variables se pone NA que luego quitaremos.
        #we use the raster with 10x10km cell size because we will 10x10 raster as input in the prediction of probability of presence for these points 
    
    #no es un data.frame, lo convertimos, 
    environment_presences<-data.frame(environment_presences) 
        #este data.frame tiene una columna por variable, una fila por punto de presencia, y esa tabla la vamos a unir con la table de occurrences (coordenadas pais, y todo lo demas)

    #lo unimos a las presencias
    occurrences<-data.frame(occurrences, environment_presences)
    nrow(occurrences)

    #quitamos los registros que tienen nulos en los valores de variables ambientales
    occurrences<-occurrences[!(is.na(occurrences$environment_presences)), ] 
    nrow(occurrences)
        #in this way we drop the points without environmental data (bioclim and soil)
    #stop if we lose more than 2 occurrences
    if(nrow(environment_presences)-nrow(occurrences)>2){
        stop(paste("ERROR! FALSE! WE HAVE LOST MORE THAN 2 OCCURRENCES DUE TO LACK OF ENVIRONMENTAL DATA FOR SPECIES ", species, sep=""))
    }
    rm(environment_presences)

    #create a coordinatePrecision variable if it is neccesary
    if (length(occurrences$coordinatePrecision)==0){
        occurrences$coordinatePrecision = rep(NA, nrow(occurrences))
    }


    ##data cleaning
    #Select ocurrences points inside the buffer
    outside_pa_buffer_values_of_points = extract(raster_outside_pa_buffer, occurrences[,c("longitude","latitude")]) 
        #extract the values of each ocurrence in the exsitu raster, the possibilities are two: 1 or NA. 1 inside, NA outside buffer. We use the raster instead of the because it is faster. 
    points_and_values_outside_pa_buffer = cbind(occurrences[,c("longitude","latitude")], outside_pa_buffer_values_of_points) #bind these values with the corresponding presences points
    points_outside_pa_buffer = which(!is.na(points_and_values_outside_pa_buffer$outside_pa_buffer_values_of_points)) #obtain rows in which values of the distribution buffer raster is 1, thus points outside the buffer
    points_inside_pa_buffer = which(is.na(points_and_values_outside_pa_buffer$outside_pa_buffer_values_of_points)) #obtain rows in which values of the distribution buffer raster is NA, thus points inside the buffer

    #create DF with the number of occurrences in and outside of the PA buffer before the resampling
    n_points_in_out_pa_buffer=cbind.data.frame(nrow(points_and_values_outside_pa_buffer), length(points_outside_pa_buffer), length(points_inside_pa_buffer))
    names(n_points_in_out_pa_buffer)=c("total_points", "points_outside", "points_inside")
    if(n_points_in_out_pa_buffer$points_inside > (n_points_in_out_pa_buffer$total_points*0.1)){
        stop(paste("ERROR! FALSE! WE HAVE A PROBLEM, MORE THAN 10% OF OCCURRENCES OF PERRET FALL WITHIN PA BUFFER FOR SPECIES ", species, sep=""))
    }
        #we can have occurrences inside the PA buffer as this buffer is a large area around the natural distribution of pines. It is the whole area used for modeling, including, not only presences but also pseudoabsences.

    #plot the result
    cairo_pdf(paste("./results/global_test_phylo_current/exsitu_occurrences/", species, "/", species, "_ocurrences.pdf", sep=""), width=12, height=12)
    plot(environment_var, col="gray80")
    plot(polygon_distribution_buffer, add=TRUE, lwd=0.1, border="black")
    plot(raster_outside_pa_buffer_polygon, add=TRUE, lwd=0.1, border="red")
    plot(polygon_raster_pa_buffer, add=TRUE, lwd=0.1, border="blue")
    points(points_and_values_outside_pa_buffer[points_inside_pa_buffer,]$longitude, points_and_values_outside_pa_buffer[points_inside_pa_buffer,]$latitude, cex=0.01, lwd=0.1, pch=20, col="blue")
    points(points_and_values_outside_pa_buffer[points_outside_pa_buffer,]$longitude, points_and_values_outside_pa_buffer[points_outside_pa_buffer,]$latitude, cex=0.01, lwd=0.1, pch=20, col="red")
    legend(x="topright", legend=c("Distribution buffer", "PA buffer", "Outside PA buffer"), fill=c("black", "blue", "red"), cex=0.8)
    dev.off()

    #drop the ocurrences outside the buffer in the presencia data frame
    presencia = occurrences[points_outside_pa_buffer,]

    #we stop here if there is no presence points
    if(nrow(presencia)>0){

        #LIMPIEZA DE DUPLICADOS EN LAS COORDENADAS
        #-----------------------------------------
        #buscamos registros duplicados en las coordenadas
        duplicados<-duplicated(presencia[ , c("latitude", "longitude")]) #la funcion duplicated del paqeute dismo busca registros duplicados para unas columnas que yo le diga. 
        #¿cuantos duplicados hay?
        if(length(duplicados[duplicados==TRUE]) > (nrow(presencia)*0.40)){
            stop(paste("ERROR! FALSE! WE HAVE LOST MORE THAN 40% OF OCCURRENCES DUE TO DUPLICATION FOR SPECIES ", species, sep=""))
        }    
        #selecciona los no duplicados (el símbolo ! significa "todo menos los duplicados")
        presencia<-presencia[!duplicados, ]
        #comprobamos que no han quedado duplicados (solo para cultivar vuestra fé)
        duplicados_2<-duplicated(presencia[ , c("latitude", "longitude")])
        if(sum(duplicados_2) != 0){
            stop("ERROR! FALSE! WE HAVE A PROBLEM REMOVING DUPLICATE OCCURRENCES")
        }


        ##Resampling ocurrences points##########
        #we are not using GBIF data here, so we do not have information about coordinate precision, so we should not use altitudinal resampling, just random sampling ensuring we have 3 occurrences per cell 
        #extrac values of elevation for presences
        elevation_of_presences = extract(elev, presencia[,c("longitude","latitude")])

        #Calculate the cell in which each point is located
        #we create a raster with the number of cell as values
        index_raster <- raster_outside_pa_buffer #the raster is created from raster_outside_pa_buffer
        index_raster[] <- 1:ncell(raster_outside_pa_buffer) #the new raster take as value the number of cells of raster_outside_pa_buffer
        index_raster <- mask(index_raster, raster_outside_pa_buffer) #remove now all cells that are not included the outside PA buffer (cells inside the PA buffer or water bodies)
        #Calculate the cell from halepensis map (buffer+bianca) in which each point is located
        cell_id_presences= extract(index_raster, presencia[,c("longitude","latitude")])
        #create a variable as id_presencia
        id_presence=1:nrow(presencia[ ,c("longitude","latitude")])
        #join data
        table_stratified_sample = cbind(cell_id_presences, id_presence, presencia[,c("longitude","latitude", "coordinatePrecision")], elevation_of_presences)

        #check 0
        #check that we have only one ID per occurrence
        if(length(unique(table_stratified_sample$id_presence))!=nrow(table_stratified_sample)){
            stop("ERROR! FALSE! WE HAVE A PROBLEM GIVING AN ID TO EACH OCCURRENCE DURING THE RESAMPLING")
        }

        #check 1
        #check that id of cell selected with presences has an 1 in raster buffer (inside buffer)
        if(unique(getValues(raster_outside_pa_buffer)[cell_id_presences])!=1){
            stop("ERROR! FALSE! WE HAVE A PROBLEM PREPARING OCCURRENCES OUTSIDE PA BUFFER FOR RESAMPLING")
        }

        #check 2
        #check that cell indicated for each point is correct
        test = index_raster #create a raster from index raster (cell numbers)
        test[] <- NA #set all NA
        test[cell_id_presences] <- index_raster[cell_id_presences] #select only cells with presences and put their number, i.e., the number of the cell
        #select the ids without repetition and ordered (less to more)
        cell_id_presences_order = sort(unique(cell_id_presences))
        #plot
        require(viridis)    
        colours_viridis = viridis(length(cell_id_presences_order)) #colour palette
        #open png
        png(paste("./results/global_test_phylo_current/exsitu_occurrences/", species, "/", species, "_check_cell_number.png", sep=""), width=1100, height=800, pointsize=30)
        #plot raster with cell numbers of cells with presences
        plot(test)
        #for each id cell with presence
        #k=1
        for(k in 1:length(cell_id_presences_order)){
            #select presences of the [k] cell
            selected_cell = table_stratified_sample[which(table_stratified_sample$cell_id_presences == cell_id_presences_order[k]),]
            #add these points to the plot 
            points(selected_cell[,"longitude"], selected_cell[,"latitude"], col=colours_viridis[k])
        }
        #add legend
        legend("topright", legend=cell_id_presences_order, fill=colours_viridis, cex=0.3)
        dev.off() 
            #you have to check that each color hf the first legend (inside plot) correspond to the points included in only ONE cell. Then, you have to cheack the the color of the cell correspond in second legend with the same number than in the second: Por ejemplo: La primera celda tiene solo PUNTOS de color muuuuy violetas, todos están ahi métidos, según la primera leyenda eso corresponde XXXXX. Ahora miras el color de la CELDA, y te vas a la segunda leyenda, su color debe corresponde con el número XXXXX ó cercano. Comprobado con canariensis y todo correcto. Además los puntos seleccionados al final son 3 por celdas en varias especies que he chequeado.
            #the point here is to visually check for a few species that the value of the cell is indeed the index of the cell, i.e., the first cell has value 1 and all points inside that cell have a value of 1 for the cell value index value. Cell 2 has value of 2....
            #so we are correctly selecting the number of each cell


        ##Resampling
        #create a categorical variable of precision
        #we consider as high precision ocurrences those with a coordinatePrecision between 1 and 25, the rest of ocurrences are considered as low precision.  

        #create the empty variable for precision
        precision_weight = NULL

        #run loop across presences
        if(nrow(table_stratified_sample)>0){
            #k=1
            for (k in 1:nrow(table_stratified_sample)){ #for each row of table_stratified_sample
                subset=table_stratified_sample[k,] #subset this row
                if (is.na(subset$coordinatePrecision)){ #if the point have cP=Na
                    precision_weight = append(precision_weight, 0.5) #precision_weight = 0.5
                } else { #if not
                    if (subset$coordinatePrecision>=1 && subset$coordinatePrecision<=25){ #if the point has a cP between 1 and 25, inclugin both extremes. 
                        precision_weight = append(precision_weight, 1)
                    } else {
                        precision_weight = append(precision_weight, 0.5)
                    } 
                }   
            }
        } 

        #check
        if(length(precision_weight)!=nrow(table_stratified_sample)){
            stop("ERROR! FALSE! WE HAVE NOT CALCULATED THE PRECISION WEIGHT FOR ALL OCCURRENCES")
        }
        if(length(unique(precision_weight))!=1){
            stop("ERROR! FALSE! WE SHOULD HAVE ALL OCCURRENCES AS PRECION =0.5, BUT THIS IS NOT THE CASE")
        } else {
            if(unique(precision_weight)!=0.5){
                stop("ERROR! FALSE! WE SHOULD HAVE ALL OCCURRENCES AS PRECION =0.5, BUT THIS IS NOT THE CASE")
            }  
        }

        #bind precision_weight with table_stratified_sample
        table_stratified_sample = cbind(table_stratified_sample, precision_weight)    

        #check the weights have been correctly calculated
        subset_1=table_stratified_sample[is.na(table_stratified_sample$coordinatePrecision),]$precision_weight
        subset_2=table_stratified_sample[!is.na(table_stratified_sample$coordinatePrecision) & table_stratified_sample$coordinatePrecision<1 | !is.na(table_stratified_sample$coordinatePrecision) & table_stratified_sample$coordinatePrecision>25,]$precision_weight
        subset_3=table_stratified_sample[!is.na(table_stratified_sample$coordinatePrecision) & table_stratified_sample$coordinatePrecision>=1 & table_stratified_sample$coordinatePrecision<=25,]$precision_weight
        checks=list()
        if(length(subset_1)>0){
            if(length(unique(subset_1))==1){
                checks[[1]]=unique(subset_1)==0.5
            } else { 
                checks[[1]]=FALSE
            }
        }
        if(length(subset_2)>0){
            checks[[2]]=FALSE
        }
        if(length(subset_3)>0){
            checks[[3]]=FALSE
        }
        if(FALSE %in% checks){
            stop("ERROR! FALSE! WE HAVE A PROBLEM WHEN CALCULATING THE PRECISION WEIGHT FOR OCCURRENCES") 
        }

        #create the variable for final presences of gbif
        final.presences = NULL

        #make the resampling
        if (nrow(table_stratified_sample)>0){ #if there is gbif points with a high precision
            #cell_ids_several_high_precision_points=table_stratified_sample[duplicated(table_stratified_sample$cell_id_presences) & table_stratified_sample$precision_weight==1,]$cell_id_presences
                #all cells with several high precision points
            #for(l in cell_ids_several_high_precision_points){if(length(which(cell_ids_several_high_precision_points==l))==3){k=l}}
                #cell with 3 high precision points
            #for(l in cell_ids_several_high_precision_points){if(length(which(cell_ids_several_high_precision_points==l))<2){k=l}}
                #cell with less than 3 high precision points. You have to use 2, because less than 3 means you can have at least 2 points in cell_ids_several_high_precision_points, but remember that this variable includes those cells with several presences, i.e., duplicated cell ID, but the first appearance of the cell ID is not included
            #k=unique(table_stratified_sample$cell_id_presences)[1]
            for (k in unique(table_stratified_sample$cell_id_presences)){ #for each cell of the distribution buffer
                subset=table_stratified_sample[table_stratified_sample$cell_id_presences==k,] #select the points inside the [i] cell. 
                subset_high_precision = subset[subset$precision_weight==1,] #subset the high precision points
                subset_high_precision_with_elevation = subset_high_precision[!is.na(subset_high_precision$elevation_of_presences),] #select the points with high precision and elevation data because of the elevation resampling
                subset_low_precision = subset[subset$precision_weight==0.5,] #subset the low precision points
                if (nrow(subset_high_precision_with_elevation)>3){ #if there is more than 3 high precision points
                    elevations_in_cell=subset_high_precision_with_elevation$elevation_of_presences 
                        #select the elevation of these points
                    p10=as.vector(quantile(elevations_in_cell, 0.10, na.rm=T))
                    p50=as.vector(quantile(elevations_in_cell, 0.50, na.rm=T))
                    p90=as.vector(quantile(elevations_in_cell, 0.90, na.rm=T))
                        #calculate the percentile 10, 50 and 90 of the elevation of all points inside the cell 
                    closest_p10 = which(abs(elevations_in_cell-p10)==min(abs(elevations_in_cell-p10))) 
                    closest_p50 = which(abs(elevations_in_cell-p50)==min(abs(elevations_in_cell-p50)))
                    closest_p90 = which(abs(elevations_in_cell-p90)==min(abs(elevations_in_cell-p90))) 
                        #which elevation is closest to each percentile
                    final.presences = rbind(
                    final.presences,
                    subset_high_precision[c(closest_p10[sample(1:length(closest_p10),1)],closest_p50[sample(1:length(closest_p50),1)], closest_p90[sample(1:length(closest_p90),1)]),]) 
                        #select the the points closest to each percentile, the point is selected randomly from the number of points closest to the percentile.
                        #it would have been better to check whether we have the same point closes to two percentiles, but it is very unlikely this happens:
                            #first, a point has to be the closest to TWO percentiles of elevation
                            #second, the point is randomly selected for the first percentile and the second
                            #alternatively, there is no more points closest to these percentiles so it is selected in both cases again
                            #also note the lack of high precision occurrences we have in the data, therefore this seems to be very unlikely to happen. Just a small bug that could produce a bit more correlation of occurrences.
                } else { #if not 
                    if (nrow(subset_high_precision)<=3 && nrow(subset_high_precision)>=1 && nrow(subset_low_precision)==0){ #if the number of high precision points is equal or lower than  3 and equal or higher than 1 and the number of low precision points is equal to 0. From here, we don`t need use elevation because elevation resampling can be made only with more than 3 high precision points and elevation data. In this way we don't lose high precision points without elevation in cells with 3 or less high precision points.
                        final.presences = rbind(
                        final.presences,
                        subset_high_precision) #Select all the points
                    } else { #if not 
                        if (nrow(subset_high_precision)<=3 && nrow(subset_high_precision)>=1 && nrow(subset_low_precision)>0){ #if the number of precision points is equal, lower than and equal or higher than 1 and the number of low precision points is higher to 0.

                            if (nrow(subset_low_precision)>3-nrow(subset_high_precision)){ #if the number of low precision points is higher than the number of point that we need to reach 3 (3 less the number of high precision points)
                                final.presences = rbind(final.presences, 
                                subset_low_precision[sample(x=1:nrow(subset_low_precision), size=3-nrow(subset_high_precision)),], 
                                subset_high_precision) #select all the high precision points and select randomly 3-nrow(subset_high_precision) points form all low precision points
                            } else { #if not
                                final.presences = rbind(final.presences, subset_high_precision, 
                                    subset_low_precision) #select all points
                            }
                        } else{ #if not 
                            if(nrow(subset_high_precision)==0 && nrow(subset_low_precision)>0){ #if the nubmer of high precision points is 0 and the number of low precision ir higher than 0

                                if (nrow(subset_low_precision)>3){ #if the number of low precision points is higher tna 3 
                                    final.presences = rbind(final.presences, subset_low_precision[sample(x=1:nrow(subset_low_precision), size=3),]) #select randomly three low precision points
                                } else { #if not 
                                    final.presences = rbind(final.presences, subset_low_precision) #select all low precision points. 
                                }
                            } 
                        }
                    }
                }   
            }
        } else { #If there is not gbif points 
            final.presences = NULL #create a empty vector
        }

        #stop if
        if(nrow(final.presences) < (nrow(table_stratified_sample)*0.5)){
            stop(paste("ERROR! FALSE! WE HAVE LOST MORE THAN 50% OF OCCURRNCES DUE TO THE RESAMPLING ", species, sep=""))
        }

        #checks
        if(nrow(table_stratified_sample)>0 && nrow(final.presences)>0){
            
            #only 3 points for each cell
            test = NULL
            for (k in unique(final.presences$cell_id_presences)){
                test = append(test, nrow(final.presences[final.presences$cell_id_presences==k,])<=3)
            }
            if(sum(test)!=length(unique(final.presences$cell_id_presences))){
                stop("ERROR! FALSE! WE HAVE A PROBLE, IN RESAMPLING, MORE THAN 3 POINTS IN SOME CELLS")
            }

            #check whether all cell IDs obtained meet one of the conditions to be included regarding number and precision of the points. We take again points using the conditions from the original data.frame
            test_cells = data.frame(id_cell=NA, points=NA, high_precision=NA, high_p_elev=NA, low_precision=NA)
            #k=unique(table_stratified_sample$cell_id_presences)[1]
            for (k in unique(table_stratified_sample$cell_id_presences)){
                subset=table_stratified_sample[table_stratified_sample$cell_id_presences== k,] 
                subset_high_precision = subset[subset$precision_weight==1,] 
                subset_high_precision_with_elevation = subset_high_precision[!is.na(subset_high_precision$elevation_of_presences),] 
                subset_low_precision = subset[subset$precision_weight==0.5,]
                test_cells = rbind(test_cells, c(k, nrow(subset), nrow(subset_high_precision), nrow(subset_high_precision_with_elevation), nrow(subset_low_precision))) #for each cell calculate the different subsets, extract the number of rows (number of points) and then save in a data.frame.
            }
            test_cells=test_cells[-1,] #remove row with NAs
            #check whether the total number of presences inside a cell is the same than the sum of high and low precision presences
            if(sum(test_cells$points == test_cells$low_precision)!=nrow(test_cells)){
                stop("ERROR! FALSE! WE SHOULD HAVE ALL POINTS LOW-PRECISION AFTER RESAMPLING, NOT THE CASE")
            }
            #Extract the cell that accomplish each condition
            test_cells[test_cells$high_p_elev>3,]$id_cell
            test_cells[test_cells$high_precision<=3 & test_cells$high_precision>=1 & test_cells$low_precision==0,]$id_cell
            test_cells[test_cells$high_precision<=3 & test_cells$high_precision>=1 & test_cells$low_precision>0 & test_cells$low_precision>3-test_cells$high_precision,]$id_cell
            test_cells[test_cells$high_precision<=3 & test_cells$high_precision>=1 & test_cells$low_precision>0 & test_cells$low_precision<=3-test_cells$high_precision,]$id_cell  
            test_cells[test_cells$high_precision==0 & test_cells$low_precision>0 &  test_cells$low_precision>3,]$id_cell
            test_cells[test_cells$high_precision==0 & test_cells$low_precision>0 & test_cells$low_precision<=3,]$id_cell

            #test if the different conditions cover all cells
            eso = c(test_cells[test_cells$high_p_elev>3,]$id_cell,
                test_cells[test_cells$high_precision<=3 & test_cells$high_precision>=1 & test_cells$low_precision==0,]$id_cell,
                test_cells[test_cells$high_precision<=3 & test_cells$high_precision>=1 & test_cells$low_precision>0 & test_cells$low_precision>3-test_cells$high_precision,]$id_cell,
                test_cells[test_cells$high_precision<=3 & test_cells$high_precision>=1 & test_cells$low_precision>0 & test_cells$low_precision<=3-test_cells$high_precision,]$id_cell,
                test_cells[test_cells$high_precision==0 & test_cells$low_precision>3,]$id_cell,
                test_cells[test_cells$high_precision==0 & test_cells$low_precision>0 & test_cells$low_precision<=3,]$id_cell)
            if(!identical(sort(eso), sort(unique(final.presences$cell_id_presences)))){
                stop("ERROR! FALSE! WE HAVE A PROBLEM CHECKING THE RESAMPLE")
            }
        }

        #Do the last steps only if we have at least one presence. In some cases like amaniana, there are no occurrences outside of the PA buffer
        if(!is.null(final.presences)){
            
            #if there are high precision points final selected plot them into altitudinal plot to see the performance of the altitudinal sampling.
            if(nrow(final.presences[which(final.presences$precision_weight==1),]) > 1){
                
                #open png
                png(paste("./results/global_test_phylo_current/exsitu_occurrences/", species, "/", species, "_check_altitudinal_sampling.png", sep=""), width=2200, height=1600, res=300)

                #plot elevation
                if(species=="radiata"){
                    #if radiata, plot the canary islands where it has been introduced and where we can see selection of high precision points across the Teide
                    canariensis_buffer=raster(paste("results/ocurrences/canariensis_distribution_buffer", ".asc", sep=""))
                    canariensis_buffer_polygon=rasterToPolygons(canariensis_buffer, fun=function(x){x==1}, dissolve=TRUE)
                    plot(crop(elev, canariensis_buffer_polygon), main="Altitudinal sampling of high precision points")
                    plot(raster_outside_pa_buffer_polygon_non_dissolved, add=TRUE, border="blue", lwd=0.1)
                } else {
                    plot(elev, main="Altitudinal sampling of high precision points")
                    plot(raster_outside_pa_buffer_polygon, add=TRUE, border="blue", lwd=0.02)
                }

                #add all points high precision points before the altitudinal sampling
                points(table_stratified_sample[which(table_stratified_sample$precision_weight==1),]$longitude, table_stratified_sample[which(table_stratified_sample$precision_weight==1),]$latitude, pch=20, cex=0.01, col="black", lwd=0.1)
                
                #add high precision points after altitudinal sampling
                points(final.presences[which(final.presences$precision_weight==1),]$longitude, final.presences[which(final.presences$precision_weight==1),]$latitude, pch=20, cex=0.01, col="red", lwd=0.1)
                #add legend
                legend("topright", legend=c("outside PA buffer", "all points high", "high selected"), fill=c("blue", "black", "red"), cex=0.5)
                dev.off()
            }
                #The resampling works great, the only little problem i see is in the case of cells with a lot of sea. In those cases can be selected 3 points close between them, but i don't think this could have great impact on the models. Only more correlation between very few points (only from buffer cells with a lot of sea).

            #visualize high and low precision occurrences
            png(paste("./results/global_test_phylo_current/exsitu_occurrences/", species, "/", species, "_check_high_low_occurrences.png", sep=""), width=2200, height=1600, res=300)

            #plot elevation
            if(species=="radiata"){
                #if radiata, plot only south east australia to see in more detail
                plot(crop(elev, extent(125, 155, -40, -20)), main="Altitudinal sampling of high precision points")
                plot(raster_outside_pa_buffer_polygon_non_dissolved, add=TRUE, border="blue", lwd=0.1)
            } else {
                plot(elev, main="Altitudinal sampling of high precision points")
                plot(raster_outside_pa_buffer_polygon, add=TRUE, border="blue", lwd=0.02)
            }

            #add high precision points after altitudinal sampling
            points(final.presences[which(final.presences$precision_weight==1),]$longitude, final.presences[which(final.presences$precision_weight==1),]$latitude, pch=20, cex=0.01, col="red", lwd=0.1)

            #add low precision points
            points(final.presences[which(final.presences$precision_weight==0.5),]$longitude, final.presences[which(final.presences$precision_weight==0.5),]$latitude, pch=20, cex=0.01, col="blue", lwd=0.1)

            #add legend
            legend("topright", legend=c("high precision points", "low precision points"), fill=c("red", "blue"), cex=0.5)
            dev.off()
                #checked for pinus radiata in canary islands, 3 points per cell

            #plot final occurrences
            pdf(paste("./results/global_test_phylo_current/exsitu_occurrences/", species, "/", species, "_final_presences.pdf", sep=""), width=12, height=12)
            plot(environment_var, col="gray80")
            points(final.presences$longitude, final.presences$latitude,  cex=0.1, lwd=0.1, pch=20, col="red")
            dev.off()

            #select the columns of occurrences we are interested in
            ultimate_ocurrences = final.presences[, c("longitude", "latitude", "precision_weight")]

            #create a variable of final.presences
            ultimate_ocurrences$presence = rep(x=1, time=nrow(ultimate_ocurrences))

            #last check weight
            if(length(unique(ultimate_ocurrences$precision_weight))==1){
                if(unique(ultimate_ocurrences$precision_weight)!=0.5){
                    stop(paste("ERROR! FALSE! WE HAVE A PROBLEM, STILL HIGH PRECISION WEIGHTS FOR SPECIES ", species, sep=""))
                }
            } else {
                stop(paste("ERROR! FALSE! WE HAVE A PROBLEM, STILL HIGH PRECISION WEIGHTS FOR SPECIES ", species, sep=""))
            }

            #write the resulting data.frame in a csv. 
            write.table(ultimate_ocurrences, paste("./results/global_test_phylo_current/exsitu_occurrences/", species, "/", species, "_final_presences.tsv", sep=""), sep="\t", row.names=FALSE)
        }
    }

    #return the number of points before resampling
    return(n_points_in_out_pa_buffer)
}

#run it for one species
#exsitu_occurrences("radiata")




##############################################
##### PREDICT AND EVALUATE WITHOUT PHYLO #####
##############################################

#packages
require(randomForest)
require(gam)

#species="radiata"
predict_eval_no_phylo = function(species){

    #open folder
    system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/", sep=""))

    #check if we have output from the previous step
    presence_file_exist=file.exists(paste("./results/global_test_phylo_current/exsitu_occurrences/", species, "/", species, "_final_presences.tsv", sep=""))

    #do stuff if we have the presence file
    if(presence_file_exist){

        ###obtain environmental predictors
        ##obtain the list of variables for the selected species
        #load the group of species according to cluster
        group_species = read.csv("./code/variables/variable_selection/species_clustering/_tables/complete_2_10_g_2.csv", header=TRUE)
            #generated in "/home/dftortosa/diego_docs/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/species_clustering_v2.R"

        #select the group of the corresponding species
        variables_cluster = group_species[which(group_species$species == species),]$groups   

        #load names of selected variables  
        load("./results/final_variables/list_selected_variables.rda") 
            #generated in "/home/dftortosa/diego_docs/science/phd/nicho_pinus/code/variables/variable_selection/variable_selection_inside_clusters/variable_selection_inside_clusters_v2.R"

        #take selected variables from the list of variables
        selected_variables = ultimate_variables[[variables_cluster]]

        #calculate the number of ocurrences for each species, this is the number of occurrences used for modelling. We need to select the same variables that the ones used for modeling, so we have to consider the number of occurences inside the natural range to get the same variable selection obtained to model the current range
        number_ocurrences = read.csv("./results/ocurrences/ocurrences_per_species.csv", header=TRUE)
            #generated from "/home/dftortosa/diego_docs/science/phd/nicho_pinus/code/presences/ocurrences_resample/number_ocurrences_by_species.R"
        n_ocurrence = number_ocurrences[number_ocurrences$species==species,]$number_ocurrences

        #update the selected variables if the number of occurrences is too low
        if (n_ocurrence<length(selected_variables)*10){ #if there is not 10 ocurrences for each variable 
            
            #load the names of selected variables of low number ocurrences species and save
            load("./datos/finals/final_variables_low_number_ocurrence_species.rda")
                #generated in "/home/dftortosa/diego_docs/science/phd/nicho_pinus/code/variables/variable_selection/variable_selection_inside_clusters/variable_selection_inside_clusters_v2.R"
            selected_variables=final_variables_low_number_ocurrence_species_new[[species]]
        }


        ##load the selected variables
        #prepare paths
        path_variables=paste("./datos/finals/", selected_variables, ".asc", sep="")
            #The rasters stored in "datos/finals" were used as input in the species clustering and the selection of variables inside each cluster, so I assume these are the rasters used in the manuscript
                #generated in "/home/dftortosa/diego_docs/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/species_clustering_v2.R"
                #generated in "/home/dftortosa/diego_docs/science/phd/nicho_pinus/code/variables/variable_selection/variable_selection_inside_clusters/variable_selection_inside_clusters_v2.R"
        
        #load
        variables_stack = stack(path_variables)    
        
        #check
        if(sum(names(variables_stack) == selected_variables)!=length(selected_variables)){
            stop("ERROR! FALSE! WE HAVE A PROBLEM WHEN STACKING VARIABLES")
        }


        ##remove the areas inside the PA buffer
        #we have selected naturalized occurrences outside the PA buffer, see occurrences script for further details
        raster_pa_buffer = raster(paste("./results/pseudo_absences/", species, "_PA_buffer.asc", sep=""))

        #convert to polygon 
        polygon_raster_pa_buffer = rasterToPolygons(raster_pa_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon

        #IMPORTANT STEP: do an inverse mask of the environmental variable using the PA buffer as mask. This will give us all areas outside of the PA buffer.
        variables_stack_masked=mask(variables_stack, polygon_raster_pa_buffer, inverse=TRUE)
            #we do not need the area inside the PA buffer, as we have removed all occurrences occurring there. We do not want to use regions considered during the modeling.
            #also, the boyce index (evaluation metric for presence-only data) uses the whole raster of predictions, considering as absences everything that does not have a presence. It would consider as absences the PA buffer if we do not remove this area because we have removed the occurrences inside that buffer in these exsitu analyses.


        ##Extract values of variables 
        presences = read.table(paste("./results/global_test_phylo_current/exsitu_occurrences/", species, "/", species, "_final_presences.tsv", sep=""), sep="\t", header=TRUE) #read the data final data with presences and pseudoabsences
        variables = extract(variables_stack_masked, presences[, c("longitude", "latitude")]) #extract the value of the variable in these points
        data = cbind(presences, variables) #bind the presence data and the variable data in one data frame 

        #change the names of variables if the number of selected variables is 1 (it is to say, only 5 columns in data)
        if (ncol(data)<6){
            colnames(data)[5] = names(variables_stack_masked)
        }

        #check
        if(!identical(presences, data[, names(presences)])){
            stop("ERROR! FALSE! WE DO NOT HAVE THE CORRECT PRESENCES IN THE FINAL DATASET")
        }

        #write the final.data file
        write.table(data, gzfile(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/", species, "_env_predictors.tsv.gz", sep="")), sep="\t", row.names=FALSE) 
   


        ###predict suitability and evaluate
        ##load models
        #open folder to save the decompressed files
        system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/decompressed_models", sep=""))

        #decompress into the folder of the species
        system(paste("unzip -o ./results/models/models_", species, ".zip -d ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/decompressed_models/", sep=""))
            #the models have been obtained from Rafa PRO: "/Users/dsalazar/nicho_pinus/results/final_analyses/models". This is the path indicated in "/home/dftortosa/diego_docs/science/phd/nicho_pinus/code/models/fit_predict_current_suit/fit_predict_current_suit_v5.R"

        #load the RDA files of the models
        load(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/decompressed_models/", species, "_glm_model.rda", sep=""))
        load(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/decompressed_models/", species, "_gam_model.rda", sep=""))
        load(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/decompressed_models/", species, "_rf_model.rda", sep=""))


        ##extract thresholds
        #These thresholds were previously calculated using the evaluation data of each partition hence selecting the threshold maximizing TSS for each data partition.
        #decompress thresholds
        system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/thresholds/", sep=""))
        system(paste("unzip -o ./results/threshold/threshold_", species, ".zip -d ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/thresholds/", sep=""))
            #the thresholds have been obtained from Rafa PRO: "/Users/dsalazar/nicho_pinus/results/final_analyses/thresholds". This is the path indicated in "/home/dftortosa/diego_docs/science/phd/nicho_pinus/code/models/binarize_ensamble_current_suit/binarize_ensamble_current_suit_v5.R"

        #load the thresholds
        load(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/thresholds/", species, "_glm_threshold.rda", sep=""))
        load(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/thresholds/", species, "_gam_threshold.rda", sep=""))
        

        ##prepare folders to save outputs
        #to save prediction rasters
        system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/prediction_rasters/", sep=""))

        #open folder to save boyce index results
        system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/boyce_partitions/", sep=""))
        system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/plots/", sep=""))


        ##predict and evaluate
        #check
        if(length(glm_resample)!=12 | length(gam_resample)!=12 | length(rf_resample)!=12 ){
            stop("ERROR! FALSE! WE DO NOT HAVE 12 PARITIONS OF THE DATA, I.E., WE DO NOT HAVE 12 MODELS FOR EACH ALGORITHM TYPE")
        }
        
        #open empty lists
        glm_predict=list()
        gam_predict=list()
        rf_predict=list()
        glm_bin_predict=list()
        gam_bin_predict=list()
        rf_bin_predict=list()
        glm_eval=list()
        gam_eval=list()
        rf_eval=list()
        glm_eval_no_dup=list()
        gam_eval_no_dup=list()
        rf_eval_no_dup=list()

        #run a loop using the models trained on the 12 data partitions
        #k=1
        for(k in 1:length(glm_resample)){

            ##predict obtaining a continuous probability of presence for each cell outside of the PA_buffer
            glm_predict[[k]] = predict(variables_stack_masked, glm_resample[[k]], type="response")
            gam_predict[[k]] = predict(variables_stack_masked, gam_resample[[k]], type="response")
                #type="response" for both the predict function of GLM (stats::predict.glm) and GAM (gam::predict.Gam)
                    #The default is on the scale of the linear predictors (link); the alternative "response" is on the scale of the response variable.
                        #This is what we want, predicted values between 0 (absence) and 1 (presence)
                        #using the link option we get weird results, with values between 0 and -1500
                    #Nick uses type="response" for both type of models in the book. You can look for the term "type=“response”" in the PDF, and you will see.
            rf_predict[[k]] = predict(variables_stack_masked, rf_resample[[k]], type="prob", index=2)
                #Special case of prediction for RF. In order to get a continuous probability of presence for each cell instead of 0/1, we need to use type="prob", index=2.
                    #We modeled presences with RF considering them as a factor (as.factor(presence)), just like Nick did in his book in page 284. We model presences/absences as binomial using a classification tree as we wanted to have presence/absence to make the ensemble.
                    #In that page, you can see how he modeled the presence of vulpes as a factor using random forest, and then predict using type="prob" and selecting the second column of the output. This makes all the sense.
                    #You will notice that raster::predict has some arguments that can control output. 
                        #The fun argument lets you pass a custom predict function, superseding any existing predict method for the model object. 
                        #The index argument lets you define the column of a multi-column data.frame or matrix that is returned from a given predict method. In other words, if the output of the predict method of the model is a data.frame with several columns, index will let you to select the column you want to use as value for the cells of the output raster 
                    #this is the case for randomForest:::predict.randomForest. In a classification model, if type="prob" or "votes" a data.frame is returned, with n columns, representing each class. 
                        #For example, for binomial (0/1), you will get the probability of each observation to be 0 and to be 1, i.e., you get 2 columns.
                        #if you want to use probability of each cell (we are predicting across the whole raster, so each 10x10 cell is an observation) of having a presence (1), you have to use index=2. This is what we want.
                        #if you use index 1, you get the probability of a cell of having absences, which is the opposite to the previous plot.
                    #https://gis.stackexchange.com/a/339856
                    #https://rdrr.io/cran/raster/man/predict.html
                    #https://www.rdocumentation.org/packages/randomForest/versions/4.7-1.1/topics/predict.randomForest

            #NOTE about using presence as continuous (not factor) for glm and gam
                #this was the approach explained in the book for both GLM (126-128) and GAM (213). They use presence (0/1) without converting it to a factor.
                #I have also check with radiata in gam that setting the presence variable as a factor does not change a thing!

            #NOTE about the weights based on precision
                #the weights were not considered in the evaluation of the original models, they were only used during fitting. These weights were used for giving more importance to some occurrences in the fitting, i.e., to their environmental conditions. As we are not fitting the models, just predicting and evaluating, we do not need the weights.


            ##calculate the binary predictions for ensemble
            #convert each projection in binary using the threshold calculated with fit_eval_models (TSS)
            #glm
            glm_bin_predict[[k]] = glm_predict[[k]] #copy the raster 
            glm_bin_predict[[k]][glm_bin_predict[[k]]>=(glm_threshold[[k]][2,2]),] <- 1 #give 1 to the pixels with a predicted value higher or equal than the threshold
            glm_bin_predict[[k]][glm_bin_predict[[k]]<(glm_threshold[[k]][2,2]),] <- 0 #give 0 to the pixels with a predicted value lower than the threshold

            #gam
            gam_bin_predict[[k]] = gam_predict[[k]] #copy the raster 
            gam_bin_predict[[k]][gam_bin_predict[[k]]>=(gam_threshold[[k]][2,2]),] <- 1 #give 1 to the pixels with a predicted value higher or equal than the threshold
            gam_bin_predict[[k]][gam_bin_predict[[k]]<(gam_threshold[[k]][2,2]),] <- 0 #give 0 to the pixels with a predicted value lower than the threshold
            
            #RF
            #already binarized as we used a classification tree
            rf_bin_predict[[k]] = predict(variables_stack_masked, rf_resample[[k]], type="response")
                #calculate also predictions using the binary results of RF so we can use it later when ensembling predictions


            ##calculate the boyce index
            #load packages
            require(modEvA)
                #for boyce index. I am using this instead of ecospat.boyce because the latter gives me errors when loading the predicted raster
            require(terra)
                #to convert RasterLayers to SpatRasters, which is the required input format of the boyce index function

            #run the function removing or not duplicates (see below)
            jpeg(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/plots/", species, "_part_", k, "_boyce_index_plot.jpeg", sep=""), width=960, height=960, pointsize=24)
            par(mfcol=c(3,2))
            glm_boyce=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(glm_predict[[k]]), n.bins=NA, bin.width=0.1, res=100, method="spearman", rm.dup.classes=FALSE, rm.dup.points=FALSE, na.rm=TRUE, plot=TRUE, main=paste("GLM", sep=""))
            mtext(paste("radiata - part ", k, sep=""), side=3, line=-1.9, outer=TRUE, cex=1, font=2)
            gam_boyce=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(gam_predict[[k]]), n.bins=NA, bin.width=0.1, res=100, method="spearman", rm.dup.classes=FALSE, rm.dup.points=FALSE, na.rm=TRUE, plot=TRUE, main=paste("GAM", sep=""))
            rf_boyce=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(rf_predict[[k]]), n.bins=NA, bin.width=0.1, res=100, method="spearman", rm.dup.classes=FALSE, rm.dup.points=FALSE, na.rm=TRUE, plot=TRUE, main=paste("RF", sep=""))
            glm_boyce_no_dup=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(glm_predict[[k]]), n.bins=NA, bin.width=0.1, res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=TRUE, na.rm=TRUE, plot=TRUE, main=paste("GLM - no duplicates", sep=""))
            gam_boyce_no_dup=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(gam_predict[[k]]), n.bins=NA, bin.width=0.1, res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=TRUE, na.rm=TRUE, plot=TRUE, main=paste("GAM - no duplicates", sep=""))
            rf_boyce_no_dup=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(rf_predict[[k]]), n.bins=NA, bin.width=0.1, res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=TRUE, na.rm=TRUE, plot=TRUE, main=paste("RF - no duplicates", sep=""))
            dev.off()
                #this is calculating the boyce index, which is well suited to evaluate presence-only models. In particular, it has been used to evaluate the ability predict species invasions (Petitpierre et al., 2012), because absences are not reliable when investigating colonizing species. Therefore it suits perfectly our case here.
                    #Initially, the boyce index was implemented by making bins (usually 10) of the predicted probability of habitat suitability. For example, from 0 to 0.1, from 0.1 to 0.2, and so on.... Of course, you are classifying the different cells into these bins, if the cell has a predicted probability of 0.05, it will go into the first bin.
                    #In each bin, 
                        #we count the number of presences, and divide it by the total number of presences we have. In other words, we calculated the proportion of presences in that bin respect to the total of presences.
                        #we also count the total number of cells falling inside the bin, and divide by the number of cell in the whole study area. If the presences were randomly distributed, this would be the proportion of presences to be found. 
                    #As more cells are in bin, more presences should be there, independently of the suitability of the bin, if the presences are randomly distributed.
                    #In contrast, if presences tend to accumulate in high suitability areas (i.e., the model predict suitability correctly), we should find more and more proportion of presences as we move along the bins, i.e., higher suitability bins. 
                    #Therefore, we can calculate the correlation between the median suitability of the bins and the ratio between the proportion of presences observed and the proportion expected under random distribution. If the correlation is positive, i.e., more presences observed as suitability increases, the model works correctly. If the correlation is zero, then works just like chance, there is no pattern between presences and the suitability of the bins. If the correlation is negative, then presences tend to accumulate in low-suitable areas.
                    #this approach has the problem of selecting arbitrary bins, but it was solved by creating the continuous boyce
                        #define a fixed-width window, for example 0.1 and start moving the window.
                        #take all cells with probability between 0 and 0.1 and calculate again the ratio between proportion of presences observed and expected by random distribution.
                        #move the window 0.05 and repeat, taking cells with a probability between 0.05 and 0.15, and so on... until 100.
                        #you have 100 observed/expected ratio and median suitability for 100 windows. calculate the correlation and you have the boyce index without selecting a number of bins. 
                        #this is the boyce index we are using here.
                #The question about duplicated presences:
                    #the function is saying we have duplicated presences. I am sure the selected presences do not have the same coordinates and there no more than 3 occurrences within each 50x50km cells.
                    #Note, however, we can have several occurrences within a 10x10 cell, and this is what the function it looking at, as the resolution of the predictions is 10x10km.
                    #the question here is if these occurrences inside the same cell are statistically independent and if they bias the boyce index in any way.
                    #following our definition, these points are independent as they fulfill the requeriment of 3 points per 50x50km cell. We consider it as a good balance between covering the environmental conditions of the cell and avoiding pseudo-replication. Within 10x10km, you can have very distant forest patches that are independent, it is biologically plausible.
                    #how this affect boyce?
                        #imagine you have only 4 cells in a bin, with 2 presences being inside the same 10x10km cell. 4 is also the total of presences in the study and the total number of cells in the study is 8
                            #2/4=0.5 is the proportion of presences in the bin
                            #4/8=0.5 is the proportion of cells/pixels in the bin
                            #0.5/0.5=1
                        #now imagine that we remove the duplicates, so only 1 of the two presences remains
                            #1/3=0.33 is the proportion of presences in the bin
                            #4/8=0.5 is the proportion of cells/pixels in the bin
                            #0.33/0.5=0.66
                        #in the proportion of presences, the denominator is going to be equal or higher than the numerator because, of course, the total number of presences is going to be higher than the number of presences in a given bin. Therefore, if you remove a few duplicated presences from both numerator and denominator (if we remove a duplicated presence in a bin that presence should be removed from the total count also), the impact is going to be higher in the numerator, which is a smaller number.
                            #if you have 2/5 and remove 1 presence, it makes more impact move from 2 to 1 than from 5 to 4.
                            #larger reduction of the numerator makes the total smaller
                            #therefore, the proportion of presences is going to be smaller
                        #a lower proportion of presences observed makes the observed/expected ratio lower because the expected proportion does not change, as we have the same number of pixels.
                        #Therefore, removing duplicates reduces the ratio. That could decrease or increase the boyce index depending on where the presences are removed, in the lower or the higher bins.
                    #but what happen if would have more presences than pixels? this is the worst case scenario
                        #imagine a small study area with three 10x10km cells. The predicted suitability of each cell is different enough to be each one in a different bin, we have 3 bins. Each cell has 2 presences, making a total of 6 presences.
                        #in a given cell, the proportion of presences is 2/6=0.33, while the proportion of cells is 1/3=0.33.
                        #YES, even having double of presences, we still see the same proportion of presences in each cell than the expected under a random distribution.
                        #thanks to the fact that observed and predicted are both proportions respect to the total number of observed and predicted, respectively, having more presences does not jeopardize the calculation.
                        #you have your presences distributed with the same PROPORTION of the cells in the bins? then the ratio is 1 and does contribute to a correlation closer to zero between the ratio and bin value, i.e., boyce close to zero. 
                        #It does not matter if the absolute value of presences is higher than cells. The important thing is the proportion with respect to the total of presences.
                    #The remaining point then, is whether we want to consider the impact of these duplicated presences. It could increase or reduce the boyce index depending on where they are present. In high suitability bins it will increase boyce index, but in low bins it will decrease it.
                    #again, I think everything is reduced to decide whether we consider as independent points two occurrences inside the same 10x10km.
                        #yes, we could have two points close after the resampling, but this could also happen in the resampling within the natural range.
                        #there is any reason for this to be a problem now?
                        #Now, we want to check whether the regions considered suitable outside of the current range, tend to have more presences than expected by chance. 
                        #if we more presences accumulated in an area with high suitability, we should consider that, and the other way around, if we have occurrences accumulated in a low suitability area, we should take that into account?
                        #Given there is not inherent problem to have several occurrences in the same cell for boyce index, provided that these occurrences are independent, then we should keep them. These are independent according to our definition in the manuscript, and we should be consistent with that. Do this validation as close as possible to the original analyses, only changing stuff required for the particular characteristics of invasive analyses.
                #modEvA::Boyce is mostly based on ecospat.boyce which is the function used in Nick's book (page 268).
                    #obs: 
                        #a set of presence point coordinates. 
                        #In other words, a two-column matrix or data frame containing, respectively, the x (longitude) and y (latitude) coordinates of the presence points, in which case the 'obs' vector will be extracted with 'ptsrast2obspred'.
                    #pred: 
                        #a SpatRaster map with the predicted values for the entire model evaluation area, in which case the 'pred' vector will be extracted with 'ptsrast2obspred'
                            #this has to be "SpatRaster", not "RasterLayer"
                    #n.bins
                        #number of classes or bins (e.g. 10) in which to group the 'pred' values, or a vector with the bin thresholds (e.g., c(0.1., 0.4, 0.6...). If 'n.bins = NA' (the default), a moving window is used (see next parameters), so as to compute the "continuous Boyce index" (Hirzel et al. 2006).
                    #bin.width
                        #width of the moving window (if n.bins = NA), in the units of 'pred'. By default, it is 1/10th of the 'pred' range). Therefore, by default, the window always has a width of 0.1.
                        #Nick says in his book that we should select the smallest window possible. 
                        #I have checked several values for halepensis
                            #increasing it 7 times (e.g., 0.7) makes the index much higher (around 0.9), but the range of predicted suitability values considered is much lower. The first bin includes suitability values between 0 and 0.7 (median of 0.35). Therefore, we are starting to look at suitability 0.35 in the correlation, ending at median suitability of 0.647. Therefore, we are only looking at a small portion of the predicted suitability values.
                            #decreasing it 7 times (e.g., 0.04) decreases the boyce index a lot (0.11), but we cover the whole range of suitability values, as we start with a bin from 0 to 0.014, i.e., median 0.007. The last bin has a median of 0.98. So we cover the whole range of suitability in very detail.
                            #the default (0.1) gives an index of 0.38, while still covering the whole range of suitability. The first bin starts at 0 and ends at 0.1, i.e., median value of 0.050. The last bin has a median of 0.941.
                        #1/10 (0.1) is also the default in ecospat.boyce, being used by Nick in the book
                        #Importantly, herzel et al 2006 showed 
                            #that a window size of 0.1-0.2 makes the boyce index very correlated with AUC measured on the same presences but also adding true absences.
                            #they also say that increase to much the number of classes (i.e., lower width of the window) increases the variance among cross-validation partitions
                            #In their opinion, 0.1 is a good optimum. So we are going to use this
                            #Hirzel, A.H., G. Le Lay, V. Helfer, C. Randin & A. Guisan (2006) Evaluating the ability of habitat suitability models to predict species presences. Ecological Modelling 199: 142-152 
                    #res
                        #resolution of the moving window (if n.bins = NA). By default it is 100 focals, providing 100 moving bins. In other words, we move the window 100 times. If you increase the number of windows, each steps will be smaller in order to fit the larger number of windows within the same range of suitability values
                        #in the guisan paper, they just say the window is moved a small amount. The default of the function is 100, and the book of Nick also used 100, being the default also for ecospat.boyce.
                    #method
                        #argument to be passed to 'cor' indicating which correlation coefficient to use. The default is ''spearman'' as per Boyce et al. (2002), but ''pearson'' and ''kendall'' can also be used.
                        #the Guisan paper also uses Spearman
                    #rm.dup.classes
                        #rm.dup.classes: if 'TRUE' (as in 'ecospat::ecospat.boyce') and if there are different bins with the same predicted/expected ratio, only one of each is used to compute the correlation.
                        #ONLY removes contiguous bins, i.e., bins that are together (e.g., 0.7-0.8 and 0.71-0.81) and have the same P/E ratio.
                        #default here is FALSE, while in ecospat::ecospat.boyce is TRUE. Nick used the default, which is True.
                        #I am not sure it is a good idea to remove bins with the same P/E ratio. If two close bins have a similar value, that should be considered to strength the correlation between presence and suitability. These are legit values.... If they are in very different suitability bins, then this will destroy the correlation and that should be the case, because different suitability values have the same P/E ratio.
                        #using the default of the modEVa function
                        #IMPORTANT, DO NOT SET THIS AS TRUE
                            #Look first the boyce plots for radiata. You can see cases like in GAM, where there are only two P/E ratios above zero, being the rest zero. If you remove bins with the same value (i.e., 0), you indeed get a correlation, but it is completely artificial. Most of the range of suitability, there is no association with the P/E ratio.
                    #na.rm: 
                        #Logical value indicating if missing values should be removed from computations. The default is TRUE.
                        #TRUE, because we want to avoid cells with suitability values. There is nothing to do there.
                    #rm.dup.points: if 'TRUE' and if 'pred' is a SpatRaster and if there are repeated points within the same pixel, a maximum of one point per pixel is used to compute the presences. See examples in 'ptsrast2obspred'. The default is FALSE
                        #note that the function does not count duplicated presences for calculating the expected proportion. In other words, the expected proportion is calculated by just dividing the number of PIXELS within a given bin respect to the total number of PIXELS in the study area. These presences are used to calculate the proportion of presences in a bin respect to the total number of presences.
                #ptsrast2obspred
                    #This function is internally used by our boyce funtion. 
                    #It takes presence points or coordinates and a raster map of model predictions, and it returns a data frame with two columns containing, respectively, the observed (presence or no presence) and the predicted value for each pixel. 
                    #Duplicate points (i.e., points falling in the same pixel, WHETHER OR NOT THEY HAVE THE EXACT SAME COORDINATES) can be kept or removed.
                        #pts:
                            #two-column matrix or data frame containing their x (longitude) and y (latitude) coordinates
                        #rst:
                            #a one-layer 'SpatRaster' map of the model predictions, in the same CRS as 'pts'. If you have a raster map in another format, you can try to convert it with 'terra::rast()'
                        #rm.dup
                            #logical, whether repeated points within the same pixel should be removed. Default is FALSE
                        #na.rm
                            #logical, whether presence points with missing or non-finite values of 'rst' should be excluded from the output. The default is FALSE.
                    #output:
                        #outputs a data frame with one column containing the observed (1 for presence, 0 for absence) and another column containing the corresponding predicted values from 'rst'.
                        #In other words, you get a value per pixel and also for each additional presence included in the same pixel with other presence. Each row has a value for presence/absence.
  

            ##check ptsrast2obspred, which is internally used by boyce
            #run the function to calculate DF with presence/absence and predicted probability
            obspred_dup=ptsrast2obspred(rst=terra::rast(rf_predict[[k]]), pts=presences[,c("longitude", "latitude")], rm.dup=FALSE, na.rm=TRUE, verbosity=2)
            obspred_no_dup=ptsrast2obspred(rst=terra::rast(rf_predict[[k]]), pts=presences[,c("longitude", "latitude")], rm.dup=TRUE, na.rm=TRUE, verbosity=2)

            #the number of presences/absences after removing dups (second point in the same cell) should be the same than the total number of cell in the raster of predictions. If the cell has presence 1, if the cell has no presence then 0.
            if(nrow(obspred_no_dup) != length(which(!is.na(getValues(rf_predict[[k]]))))){
                stop(paste("ERROR! FALSE! WE ARE NOT CORRECTLY CALCULATING THE PRESENCES/ABSENCES FOR EVALUATION IN SPECIES ", species, sep=""))
            }
                #the additional presences in the same cell are included as additional rows in the data.frame
                #these will be also considered when counting the number of presences in each suitability bin and when counting the total number of presences, but this is ok, because these are legit and independent observations following our definition (see above).

            #check if we have too many points within the same 10x10 cells
            if((nrow(obspred_dup)-nrow(obspred_no_dup))>(nrow(presences)*0.2)){
                stop(paste("ERROR! FALSE! WE HAVE A PROBLEM WITH THE NUMBER OF PRESENCES WITIN THE SAME 10x10km CELL. MORE THAN 20% OF THE PRESENCES ARE PRESENT IN THE SAME CELL WITH OTHER PRESENCES FOR SPECIES ", species, sep=""))
            }
                #the difference in points between obspred_dup and obspred_no_dup is the number of presences being in 10x10km cells with other presences, as the latter was created by removing these presences, in contrast with the former.


            ##check we do not have a low number of pixels in any bin
            #They say the following: "In bins with overly small sample sizes, the comparison between median prediction and random expectation may not be meaningful, although these bins will equally contribute to the overall Boyce index. When there are bins with less than 30 values, a warning is emitted and their points are plotted in red, but mind that 30 is a largely arbitrary number. See the $bins$bin.N section of the console output, and use the 'bin.width' argument to enlarge the bins."
            #we have a veeery large study area, so the number of cells in each bin is always relatively high. We are going to check that just in case
            min_bin_n_glm=min(glm_boyce[[1]]$bin.N) 
            min_bin_n_gam=min(gam_boyce[[1]]$bin.N) 
            min_bin_n_rf=min(rf_boyce[[1]]$bin.N) 
                #calculate the min number of points across bins in each model. You can have different sample size between models because these are bins based on suitability, and the suitability is differently predicted by each model
            if((min_bin_n_glm<100) | (min_bin_n_gam<100) | (min_bin_n_rf<100)){
                stop(paste("ERROR! FALSE! WE HAVE AT LEAST 1 SUITABILITY BIN WITH LESS THAN 100 DATA.POINTS FOR SPECIES ", species, ". THIS CAN MAKE THAT THE COMPARISON BETWEEN MEDIAN PREDICTION AND RANDOM EXPECTATION MAY NOT BE MEANINGFUL", sep=""))
            }


            ##FINAL NOTES about modEvA boyce
                #The modEVA boyce function says (note at the end of the documentation) that we should only use this for habitat suitability predictions, not presence predictions
                    #These models (glm, gam....) predict presence probability, which (unless presences and absences are given different weights) incorporates the prevalence (proportion of presences) of the species in the modelled sample. So, predictions for restricted species are always generally low, while predictions for widespread species are always generally higher, regardless of the actual environmental quality.
                    #they specifically say that "This index is designed for evaluating predictions of habitat suitability, not presence probability (which also depends on the species' presence/absence ratio: rare species do not usually show high proportions of presences, even in highly suitable areas)."
                    #I guess what they are saying here is that if you have a species with a low number of presences, it is difficult for the model to detect its suitable areas as you do not have enough points of presence in these suitable areas, so the prediction are going to be always low, just because its low prevalence, but this is not a problem for us.
                    #But we have used different weights for presences and absences, we have also ensure the same ratio of PA/presences across species independently of their range size, and we have also used the whole known range of each species to create occurrences reducing the number if there were many occurrences together. Therefore, even small-range species have a high proportion of presences inside their range, where they, of course, have high suitability. Indeed, you can check the ensemble suitability maps for small-range species like canariensis, amamiana or culminicola, having all high certainty of suitability inside their natural range.
                    #Also, I have checked the predictions across the globe of canariensis, and it has different predicted suitable areas in the planet.
                    #I think we are good here.


            ##save the boyce output (with and without dups) in a list
            glm_eval[[k]]=glm_boyce
            gam_eval[[k]]=gam_boyce
            rf_eval[[k]]=rf_boyce
            glm_eval_no_dup[[k]]=glm_boyce_no_dup
            gam_eval_no_dup[[k]]=gam_boyce_no_dup
            rf_eval_no_dup[[k]]=rf_boyce_no_dup
        }

        #stop if we do not have 12 elements in each list
        if(length(glm_predict)!=12 | length(gam_predict)!=12 | length(rf_predict)!=12 | length(glm_bin_predict)!=12 | length(gam_bin_predict)!=12 | length(rf_bin_predict)!=12 | length(glm_eval)!=12 | length(gam_eval)!=12 | length(rf_eval)!=12 | length(glm_eval_no_dup)!=12 | length(gam_eval_no_dup)!=12 | length(rf_eval_no_dup)!=12){
            stop("ERROR! FALSE! WE HAVE NOT ANALYZED ALL PARTITIONS")
        }


        ##calculate a table with boyce index across models and partitions
        #add the boyce index of each combination as a row
        boyce_table=cbind.data.frame(
            1:length(glm_eval),
            cbind.data.frame(sapply(glm_eval, function(x){return(x[[2]])})), 
            cbind.data.frame(sapply(gam_eval, function(x){return(x[[2]])})),
            cbind.data.frame(sapply(rf_eval, function(x){return(x[[2]])})), 
            cbind.data.frame(sapply(glm_eval_no_dup, function(x){return(x[[2]])})), 
            cbind.data.frame(sapply(gam_eval_no_dup, function(x){return(x[[2]])})),
            cbind.data.frame(sapply(rf_eval_no_dup, function(x){return(x[[2]])})))
        names(boyce_table)=c("partition", "glm_eval", "gam_eval", "rf_eval", "glm_eval_no_dup", "gam_eval_no_dup", "rf_eval_no_dup")
            #for each partition, get the second element which is the boyce index, obtaining a vector with all of them. 
            #Then, convert to a column of DF. 
            #repeate for all models and combine the different columns in a DF along with the number of the partition.

        #add the confidence interval as a new row
        #calculate the different percentiles per column
        low_interval_boyce=apply(X=boyce_table[,which(colnames(boyce_table) != "partition")], MARGIN=2, FUN=quantile, probs=0.025)
        median_boyce=apply(X=boyce_table[,which(colnames(boyce_table) != "partition")], MARGIN=2, FUN=quantile, probs=0.5)
        high_interval_boyce=apply(X=boyce_table[,which(colnames(boyce_table) != "partition")], MARGIN=2, FUN=quantile, probs=0.975)
            #calculate the median of each column (excluding the first one with the number of the partition)
            #margin lets you select if you want to apply the function across columns (2) or rows (1)

        #add them as new rows
        boyce_table=rbind.data.frame(
            boyce_table, 
            c("percentile_2.5", low_interval_boyce),
            c("percentile_50", median_boyce),
            c("percentile_97.5", high_interval_boyce))
            #add the percentiles with their names as new rows


        ##make stacks
        predictions_glm=stack(glm_predict)
        predictions_gam=stack(gam_predict)
        predictions_rf =stack(rf_predict)
        predictions_bin_glm =stack(glm_bin_predict)   
        predictions_bin_gam =stack(gam_bin_predict)    
        predictions_bin_rf =stack(rf_bin_predict)   


        ##prepare ensemble
        #stack all binary predictions 
        binary_predictions = stack(predictions_bin_glm, predictions_bin_gam, predictions_bin_rf)   

        #calculate the percentage of models for which a pixel is suitable
        ensamble_predictions_bin = calc(binary_predictions, function(x) (sum(x)*100)/nlayers(binary_predictions))


        ##save the results
        #save predictions
        writeRaster(predictions_glm, filename=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/prediction_rasters/", species, "_predictions_glm", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        writeRaster(predictions_gam, filename=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/prediction_rasters/", species, "_predictions_gam", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        writeRaster(predictions_rf, filename=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/prediction_rasters/", species, "_predictions_rf", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        writeRaster(predictions_bin_glm, filename=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/prediction_rasters/", species, "_predictions_bin_glm", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        writeRaster(predictions_bin_gam, filename=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/prediction_rasters/", species, "_predictions_bin_gam", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        writeRaster(predictions_bin_rf, filename=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/prediction_rasters/", species, "_predictions_bin_rf", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        writeRaster(ensamble_predictions_bin, filename=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/prediction_rasters/", species, "_ensemble", sep=""), options="COMPRESS=LZW", overwrite=TRUE)

        #save evaluations
        save(glm_eval, file=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/boyce_partitions/", species, "_glm_boyce.rda", sep=""))
        save(gam_eval, file=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/boyce_partitions/", species, "_gam_boyce.rda", sep=""))
        save(rf_eval, file=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/boyce_partitions/", species, "_rf_boyce.rda", sep=""))
        save(glm_eval_no_dup, file=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/boyce_partitions/", species, "_glm_boyce_no_dup.rda", sep=""))
        save(gam_eval_no_dup, file=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/boyce_partitions/", species, "_gam_boyce_no_dup.rda", sep=""))
        save(rf_eval_no_dup, file=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/boyce_partitions/", species, "_rf_boyce_no_dup.rda", sep=""))

        #save median boyce
        write.table(boyce_table, gzfile(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/", species, "_boyce_table.tsv.gz", sep="")), sep="\t", col.names=TRUE, row.names=FALSE)

        #compress and remove files
        system(paste("rm -rf ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/decompressed_models/", sep=""))
        system(paste("rm -rf ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/thresholds/", sep=""))
        system(paste("
            cd ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/; 
            tar -zcvf ", species, "_prediction_rasters.tar.gz ./prediction_rasters; 
            rm -rf ./prediction_rasters", sep=""))
            #https://unix.stackexchange.com/a/93158
    }
}

#run it for one species
#predict_eval_no_phylo("radiata")




###########################################
##### PREDICT AND EVALUATE WITH PHYLO #####
###########################################

#species="radiata"
predict_eval_phylo = function(species){

    #open folder
    system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_phylo/", species, "/", sep=""))

    #check the outputs of the previous step exists
    pred_raster_check=file.exists(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/", species, "_prediction_rasters.tar.gz", sep=""))
    boyce_table_check=file.exists(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/", species, "_boyce_table.tsv.gz", sep=""))

    #if the outputs exists, do operations
    if(pred_raster_check & boyce_table_check){
        
        ##perform the phylogenetic correction
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

        #load bio4 and bio17
        bio4_raw = raster("./datos/finals/bio4.asc")
        bio17_raw = raster("./datos/finals/bio17.asc")

        #load ancestral reconstruction
        final_anc_ou_bio4 = read.table("./results/phylo_reconstruction/final_recons/final_anc_ou_bio4.csv", sep=",", header=TRUE)
        final_anc_ou_bio17 = read.table( "./results/phylo_reconstruction/final_recons/final_anc_ou_bio17.csv", sep=",", header=TRUE)
        final_anc_bm_bio4 = read.table("./results/phylo_reconstruction/final_recons/final_anc_bm_bio4.csv", sep=",", header=TRUE)
        final_anc_bm_bio17 = read.table( "./results/phylo_reconstruction/final_recons/final_anc_bm_bio17.csv", sep=",", header=TRUE)

        #remove the areas inside the PA buffer
        #we have selected naturalized occurrences outside the PA buffer, see occurrences script for further details
        raster_pa_buffer = raster(paste("./results/pseudo_absences/", species, "_PA_buffer.asc", sep=""))

        #convert to polygon 
        polygon_raster_pa_buffer = rasterToPolygons(raster_pa_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon

        #IMPORTANT STEP: do an inverse mask of the environmental variable using the PA buffer as mask. This will give us all areas outside of the PA buffer.
        bio4=mask(bio4_raw, polygon_raster_pa_buffer, inverse=TRUE)
        bio17=mask(bio17_raw, polygon_raster_pa_buffer, inverse=TRUE)
            #we do not need the area inside the PA buffer, as we have removed all occurrences occurring there. We do not want to use regions considered during the modeling.
            #also, the boyce index (evaluation metric for presence-only data) uses the whole raster of predictions, considering as absences everything that does not have a presence. It would consider as absences the PA buffer if we do not remove this area because we have removed the occurrences inside that buffer in these exsitu analyses.

        #load ensamble of binary suitability
        system(paste("tar --directory ./results/global_test_phylo_current/predict_eval_phylo/", species, " -zxvf ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/", species, "_prediction_rasters.tar.gz ./prediction_rasters/nigra_ensemble.grd ./prediction_rasters/nigra_ensemble.gri", sep=""))
            #--directory has to go first

        ensamble_suitability = raster(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters/", species, "_ensemble.grd", sep=""))



        ##you could apply the correction to 25-75 and to the whole area, and compare the boyce index, we are already comparing with random!!
        ##we can also compare phylo scaling or not

    }

    #load the phylogenetic range for the species
    #calculate the probability for both variables in the whole world except for PA buffer and combine to get a phylo probability
    #you have two options to combine the predicted probability of a model (e.g., glm) and the phylo probability, being both from 0 to 1
        #make an average of the two probability predictions as Nick in page 284 of the book
        #manually calculate the boyce index adding new presences that fall within the bin of probability of the phylo-map
            #you could calculate the "obs" object, i.e., presences- absences, and the "pred" object manually using "ptsrast2obspred", see boyce function docs.
            #if you do this, be very careful with the calculation. the expected ratio is calculated just with the number of pixels in the bin divided by the total number of pixels of the study area. It does not count repeated presences in this proportion, these are considered just in the predicted ratio.
    #calculate the boyce index with the new probability
    #compare the average boyce index with and without phylo and then check cases with important differences

    #make ensemble using a phylo threshold the same used in the manuscript, 0.75 I think.
}

#took old version of _check_altitudinal_sampling.png in results/occurrences folders for radiata

#check the thing about glm non binary!
    #mira pagina 284 Nivk's book, he seems to do as we



##########################
##### RUN EVERYTHING #####
##########################

#run all the functions
#better this way, so if a species is too slow it will not get the other species stuck in the first step, as the other ones will continue running
master_processor=function(species){

    #occurrences preparation
    n_points_before_resampling=exsitu_occurrences(species)
    
    #predict and evaluate without phylo
    predict_eval_no_phylo(species)

    #predict and evaluate with phylo
    predict_eval_phylo(species)

    #return the number of points before resampling
    return(n_points_before_resampling)
}


#run in parallel across species with naturalized occurrences
unique(naturalized_occurrences$species)

n_points_before_resampling=master_processor


#create table across species with boyce


#check how many naturalized presences inside PA buffer across species in "n_points_before_resampling". Use this to answer comment 8 of review about buffer size too big.

######################
##### NEXT STEPS #####
######################
##STEPS

    #evaluate again, but this time considering also as suitable regions fully within the phylogenetic range. 
        #evaluation with presence only!!!
            #This refined Boyce evaluation approach was, for instance, used in a study assessing our ability to project niche models into new territories to predict species invasions (Petitpierre et al., 2012), because absences are not reliable when investi- gating colonizing species.
                #Nick's book

    #we can do this for all species, then compare the median TSS with and without phylo, also look for species where phylo improves a lot

#maybe separate NA, EU and AUS? maybe phylo works different depending on the breath of climatic conditions of the continent?
    #Models were calibrated with 70% of the data in the native range (training dataset) and then projected onto the remaining 30% of the native range (evaluation dataset), the reciprocal Holarctic invaded range (EU or NA) and the Australian range (AU). For each technique, presences and pseudo-absences used to calibrate the model were weighted such as to ensure neutral (0.5) prevalence. The procedure was replicated 10 times, with random training and evaluation datasets, such that we obtained 40 models (10 replicates x 4 techniques) per species
    #We used the Boyce index (51, 52) to assess model performance. Contrary to common evaluation measures for presence-absence data, like AUC (53) or TSS (54) which present problems for presence-only data or for comparisons between species when the extents of calibration area differ (see 49, 55, 56), the Boyce index only requires presences and measures how much model predictions differ from random distribution of the observed presences across the prediction gradients (51). It is thus the most appropriate metric in the case of presence-only models. It is continuous and varies between -1 and +1. Positive values indicate a model which present predictions are consistent with the distribution of presences in the evaluation dataset, values close to zero mean that the model is not different from a random model, negative values indicate counter predictions, i.e., predicting poor quality areas where presences are more frequent (52). We calculated the Boyce index for each of the 40 models per species and averaged the values to obtain a final evaluation for the 30% left-out data in the native range (Bnat), in the invaded range in Holarctic (i.e. EU or NA, B inv ) and in AU when data were available (Bau
    #Climatic Niche Shifts Are Rare Among Terrestrial Plant Invaders
        #from the suple
