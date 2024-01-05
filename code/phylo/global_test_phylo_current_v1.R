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

#pre-defined functions
plot_sin=function(input){
    jpeg("./singularity_plot.jpeg", height=2000, width=2000, res=300)
    plot(input)
    dev.off()
}

#require packages
require(raster)
require(sf)

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

#NOTE: In general, I have used extract to obtain any information about the points: elevation, environmental data and number of cell. Extract only consider that a point falls inside a cell if its center is inside that cell. It is important consider this. 




###################################################
##### SELECT OCCURRENCES OUTSIDE DISTRIBUTION #####
###################################################

#species="radiata"
exsitu_occurrences=function(species){

    ##open folders
    system(paste("mkdir -p ./results/global_test_phylo_current/exsitu_occurrences/", species, sep=""))

    ##load the distribution buffer
    #read the raster
    distribution_buffer = raster(paste("results/ocurrences/", species, "_distribution_buffer.asc", sep=""))

    #convert to polygon
    polygon_distribution_buffer = rasterToPolygons(distribution_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon


    ##load the buffer used to create PAs.
    #This is out area of known absences due to edaphoclimatic conditions, i.e., no migration, etc... so we will consider as exsitu all occurrences outside this buffer
    #Note that we have to include EVERYTHING outside of the PA buffer, even if it is outside of the global distribution of pines. Pinus radiata is present in Australia...
    raster_pa_buffer = raster(paste("./results/pseudo_absences/", species, "_PA_buffer.asc", sep=""))  
    
    #convert to polygon 
    polygon_raster_pa_buffer = rasterToPolygons(raster_pa_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon


    ##obtain a polygon for all the world outside of the PA buffer
    #we will use the environmental raster, as it already includes areas with environmental data (both climatic and edaphic) across the globe and exclude water bodies
    #do an inverse mask of the environmental variable using the PA buffer as mask. This will give us all areas outside of the PA buffer.
    raster_outside_pa_buffer=mask(environment_var_low_res, polygon_raster_pa_buffer, inverse=TRUE)
        #we are masking using the high resolution climatic variable as input because if we use a coarser resolution (0.5), some areas inside the PA buffer are left in the coast for Pinus radiata.
        #it is ok, because we will then convert to low resolution to match the PA buffer we used to create occurrences

    #convert to 1 all cells with no NA, i.e., with edaphoclimatic data
    raster_outside_pa_buffer[which(!is.na(getValues(raster_outside_pa_buffer)))] = 1

    #create a polygon with all rows having edaphoclimatic data
    raster_outside_pa_buffer_polygon = rasterToPolygons(raster_outside_pa_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE)

    #create the same polygon but with cell borders inside just for checking occurrences resampling in the plots
    raster_outside_pa_buffer_polygon_non_dissolved = rasterToPolygons(raster_outside_pa_buffer, fun=function(x){x==1}, n=16, dissolve = FALSE)


    ##load occurrence data
    #read the table
    occurrences=read.csv(paste("./datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus_", species, ".csv", sep=""), header=TRUE, fill=TRUE, check.names=TRUE, stringsAsFactors=FALSE)
        #stringsAsFactors=FALSE sirve para procesar mejor la columna coordinatePrecision, que lleva mezclados números y caracteres

    #Cambia el nombre de lat a latitud y lon a longitud
    n.column.lon = which(colnames(occurrences)=="lon") 
    n.column.lat = which(colnames(occurrences)=="lat") 
        #look for the number of the column with longitude and latitude
    colnames(occurrences)[n.column.lon] = "longitude" 
    colnames(occurrences)[n.column.lat] = "latitude" #Use these number as index for select the correct column and change the name. 
    
    #drop rows with NA for longitude and latitude
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
    rm(environment_presences)

    #quitamos los registros que tienen nulos en los valores de variables ambientales
    occurrences<-occurrences[!(is.na(occurrences$environment_presences)), ] 
    nrow(occurrences)
        #in this way we drop the points without environmental data (bioclim and soil)

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

    #Pedimos nombres de especies (ya realizada limpaiza taxonomicoa, solo lo dejamos por si aca)
    #------------------------------------
    print(unique(presencia$specificEpithet))
    print(unique(presencia$scientificName))
    print(unique(presencia$speciesKey))
    for (k in unique(presencia$speciesKey)){
        print(unique(presencia[presencia$speciesKey==k,]$scientificName))
    } #loop for taking a look of scientificName of each speciesKey


    #QUITAMOS LOS REGISTROS FÓSILES (SI LOS HAY)
    #-------------------------------------------
    #barplot(table(presencia$basisOfRecord)) #Las observaciones humanas y observaciones nos interesan, y tambien los pliegos de herbario (dentro del area de distribucion), pero si hay registros fosiles entonces hay que echar un ojo a ver si cuadranc on los registros actuales (mira codgio blas), no es nuestro caso. 
    #vemos la distribución de los especímenes fósiles y preservados. Para ello: 
    #plot(elev, col="gray80") #ploteo una variable
    #points(presencia$longitude,presencia$latitude, cex=0.6) #pongo los puntos de presencia
    #points(presencia[presencia$basisOfRecord == "FOSSIL_SPECIMEN", ]$longitude,presencia[presencia$basisOfRecord == "FOSSIL_SPECIMEN", ]$latitude, col="red", cex=1.2) #pongo encima los puntos fosiles. Poneindo la longitudeigtid y latitudeitis de unicamente las filas de regristors fosiles. 
    #nos quedamos con todos los registros que no son FOSSIL_SPECIMEN o UNKNOWN (tipo de registro desconocido)
    presencia<-presencia[which(!presencia$basisOfRecord %in% c("FOSSIL_SPECIMEN")), ]
    #barplot(table(presencia$basisOfRecord)) #Ya no tenemos datos fosiles, hay observaciones, datos de especimenes preservados y datos desconocisos que son datos en los que no se han rellenado ese campo. 
    if(sum(presencia$basisOfRecord=="FOSSIL_SPECIMEN") != 0){
        stop("ERROR! FALSE! WE HAVE A PROBLEM REMOVING FOSSILS")
    }

    #LIMPIEZA DE DUPLICADOS EN LAS COORDENADAS
    #-----------------------------------------
    #buscamos registros duplicados en las coordenadas
    duplicados<-duplicated(presencia[ , c("latitude", "longitude")]) #la funcion duplicated del paqeute dismo busca registros duplicados para unas columnas que yo le diga. 
    #¿cuantos duplicados hay?
    length(duplicados[duplicados==TRUE])
    #selecciona los no duplicados (el símbolo ! significa "todo menos los duplicados")
    presencia<-presencia[!duplicados, ]
    #comprobamos que no han quedado duplicados (solo para cultivar vuestra fé)
    duplicados<-duplicated(presencia[ , c("latitude", "longitude")])
    if(sum(duplicados) != 0){
        stop("ERROR! FALSE! WE HAVE A PROBLEM REMOVING DUPLICATE OCCURRENCES")
    }


    ##Create a variable of precision weight
    #Hasta abril de 2016 había solo una variable de accuracy of coordinates, pero gbif cambio en favor de los terminos de Darwin Core. Ahora, o tienes dato de coordinateUncertaintyInMeters ó de coordinatePrecision, pero nunca de los dos. 

    #coordinatePrecision
    #Hasta abril de 2016 este termino era el radio de error que tiene cada coordenada. La unidad es en metros, por tanto un valor de 25 indica un radio de 25x25 metros, y eso es lo que tenemos nosotros. 
    ###Desde abril de 2016 ha cambiado. Una representación lineal de la precisión de las coordenadas dada en latitud y longitud decimal. Debe estar entre 0 y 1. Nosotros NO tenemos esto, por que la peña sigue subiendo los datos para esta variable como metros. 
    #for further details: 
        #http://lists.gbif.org/pipermail/api-users/2016-April/000300.html 
        #http://www.gbif.org/publishing-data/quality
        #http://tdwg.github.io/dwc/terms/index.htm#coordinateUncertaintyInMeters
    #vemos el tipo de los datos
    #str(presencia$coordinatePrecision)  
    #vemos los valores de los datos
    #unique(presencia$coordinatePrecision) 
    #vemos su distribución
    #barplot(table(presencia$coordinatePrecision)) #Mis variables tenian una resolucion de 4x4 km  y aqui tenemos valores por debajo (1,10, 25 = 1.85 km) y algunos por encima (0 =  18 km), los calculos están explicados en niche_questions.txt. Points with 1,10,25 will receive a weight of 1, whilst points with 0 will receive a weight of 0.5. 

    #coordinateUncertaintyInMeters
    #Esta variable indica la distancia horizontal en metros desde la latitud decimal dada y la longitud decimal ada describrined the smallest circle containing the whole of the location. Se deja empty si la incertidumbre es desconocida, no puede ser estimada o no es aplicable (porque no hay coordenadas). Zero no es valor valido para este termino. Debe ser mayor y diferente de cero además de menor de 5000000 (5000 km). 
    #vemos el tipo de los datos
    #str(presencia$coordinateUncertaintyInMeters) 
    #vemos los valores de los datos
    #unique(presencia$coordinateUncertaintyInMeters) 
    #vemos su distribución
    #barplot(table(presencia$coordinateUncertaintyInMeters)) #points with precision lower than 4 km will receive a weight of 1, whilts point with coarser precision will receive a weight of 0.5.
    #It will not be used! becasue we have no evidence that this variable explcain anything about precision (see notebook and example sylvestris).


    ##Resampling ocurrences points##########    
    #extrac values of elevation for presences
    elevation_of_presences = extract(elev, presencia[,c("longitude","latitude")])

    #Calculate the cell in which each point is located
    #we create a raster with the number of cell as values
    index_raster <- raster_outside_pa_buffer #the raster is created from raster_outside_pa_buffer
    index_raster[] <- 1:ncell(raster_outside_pa_buffer) #the new raster take as value the number of cells of raster_outside_pa_buffer
    index_raster <- mask(index_raster, raster_outside_pa_buffer) #remove now all cells that are not included the outside PA buffer (cells inside the PA buffer or water bodies)
    #Calculate the cell from halepensis map (buffer+bianca) in which each point is located
    cell_id_presences= extract(index_raster, presencia[,c("longitude","latitude")])
    #join data
    if (nrow(presencia)>0){
        id_presence=1:nrow(presencia[ ,c("longitude","latitude")]) #create a variable as id_presencia
    } else {
        id_presence = NULL
    }
    table_stratified_sample = cbind(cell_id_presences, id_presence, presencia[,c("longitude","latitude", "coordinatePrecision")], elevation_of_presences)
    
    #check 0
    #check that we have only one ID per occurrence
    if(length(unique(table_stratified_sample$id_presence))!=nrow(table_stratified_sample)){
        stop("ERROR! FALSE! WE HAVE A PROBLEM GIVING AN ID TO EACH OCCURRENCE DURING THE RESAMPLING")
    }

    #check 1
    #check that id of cell selected with presences has an 1 in raster buffer (inside buffer) and the rest cells have NA (outside of the buffer)
    if((unique(raster_outside_pa_buffer[cell_id_presences])!=1) | (!is.na(unique(raster_outside_pa_buffer[-cell_id_presences])))){
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


    ##Resampling gbif points according to precision######
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
        if(length(unique(subset_2))==1){
            checks[[2]]=unique(subset_2)==0.5
        } else { 
            checks[[2]]=FALSE
        }
    }
    if(length(subset_3)>0){
        if(length(unique(subset_3))==1){
            checks[[3]]=unique(subset_3)==1
        } else { 
            checks[[3]]=FALSE
        }
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
                                final.presences = rbind(final.presences, subset_low_precision[sample(x=1:nrow(subset_low_precision), size=3),]) #selec randomly three low precision points
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
        names(test_cells) <- c("id_cell", "points", "high_precision", "high_p_elev", "low_precision") #change the name of the data.frame
        test_cells=test_cells[-1,] #remove row with NAs
        #check whether the total number of presences inside a cell is the same than the sum of high and low precision presences
        if(sum(test_cells$points == test_cells$high_precision+test_cells$low_precision)!=nrow(test_cells)){
            stop("ERROR! FALSE! WE HAVE A PROBLEM CHECKING THE RESAMPLE")
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
        dev.off() #you have to check that the three high preicion points selected (red) in each cell are close to the maximum, mediam and minimun altitude considering altitude of ALL high precision points (white). Checked in pinus canriensis, it is ok, pay attention to Teide, it is very clear. 
    }
        #The resampling works great, the only little problem i see is in the case of cells with a lot of sea. In those cases can be selected 3 points close between them, but i don't think this could have great impact on the models. Only more correlation between very few points (only from buffer cells with a lot of sea).

    #visualize high and low precision occurrences
    png(paste("./results/global_test_phylo_current/exsitu_occurrences/", species, "/", species, "_check_high_low_occurrences.png", sep=""), width=2200, height=1600, res=300)

    #plot elevation
    if(species=="radiata"){
        #if radiata, plot the canary islands where it has been introduced and where we can see selection of points in more detail
        canariensis_buffer=raster(paste("results/ocurrences/canariensis_distribution_buffer", ".asc", sep=""))
        canariensis_buffer_polygon=rasterToPolygons(canariensis_buffer, fun=function(x){x==1}, dissolve=TRUE)
        plot(crop(elev, canariensis_buffer_polygon), main="Altitudinal sampling of high precision points")
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

    #write the resulting data.frame in a csv. 
    write.table(ultimate_ocurrences, paste("./results/global_test_phylo_current/exsitu_occurrences/", species, "/", species, "_final_presences.tsv", sep=""), sep="\t", row.names=FALSE)
}

#run it for one species
#exsitu_occurrences("halepensis")




########################################
##### EXTRACT ENVIRONMENTAL VALUES #####
########################################


##RUN only ONE time???

##CHECK WE ARE USING THE CORRECT INPUT FILES HERE

#species="radiata"
env_predictors = function(species) {

    #load the group of species according to cluster
    group_species = read.csv("./code/variables/variable_selection/species_clustering/_tables/complete_2_10_g_2.csv", header=TRUE)  
    
    #select the group of the corresponding species
    variables_cluster = group_species[which(group_species$species == species),]$groups   
    #select the number cluster of this species 
    rasters_list = list.files("./datos/finals", pattern=".asc", full.names=TRUE) #list the corresponding group of variables
    rasters_names = list.files("./datos/finals", pattern=".asc", full.names=FALSE) #list names of raster

    #loop for separate raster names from extension ".asc"
    names_variables = NULL
    #i=rasters_names[1]
    for (i in rasters_names){
        names_variables = append(names_variables, strsplit(i, split=".asc")[[1]])
    }
    variables_stack = stack(rasters_list) #stack them 
    names(variables_stack) = names_variables #give the names to layers of the stack

    #load names of selected variables  
    load("/Users/dsalazar/nicho_pinus/data/climate/list_selected_variables.rda") 

    #take selected variables from the list of variables
    selected_variables = ultimate_variables[[variables_cluster]]

    #calculate the number of ocurrences for each species 
    number_ocurrences = read.csv("/Users/dsalazar/nicho_pinus/data/ocurrences/ocurrences_per_species.csv", header=TRUE)
    n_ocurrence = number_ocurrences[number_ocurrences$species==species,]$number_ocurrences

    #load names of selected variables 
    if (n_ocurrence<length(selected_variables)*10){ #if there is not 10 ocurrences for each variable 
        
        #load the names of selected variables of low number ocurrences species
        load("./datos/finals/final_variables_low_number_ocurrence_species.rda")

        #select the selecte variables for low number ocurrences species
        variables_stack = variables_stack[[final_variables_low_number_ocurrence_species_new[[species]]]]    
    
    } else {

        variables_stack = variables_stack[[selected_variables]]
    
    }     
 
    ##Extract values of variables 
    presences = read.csv(paste("/Users/dsalazar/nicho_pinus/data/pseudo_absences", paste(species, "complete.presences.csv", sep="_"), sep="/")) #read the data final data with presences and pseudoabsences
    variables = extract(variables_stack, presences[, c("longitude", "latitude")]) #extract the value of the variable in these points
    data = cbind(presences, variables) #bind the presence data and the variable data in one data frame 

    #change the names of variables if the number of selected variables is 1 (it is to say, only 5 columns in data)
    if (ncol(data)<6){
        colnames(data)[5] = names(variables_stack)
    }

    #write the final.data file
    write.csv(data, paste("/Users/dsalazar/nicho_pinus/data/data_prepare_modelling", paste(species, "csv", sep="."), sep="/"), row.names=FALSE) 

    #name of species
    print(paste(species, "ended"))   
}




#took old version of _check_altitudinal_sampling.png in results/occurrences folders for radiata



#run all the functions
master_processor=function(species){

    #occurrences preparation
    exsitu_occurrences(species)

    #

}


##STEPS

    #once we have 3 occurrences per 50x50 cell (as we did for training), we can select the environmental variables used as predictors in the corresponding species, and then extract their values for the occurrences outside the PA buffer

    #we can then predict the probability of occurrence based on env variables and then check how well predict true presences, TSS....
        #high vs low precision was not considered in the evaluation of the models we did, so we could skip that for evaluation here

    #evaluate again, but this time considering also as suitable regions fully within the phylogenetic range. 
        #evaluation with presence only!!!
            #This refined Boyce evaluation approach was, for instance, used in a study assessing our ability to project niche models into new territories to predict species invasions (Petitpierre et al., 2012), because absences are not reliable when investi- gating colonizing species.
                #Nick's book

    #we can do this for all species, then compare the median TSS with and without phylo, also look for species where phylo improves a lot



####################################################
########### Loop preparing data #########################
####################################################


####################################################
########### Loop of fitting ########################
####################################################

#In each step we will create a dataset for training and evaluation and make both process, then we will repeat with other partition of training and evaluation. 
#We have to create a loop with for each species and before create empty lists (glm_species...) for include all the models for each species. In these lists we will sabe XX_resample lists


#function for fitting models of current habitat suitability
fit_eval_models = function(species){ #for the corresponding species 
    
    #required libraries
    #library(raster)
    #library(dismo)
    #library(gam) #OTRA LIBRERIA PARA GAM
    #library(randomForest) #RANDOM FOREST

    #begin the species
    print(paste("begin",species))

    #load cluster number for each species
    group_species = read.csv("/Users/dsalazar/nicho_pinus/data/climate/complete_2_10_g_2.csv", header=TRUE)[,c("species", "groups")] 

    #select the group of the corresponding species
    variables_cluster = group_species[which(group_species$species == species),]$groups 
    #Load data with presences and values of environmental variables
    data = read.csv(paste("/Users/dsalazar/nicho_pinus/data/data_prepare_modelling", paste(species, "csv", sep="."), sep="/"), header=TRUE)

    #select the number cluster of this species 
    rasters_list = list.files("/Users/dsalazar/nicho_pinus/data/climate/finals", pattern=".asc", full.names=TRUE) #list the corresponding group of variables
    rasters_names = list.files("/Users/dsalazar/nicho_pinus/data/climate/finals", pattern=".asc", full.names=FALSE) #list names of raster
    names_variables = NULL #loop for separate raster names from extension ".asc"
    for (i in rasters_names){
        names_variables = append(names_variables, strsplit(i, split=".asc")[[1]])
    }
    variables_stack = stack(rasters_list) #stack them 
    names(variables_stack) = names_variables #give the names to layers of the stack

    #load names of selected variables  
    load("/Users/dsalazar/nicho_pinus/data/climate/list_selected_variables.rda") 

    #take selected variables from the list of variables
    selected_variables = ultimate_variables[[variables_cluster]]

    #calculate the number of ocurrences for each species 
    number_ocurrences = read.csv("/Users/dsalazar/nicho_pinus/data/ocurrences/ocurrences_per_species.csv", header=TRUE)
    n_ocurrence = number_ocurrences[number_ocurrences$species==species,]$number_ocurrences

    if (n_ocurrence<length(selected_variables)*10){ #if there is not 10 ocurrences for each variable 
        
        #load the names of selected variables for los number ocurrences species
        load("/Users/dsalazar/nicho_pinus/data/climate/finals/final_variables_low_number_ocurrence_species.rda")

        #select the selecte variables for low number ocurrences species
        variables_stack = variables_stack[[final_variables_low_number_ocurrence_species_new[[species]]]]    
    
    } else {

        variables_stack = variables_stack[[selected_variables]]
    
    }       
 
    #load raster of PA buffer as a raster
    raster_PA_buffer = raster(paste("/Users/dsalazar/nicho_pinus/data/pa_buffers", paste(species, "PA_buffer.asc", sep="_"), sep="/"))

    #convert to polygon
    polygon_PA_buffer = rasterToPolygons(raster_PA_buffer, fun=function(x){x==1}, dissolve=TRUE)

    #crop variables using PA buffer to reduce the extension of the map (reduce running time)
    variables_stack = crop(variables_stack, polygon_PA_buffer)

    #mask variables using PA buffer to reduce all areas outside the buffer (we use crop before because mask not reduce the extension of the map)
    variables_stack = mask(variables_stack, polygon_PA_buffer)#no problem with sea areas inside the buffer because presence/absence points were only selected from areas with soil/climate data.

    #lists for saving models 
    glm_resample = list() 
    gam_resample = list() 
    rf_resample = list()

    #lists for saving predictions 
    glm_predict = list() 
    gam_predict = list() 
    rf_predict = list()

    #lists for saving predicted values on evaluation points 
    glm_evaluation_predict = list() 
    gam_evaluation_predict = list() 
    rf_evaluation_predict = list()

    #lists for saving evaluation of models 
    glm_evaluation = list()
    gam_evaluation = list()
    rf_evaluation = list()

    #list for save training data
    training_data = list()
    
    #set the seed for reproducibility
    #set.seed(56756)

    #For 12 times 
    for (k in 1:12){

        ######################################################
        ######separate data for training and evaluation#######
        ######################################################        
        presences = data[data$presence==1,] #subset the presences of the final data
        presences = presences[sample(1:nrow(presences), nrow(presences)),] #random sort of the presences for avoid that some aggregations of low or high precision occurrences could induce bias
        pseudo_absences = data[data$presence==0,] #subset the psuedo-absences
        pseudo_absences = pseudo_absences[sample(1:nrow(pseudo_absences), nrow(pseudo_absences)),] #random reorder of the pseudoabsences for following a similar approach than in occurrences. 
        index_evaluating_presences = sample(1:nrow(presences), round((nrow(presences)*30)/100)) #create index for selecting randmoly the 30% of presences for evaluation
        index_evaluating_absences = sample(1:nrow(pseudo_absences), round((nrow(pseudo_absences)*30)/100)) #create index for selecting randmoly the 70% of presences for training
        evaluate_presen = presences[index_evaluating_presences, ] #subset presences for evaluation
        evaluate_pseudo_absences = pseudo_absences[index_evaluating_absences,] #subset PAs for evaluation
        train_presen = presences[-index_evaluating_presences, ]
        train_pseudo_absences = pseudo_absences[-index_evaluating_absences,] #the same but wit the rest of ocurrences for obtaining the training data set. 
        training = rbind(train_presen, train_pseudo_absences) #bind presences and PAs for training
        evaluation = rbind(evaluate_presen, evaluate_pseudo_absences) #bind presences and PAs for evaluation

        ####changes in precision weight

        #calculate again the weight of PAs in training data
        correct_PA_weight = sum(training[which(training$presence==1),]$precision_weight)/nrow(training[which(training$presence==0),])

        #set the new weight for PAs according to the final number of PAs in the training data set
        training[which(training$presence==0),]$precision_weight <- correct_PA_weight        
        #if there are not presences with high precision (i.e. precision weight=1) we will change the precision weight of low precision points to 1 for glm and gam. In RF will be done in the sample size vector
        if(!1 %in% unique(training$precision_weight)){

            #create a new variable with precision weights
            training$weight_2_times_for_model <- training$precision_weight

            #change the precision weight of all presences (i.e. low precision presences) to 1
            training[which(training$presence==1 & training$weight_2_times_for_model==0.5),]$weight_2_times_for_model <- 1

            #calculate again the weight of PAs in training data
            correct_PA_weight_for_model = sum(training[which(training$presence==1),]$weight_2_times_for_model)/nrow(training[which(training$presence==0),])

            #set the new weight for PAs according to the final number of PAs in the training data set
            training[which(training$presence==0),]$weight_2_times_for_model <- correct_PA_weight_for_model  

            #weights for glm and gam
            weigth_glm_gam = training$weight_2_times_for_model
                #In that way, low precision points are now considered always at the maximum probability. The models retaining weight=0.5 for these points have exactly the same coefficients in glm and mostly in gam (so predictions should be similar), BUT the deviance values and AIC was different, so maybe this could affect to the stepwise regression. Indeed, in Pinus arizonica, selection of preditors in the stepwise regression changed. Because of this, I have amended the wight of low precision occurrences here.
                    #Now, the weight of PAs when only low precision presences are present are the same than in RAndom forest. For random forest we considered that the same number of PAs were selected relative to the number of low precision points when high precision points are not present.

            #check
            #summary(round(training$weight_2_times_for_model,4) == round(training$precision_weight*2,4))           
        } else{#if high precision points are present, then we will use the initial precision weight, with 1, 0.5...

            #weights for glm and gam
            weigth_glm_gam = training$precision_weight        

        }

        #extract the number of the columns with predictors
        col_numbers_predictors = which(colnames(training) %in% names(variables_stack))

        ###################################
        ############Fit models ############
        ###################################

        ################################
        #GENERALIZED LINEAR MODELS (GLM)
        ################################
        formula.regresion.poly = as.formula(paste("presence ~ poly(", paste(names(training)[col_numbers_predictors], collapse=", 2) + poly("), ", 2)", collapse=""))
        glm_pseudo_absence<-glm(formula.regresion.poly, family=binomial(link=logit), weights=weigth_glm_gam, data=training) #usa presencia como variable respuesta y como predictores solo dos variable. La familia e sbinomial porque nuestros datos de presencia tienen 1 y 0, y la curva que voy a usar es logística (polínomo de primer grado). Scamos los datos de la tabla presencia.pseudoausencia.entrenamiento, que tiene todas las presencias menos las de evaluacion, y todas las pseudoausencias. "quasibinomial" family bacause of the presences are weighted by precision_weight, and thus there is numbers with decimals, not only 0 and 1. We will use polynomial the level 2 with the purpose that the response curve can differ from the lineality, but no more than 2, because we want a simple model. 
        #The  warning "fitted probabilities numerically 0 or 1 occurred" indicates that the predicted probaiblities for presences/absences used in fittinh tienen valores entre 0 y 1, not only 0 and 1, this is caused by the weight. With quasibinominal the warning disappear, but the results are EXACTLY the same.

        #make stepwise
        glm_pseudo_absence = step(glm_pseudo_absence, direction="both", trace=0) #direction: the mode the mode of stepwise search, can be one of ‘"both"’, ‘"backward"’, or ‘"forward"’, with a default of ‘"both"’. 
            #trace: if positive, information is printed during the running of ‘step’. Larger values may give more detailed information. 

        #####################################
        #MODELOS ADITIVOS GENERALIZADOS (GAM)
        #####################################
        formula.gam<-as.formula(paste("presence ~ s(", paste(names(training)[col_numbers_predictors], collapse=",4) + s("), ",4)", collapse="")) #cada varaible está envuelta en una funcion que suaviza (smooth), gam traaja sobre valore suavizados de las variables. Por eso es s(). 
        if (is.element("package:mgcv", search())) detach("package:mgcv", force=TRUE) #make sure the mgcv package is not loaded to avoid conflicts with gam package.   
        gam_pseudo_absence<-gam(formula.gam, family=binomial(link=logit), data=training, weights=weigth_glm_gam)

        #make stepwise
        #genereate a scope for step.gam
        #Given a data.frame as an argument, generate a scope list for use in step.gam, each element of which gives the candidates for that term.
        gam_scope = gam.scope(frame=training[,c(which(colnames(training)=="presence"), col_numbers_predictors)], response=1, smoother="s", arg="df=4") #frame is the data: response variable and explicativa variables; response indicate the column of the response; smoother is the type of smooth that we use; arg: additional argument, like for example the freedom degree of smooth.

        #stepwise regression in GAM
        #In the case of amamiana, no one model is included in the output, maybe the reason is that if there is only one variable and it is not significant. We include the only variable selected for this species (we have few points for this species, so models will not be very good)
        if(species == "amamiana"){
            gam_pseudo_absence = gam_pseudo_absence
        } else {
            gam_pseudo_absence = step.Gam(gam_pseudo_absence, scope=gam_scope, direction="both", trace=0, data=training) #we include the list of terms in scope.             
        }    
        
        ##############
        #RANDOM FOREST
        ##############

        ###WEIGHT the random forest analysis###
        ##Strata variable
        #create a strata variable with diferent value for each precision_weight. This variables let us increase the probability of sample of data in relation to precision weight. Points with higher precision weight were more proably sampled than those with lower precision weight. We will give a letter to each point according to the precision weight (A = high precision presences; B = low precision presences; C = pseudoabsences).  
        if(nrow(training[training$precision_weight==1,])>0){ #if there are points with high precision (we have to consider the posibility that some species don't have high precision points)
            if(nrow(training[training$precision_weight==0.5,])>0){ #if there is low precision ocurrences

                training$strata = factor(NA, levels=c("A", "B", "C")) #create an empty factor 
                training[training$precision_weight==1,]$strata <- "A" #if precision_weight==1 strata is equal to "A"
                training[training$precision_weight==0.5,]$strata <- "B" #if precision_weight==0.5 strata is equal to "B"
                training[!(training$precision_weight==1 | training$precision_weight==0.5),]$strata <- "C" #if precision_weight is different of 1 and 0.5, strata is equal to "C"
            } else { #if there is not low precision ocurrences
                training$strata = factor(NA, levels=c("A", "C")) #create an empty factor 
                training[training$precision_weight==1,]$strata <- "A" #if precision_weight==0.5 strata is equal to "B"
                training[!(training$precision_weight==1 | training$precision_weight==0.5),]$strata <- "C" #if precision_weight is different of 1 and 0.5, strata is equal to "C"
            }
        } else { #if there is not high precision points 
            training$strata = factor(NA, levels=c("B", "C")) #create an empty factor 
            training[training$precision_weight==0.5,]$strata <- "B" #if precision_weight==0.5 strata is equal to "B"
            training[!(training$precision_weight==1 | training$precision_weight==0.5),]$strata <- "C" #if precision_weight is different of 1 and 0.5, strata is equal to "C"
        }#these lines are prepared to run with the initial precision_weight variable. If a species only have low precision presences, all of the will be selected and the same number of PAs will be also selected. In the following lines you can see as the B points are all selected when no high precision points are present. This is the same than than give a weight of 1 to low precision occurrences when no high precision points are available (as done in weight_2_times_for_model)

        #comprobation
        #table(training[training$strata=="A",]$precision_weight)
        #table(training[training$strata=="B",]$precision_weight)
        #table(training[training$strata=="C",]$precision_weight)

        #save the new factor as a vector 
        strata_vector = training$strata 

        #delete the factor from the data.frame         
        training = training[,-which(names(training)=="strata")] 

        ##N sample size
        #create the N sample size
        #this vector have the same length than unique(strata_vector), and indicate the number of rows that will be taken for each strata (category=high precision, low precision and PAs).
        #we need that the probaiblity of a high precision point to be take is 1 (take all of them), the probability for a low precision point the half (take half of low precision points), and the probability for a PA probability of high + low divided by total number of PAs. If only high or low precision points exist, all of them are taken (high and low for each case) and the same number for PAs. The occurrence type with the highest precision is taken in a probability of 1, and the PAs are taken with a lower probaiblity (sum of occurrence weight/number PAs)
 
        if (nrow(training[training$precision_weight==1,])>0){ #if there are points with high precision (A higher than 0):
            if(nrow(training[training$precision_weight==0.5,])>0){ #if there is low precision ocurrences

                #select all high precision points (probability  of a high precision point to be taken = 1)
                a = nrow(training[training$precision_weight==1,])

                #select half of low precision points (probability  of a high precision point to be taken = 0.5)
                b = nrow(training[training$precision_weight==0.5,])/2

                #select the same number of PAs than total occurrences selected (probability  of a PA to be taken = occurrence_weight/number total PAs)
                c = a + b 
                #c/nrow(training[which(training$presence==0),]) == unique(training[which(training$presence==0),]$precision_weight) #as internal checking, the number of PAs selected divided by the total number of PAs is equal to PA weight previously established for PA in training data set.  

                #bind all numbers in a unique vector of length equal to 3. 
                sampsize_vector = c("A"=as.vector(a), "B"=as.vector(b), "C"=as.vector(c)) 

                #Not rounded. Results are similar without rounding, rounding to floor or ceiling, so we use not rounded value to obtain exactly the same proportion than used in GLM and GAM (I have checked this for P. canariensis).

            } else { #if there is not low precision ocurrences. We take all high precision points (P(A) = 1) and the same number of PAs (P = P(A)/nº PAs)
                
                #select all high precision points (probability  of a high precision point to be taken = 1)
                a = nrow(training[training$precision_weight==1,])

                #select the same number of PAs than high precision occurrences selected (probability  of a PA to be taken = high precision occurrence weight/total PAs)
                c = a 
                #c/nrow(training[which(training$presence==0),]) == unique(training[which(training$presence==0),]$precision_weight) #as internal checking, the number of PAs selected divided by the total number of PAs is equal to PA weight previously established for PA in training data set.

                sampsize_vector = c("A"=as.vector(a), "C"=as.vector(c)) #bind all numbers in a unique vector of length equal to 2 (there is no ocurrences with high precision, A=0)  
            }
        } else { #if not, it is to say, there is not high precision ocurrences (A=0), We take all low precision points (P(B) = 1) and the same number of PAs (P = P(B)/nº PAs)

            #select all low precision points (probability  of a low precision point to be taken = 1)
            b = nrow(training[training$precision_weight==0.5,])

            #select the same number of PAs than low precision occurrences selected (probability of a PA point to be taken = low precision occurrence weight/total PAs)
            c = b
            #c/nrow(training[which(training$presence==0),]) == unique(training[which(training$presence==0),]$weight_2_times_for_model) #as internal checking, the number of PAs selected divided by the total number of PAs is equal to PA weight previously established for PA in training data set.
                #weight_2_times_for_model is the correct weight in this case, as this weight was calcualte when no precision points were avaiables so we change the weight of 0.5 to 1.

            sampsize_vector = c("B"=as.vector(b), "C"=as.vector(c)) #bind all numbers in a unique vector of length equal to 2 (there is no ocurrences with high precision, A=0)     
        }

        #write regresion formula
        formula.regresion<-as.formula(paste("as.factor(presence) ~ ", paste(names(training)[col_numbers_predictors], collapse="+"), collapse=""))

        #run the model
        rf_pseudo_absence<-randomForest(formula.regresion, data=training, importance=TRUE, ntree=500, strata=strata_vector, sampsize=sampsize_vector) #It is a classification tree, we select this becasue our variable is a factor of two levels. See "http://www.simafore.com/blog/bid/62482/2-main-differences-between-classification-and-regression-trees" for further information. 

        ###################################
        ############Save training data ####
        ###################################
        #save training data and sample size used for stratified random forest
        training_data[[k]]=list(sampsize_vector, training)
        names(training_data[[k]]) <- c("sampsize_rf", "training_data")

        ###################################
        ############Save models ###########
        ###################################
        glm_resample[[k]] = glm_pseudo_absence
        gam_resample[[k]] = gam_pseudo_absence
        rf_resample[[k]] = rf_pseudo_absence

        ###################################
        ############Predictions############
        ################################### 
        glm_predict[[k]] = predict(variables_stack, glm_resample[[k]], type="response")
        gam_predict[[k]] = predict(variables_stack, gam_resample[[k]], type="response")
        rf_predict[[k]] = predict(variables_stack, rf_resample[[k]], type="response") #the prediction of rf is binary because we are using classification forest. therefore, we don't have to binarize the predictions of this model. 


        ###################################
        ############Evaluation ############
        ################################### 

        #extract presence probability on evaluation points
        predict_glm = extract(glm_predict[[k]], evaluation[,c("longitude", "latitude")])
        predict_gam = extract(gam_predict[[k]], evaluation[,c("longitude", "latitude")])
        predict_rf = extract(rf_predict[[k]], evaluation[,c("longitude", "latitude")])


        #bind predictions and cordinates of evaluation points
        glm_evaluation_predict[[k]] = cbind(evaluation[,c("longitude", "latitude", "presence")], predict_glm)
        gam_evaluation_predict[[k]] = cbind(evaluation[,c("longitude", "latitude", "presence")], predict_gam)
        rf_evaluation_predict[[k]] = cbind(evaluation[,c("longitude", "latitude", "presence")], predict_rf)     

        #make the evaluation
        glm_evaluation[[k]] = evaluate(p=glm_evaluation_predict[[k]][glm_evaluation_predict[[k]]$presence==1, "predict_glm"], a=glm_evaluation_predict[[k]][glm_evaluation_predict[[k]]$presence==0, "predict_glm"])
        gam_evaluation[[k]] = evaluate(p=gam_evaluation_predict[[k]][gam_evaluation_predict[[k]]$presence==1, "predict_gam"], a=gam_evaluation_predict[[k]][gam_evaluation_predict[[k]]$presence==0, "predict_gam"])
        rf_evaluation[[k]] = evaluate(p=rf_evaluation_predict[[k]][rf_evaluation_predict[[k]]$presence==1, "predict_rf"], a=rf_evaluation_predict[[k]][rf_evaluation_predict[[k]]$presence==0, "predict_rf"])
    }

    #stack all continuous predictions 
    continuous_predictions_glm= stack(glm_predict)
    continuous_predictions_gam= stack(gam_predict)

    #strack predictions of random forest (binary)
    binary_predictions_rf = stack(rf_predict)   

    #save the training data
    save(training_data, file=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/models", paste(species, "training_data.rda", sep="_"), sep="/"))

    #save models 
    save(glm_resample, file=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/models", paste(species, "glm_model.rda", sep="_"), sep="/"))
    save(gam_resample, file=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/models", paste(species, "gam_model.rda", sep="_"), sep="/"))
    save(rf_resample, file=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/models", paste(species, "rf_model.rda", sep="_"), sep="/"))  

    #save the continuous predictions
    writeRaster(continuous_predictions_glm, filename=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/continuous_predictions/continuous_predictions_glm", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
    writeRaster(continuous_predictions_gam, filename=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/continuous_predictions/continuous_predictions_gam", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)

    #save the predictions of random forest (binary)
    writeRaster(binary_predictions_rf, filename=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_predictions/binary_predictions_rf", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)

    #save data for evaluation
    save(glm_evaluation_predict, file=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", paste(species, "glm_evaluation_data.rda", sep="_"), sep="/"))
    save(gam_evaluation_predict, file=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", paste(species, "gam_evaluation_data.rda", sep="_"), sep="/"))
    save(rf_evaluation_predict, file=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", paste(species, "rf_evaluation_data.rda", sep="_"), sep="/")) 

    #save evaluations
    save(glm_evaluation, file=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", paste(species, "glm_evaluation.rda", sep="_"), sep="/"))
    save(gam_evaluation, file=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", paste(species, "gam_evaluation.rda", sep="_"), sep="/"))
    save(rf_evaluation, file=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", paste(species, "rf_evaluation.rda", sep="_"), sep="/"))

    #name of species
    print(paste(species, "ended"))  
}


########Paralelize the process######
require(foreach)
require(doParallel) #for parallel


#load data about problem with weight of PAs
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


# set up cluster
clust <- makeCluster(12) 
registerDoParallel(clust)

#####################
###PREPARING DATA####
#####################
#preparing data for non-problematic species
foreach(i = epithet_species_list, .packages=c("raster", "dismo")) %dopar% { 
    prepar_data(species = i)
} 

###########################
###FITTING & EVALUATION####
###########################
#fitting and evaluation for non-problematic species
foreach(i = epithet_species_list, .packages=c("raster", "dismo", "gam", "randomForest")) %dopar% { 
    fit_eval_models(species = i)
} 


#stop the cluster 
stopCluster(clust)