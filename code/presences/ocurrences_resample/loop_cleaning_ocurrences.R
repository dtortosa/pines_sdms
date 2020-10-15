#loop that cleans ocurrences and save them for modelling

###definimos el directorio de trabajo
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus")

#Librerias
require(raster) #for work with rasters
require(rgeos) #for creating the buffer and the centroids of the cells without gbif points

#load list of species
list_species = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/species.txt", sep="\t", header=T)

#extract epithet from species list
epithet_species_list = NULL
for(i in 1:nrow(list_species)){

    #selected species
    selected_species = as.vector(list_species[i,])

    #extract epithet
    epithet_species_list = append(epithet_species_list, strsplit(selected_species, split=" ")[[1]][2])
}
summary(is.na(epithet_species_list)) #all false

##load elev for stratifies sampling
#elev_high_resolution = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/topografia/alt_wls_instra/alt/w001001.adf")

#agregate cells of elevation to reach the same resolution of bioclim variables (see notebook_pine_distribution_models.tex for further information). 
#elev_low_resolution = aggregate(elev_high_resolution, fact = 10, fun = mean)
#writeRaster(elev_low_resolution, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/topografia/elev_low_resolution.asc", format="ascii", overwrite=TRUE)

#load elevation at low resolution
elev_low_resolution = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/topografia/elev_low_resolution.asc")
#res(elev_high_resolution)
#res(elev_low_resolution)[1] #Now, it has a similar resolution of moisture bioclimvariables
elev = elev_low_resolution 

#Load a one bioclim and one soil variable 
bio1 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio1.asc")
clay = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/clay.asc")
environment = clay*bio1 #Multiply both variables for obtaining a raster with all NAs. We will use them to ensure that all presences have environmental data, and specific we want to ensure that 0.5 fall in areas with value of environmental variables. We use a soil variable because there is bioclim data for some water bodies inside the continents, and we don't want create low precision points in these areas. 

#NOTE: In general, I have used extract to obtain any information about the points: elevation, environmental data and number of cell. Extract only consider that a point falls inside a cell if its center is inside that cell. It is important consider this. 

#make the loop for make cleaning of ocurrences for each species
for (i in epithet_species_list){

    #Cargamos area distribucion
    distribution = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(i, "01.img", sep="_") ,sep="_")) #select the path of the distribution file of the corresponding species

    #create a polygon of the dsitribution
    polygon = rasterToPolygons(distribution, fun=function(x){x==1}) #convertimos en poligono the raster using the cell with values=1

    #create a buffer
    polygon_buffer = gBuffer(polygon, byid=FALSE, width=1) #aumentamos el area del poligono en 2 celdas (width=1), según nos indico el test de las especies con datos de euforgen. byid determining if the function should be applied across subgeometries (TRUE) or the entire object (FALSE).  

    ###create a raster from buffer+bianca area converting the buffer+bianca polygon in a raster (rasterize) 
    #We will use the cells of this raster to know the number of points inside a cell, it will be our grid of distribution. Because we will obtaine the centroid of this cells (points of atlas without gbif ocurrence), thus we have to calculate cells with more than three points using this grid, not the bioclim grid, with the purpose of that the resampling of points is made in a same resolution than the creation of centroids without gbif points. 
    raster_buffer = raster() #create a empty raster
    extent(raster_buffer) = extent(distribution) #give it a extention similar to distribution area of pinus halepensis polygon
    res(raster_buffer) = res(distribution) #give it a resolution similar to he biancas maps because the polygon tha we will rasterize was obtained from this maps (I can`t extract resolution from polygon_buffer). In this way, we maintenain the same resolution. 
    raster_buffer  = rasterize(polygon_buffer,raster_buffer) #convert the polygon of the halepensis buffer in a raster

    #Drop marine areas of the buffer
    distribution #take a raster of halepensis distribution (without buffer) for obtaining coast limits
    raster_buffer = distribution*raster_buffer #multiply both raster for creating a raster with NAs as much as these rasters, it is to say, we want a raster with NA in pixel where some of these raster have a NA. The zones without cells in distribution are places of sea only, if we would have used a environmental variable, terrestrial cells without data of the variable also would be computed as sea cells. 
    raster_buffer[!is.na(raster_buffer)] <- 1 #give value of 1 to all cells without NA, in this way include cell of distribution raster with 1 and cell outside the distribution raster but inside the raster buffer (with 0).  

    #create the new polygon halepensis distribution without sea areas
    polygon_buffer = rasterToPolygons(raster_buffer, fun=function(x){x==1})

    #save raster of distribution buffer
    writeRaster(raster_buffer, paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences", paste(i, "distribution_buffer.asc", sep="_"), sep="/"), format="ascii", overwrite=TRUE)

    #crop the raster of elevation to the distribution area using the buffer
    elev_plots = crop(elev, polygon_buffer) #will be used for elevation sampling

    #crop bio1 to distribution area of the species. 
    environment_crop = crop(environment, polygon_buffer) #will be used for the plots and for ensure that low precision points fall in areas with values of environmental variables. 

    ################################################################
    ################################################################
    #PREPARACION DE LAS TABLAS DE DATOS PARA HACER LOS MODELOS
    ################################################################
    ################################################################

    #IMPORTA REGISTROS DE PRESENCIA
    #------------------------------
    #importa la tabla
    presencia.completa<-read.csv(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus", paste(i, "csv", sep="."), sep="_"), header=T, fill=TRUE, check.names=TRUE, stringsAsFactors=FALSE) #stringsAsFactors=FALSE sirve para procesar mejor la columna coordinatePrecision, que lleva mezclados números y caracteres

    #hacemos una copia sobre la que haremos los cambios
    presencia<-presencia.completa

    #Cambia el nombre de lat a latitud y lon a longitud
    n.column.lon = which(colnames(presencia)=="lon") 
    n.column.lat = which(colnames(presencia)=="lat") #look for the number of the column with longitude and latitude
    colnames(presencia)[n.column.lon] = "longitude" 
    colnames(presencia)[n.column.lat] = "latitude" #Use these number as index for select the correct column and change the name. 
    
    #drop rows with NA for longitude and latitude
    presencia<-presencia[!(is.na(presencia$longitude)),] 
    presencia<-presencia[!(is.na(presencia$latitude)),] 

    #drop points without environmental variables.
    environment_presences = extract(x=environment_crop, y=presencia[c("longitude", "latitude")]) #extract coge los puntos d epresencia, coge las varialbes, y para cada punto de presencai coge los valoers de las variables. Los puntos que quedan fuera del area de las variables se pone NA que luego quitaremos. 
    is.data.frame(environment_presences)
    #no es un data.frame, lo convertimos, 
    environment_presences<-data.frame(environment_presences) 
    is.data.frame(environment_presences) #este data.frame tiene una columna por variable, una fila por punto de presencia, y esa tabla la vamos a unir con la table de presencia (coordenadas pais, y todo lo demas)

    #lo unimos a las presencias
    presencia<-data.frame(presencia, environment_presences)
    nrow(presencia)
    rm(environment_presences)

    #quitamos los registros que tienen nulos en los valores de variables ambientales
    presencia<-presencia[!(is.na(presencia$environment_presences)), ] 
    nrow(presencia) #in this way we drop the points with environmental data (bioclim and soil)

    #create a coordinatePrecision variable if it is neccesary
    if (length(presencia$coordinatePrecision)==0){
        presencia$coordinatePrecision = rep(NA, nrow(presencia))            
    }

    #LIMPIEZA DE LA TABLA
    #####################
    ###Select ocurrences points inside the buffer
    distribution_buffer_values_of_points = extract(raster_buffer, presencia[,c("longitude","latitude")])  #extract the values of each ocurrence in the raster distribution buffer, the possibilities are two: 1 or NA. 1 inside distribution buffer, NA outside distribution buffer. We use the raster instead of the buffer because it is faster. 
    points_and_values_distribution_buffer = cbind(presencia[,c("longitude","latitude")], distribution_buffer_values_of_points) #bind these values with the corresponding presences points
    points_inside_distribution_buffer = which(!is.na(points_and_values_distribution_buffer$distribution_buffer_values_of_points)) #obtain rows in which values of the distribution buffer raster is 1, thus points outside the buffer
    points_outside_distribution_buffer = which(is.na(points_and_values_distribution_buffer$distribution_buffer_values_of_points)) #obtain rows in which values of the distribution buffer raster is NA, thus points inside the buffer
    
    #plot the result
    png(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences", paste(i, "ocurrences_inside_buffer.png", sep="_"), sep="/"), width=1500, height=800, pointsize=30)
    plot(environment_crop, col="gray80")
    plot(polygon_buffer, add=T)
    plot(polygon, col="green", add=T)
    points(points_and_values_distribution_buffer[points_outside_distribution_buffer,]$longitude, points_and_values_distribution_buffer[points_outside_distribution_buffer,]$latitude, cex=0.2, col="black")
    points(points_and_values_distribution_buffer[points_inside_distribution_buffer,]$longitude, points_and_values_distribution_buffer[points_inside_distribution_buffer,]$latitude, cex=0.2, col="red")
    legend(x="topright", legend=c("Bianca", "Buffer", "Points outside buffer", "Points inside buffer"), fill=c("green", "white", "black", "red"), cex=0.8)
    dev.off()

    #drop the ocurrences outside the buffer in the presencia data frame
    presencia = presencia[points_inside_distribution_buffer,]

    #Pedimos nombres de especies (ya realizada limpaiza taxonomicoa, solo lo dejamos por si aca)
    #------------------------------------
    unique(presencia$specificEpithet) 
    unique(presencia$scientificName) 
    unique(presencia$speciesKey) 
    for (k in unique(presencia$speciesKey)){
        print(unique(presencia[presencia$speciesKey==k,]$scientificName))
    } #loop for taking a look of  scientificName of each speciesKey


    #QUITAMOS LOS REGISTROS FÓSILES (SI LOS HAY)
    #-------------------------------------------
    #barplot(table(presencia$basisOfRecord)) #Las observaciones humanas y observaciones nos interesan, y tambien los pliegos de herbario (dentro del area de distribucion), pero si hay registros fosiles entonces hay que echar un ojo a ver si cuadranc on los registros actuales (mira codgio blas), no es nuestro caso. 
    #vemos la distribución de los especímenes fósiles y preservados. Para ello: 
    #plot(elev, col="gray80") #ploteo una variable
    #points(presencia$longitude,presencia$latitude, cex=0.6) #pongo los puntos de presencia
    #points(presencia[presencia$basisOfRecord == "FOSSIL_SPECIMEN", ]$longitude,presencia[presencia$basisOfRecord == "FOSSIL_SPECIMEN", ]$latitude, col="red", cex=1.2) #pongo encima los puntos fosiles. Poneindo la longitudeigtid y latitudeitis de unicamente las filas de regristors fosiles. 
    #nos quedamos con todos los registros que no son FOSSIL_SPECIMEN
    presencia<-presencia[presencia$basisOfRecord != "FOSSIL_SPECIMEN", ]
    #barplot(table(presencia$basisOfRecord)) #Ya no tenemos datos fosiles, hay observaciones, datos de especimenes preservados y datos desconocisos que son datos en los que no se han rellenado ese campo. 


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
    length(duplicados[duplicados==TRUE])

    ####################################################
    #####Create a variable of precision weight##########
    ####################################################
    ##Hasta abril de 2016 había solo una variable de accuracy of coordinates, pero gbif cambio en favor de los terminos de Darwin Core. Ahora, o tienes dato de coordinateUncertaintyInMeters ó de coordinatePrecision, pero nunca de los dos. 

    ###coordinatePrecision
    ##Hasta abril de 2016 este termino era el radio de error que tiene cada coordenada. La unidad es en metros, por tanto un valor de 25 indica un radio de 25x25 metros, y eso es lo que tenemos nosotros. 
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

    ###coordinateUncertaintyInMeters
    #Esta variable indica la distancia horizontal en metros desde la latitud decimal dada y la longitud decimal ada describrined the smallest circle containing the whole of the location. Se deja empty si la incertidumbre es desconocida, no puede ser estimada o no es aplicable (porque no hay coordenadas). Zero no es valor valido para este termino. Debe ser mayor y diferente de cero además de menor de 5000000 (5000 km). 
    #vemos el tipo de los datos
    #str(presencia$coordinateUncertaintyInMeters) 
    #vemos los valores de los datos
    #unique(presencia$coordinateUncertaintyInMeters) 
    #vemos su distribución
    #barplot(table(presencia$coordinateUncertaintyInMeters)) #points with precision lower than 4 km will receive a weight of 1, whilts point with coarser precision will receive a weight of 0.5.
    #It will not be used! becasue we have no evidence that this variable explcain anything about precision (see notebook and example sylvestris). 

    ###########################################
    #####Resampling ocurrences points##########
    ###########################################
    
    #extrac values of elevation for presences
    elevation_of_presences = extract(elev_plots, presencia[ , c("longitude","latitude")])

    ###Calculate the cell from halepensis map (buffer+bianca) in which each point is located
    #we create a raster with the number of cell as values
    index_raster <- raster_buffer #the raster is created from raster_buffer
    index_raster[] <- 1:ncell(raster_buffer) #the new raster take as value the number of cells of raster_buffer
    index_raster <- mask(index_raster, raster_buffer) #include the number of  raster_buffer cells in the new raster as values

    ###Calculate the cell from halepensis map (buffer+bianca) in which each point is located
    cell_id_presences= extract(index_raster, presencia[ , c("longitude","latitude")])

    #join data
    if (nrow(presencia)>0){
        id_presence=1:nrow(presencia[ ,c("longitude","latitude")]) #create a variable as id_presencia
    } else {
        id_presence = NULL
    }
    table_stratified_sample = cbind(cell_id_presences, id_presence, presencia[ ,c("longitude","latitude", "coordinatePrecision"),], elevation_of_presences)

    ### comproboations
    ##check that id of cell selected with presences has an 1 in raster buffer (inside buffer) and the rest cells have NA (outside of the buffer)
    unique(raster_buffer[cell_id_presences]) #all 1
    unique(raster_buffer[-cell_id_presences]) # all NA.

    ##check that cell indicated for each point is correct
    test = index_raster #create a raster from index raster (cell numbers)
    test[] <- NA #set all NA
    test[cell_id_presences] <- index_raster[cell_id_presences] #select only cells with presences and put their number

    #select the ids without repetition and ordered (less to more)
    cell_id_presences_order = sort(unique(cell_id_presences))

    #plot
    require(viridis)    
    colours_viridis = viridis(length(cell_id_presences_order)) #colour palette

    #open png
    png(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/", i, "_check_cell_number.png", sep=""), width=1100, height=800, pointsize=30)

    #plot raster with cell numbers of cells with presences but cropped with polygon_buffer to better visualization
    plot(crop(test, extent(polygon_buffer)))

    #for each id cell with presenece
    for(k in 1:length(cell_id_presences_order)){

        #select points of the [k] cell
        selected_cell = table_stratified_sample[which(table_stratified_sample$cell_id_presences == cell_id_presences_order[k]),]

        #add these points to the plot 
        points(selected_cell[,"longitude"], selected_cell[,"latitude"], col=colours_viridis[k])
    }

    #add legend
    legend("topright", legend=cell_id_presences_order, fill=colours_viridis, cex=0.5)
    dev.off() #you have to check that each color hf the first legend (inside plot) correspond to the points included in only ONE cell. Then, you have to cheack the the color of the cell correspond in second legend with the same number than in the second: Por ejemplo: La primera celda tiene solo PUNTOS de color muuuuy violetas, todos están ahi métidos, según la primera leyenda eso corresponde XXXXX. Ahora miras el color de la CELDA, y te vas a la segunda leyenda, su color debe corresponde con el número XXXXX ó cercano. Comprobado con canariensis y todo correcto. Además los puntos seleccionados al final son 3 por celdas en varias especies que he chequeado. 

    #####Resampling gbif points according to precision######
    #create a categorical variable of precision
    #we consider as high precision ocurrences those with a coordinatePrecision between 1 and 25, the rest of ocurrences are considered as low precision.  

    #create the empty variable for precision
    precision_weight = NULL

    if (nrow(table_stratified_sample)>0){
        for (k in 1:nrow(table_stratified_sample)){ #for each row of table_stratified_sample
            subset=table_stratified_sample[k,] #subset this row
            if (is.na(subset$coordinatePrecision)){ #if the point have cP=Na
                precision_weight = append(precision_weight, 0.5) #precision_weight = 0.5
            } else { #if not
                if (subset$coordinatePrecision>=1 && subset$coordinatePrecision<=25){ #if the point has a cP between 1 and 25, inclugin both extremes. 
                    precision_weight = append(precision_weight, 1) #precision_weight = 1
                } else { #if not
                    precision_weight = append(precision_weight, 0.5) #precision_weight = 0.5
                } 
            }   
        }
    } 

    #bind precision_weight with table_stratified_sample
    table_stratified_sample = cbind(table_stratified_sample, precision_weight)    


    #test 
    #summary(table_stratified_sample[is.na(table_stratified_sample$coordinatePrecision),]$precision_weight) #0.5
    #summary(table_stratified_sample[!is.na(table_stratified_sample$coordinatePrecision) & table_stratified_sample$coordinatePrecision<1 | !is.na(table_stratified_sample$coordinatePrecision) & table_stratified_sample$coordinatePrecision>25,]$precision_weight) #0.5
    #summary(table_stratified_sample[!is.na(table_stratified_sample$coordinatePrecision) & table_stratified_sample$coordinatePrecision>=1 & table_stratified_sample$coordinatePrecision<=25,]$precision_weight) #1

    #create the variable for final presences of gbif
    final.presences = NULL

    #make the resampling
    if (nrow(table_stratified_sample)>0){ #if there is gbif points with a high precision
        for (k in unique(table_stratified_sample$cell_id_presences)){ #for each cell of the distribution buffer
            subset=table_stratified_sample[table_stratified_sample$cell_id_presences== k,] #select the points inside the [i] cell. 
            subset_high_precision = subset[subset$precision_weight==1,] #subset the high precision points
            subset_high_precision_with_elevation = subset_high_precision[!is.na(subset_high_precision$elevation_of_presences),] #select the points with high precision and elevation data because of the elevation resampling
            subset_low_precision = subset[subset$precision_weight==0.5,] #subset the low precision points
            if (nrow(subset_high_precision_with_elevation)>3){ #if there is more than 3 high precision points
                elevations_in_cell=subset_high_precision_with_elevation$elevation_of_presences #select the elevation of these points
                p10=as.vector(quantile(elevations_in_cell, 0.10, na.rm=T))
                p50=as.vector(quantile(elevations_in_cell, 0.50, na.rm=T))
                p90=as.vector(quantile(elevations_in_cell, 0.90, na.rm=T)) #calculate the percentile 10, 50 and 90 of the elevation of all points inside the cell 
        
                closest_p10 = which(abs(elevations_in_cell-p10)==min(abs(elevations_in_cell-p10))) 
                closest_p50 = which(abs(elevations_in_cell-p50)==min(abs(elevations_in_cell-p50)))
                closest_p90 = which(abs(elevations_in_cell-p90)==min(abs(elevations_in_cell-p90))) #which elevation is closest to each percentile

                final.presences = rbind(
                final.presences,
                subset_high_precision[c(closest_p10[sample(1:length(closest_p10),1)],closest_p50[sample(1:length(closest_p50),1)], closest_p90[sample(1:length(closest_p90),1)]),]) #select the the points closest to each percentile, the point is selected randomly frome the number of points closest to the percentile.     
    
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

    ###tests
    if(nrow(table_stratified_sample)>0 && nrow(final.presences)){
        
        #only 3 points for each cell
        test = NULL
        for (k in unique(final.presences$cell_id_presences)){
            test = append(test, nrow(final.presences[final.presences$cell_id_presences==k,])<=3)
        }
        summary(test) 

        #each condition of the cleanign code works
        test_cells = data.frame()
        for (k in unique(table_stratified_sample$cell_id_presences)){
            subset=table_stratified_sample[table_stratified_sample$cell_id_presences== k,] 
            subset_high_precision = subset[subset$precision_weight==1,] 
            subset_high_precision_with_elevation = subset_high_precision[!is.na(subset_high_precision$elevation_of_presences),] 
            subset_low_precision = subset[subset$precision_weight==0.5,] 
            test_cells = rbind(test_cells, c(k, nrow(subset), nrow(subset_high_precision), nrow(subset_high_precision_with_elevation), nrow(subset_low_precision))) #for each cell calculate the different subsets, extract the number of rows (number of points) and then save in a data.frame.
        }
        names(test_cells) <- c("id_cell", "points", "high_precision", "high_p_elev", "low_precision") #change the name of the data.frame

        #Extract the cell that accomplish each condition
        test_cells[test_cells$high_p_elev>3,]$id_cell
        test_cells[test_cells$high_precision<=3 & test_cells$high_precision>=1 & test_cells$low_precision==0,]$id_cell
        test_cells[test_cells$high_precision<=3 & test_cells$high_precision>=1 & test_cells$low_precision>0 & test_cells$low_precision>3-test_cells$high_precision,]$id_cell
        test_cells[test_cells$high_precision<=3 & test_cells$high_precision>=1 & test_cells$low_precision>0 & test_cells$low_precision<=3-test_cells$high_precision,]$id_cell  
        test_cells[test_cells$high_precision==0 & test_cells$low_precision>0 &  test_cells$low_precision>3,]$id_cell
        test_cells[test_cells$high_precision==0 & test_cells$low_precision>0 &  test_cells$low_precision<3,]$id_cell

        #test if the different conditions cover all cells
        eso = c(test_cells[test_cells$high_p_elev>3,]$id_cell,
            test_cells[test_cells$high_precision<=3 & test_cells$high_precision>=1 & test_cells$low_precision==0,]$id_cell,
            test_cells[test_cells$high_precision<=3 & test_cells$high_precision>=1 & test_cells$low_precision>0 & test_cells$low_precision>3-test_cells$high_precision,]$id_cell,
            test_cells[test_cells$high_precision<=3 & test_cells$high_precision>=1 & test_cells$low_precision>0 & test_cells$low_precision<=3-test_cells$high_precision,]$id_cell,
            test_cells[test_cells$high_precision==0 & test_cells$low_precision>3,]$id_cell,
            test_cells[test_cells$high_precision==0 & test_cells$low_precision<=3 &  test_cells$low_precision>0,]$id_cell)
        length(eso) == nrow(test_cells)
    }


    #if there are high precision points final selected plot them into altitudinal plot to see the performance of the altitudinal sampling.
    if(nrow(final.presences[which(final.presences$precision_weight==1),]) > 1){
        
        #open png
        png(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/", i, "_check_altitudinal_sampling.png", sep=""), width=1100, height=800, pointsize=30)

        #plot elevation
        plot(elev_plots, main="Altitudinal sampling of high precision points")

        #add polygon buffer to see celss
        plot(polygon_buffer, add=TRUE)

        #add all points high precision points before the altitudinal sampling
        points(table_stratified_sample[which(table_stratified_sample$precision_weight==1),]$longitude, table_stratified_sample[which(table_stratified_sample$precision_weight==1),]$latitude, cex=0.5, col="white", lwd=0.15)
        
        #add high precision points after altitudinal sampling
        points(final.presences[which(final.presences$precision_weight==1),]$longitude, final.presences[which(final.presences$precision_weight==1),]$latitude, cex=0.5, col="red")
        #add legend
        legend("topright", legend=c("all points", "high preci. selected"), fill=c("white", "red"), cex=0.5)
        dev.off() #you have to check that the three high preicion points selected (red) in each cell are close to the maximum, mediam and minimun altitude considering altitude of ALL high precision points (white). Checked in pinus canriensis, it is ok, pay attention to Teide, it is very clear. 
    }

    #The resampling works great, the only little problem i see is in the case of cells with a lot of sea. In those cases can be selected 3 points close between them, but i don't think this could have great impact on the models. Only more correlation between very few points (only from buffer cells with a lot of sea).

    #############################################################
    ######Create points in cells without gbif points#############
    #############################################################

    ########Select cells without points in the halepensis polygon divided created above

    #conditional for the case in which does not exist points with high precision
    if (length(final.presences)>0){ #If final presences have columns, and thus it is not null. 
        polygons_with_points = extract(polygon_buffer, final.presences[,c("longitude","latitude")])
        cells_with_points = unique(na.omit(polygons_with_points$poly.ID)) #Extraemos los ID de las celdas with gbif points
    } else {
        cells_with_points = NULL
    }

    all_cells = getSpPPolygonsIDSlots(polygon_buffer) #extract ID of all cell of Bianca`s polygon
    all_cells = as.numeric(all_cells) #convert in numeric

    cells_without_points = all_cells[which(!(all_cells %in% cells_with_points))]
    cells_without_points = subset(all_cells, !(all_cells %in% cells_with_points)) #subset form all_cells, the ID_cells which are not included in cells_with_points, thus the cells without points. This is get by means of "!". 
    summary(cells_with_points %in% cells_without_points) #test if cells_with_points contains cells_without_points
    summary(cells_without_points %in% cells_with_points) #test if cells_without_points contains cells_with_points
    ##Subset from the POLYGON, cells without gbif points
    polygon_no_points = polygon_buffer[cells_without_points,]
    polygon_si_points = polygon_buffer[cells_with_points,] 

    ##create one point inside each cell of the buffer withouth gbif points
    #process bio1, we will use it for ensure all low precision created fall in areas with environmental variables. Remember that you buffer has a lower precision than environmental variables used in the models, thus the cells of variables are smaller than buffer cells, this produced that in some cases, the centroid of the buffer cell could fall in a area that it is sea in the variables maps (because of their higher resolution). We will cope with this. 

    bio1_polygon = rasterToPolygons(environment_crop) #convert environment_crop in a polygon. This variable will be our mask, thanks to it we avoid the creation of points inside water bodies. I need the environmental varaible at the same resolution than we will use in the models because we want to avoid points in the sea. With low resolution points can fall in areas that a high resolution are sea. 

    points_coords = data.frame() #create a empty data.frame
    for (j in 1:length(bio1_polygon@polygons)){ #for each cell of bio1_polygon
        points = spsample(x=bio1_polygon@polygons[[j]], n=1, type="random") #create a random point
        points_coords = rbind(points_coords, points@coords) #save the coordinates in points_coords 
    }

    #create a raster from no point polygon that has id_cell as value cell
    raster_buffer_no_points = raster() 
    extent(raster_buffer_no_points) = extent(distribution) 
    res(raster_buffer_no_points) = res(distribution) #same resolution and extent than raster distribution
    raster_buffer_no_points  = rasterize(polygon_no_points, raster_buffer_no_points) #with this we have created a raster with onle areas WITHOUT gbif points

    #we create a raster with the number of cell as values
    index_raster_low_p_points <- raster_buffer_no_points #the raster is created from raster_buffer
    index_raster_low_p_points[] <- 1:ncell(raster_buffer_no_points) #the new raster take as value the number of cells of raster_buffer
    index_raster_low_p_points <- mask(index_raster_low_p_points, raster_buffer_no_points) #include the number of  raster_buffer cells in the new raster as values

    ###Calculate the cell from halepensis map (buffer+bianca) in which each random low precision point is located
    cell_id_low_p_points= extract(index_raster_low_p_points, points_coords[, c("x", "y")])
    id_cell_low_p_points = cbind(cell_id_low_p_points, points_coords[, c("x", "y")]) #bind these IDs and the coordinates of the points. 
    id_cell_low_p_points = id_cell_low_p_points[!is.na(id_cell_low_p_points$cell_id_low_p_points),]

    final.points = data.frame() #create a final data frame for low precision points
    for (m in unique(id_cell_low_p_points$cell_id_low_p_points)){ #for each cell of the buffer
        subset = id_cell_low_p_points[id_cell_low_p_points$cell_id_low_p_points==m,] #select rows of points inside that cell (m)
        selected_points = subset[sample(1:nrow(subset),1),] #select from the all points inside [m] cell, only one. 
        final.points = rbind(final.points, selected_points) #bind that point the rest of points
    }

    #test if each cell have only 1 random low precision point
    test_low_p_points = NULL
    for (k in unique(final.points$cell_id_low_p_points)){
        test_low_p_points = append(test_low_p_points, nrow(final.points[final.points$cell_id_low_p_points==k,])==1)
    }
    summary(test_low_p_points) 

    #ploteamos 
    pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences", paste(i, "bianca_cells_with_without_gbif_points.pdf", sep="_"), sep="/"), width=24, height=24)
    par(mfcol=c(2,2))
    plot(environment_crop, col="gray80", main="Cells with gbif points")
    plot(polygon_si_points, add=T, col="red")
    points(final.presences$longitude, final.presences$latitude, cex=0.5, col="yellow")
    plot(environment_crop, col="gray80", main="Cells without gbif points")
    plot(polygon_no_points, add=T, col="blue")
    points(final.presences$longitude, final.presences$latitude, cex=0.5, col="yellow")
    plot(environment_crop, col="gray80", main="Cells with and without gbif points")
    plot(polygon_no_points, add=T, col="blue")
    plot(polygon_si_points, add=T, col="red")
    points(final.presences$longitude, final.presences$latitude, cex=0.5, col="yellow")
    plot(environment_crop, col="gray80", main="Cells with and without gbif points")
    plot(polygon_si_points, col="red", add=T) 
    plot(polygon_no_points, col="blue", add=T)   
    points(final.presences$longitude, final.presences$latitude, cex=0.5, col="yellow")
    points(final.points$x, final.points$y, col="black", cex=0.5)
    legend(x="topright", legend=c("Cells with gbif points", "Cells without gbif points", "gbif points", "low precision points"), fill=c("red", "blue", "yellow", "black"))
    dev.off()  

    #Join the new atlas points with the gbif points
    if (nrow(final.points)>0){ #if the number of points created from distribution buffer and presented outside water is higher than 0
        coordinates_points_without_ocurrence = final.points[,which(colnames(final.points) == "x" | colnames(final.points) == "y")] #select the coodinates of these points. 

        ##include all variable of final.presences in this data.frame
        #change to the same names than final.presences data.frame
        colnames(coordinates_points_without_ocurrence)[1] <- "longitude"
        colnames(coordinates_points_without_ocurrence)[2] <- "latitude"

        #create atlas data.frame
        atlas_points = coordinates_points_without_ocurrence
        atlas_points$precision_weight = rep(x=0.5, times=nrow(atlas_points))#add a variable of precision_weight. All points of the atlas have 0.5

        #join atlas points with gbif points
        ultimate_ocurrences = rbind(final.presences[, c("longitude", "latitude", "precision_weight")], atlas_points)
    } else { #if not 
        ultimate_ocurrences = final.presences[, c("longitude", "latitude", "precision_weight")] #the ultimate presences are only gbif points
    }

    #create a variable of final.presences
    ultimate_ocurrences$presence = rep(x=1, time=nrow(ultimate_ocurrences))

    #plot
    pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences", paste(i, "final_presences.pdf", sep="_"), sep="/"), width=12, height=12)
    plot(environment_crop, col="gray80")
    points(ultimate_ocurrences$longitude, ultimate_ocurrences$latitude,  cex=0.5, col="red")
    dev.off()

    #write the resulting data.frame in a csv. 
    write.csv(ultimate_ocurrences, paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences", paste(i, "final.presences.csv", sep="_"), sep="/"), row.names=FALSE)

    save.image(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/rdata", paste(i, "preparation_presences.RData", sep="_"), sep="/"))
}
