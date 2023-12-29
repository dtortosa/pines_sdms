#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
    #https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
    #https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



#########################################################
####### loop that cleans ocurrences and save them for modelling ####### #########################################################



###################################################
##### DIFFERENCES RESPECT TO PREVIOUS VERSION #####
###################################################

#Respect to version 1:



#######################
##### PREPARATION #####
#######################

##extract the in line arguments
#get the arguments passed
args=commandArgs(trailingOnly=TRUE)
    #https://stackoverflow.com/questions/14980322/is-it-possible-to-specify-command-line-parameters-to-r-script-in-rstudio

#extract each specific one
#species="radiata"; exsitu="yes"
species=args[1]
exsitu=args[2]
    #species for the species to analyze
    #exsitu for including or not occurrences inside or outside of the distribution


##create folders
system(paste("mkdir -p ./results/pipeline_2024/", species, "/01_occurrences", sep=""))




########################
##### BEGIN SCRIPT #####
########################

#Librerias
require(raster) #for work with rasters
require(sf) #for creating the buffer and the centroids of the cells without gbif points

#load list of species
list_species = read.table("./code/presences/species.txt", sep="\t", header=T)

#extract epithet from species list
#i=1
epithet_species_list = NULL
for(i in 1:nrow(list_species)){

    #selected species
    selected_species = as.vector(list_species[i,])

    #extract epithet
    epithet_species_list = append(epithet_species_list, strsplit(selected_species, split=" ")[[1]][2])
}
summary(!is.na(epithet_species_list))
#check
if(FALSE){
    require(tidyverse)
    paste("Pinus", epithet_species_list, sep=" ") == str_trim(as.vector(list_species[,1])) #the seperated epithet more Pinus are equal to list_species? We use str_trim from tidyverse to remove end spaces in each element of list_species
}#it is in false because loading tidyverse lead to load several packages that have a function name "extract", and this gives problems with the extract function of raster. If you want to check run these lines manually

#remove tecunumanii, jaliscana y discolor. These species are not used for Niche paper. The two first because we should downloadad gbif data now, so we would mix gbif data form 2016 and 2019. The third one was impossible to differentiate from P. cembroides
epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("tecunumanii", "jaliscana", "discolor"))]
#check
!c("tecunumanii", "jaliscana", "discolor") %in% epithet_species_list

##load elev for stratifies sampling
elev_high_resolution = raster("./datos/topografia/alt_wls_instra/alt/w001001.adf")

#agregate cells of elevation to reach the same resolution of bioclim variables (see notebook_pine_distribution_models.tex for further information). 
if(FALSE){
    elev_low_resolution = aggregate(elev_high_resolution, fact = 10, fun = mean)
    writeRaster(elev_low_resolution, "./datos/topografia/elev_low_resolution.asc", format="ascii", overwrite=TRUE)
}

#load elevation at low resolution
elev_low_resolution = raster("./datos/topografia/elev_low_resolution.asc")#we will calculate the percentile of altitude inside each cell of 50x50 inside the buffer, so we don't need great resolution. Only 10x10 is enough (0.08333334, 0.08333334)
res(elev_high_resolution)
res(elev_low_resolution)[1] #Now, it has a similar resolution of moisture bioclimvariables
elev = elev_low_resolution 

#Load a one bioclim and one soil variable 
bio1 = raster("./datos/finals/bio1.asc")
clay = raster("./datos/finals/clay.asc")
environment = clay*bio1 #Multiply both variables for obtaining a raster with all NAs. We will use them to ensure that all presences have environmental data, and specific we want to ensure that 0.5 fall in areas with value of environmental variables. We use a soil variable because there is bioclim data for some water bodies inside the continents, and we don't want create low precision points in these areas. 

#NOTE: In general, I have used extract to obtain any information about the points: elevation, environmental data and number of cell. Extract only consider that a point falls inside a cell if its center is inside that cell. It is important consider this. 

#set the seed to have reproducible results
set.seed(5683)

#Cargamos area distribucion
distribution = raster(paste("./datos/MAPS/p_", species, "_01.img", sep="")) #select the path of the distribution file of the corresponding species
    #new notes
        #the original path was "paste("/home/dftortosa/diego_docs/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(species, "01.img", sep="_") ,sep="_")"
        #we have moved the data to nicho_pinus
            #cp -r /home/dftortosa/diego_docs/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS /home/dftortosa/diego_docs/science/phd/nicho_pinus/datos/

#create a polygon of the dsitribution
polygon = rasterToPolygons(distribution, fun=function(x){x==1}) #convertimos en poligono the raster using the cell with values=1

#create a buffer
polygon_buffer = terra::buffer(polygon, width=100000, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) 
    #original notes when using rgeos in 2016
        #aumentamos el area del poligono en 2 celdas (width=1), según nos indico el test de las especies con datos de euforgen. byid determining if the function should be applied across subgeometries (TRUE) or the entire object (FALSE).  
    #new notes
        #width
            #we are now using terra::buffer. In this function, the unit of width is meter if the input has a longitude/latitude CRS, which is our case. Therefore, width=1 is no longer 2 cells wide. 
            #given that our cells are 50x50km, thus 2 cells wide would be 50*2=100km, i.e., 100000 meters. This should be width.
            #https://rspatial.github.io/terra/reference/buffer.html
        #byid and id were set to FALSE originally, but these arguments are no longer included in terra::buffer.
        #The rest of arguments are included in terra::buffer
        #Do not worry about the differences, we are going to check the buffer is identical to the previous one calculated with rgeos in 2016


##create a raster from buffer+bianca area converting the buffer+bianca polygon in a raster (rasterize) 
#We will use the cells of this raster to know the number of points inside a cell, it will be our grid of distribution. Because we will obtaine the centroid of this cells (points of atlas without gbif ocurrence), thus we have to calculate cells with more than three points using this grid, not the bioclim grid, with the purpose of that the resampling of points is made in a same resolution than the creation of centroids without gbif points. 
raster_buffer = raster() #create a empty raster
extent(raster_buffer) = extent(distribution) #give it a extention similar to distribution area of pinus halepensis polygon
res(raster_buffer) = res(distribution) #give it a resolution similar to he biancas maps because the polygon tha we will rasterize was obtained from this maps (I can`t extract resolution from polygon_buffer). In this way, we maintenain the same resolution. 
raster_buffer  = rasterize(polygon_buffer,raster_buffer) #convert the polygon of the halepensis buffer in a raster

#Drop marine areas of the buffer
distribution #take a raster of halepensis distribution (without buffer) for obtaining coast limits
raster_buffer = distribution*raster_buffer #multiply both raster for creating a raster with NAs as much as these rasters, it is to say, we want a raster with NA in pixel where some of these raster have a NA. The zones without cells in distribution are places of sea only, if we would have used a environmental variable, terrestrial cells without data of the variable also would be computed as sea cells. 
raster_buffer[!is.na(raster_buffer)] <- 1 #give value of 1 to all cells without NA, in this way include cell of distribution raster with 1 and cell outside the distribution raster but inside the raster buffer (with 0).  

#compare with the original buffer created using rgeos in 2016
original_raster_buffer = raster(paste("./results/ocurrences/", species, "_distribution_buffer", ".asc", sep=""))
compareRaster(original_raster_buffer, raster_buffer, extent=TRUE, rowcol=TRUE, crs=TRUE, res=TRUE, orig=TRUE, rotation=TRUE, values=TRUE, tolerance=1e-9, stopiffalse=TRUE, showwarning=FALSE)
    #we check that everything is the same
    #tolerance goes from 0 to 0.5, being the default in raster (rasterOptions()) equal to 0.1. Therefore, we are using a very low tolerance by setting it at 1e-9.
    #stopiffalse=TRUE to stop the script in case they are not the same

#save raster of distribution buffer
writeRaster(raster_buffer, paste("./results/pipeline_2024/", species, "/01_occurrences/ocurrences_", species, "_distribution_buffer.asc", sep=""), format="ascii", overwrite=TRUE)

#create the new polygon without sea areas
polygon_buffer = rasterToPolygons(raster_buffer, fun=function(x){x==1})


##crop the elevation and environmental raster using the buffer
#crop the raster of elevation to the distribution area using the buffer
elev_plots = crop(elev, polygon_buffer) #will be used for elevation sampling

#crop bio1 to distribution area of the species. 
environment_crop = crop(environment, polygon_buffer) #will be used for the plots and for ensure that low precision points fall in areas with values of environmental variables. 


##import presence data
#importa la tabla
presencia.completa<-read.csv(paste("./datos/raw_ocurrences/ocurrences/tree_species/final_species/Pinus", paste(species, "csv", sep="."), sep="_"), header=T, fill=TRUE, check.names=TRUE, stringsAsFactors=FALSE) #stringsAsFactors=FALSE sirve para procesar mejor la columna coordinatePrecision, que lleva mezclados números y caracteres

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



###env is limited to distribution, we should change that for exsitu, but how we calculate PAs? just buffers around points??? maybe you should create a polygon using the GBIF data points...


#pinaster or halepensis?

plot(environment)
points(presencia$longitude, presencia$latitude, cex=0.1)


#drop points without environmental variables.
environment_presences = extract(x=environment_crop, y=presencia[c("longitude", "latitude")]) #extract coge los puntos d epresencia, coge las varialbes, y para cada punto de presencai coge los valoers de las variables. Los puntos que quedan fuera del area de las variables se pone NA que luego quitaremos. 
is.data.frame(environment_presences)
#no es un data.frame, lo convertimos, 
environment_presences<-data.frame(environment_presences) 
is.data.frame(environment_presences) #este data.frame tiene una columna por variable, una fila por punto de presencia, y esa tabla la vamos a unir con la table de presencia (coordenadas pais, y todo lo demas)

