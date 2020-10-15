#Code for making the selection of variables used for each species: We select variables for both high and low number of ocurrences species.  Selection with both ciriterias (ancient and new) is made here. Moreover here we check what variables have dramactic changes in their values.

###definimos el directorio de trabajo
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus")

#Librerias
require(raster)
require(rgdal) #for saving shapefiles

###list ocurrences
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

##Load variables
list_variables = list.files("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals", full.names=TRUE, pattern=".asc")
variables = stack(list_variables)

##Load the group variable and the ranks 
group_variables = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/_tables/complete_2_10_g_5.csv", header=TRUE)  
length(unique(group_variables$groups)) #we have 5 groups


######################################################
################Mask variables########################
######################################################
#ONLY ONE TIME. With this, we select only the distribution area + PA_buffer of all species of each cluster. We are not interested in studying correlation outside this area, because there are not ocurrences there. 

#loop
for (i in unique(group_variables$groups)){ #for each group

    #create a empty stack 
    distribution = stack() 

    #other loop inside the major loop
    for (k in group_variables[group_variables$groups==i,]$species){ #for each specie in the [i] group 

        #load the distribution of species
        distri = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(k, "01.img", sep="_"), sep="_")) 

        #load the raster of the PA buffer
        raster_PA_buffer = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences", paste(k, "PA_buffer.asc", sep="_"), sep="/"))

        #convert to a polygon
        polygon_PA_buffer = rasterToPolygons(raster_PA_buffer, fun=function(x){x==1}, dissolve=TRUE) #we change to polygon for converting again in raster but with the extent of distribution raster (global). Therefore, all raster will have the same extent and could be stacked

        #rasterize to the dimension of distribituion
        raster_buffer = raster() 
        extent(raster_buffer) = extent(distri) 
        res(raster_buffer) = res(distri) 
        raster_buffer  = rasterize(polygon_PA_buffer,raster_buffer) #create a raster from the PA buffer, with resolution and extension of distribution map
        raster_buffer[!is.na(raster_buffer)] <- 1 #areas without NA will have 1
        raster_buffer[is.na(raster_buffer)] <- 0 #areas with NA will have 0

        #save the raster in the empty stack
        distribution = stack(distribution, raster_buffer) 
    }

    sum_layers = sum(distribution>0) #sum all layers (all species group) of the stack to obtain a raster with a number in each cell that represent the presence the 1 or more PA-distribution-buffers. 
    variables_mask = rasterToPolygons(sum_layers, fun=function(x){x>0}, dissolve=TRUE) #create a polygon from the sum of layers using cells that don't have 0 (some buffer). 

    #save the polygon mask a shape file
    writeOGR(variables_mask, dsn = paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/mask_groups/mask_group", paste(i, "asc", sep="."), sep="_"), layer = 'poly', driver = "ESRI Shapefile", overwrite_layer=TRUE)
}

#read all of the masks and crop the variables
require(rgdal) #for saving shapefiles
masks = list()
variables_cropped = list()
for (i in unique(group_variables$groups)){
    masks[[i]] = readOGR(dsn=paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/mask_groups/mask_group", paste(i, "asc/poly.shp", sep="."), sep="_"), layer="poly")
    variables_cropped[[i]] = mask(variables, masks[[i]])
}

#save stacks
for (i in unique(group_variables$groups)){
    writeRaster(variables_cropped[[i]], filename=paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals_cropped/group", paste(i, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
}


#################################################
########Preparation of variables#################
#################################################

#load variables
variables_cropped = list()
for (i in unique(group_variables$groups)){
    variables_cropped[[i]] = stack(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals_cropped/group", paste(i, "tif", sep="."), sep="_"))
    names(variables_cropped[[i]]) = names(variables)
}

#copy variables_cropped
variables_clusters = variables_cropped

###All variable in each group have the same resolution and extent. If not, it had not been possible include all in a stack
for (i in unique(group_variables$groups)){
    print(res(variables_clusters[[i]]))
    print(extent(variables_clusters[[i]])) 
}

#RESOLUCION DE LAS VARIABLES
res.grados<-xres(variables_clusters[[1]])
res.minutos = res.grados*60 #Si 1 grado son 60 minutos, 0.08333333(res.grados) serán X minutos. Y salen 5 minutos, justo la resolucion de moisture index de bianca. 
#en km
res.km<-res.grados*111.19 #Si 1 grado son 111.19 km, 0.08333333(res.grados) serán X km y salen 9.265833 km. Es decir, cada celda de mi mapa mide en la realidad 9.265833 x 9.265833 km, pero en el ECUADOR!
res.km #The radious of the earth at the equator is 6378137.0 meters resulting in a circumference of 40075161.2 meters (2*pi*6378137.0; 2*pi*r). The equator is divided into 360 degrees of longitude, so each degree at the equator represents 111.32 km (40075161.2/360 divided by 1000 for kilometers). As one moves away from the equator towards a pole, however, one degree of longitude is multiplied by the cosine of the latitude (in radians), decreasing the distance, approaching zero at the pole. Therefore, if the maps is far away from the equator, the cells will be smaller. We have to multiplicate the length of the at equator (degrees * 111.19 km) by the cosine of latitude in radians. For example a cell of 10*10 at the equator, will be at 37 degrees of latitude (latitude of for example Spain):
library(NISTunits) #functon to convert degrees to radians ("https://stackoverflow.com/questions/32370485/r-convert-radians-to-degree-degree-to-radians")
cos(NISTdegTOradian(37))*res.km #7.4 km
cos(NISTdegTOradian(30))*res.km
cos(NISTdegTOradian(40))*res.km #30-40 is the range of latitude with the maximun diversity of Pines, thus most of our distribution cells have a size between 8km and 7km. See "https://gis.stackexchange.com/questions/142326/calculating-longitude-length-in-miles" and "https://en.wikipedia.org/wiki/Decimal_degrees"


#RESOLUCION DE LOS MAPAS DE DISTRIBUCION (sensu Critfield and Little)
list_distrib = list() #empty list for save species distribution maps 

#for each species
for(i in 1:length(epithet_species_list)){

    #select the [i] species
    selected_species = epithet_species_list[i]

    #create a raster with the distribution of the [i] species 
    distrib = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_", selected_species, "_01.img", sep=""))

    #save it in the list
    list_distrib[[i]] = distrib
}

#stack the rasters of the list
stack_distrib = stack(list_distrib)
nlayers(stack_distrib) == 113
xres(stack_distrib)[[1]]

#Resolucion de la distribucion de la especie
res.grados.distribution = xres(stack_distrib)[[1]] #res of 0.5 degrees in all raster of distribution. We stack in a raster stack with quick=FALSE, thus it has been checked that all raster introduced in the stack have the same resolution.
res.minutos.distribution = res.grados.distribution*60
res.km.distribution = res.grados.distribution*111.19 #Si 1 grado son 111.19 km, 0.5(res.grados) serán X km y salen 55 km. Es decir, cada celda de mi mapa mide en la realidad 55 x 5 km, pero en el ECUADOR!
res.km.distribution #The radious of the earth at the equator is 6378137.0 meters resulting in a circumference of 40075161.2 meters (2*pi*6378137.0; 2*pi*r). The equator is divided into 360 degrees of longitude, so each degree at the equator represents 111.32 km (40075161.2/360 divided by 1000 for kilometers). As one moves away from the equator towards a pole, however, one degree of longitude is multiplied by the cosine of the latitude (in radians), decreasing the distance, approaching zero at the pole. Therefore, if the maps is far away from the equator, the cells will be smaller. We have to multiplicate the length of the at equator (degrees * 111.19 km) by the cosine of latitude in radians. For example a cell of 50*50 at the equator, will be at 37 degrees of latitude (latitude of for example Spain):

library(NISTunits) #functon to convert degrees to radians ("https://stackoverflow.com/questions/32370485/r-convert-radians-to-degree-degree-to-radians")
cos(NISTdegTOradian(37))*res.km.distribution #44.5
cos(NISTdegTOradian(30))*res.km.distribution 
cos(NISTdegTOradian(40))*res.km.distribution #30-40 is the range of latitude with the maximun diversity of Pines, thus most of our distribution cells have a size between 48km and 42km. See "https://gis.stackexchange.com/questions/142326/calculating-longitude-length-in-miles" and "https://en.wikipedia.org/wiki/Decimal_degrees"

#Our variables has a resolution of 7.5x7.5 km, while distribution grid, which has been used for create points of low resolution, has 45x45 km. Thus, we have some points with a resolution coarser than variables, because of this, these point have a lower weight than high precision points.  



###############################################################################
##########Create a rank list for variables (D2) in ech cluster ################
###############################################################################
###create a list with the sum of ranks for each variable in each cluster and oredered from low to high. I will use this for select variable easier. 
results_ranks = list() #create a list

selected_var_cor_analysis = c(names(group_variables[,-which(colnames(group_variables)=="species" | colnames(group_variables)=="groups")])) #select the names of variables from the table with ranks for all variables and ranks. We will use this table for create the data frame with the sum of ranks. 

for (i in unique(group_variables$groups)){ #for each each cluster

    group = group_variables[group_variables$groups==i,-which(colnames(group_variables)=="species" | colnames(group_variables)=="groups")] #select the row if the cluster [i]

    #create empty vectors for the sum of each type of variable
    sum_ranks = NULL
    #make a loop of each type of variable
    for (k in group){ #for each temperature variable 
        sum_ranks = append(sum_ranks, sum(k))
    } #sum its ranking in all species of the cluster [i]
    sum_ranks=as.data.frame(sum_ranks) #convert to data.frame 
    sum_ranks = cbind(selected_var_cor_analysis, sum_ranks) #bind the data frame of rank_sums with the a vector with variables names
    names(sum_ranks)[1] = "Variables" 
    names(sum_ranks)[2] = "sum" #establish column names
    sum_ranks = sum_ranks[order(sum_ranks$sum, decreasing=FALSE),] #reorder the rows according to sum variable, with a increasing order (low values at first)

    results_ranks[[i]] = sum_ranks #bind all vectors (sums of each varible category) and include them as element number [i] of the previous created list. 
}

results_ranks #look at the result

################################################################
######ANÁLISIS DE CORRELACIÓN DE LAS VARIABLES PREDICTIVAS######
################################################################

###we will create a function to create variables tables and read them for each cluster. This function will be paralelized for clusters. 

#library
require(foreach) #for repeat a process several times 
require(doParallel) #for parallel

###create a function to write
write_fun = function(group){
    variables_table <-na.omit(as.data.frame(variables_clusters[[group]]))
    write.csv(variables_table, file=paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/final_variables/variables", paste("group", paste(group, "csv", sep="."), sep="_"), sep="_"), row.names=FALSE)
    rm(variables_table)
}

# set up cluster
clust <- makeCluster(3)
registerDoParallel(clust) 

# run
foreach(i = unique(group_variables$groups), .packages="raster") %dopar% { 
    write_fun(group = i)
} 

#stop the cluster 
stopCluster(clust) 


###create a function to read
read_fun = function(group){
    read.csv(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/final_variables/variables", paste("group", paste(group, "csv", sep="."), sep="_"), sep="_"), header=TRUE)
}

# set up cluster
clust <- makeCluster(3)
registerDoParallel(clust) 

# run
variables_table = foreach(i = unique(group_variables$groups)) %dopar% { 
    read_fun(group = i)
} 

#stop the cluster 
stopCluster(clust) 

#create a dendrogram with the variables of each cluster
variables_dendo = list() #create a list for save the dendogram of each cluster
for (i in unique(group_variables$groups)){ #for each cluster: 
    correlation_variables<-cor(variables_table[[i]]) #MATRIZ DE CORRELACIÓN
    variables_dist<-as.dist(abs(correlation_variables)) #MATRIZ DE DISTANCIAS ('ABS' = VALOR ABSOLUTO, PARA ELIMINAR CORRELACIONES NEGATIVAS). Nos da igual 0.9 ó -0.9, nos interesa saber entre que variables hay correlacion. Por eso usamos abs, y luego aplicamos as.dist para sacar una matriz de distancias. 
    variables_dendo[[i]]<-hclust(1-variables_dist) ##CLUSTER DE VARIABLES SEGÚN LA DISTANCIA (MENOR DISTANCIA = MAYOR CORRELACIÓN). Las varialbes con altos niveles de correlacion deben aparecer juntas, pero un dendo mide distancias, cuanto mas corta es la distnacia mas cerca aparecen en el dendograma. Pero nuestros valores de correlaione stan justo al reves, los más correlacionados tienen valores mas altos, y saldrian mas separados, por eso hacemos 1- las distancias. hclust crea el dendograma. 
}

#GRÁFICO DEL CLUSTER DE CORRELACIONES EXPORTADO A PDF
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/final_variables/dendogram.pdf", width=12, height=12, pointsize=20)
for (i in unique(group_variables$groups)){
    plot(variables_dendo[[i]])
}
dev.off()


#OBSERVAMOS EL PDF Y ELEGIMOS LAS VARIABLES CON LAS QUE QUEREMOS TRABAJAR
#El limite que establezco de correlacion es 0.4, todas aquellas varaibles correlacionadas en menos de un 0.4 y que por tanto se separan antes del 0.6 (que equivale a ese 0.4 de correlacion) las meto sin más. En el resto de casos, si la separaión se ha hecho despues, nos quedamos con solo UN representando del grupo. Esto hace que nos quitemos variables con un nivel de correlacion mayor del 0.4. ES UN NIVEL MUY ESTRICTO, POR ESO SEREMOS FLEXIBLES PARA CASOS CERCANOS AL LIMITE. In GENERAL, only 5 variables for cluster, 3 of climate (at least 1 and 1 of pp and temperature, respectively) and 2 of soil. If in a cluster there is more than 5 groups (and thus more than 5 variables), we select from the representative variable of the groups, the 5 with the highest D2. And controlling the ratio of soil and climate variables. See for more information: http://blog.minitab.com/blog/understanding-statistics/handling-multicollinearity-in-regression-analysis

variables_for_vif = list()

#group 1
results_ranks[[1]] #load sum of ranks of variables
group_variables[group_variables$groups==1,]$species #load the names of species
variables_for_vif[[1]] = c(
    "silt", #360
    "bio11", #64
    #"bio10", #231 
    "clay", #396
    #"bio18", #337
    "bio13", #462
    "bio12") #347
    #Tenemos 2 de suelo, 2 de pp y 1 de tª (3 clima), sobra 1 variable de clima. Quitamos bio10 porque junto con bio11 tiene un VIF mayor de 7, por tanto ambas parecen estar muy correlacionadas. Dejamos Bio11 porqque es la más explicativa (media del cuarto más frío). Luego quitamos bio18 porque tiene saltos muy abruptos a lo largo del espacio, y la cambiamos por bio13, bio16 sería la siguiente más explicativa de ese grupo (5), pero cuando se mete el VIF de bio 12 salta mucho. Así que teniendo en cuenta que bio12 explica más que bio13 y bio16, cogemos bio13 que no da problemas de VIF, añadimos bio13.
    plot(variables_clusters[[1]][["silt"]])
    plot(variables_clusters[[1]][["bio11"]])
    plot(variables_clusters[[1]][["bio10"]])
    plot(variables_clusters[[1]][["clay"]])
    plot(variables_clusters[[1]][["bio18"]])
    plot(variables_clusters[[1]][["bio13"]])
    plot(variables_clusters[[1]][["bio12"]]) 
    
    #selección anterior
    if(FALSE){
        #"bio3", #146
        "bio11", #64
        "clay", #396
        #"bio18", #337
        "bio12", 
        "bio13") #347
    }
        #surplus 1 variable, and I have two soil variables and 4 climate variable. Instead of drop the climatic variable with highest rank (bio12), we will drop bio3, and remain bio12, because the VIF of bio3 es 5.5. When we drop bio3 and include bio12, in any variable the VIF is higher than 5. This means that bio3 is explained by the rest of variables, but bio12 no, so is more interesting inlcude bio12, although this variable explain less variability, this variability is not explaine by the others variables. Moreover, we have to drop bio18 because has dramatic changes y vcalues of variables (see down), because of this we see thr variables of the group 5 (group of bio18) and the group 2 (bio3), groups withour variables selected. We selecte bio13 because is low correlated with the rest of variables, altougt explain less. 

#group 2
results_ranks[[2]] #load sum of ranks of variables
group_variables[group_variables$groups==2,]$species #load the names of species
variables_for_vif[[2]] = c(
    "carbon", #678
    #"bio10", #475
    "bio11", #175
    #"depth", #865
    "silt", #781  
    #"bio18", #497
    "bio13", #502
    "bio17") #315
    #Tenemos 3 varibles de suelo y 4 de clima, sobra una de suelo y una de clima. Ninguna de las de suelo tiene vif alto, así que quitamos la menos explicativa (depth). Bio10 y bio11 tienen un vif mayor de 5, ambas de pp, por lo que quitamos la menos explicativa de las dos para bajar el VIF. Quitamos ademas bio18 por los problemas de saltos bruscos, ponemos la siguiente del grupo 6, bio13. 
    plot(variables_clusters[[2]][["carbon"]])
    plot(variables_clusters[[2]][["bio10"]])
    plot(variables_clusters[[2]][["bio11"]])
    plot(variables_clusters[[2]][["depth"]])
    plot(variables_clusters[[2]][["silt"]])
    plot(variables_clusters[[2]][["bio18"]])
    plot(variables_clusters[[2]][["bio13"]]) 
    plot(variables_clusters[[2]][["bio17"]]) 

    #selección anterior
    if(FALSE){
        #"depth", #865
        "bio17", #315
        #"bio18", #497
        "bio13",
        "silt", #781
        "carbon", #678
        "bio11") #175
    }
    #surplus in 1 variable. I have 3 climatic variables (1 of temperature and 2 of pp) and 3 of soil. We drop the soil variable with highest rank to accomplish 3 climatic and 2 soil. Moreover, we drop bio18 becuase dramatic changes on vlaues, and include the next most explicativa variable of its cluster (bio13).

#group 3
results_ranks[[3]] #load sum of ranks of variables
group_variables[group_variables$groups==3,]$species #load the names of species
variables_for_vif[[3]] = c(
    "bio12", #71
    "silt", #120
    #"carbon", #157
    "bio4", #65
    #"bio5", #115
    "clay", #84
    "bio13") #48
    #Tenemos 3 variables de suelo y 4 de clima. Sobran 2, una de suelo y una de clima. Quitamos la variable de suelo menos explicativa, ya que todas tienen VIF bajo, carbon. Quitamos la variable de clima menos explicativa (bio5) y con eso ya el VIF queda por debajo de 5 en todos los casos.
    plot(variables_clusters[[3]][["bio12"]])
    plot(variables_clusters[[3]][["silt"]])
    plot(variables_clusters[[3]][["carbon"]])
    plot(variables_clusters[[3]][["bio4"]])
    plot(variables_clusters[[3]][["bio5"]])
    plot(variables_clusters[[3]][["clay"]])
    plot(variables_clusters[[3]][["bio13"]])

    #selección anterior
    if(FALSE){
        #"carbon", #157
        "bio4", #65
        "clay", #84
        "silt", #120
        "bio13", #48
        "bio12") #71
    }
    #surplus in 1 variable. I have 3 climatic variables (1 of temperature and 2 of pp) and 3 of soil. We drop the soil variable with highest rank to accomplish 3 climatic and 2 soil. 

#group 4
results_ranks[[4]] #load sum of ranks of variables
group_variables[group_variables$groups==4,]$species #load the names of species
variables_for_vif[[4]] = c(
    #"depth", #434
    #"ph", #490
    #"bio19", #291
    "bio17", #348
    "sand", #128
    #"carbon", #458
    "bio4", #114
    "clay", #43  
    "bio12") #503
    #Tenemos 5 variables de suelo y de 2 clima. Hay que quitar 3 de suelo y añadir una de clima. Quitamos las tres variables de suelo menos explicativas (depth, ph y carobn). Y añadimos la variable de clima más explicativa en cuyo grupo no haya ninunga de las 2 variables de clima seleccionadas para evitar así que aumente el VIF, cogemos bio12. Quitamos además bio19 porque tiene cambiso bruscos en la zona de brasil, lo cambiamos por la siguiente variable en ese grupo (2), que es bio17.
    plot(variables_clusters[[4]][["depth"]])
    plot(variables_clusters[[4]][["ph"]])
    plot(variables_clusters[[4]][["bio17"]])
    plot(variables_clusters[[4]][["sand"]])
    plot(variables_clusters[[4]][["carbon"]])
    plot(variables_clusters[[4]][["bio4"]])
    plot(variables_clusters[[4]][["clay"]])
    plot(variables_clusters[[4]][["bio12"]])

    #selección anterior
    if(FALSE){
        #"depth", #434
        #"ph", #490
        #"bio19", #291
        "sand", #128
        #"carbon", #458
        "bio4", #114
        "clay", #43
        "bio12",
        "bio17")
    }
    #Surplus 2 variables. I have 2 climatic variables (1 of temperature and 1 of pp) and 5 of soil. We drop the three soil variables with highest rank, and select the climatic variable less correlated and with lowest rank of all variables don't selected. We select the most explicative variable of the group of pH (variable dropped), bio12. Moreover, we drop bio19 becuase dramatic changes on vlaues, and include the next most explicativa variable of its cluster (bio17)

#group 5
results_ranks[[5]] #load sum of ranks of variables
group_variables[group_variables$groups==5,]$species #load the names of species
variables_for_vif[[5]] = c(
    #"carbon", #217
    "bio4", #37
    "silt", #128
    "bio17", #82
    #"bio18", #110
    "bio16", #138
    "depth") #117
    #Tenemos 3 variables de suelo y 3 de clima. Hay que quitar una de suelo, quitamos la menos explicativa, carbon. Quitamos bio18 por los cambiso bruscos y ponemos la siguiente más explicativa de ese grupo (5), que es bio16.
    plot(variables_clusters[[5]][["carbon"]])
    plot(variables_clusters[[5]][["bio4"]])
    plot(variables_clusters[[5]][["silt"]])
    plot(variables_clusters[[5]][["bio17"]])
    plot(variables_clusters[[5]][["bio18"]])
    plot(variables_clusters[[5]][["bio16"]])    
    plot(variables_clusters[[5]][["depth"]])
    
    #selección anterior
    if(FALSE){
        "depth", #114
        "bio17", #82
        #"bio18", #110
        "bio16",
        "silt", #128
        #"carbon", #217
        #"bio5", #238
        "bio4") #37
    }     
    #surplus 2 variables. There are 4 climatic variables (2 of pp and 2 of tª), whilst there are 3 soil variables. We will drop 1 climatic and 1 soil variables, those with highest rank (lowest D2). Moreover, we drop bio18 becuase dramatic changes on vlaues, and include bio16, the most explicative variable in the group of bio18. 

#HACEMOS UNA NUEVA TABLA SOLO CON ESAS VARIABLES
variables_table_2 = list()
for (i in unique(group_variables$groups)){
    variables_table_2[[i]]<-variables_table[[i]][, variables_for_vif[[i]]]
}

#check that the correct variables has been selected
for (i in unique(group_variables$groups)){
    print(ncol(variables_table_2[[i]]) == length(variables_for_vif[[i]]))
}

for (i in unique(group_variables$groups)){
    print(colnames(variables_table_2[[i]]) == variables_for_vif[[i]])
}


#PERO PUEDE HABER VARIABLES QUE SON COMBINACIÓN LINEAL DE OTRAS VARIABLES...
#CALCULAMOS EL VARIANCE INFLATION FACTOR (explicación en la diapositiva)
require(HH)
#VIF = 1/(1-Ri^2), where Ri^2 is the R squared of the linear regresion of variable  (response) i in front of the rest of variables (predictors). Therefore, how much more high is Ri^2, lower is (1-Ri^2) y thus higher is VIF. A variable with high vif has a lot of its variability explained by the other variables. It is to say, the information contained in this variables is also contained in the rest of varibles. This variable is a lineal combination of the rest of variables. 
#Cualquier variable con un VIF mayor que 5 es una combinacion lineal de otras variables. But we don't drop that variable, we look for the most correlated varaibles (according to dendogram) that have been selected, and drop the variable with the lowest explicative power. See cluster two as an example: 

#cluster 1
results_ranks[[1]] #load sum of ranks of variables
group_variables[group_variables$groups==1,]$species #load the names of species
final_variables_cluster_1 = vif(variables_table_2[[1]]) 
final_variables_cluster_1 #No problem

#cluster 2
results_ranks[[2]] #load sum of ranks of variables
group_variables[group_variables$groups==2,]$species #load the names of species
final_variables_cluster_2 = vif(variables_table_2[[2]]) 
final_variables_cluster_2 #No problem

#cluster 3
results_ranks[[3]] #load sum of ranks of variables
group_variables[group_variables$groups==3,]$species #load the names of species
final_variables_cluster_3 = vif(variables_table_2[[3]]) 
final_variables_cluster_3 #No problem 

#cluster 4
results_ranks[[4]] #load sum of ranks of variables
group_variables[group_variables$groups==4,]$species #load the names of species
final_variables_cluster_4 = vif(variables_table_2[[4]]) 
final_variables_cluster_4 #No problem 

#cluster 5
results_ranks[[5]] #load sum of ranks of variables
group_variables[group_variables$groups==5,]$species #load the names of species
final_variables_cluster_5 = vif(variables_table_2[[5]]) 
final_variables_cluster_5 #No problem

#sacamos los nombres del resultado para usarlos como lista de variables seleccionadas
ultimate_variables = list(names(final_variables_cluster_1), names(final_variables_cluster_2), names(final_variables_cluster_3), names(final_variables_cluster_4), names(final_variables_cluster_5))

#write the selected variables as names only 
save(ultimate_variables, file="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/final_variables/list_selected_variables.rda")


##########################################################################
############### Dramatic changes in variable values ######################
##########################################################################
###Test if some of the variables selected have extremely changes in values (big change in a small area)
bio_vars = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

#loop for see what bio variables are included in the final selection
for (k in bio_vars){
    print("******************")
    print(k)
    for (i in 1:length(ultimate_variables)){
        print(k %in% ultimate_variables[[i]])
    }
}

#bio4: 
bio4 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio4.asc")
plot(bio4) #no problem

#bio11: 
bio11 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio11.asc")
plot(bio11) #no problem

#bio12: 
bio12 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio12.asc")
plot(bio12) #no problem

#bio13: 
bio13 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio13.asc")
plot(bio13) #no problem

#bio16: 
bio16 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio16.asc")
plot(bio16) #no problem

#bio17: 
bio17 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio17.asc")
plot(bio17) #no problem

#bio18: 
bio18 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio18.asc")
plot(bio18) #problems

#bio19: 
bio19 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio19.asc")
plot(bio19) #problems
#bio18 is in cluster 1, 2 and 5. The areas with problems are included in distributio fo species to these clusters. 
#bio19 is in cluster 4. The areas with problems are included in distribution of species of this cluster. 

##########################################################
#########Variable selection for small species#############
##########################################################

###make the selection of variable for the selected species according the two criterias of variables selection (new and ancient). The difference between these two criterias is the number of variables selectec, in the new criteria we have a fixed number of variables (5) and ratio soil/climate (2/3). In the ancient, can be selected more than 5 variables, and the ratio soil/climate  is variable (bias potential response to climate of species)


###definimos el directorio de trabajo
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus")

#select variables according to the ancient criteria (more than 5 variables)
load("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/variables/final_variables_MAL/list_selected_variables_ancient_criteria.rda") 
selected_variables_ancient = ultimate_variables 

#select variables according to the new criteria (5 variables: 3 climatic and 2 of soil)
load("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/final_variables/list_selected_variables.rda") #this change the ultimate variable object
selected_variables_new = ultimate_variables

###list ocurrences
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

#load cluster number for each species
##Load the group variable and the ranks 
group_species = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/_tables/complete_2_10_g_5.csv", header=TRUE)[,c("species", "groups")] 

##Load the group variable and the ranks 
group_variables = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/_tables/complete_2_10_g_5.csv", header=TRUE)  #all the info

#load the number of ocurrences for each species 
number_ocurrences = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/ocurrences_per_species.csv", header=TRUE)

#Create a rank for variables (D2) in ech cluster 
###create a list with the sum of ranks for each variable in each cluster and oredered from low to high. 
results_ranks = list() #create a list
selected_var_cor_analysis = c(names(group_variables[,-which(colnames(group_variables)=="species" | colnames(group_variables)=="groups")])) #select the names of variables from the table with ranks for all variables and ranks. We will use this table for create the data frame with the sum of ranks. 
for (i in unique(group_variables$groups)){ #for each each cluster

    group = group_variables[group_variables$groups==i,-which(colnames(group_variables)=="species" | colnames(group_variables)=="groups")] #select the row if the cluster [i]

    #create empty vectors for the sum of each type of variable
    sum_ranks = NULL
    #make a loop of each type of variable
    for (k in group){ #for each temperature variable 
        sum_ranks = append(sum_ranks, sum(k))
    } #sum its ranking in all species of the cluster [i]
    sum_ranks=as.data.frame(sum_ranks) #convert to data.frame 
    sum_ranks = cbind(selected_var_cor_analysis, sum_ranks) #bind the data frame of rank_sums with the a vector with variables names
    names(sum_ranks)[1] = "Variables" 
    names(sum_ranks)[2] = "sum" #establish column names
    sum_ranks = sum_ranks[order(sum_ranks$sum, decreasing=FALSE),] #reorder the rows according to sum variable, with a increasing order (low values at first)

    results_ranks[[i]] = sum_ranks #bind all vectors (sums of each varible category) and include them as element number [i] of the previous created list. 
}

 
############## New criteria #################

#make a loop for calculate species with not enough ocurrences 
number_variables_new = data.frame()
for (i in epithet_species_list){

    #select the group (cluster) of the species [i]
    variables_cluster = group_species[which(group_species$species == i),]$groups 

    #select the number of ocurrences of the species [i]
    n_presences = number_ocurrences[number_ocurrences$species==i,]$number_ocurrences

    #select the selected variables of the corresponding cluster
    selected_variables = selected_variables_new[[variables_cluster]]

    #if - else    
    if (n_presences<length(selected_variables)*10){ #if there is not 10 ocurrences for each variable 
        n_variables = floor(n_presences/10)

    } else { 

        n_variables = length(selected_variables)

    }

    number_variables_new = rbind(number_variables_new, c(length(selected_variables), n_presences, n_variables))
}

#bind species names
data_number_variables_new = cbind(epithet_species_list, number_variables_new)
names(data_number_variables_new) = c("species", "n_variables_stack", "n_ocurrences", "final_n_variable") #stablish the names of the columns
str(data_number_variables_new)
summary(data_number_variables_new)

#subset species in which number of ocurrences is lower
subset_test_new = data_number_variables_new[data_number_variables_new$n_ocurrences<data_number_variables_new$n_variables_stack*10,]

#test in this species if the number of ocurrences is higher or equal than the final number of variables by 10, thus we test if there is 10 ocurrences for each variable
subset_test_new$n_ocurrences>=subset_test_new$final_n_variable*10 #YES

#write the resulting file
write.csv(data_number_variables_new, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/variable_selection_inside_clusters/data_number_variables.csv", row.names=FALSE)


#create a list with the ranks of the selected variables and the final number of variables for each species. This is for help in manual variable selection
ranks_variables_cluster_new = list() #empty list
for (i in subset_test_new$species){ #for each species with low number of ocurrences

    #select the cluster of the [i] species
    variables_cluster = group_species[group_species$species==i,]$groups 
    
    #select the selected variables of that cluster
    selected_variables = selected_variables_new[[variables_cluster]] 

    #select the ranks of that cluster
    ranks = results_ranks[[variables_cluster]] 

    #select the final number of variables for [i] species
    final_n_variable = subset_test_new[subset_test_new$species==i,]$final_n_variable

    #bind the ranks and the final number of variables
    ranks = cbind(ranks, final_n_variable)

    #select only the ranks of variables selected and save in the list
    ranks_variables_cluster_new[[i]] = ranks[ranks$Variables %in% selected_variables,]
}
length(ranks_variables_cluster_new)

##selection by hand for each species 
#If: 
    #Final variables = 4, we select the three most explicative climatic variables and the most explicative soil variable. 
    #Final variables = 3, we select the two most explicative climatic variables and the most explicative soil variable
    #Final variable = 2, we select the most explicative climatic variable and the most explicative soil variables
    #Final variable = 1, we select the most explicative climatic variable 

final_variables_low_number_ocurrence_species_new = list() #empty list

#amamiana
ranks_variables_cluster_new$amamiana
final_variables_low_number_ocurrence_species_new$amamiana <- "bio13"
#canariensis
ranks_variables_cluster_new$canariensis
final_variables_low_number_ocurrence_species_new$canariensis <- c("bio13", "bio4", "clay")
#cubensis
ranks_variables_cluster_new$cubensis
final_variables_low_number_ocurrence_species_new$cubensis <- c("clay", "bio4")
#culminicola
ranks_variables_cluster_new$culminicola
final_variables_low_number_ocurrence_species_new$culminicola <- c("bio4", "depth")
#dalatensis
ranks_variables_cluster_new$dalatensis
final_variables_low_number_ocurrence_species_new$dalatensis <- c("bio4", "bio17", "depth")
#fragilissima
ranks_variables_cluster_new$fragilissima
final_variables_low_number_ocurrence_species_new$fragilissima <- c("bio13", "clay")
#krempfii
ranks_variables_cluster_new$krempfii
final_variables_low_number_ocurrence_species_new$krempfii <- c("bio4", "bio17", "depth")
#kwangtungensis
ranks_variables_cluster_new$kwangtungensis
final_variables_low_number_ocurrence_species_new$kwangtungensis <- c("bio11", "silt")
#luchuensis
ranks_variables_cluster_new$luchuensis
final_variables_low_number_ocurrence_species_new$luchuensis <- c("bio11")
#maestrensis
ranks_variables_cluster_new$maestrensis
final_variables_low_number_ocurrence_species_new$maestrensis <- c("clay", "bio4")
#maximartinezii
ranks_variables_cluster_new$maximartinezii
final_variables_low_number_ocurrence_species_new$maximartinezii <- c("clay", "bio4")
#morrisonicola
ranks_variables_cluster_new$morrisonicola
final_variables_low_number_ocurrence_species_new$morrisonicola <- c("bio13", "bio4", "clay")
#squamata
ranks_variables_cluster_new$squamata
final_variables_low_number_ocurrence_species_new$squamata <- c("bio13", "clay")
#taiwanensis
ranks_variables_cluster_new$taiwanensis
final_variables_low_number_ocurrence_species_new$taiwanensis <- c("bio13", "bio4", "bio12", "clay")
#torreyana
ranks_variables_cluster_new$torreyana
final_variables_low_number_ocurrence_species_new$torreyana <- c("bio11", "bio13", "bio12", "silt")
#tropicalis
ranks_variables_cluster_new$tropicalis
final_variables_low_number_ocurrence_species_new$tropicalis <- c("bio11", "silt")

#length of the list must be equal to the number of rows of ancient subset of species with a low number of ocurrences
length(final_variables_low_number_ocurrence_species_new) == nrow(subset_test_new)

#save the result
save(final_variables_low_number_ocurrence_species_new, file="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/final_variables_low_number_ocurrence_species.rda")

################# ancient criteria ###################
#make a loop for calculate species with not enough ocurrences
number_variables_ancient = data.frame()
for (i in epithet_species_list){

    #select the group (cluster) of the species [i]
    variables_cluster = group_species[which(group_species$species == i),]$groups 

    #select the number of ocurrences of the species [i]
    n_presences = number_ocurrences[number_ocurrences$species==i,]$number_ocurrences

    #select the selected variables of the corresponding cluster
    selected_variables = selected_variables_ancient[[variables_cluster]]

    #if - else    
    if (n_presences<length(selected_variables)*10){ #if there is not 10 ocurrences for each variable 
        n_variables = floor(n_presences/10)

    } else { 

        n_variables = length(selected_variables)

    }

    number_variables_ancient = rbind(number_variables_ancient, c(length(selected_variables), n_presences, n_variables))
}

#bind species names
data_number_variables_ancient = cbind(epithet_species_list, number_variables_ancient)
names(data_number_variables_ancient) = c("species", "n_variables_stack", "n_ocurrences", "final_n_variable") #stablish the names of the columns
str(data_number_variables_ancient)
summary(data_number_variables_ancient)

#subset species in which number of ocurrences is lower
subset_test_ancient = data_number_variables_ancient[data_number_variables_ancient$n_ocurrences<data_number_variables_ancient$n_variables_stack*10,]

#test in this species if the number of ocurrences is higher or equal than the final number of variables by 10, thus we test if there is 10 ocurrences for each variable
subset_test_ancient$n_ocurrences>=subset_test_ancient$final_n_variable*10 #YES

#write the resulting file
write.csv(data_number_variables_ancient, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/variable_selection_inside_clusters/comparison_selection_criterias/data_number_variables_ancient.csv", row.names=FALSE)

#select variables automaticallty in basis on D2 rank, because this is ancient criteria
final_variables_low_number_ocurrence_species_ancient = list() #empty list

for (i in subset_test_ancient$species){ #for each species with low occurrences according to ancient criteria (lower than 70 - 7 variables)

    #select the cluster of the [i] species
    variables_cluster = group_species[group_species$species==i,]$groups
    
    #extract the variables of the corresponding cluster
    selected_variables = selected_variables_ancient[[variables_cluster]]

    #select sum ranks of variables in this cluster
    ranks = results_ranks[[variables_cluster]]

    #final number of variables according to ancient criteria
    total_number_variables = subset_test_ancient[subset_test_ancient$species==i,]$final_n_variable

    #ranks of the variables used of this species (ordenado de menor a mayor)
    ranks_variables_cluster = ranks[ranks$Variables %in% selected_variables,]

    #select the variables with lower sum rank, the number depends on total_number_variables
    final_variables_low_number_ocurrence_species_ancient[[i]] = ranks_variables_cluster[c(1:total_number_variables), "Variables"]
}


#length of the list must be equal to the number of rows of ancient subset of species with a low number of ocurrences
length(final_variables_low_number_ocurrence_species_ancient) == nrow(subset_test_ancient)

#other test 
test = NULL
for (i in subset_test_ancient$species){ #for each species with low occurrences according to ancient criteria (lower than 70 - 7 variables)

    #test if the number of variables selected is the same the final number calculated for each one
    test_result = length(final_variables_low_number_ocurrence_species_ancient[[i]]) == subset_test_ancient[subset_test_ancient$species==i,]$final_n_variable
    
    #save it
    test = append(test, test_result)
}
summary(test) #all speceis have the number of variables that correspong in basis on the number of ocurrences

#save the result
save(final_variables_low_number_ocurrence_species_ancient, file="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/final_variables_low_number_ocurrence_species_ancient.rda")