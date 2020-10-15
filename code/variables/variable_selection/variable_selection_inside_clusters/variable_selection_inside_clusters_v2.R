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
#check
if(FALSE){
    require(tidyverse)
    paste("Pinus", epithet_species_list, sep=" ") == str_trim(as.vector(list_species[,1])) #the seperated epithet more Pinus are equal to list_species? We use str_trim from tidyverse to remove end spaces in each element of list_species
}#it is in false because loading tidyverse lead to load several packages that have a function name "extract", and this gives problems with the extract function of raster. If you want to check run these lines manually

#remove tecunumanii, jaliscana y discolor. These species are not used for Niche paper. The two first because we should downloadad gbif data now, so we would mix gbif data form 2016 and 2019. The third one was impossible to differentiate from P. cembroides
epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("tecunumanii", "jaliscana", "discolor"))]
#check
c("tecunumanii", "jaliscana", "discolor") %in% epithet_species_list

##Load variables
list_variables = list.files("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals", full.names=TRUE, pattern=".asc")
variables = stack(list_variables)

##Load the group variable and the ranks 
group_variables = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/_tables/complete_2_10_g_2.csv", header=TRUE)  
length(unique(group_variables$groups)) #we have 2 groups


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
cos(NISTdegTOradian(30))*res.km #8 km
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
nlayers(stack_distrib) == 112
xres(stack_distrib)[[1]]

#Resolucion de la distribucion de la especie
res.grados.distribution = xres(stack_distrib)[[1]] #res of 0.5 degrees in all raster of distribution. We stack in a raster stack with quick=FALSE, thus it has been checked that all raster introduced in the stack have the same resolution.
res.minutos.distribution = res.grados.distribution*60
res.km.distribution = res.grados.distribution*111.19 #Si 1 grado son 111.19 km, 0.5(res.grados) serán X km y salen 55 km. Es decir, cada celda de mi mapa mide en la realidad 55 x 5 km, pero en el ECUADOR!
res.km.distribution #The radious of the earth at the equator is 6378137.0 meters resulting in a circumference of 40075161.2 meters (2*pi*6378137.0; 2*pi*r). The equator is divided into 360 degrees of longitude, so each degree at the equator represents 111.32 km (40075161.2/360 divided by 1000 for kilometers). As one moves away from the equator towards a pole, however, one degree of longitude is multiplied by the cosine of the latitude (in radians), decreasing the distance, approaching zero at the pole. Therefore, if the maps is far away from the equator, the cells will be smaller. We have to multiplicate the length of the at equator (degrees * 111.19 km) by the cosine of latitude in radians. For example a cell of 50*50 at the equator, will be at 37 degrees of latitude (latitude of for example Spain):

library(NISTunits) #functon to convert degrees to radians ("https://stackoverflow.com/questions/32370485/r-convert-radians-to-degree-degree-to-radians")
cos(NISTdegTOradian(37))*res.km.distribution #44.5
cos(NISTdegTOradian(30))*res.km.distribution #48
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
    col_name = NULL
    #make a loop of each type of variable
    for (k in 1:ncol(group)){ #for each column of group (i.e. for each predictor)

        #select the [k] column
        selected_column = colnames(group)[k]

        #sum the ranks of the [k] column (variable)
        sum_of_column = sum(group[,which(colnames(group) == selected_column)])

        #save the sum
        sum_ranks = append(sum_ranks, sum_of_column)

        #save the column name
        col_name = append(col_name, selected_column)
    } #sum its ranking in all species of the cluster [i]

    #convert to data.frame     
    sum_ranks=as.data.frame(sum_ranks)

    #check colnames in "selected_var_cor_analysis" are correct
    summary(col_name == selected_var_cor_analysis)

    #bind the data frame of rank_sums with the a vector with variables names
    sum_ranks = cbind(selected_var_cor_analysis, sum_ranks)
    names(sum_ranks)[1] = "Variables" 
    names(sum_ranks)[2] = "sum" #establish column names

    #reorder the rows according to sum variable, with a increasing order (low values at first)    
    sum_ranks = sum_ranks[order(sum_ranks$sum, decreasing=FALSE),]

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
clust <- makeCluster(2)
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
clust <- makeCluster(2)
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
    plot(variables_dendo[[i]], main=paste("group ", i, sep=""))
}
dev.off()


#OBSERVAMOS EL PDF Y ELEGIMOS LAS VARIABLES CON LAS QUE QUEREMOS TRABAJAR
#El limite que establezco de correlacion es 0.4, todas aquellas varaibles correlacionadas en menos de un 0.4 y que por tanto se separan antes del 0.6 (que equivale a ese 0.4 de correlacion) las meto sin más. En el resto de casos, si la separaión se ha hecho despues, nos quedamos con solo UN representando del grupo. Esto hace que nos quitemos variables con un nivel de correlacion mayor del 0.4. ES UN NIVEL MUY ESTRICTO, POR ESO SEREMOS FLEXIBLES PARA CASOS CERCANOS AL LIMITE. In GENERAL, only 5 variables for cluster, 3 of climate (at least 1 and 1 of pp and temperature, respectively) and 2 of soil. If in a cluster there is more than 5 groups (and thus more than 5 variables), we select from the representative variable of the groups, the 5 with the highest D2. And controlling the ratio of soil and climate variables. See for more information: http://blog.minitab.com/blog/understanding-statistics/handling-multicollinearity-in-regression-analysis

variables_for_vif = list()

#group 1
results_ranks[[1]] #load sum of ranks of variables
group_variables[group_variables$groups==1,]$species #load the names of species
variables_for_vif[[1]] = c(
    "clay", #clay=440
    "carbon", #carbon=424: cec=491
    "bio4", #bio4=161: bio7=289; bio3=245; bio8=303; bio1=164; bio11=171;  bio6=255;  
    #"silt", #silt=444: depth=500; sand=487
    "bio5", #bio5=168: bio19=360; bio9=214; bio10=198
    "bio12") #bio12=277; bio2=489; ph=365; bio18=332; bio14=419; bio17=377
    #"bio16") #bio16=387: bio15=402; bio13=413
    
    #Tenemos 2 variables de temperatura, 2 de precipitation, y 3 de suelo, haciendo un total de 7 variables, sobran 2. Sobra 1 de suelo y otra de clima, eliminiamos las que tienen el ranking mas alto de esos dos grupos. Son silt y bio16

    plot(variables_clusters[[1]][["clay"]])
    plot(variables_clusters[[1]][["carbon"]])
    plot(variables_clusters[[1]][["bio4"]])
    plot(variables_clusters[[1]][["silt"]])
    plot(variables_clusters[[1]][["bio5"]])
    plot(variables_clusters[[1]][["bio12"]])
    plot(variables_clusters[[1]][["bio16"]]) 
    
 

#group 2
results_ranks[[2]] #load sum of ranks of variables
group_variables[group_variables$groups==2,]$species #load the names of species
variables_for_vif[[2]] = c(
    "clay", #clay=978
    "bio4", #bio4=432: bio8=745; bio1=1027; bio11=569; bio6=638; bio3=564; bio7=456;
    "bio15", #bio15=1509; bio13=1604; bio16=1607  
    #"carbon", #carbon=1671; bio9=1721; bio10=1688; bio5=1736; cec=1767  
    "silt", #silt=1489; sand=1386; depth=1711;
    #"bio19", #bio19=1024
    "bio17") #bio17=732; bio2=1415; ph=1284; bio18=848; bio12=1194; bio14=742 
    
    #Tenemos 3 variables de suelo, 3 variables de precipitacion y 1 de temperatura, haciendo un total de 7, sobran 2. Quitamos una de pp y una de suelo, las que tengan el ranking más alto de esa categoria: Carbon y bio15. 

    #Hemos cambiado sand por silt, porque sand y clay tienen un VIF muy alto cuando estan juntas, así que hemos cogido del grupo de silt/sand/depth un sustituto para sand Todas esas variables tienen más ranking que clay, así que es mejor conservar clay. 

    #bio19 tiene cambios bruscos, que en los datos actuales no afectan a zonas con pinos, pero en las proyecciones del futuro se acercan más, así que la hemos quitado, y añadido la variable de pp siguiente, bio15. No tiene sentido buscar la siguiente más cercana en correlacion a bio19, porque sería del grupo de bio17 y por tanto estaría muy correlacionada con esa variable.
        #bio19 is in cluster 2. The areas with problems and they are not included in the in the PA buffer of species in cluster 2. However, the future projections have problems in areas closer distribution of pines, so to avoid risks, we will change it for the next most explicative pp variable.

    plot(variables_clusters[[2]][["clay"]])
    plot(variables_clusters[[2]][["bio4"]])
    plot(variables_clusters[[2]][["bio15"]])
    plot(variables_clusters[[2]][["carbon"]])
    plot(variables_clusters[[2]][["silt"]])
    plot(variables_clusters[[2]][["bio19"]])
    plot(variables_clusters[[2]][["bio17"]]) 



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


#sacamos los nombres del resultado para usarlos como lista de variables seleccionadas
ultimate_variables = list(names(final_variables_cluster_1), names(final_variables_cluster_2))

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

#bio5: 
bio5 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio5.asc")
plot(bio5) #no problem

#bio12: 
bio12 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio12.asc")
plot(bio12) #no problem

#bio15:
bio15 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio15.asc")
plot(bio15) #no problem

#bio17: 
bio17 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio17.asc")
plot(bio17) #no problem
    #We don't have any variable with quarters, so are we are fine. No strong changes..

##########################################################
#########Variable selection for small species#############
##########################################################

###make the selection of variable for the selected species according the two criterias of variables selection (new and ancient). The difference between these two criterias is the number of variables selected, in the new criteria we have a fixed the ratio soil/climate (2/3). In the ancient, the ratio soil/climate  is variable (bias potential response to climate of species), because the variables are selected automatically in basis on the predictive power (without considering the ratio soil/climate)


###definimos el directorio de trabajo
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus")

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
#check
if(FALSE){
    require(tidyverse)
    paste("Pinus", epithet_species_list, sep=" ") == str_trim(as.vector(list_species[,1])) #the seperated epithet more Pinus are equal to list_species? We use str_trim from tidyverse to remove end spaces in each element of list_species
}#it is in false because loading tidyverse lead to load several packages that have a function name "extract", and this gives problems with the extract function of raster. If you want to check run these lines manually

#remove tecunumanii, jaliscana y discolor. These species are not used for Niche paper. The two first because we should downloadad gbif data now, so we would mix gbif data form 2016 and 2019. The third one was impossible to differentiate from P. cembroides
epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("tecunumanii", "jaliscana", "discolor"))]
#check
c("tecunumanii", "jaliscana", "discolor") %in% epithet_species_list

#load cluster number for each species
##Load the group variable and the ranks 
group_species = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/_tables/complete_2_10_g_2.csv", header=TRUE)[,c("species", "groups")] 

##Load the group variable and the ranks 
group_variables = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/_tables/complete_2_10_g_2.csv", header=TRUE)  #all the info

##check
group_variables$species == group_species$species

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
names(ranks_variables_cluster_new)

##selection by hand for each species 
#If: 
    #Final variables = 4, we select the three most explicative climatic variables and the most explicative soil variable. 
    #Final variables = 3, we select the two most explicative climatic variables and the most explicative soil variable
    #Final variable = 2, we select the most explicative climatic variable and the most explicative soil variables
    #Final variable = 1, we select the most explicative climatic variable 

final_variables_low_number_ocurrence_species_new = list() #empty list

#amamiana
ranks_variables_cluster_new$amamiana
final_variables_low_number_ocurrence_species_new$amamiana <- "bio4"
#canariensis
ranks_variables_cluster_new$canariensis
final_variables_low_number_ocurrence_species_new$canariensis <- c("bio4", "bio5", "carbon")
#cubensis
ranks_variables_cluster_new$cubensis
final_variables_low_number_ocurrence_species_new$cubensis <- c("bio4", "clay")
#culminicola
ranks_variables_cluster_new$culminicola
final_variables_low_number_ocurrence_species_new$culminicola <- c("bio4", "clay")
#fragilissima
ranks_variables_cluster_new$fragilissima
final_variables_low_number_ocurrence_species_new$fragilissima <- c("bio4", "clay")
#krempfii
ranks_variables_cluster_new$krempfii
final_variables_low_number_ocurrence_species_new$krempfii <- c("bio4", "bio17", "clay")
#luchuensis
ranks_variables_cluster_new$luchuensis
final_variables_low_number_ocurrence_species_new$luchuensis <- c("bio4", "clay")
#maestrensis
ranks_variables_cluster_new$maestrensis
final_variables_low_number_ocurrence_species_new$maestrensis <- c("bio4", "clay")
#maximartinezii
ranks_variables_cluster_new$maximartinezii
final_variables_low_number_ocurrence_species_new$maximartinezii <- c("bio4", "bio17", "bio15", "clay")
#morrisonicola
ranks_variables_cluster_new$morrisonicola
final_variables_low_number_ocurrence_species_new$morrisonicola <- c("bio4", "bio17", "clay")
#squamata
ranks_variables_cluster_new$squamata
final_variables_low_number_ocurrence_species_new$squamata <- c("bio4", "bio17", "bio15", "clay")
#taiwanensis
ranks_variables_cluster_new$taiwanensis
final_variables_low_number_ocurrence_species_new$taiwanensis <- c("bio4", "bio17", "bio15", "clay")
#torreyana
ranks_variables_cluster_new$torreyana
final_variables_low_number_ocurrence_species_new$torreyana <- c("bio4", "bio17", "bio15", "clay")
#tropicalis
ranks_variables_cluster_new$tropicalis
final_variables_low_number_ocurrence_species_new$tropicalis <- c("bio4", "carbon")

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
    selected_variables = selected_variables_new[[variables_cluster]]

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


#select variables automaticallty in basis on D2 rank, because this is ancient criteria
final_variables_low_number_ocurrence_species_ancient = list() #empty list

for (i in subset_test_ancient$species){ #for each species with low occurrences according to ancient criteria (lower than 70 - 7 variables)

    #select the cluster of the [i] species
    variables_cluster = group_species[group_species$species==i,]$groups
    
    #extract the variables of the corresponding cluster
    selected_variables = selected_variables_new[[variables_cluster]]

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

#compare both criterias
for(i in subset_test_ancient$species){

    print("")
    print("###########")
    print(i)

    print("old criteria (utomaitic in basis on predicitive power")
    print(final_variables_low_number_ocurrence_species_ancient[[which(names(final_variables_low_number_ocurrence_species_ancient) == i)]])

    print("new criteria considering the ratio climate/soil (manual)")
    print(final_variables_low_number_ocurrence_species_new[[which(names(final_variables_low_number_ocurrence_species_new) == i)]])
}

#save the result
save(final_variables_low_number_ocurrence_species_ancient, file="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/final_variables_low_number_ocurrence_species_ancient.rda")