#code for downloading ocurrences data 
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus")

#########################
########OCURRENCIAS######
#########################

#Descargamos las ocurriencias de gbif, pero hay que tener muuucho cuidado y limpiar muy bien. Varios puntos a destacar: 

    #Hay datos de jardines botánicos que no nos valen una mierda, por ejemplo una presencia de P. halepensis en sudafrica por un jardín botánico que lo tiene guardado. 
        #Revisión 2018: Por lo que eh revisado en mis apuntes del master de la asignatura de datos, las coordenadas que se le ponen a un pliego deben ser las del sitio donde se cogió, no dle herbario, de hecho esos datos se usan luego en la identificación. Por ejemplo, hay halepensis en sudáfrica según datos de especímen de herbario, pero cogidos de la propia sudafrica ("https://www.gbif.org/occurrence/685237824"). Eso no es raro porque hay observaciones directas también de esta especie en sudáfrica. 
    #Tambien mucho ojo con las invasiones, por ejemplo sylvestris estas en muchos tiios feura de españa, esos puntos no son vale porque probablmenet las condicioens en las que se desarrollara fuera de europa no puede desarrollarse dentro de eruopa, porque aquí habrá otras restricciones como interacciones bióticas. Entonces Tienes que mirar el mapa de distribucion que Bianca te dio y qutiar todos los ptnuiso que queden fuera de él. 
    #Sesgos en torno a carreteras, universidades, y en general en determinados continenetes (Europa y USA muchisiomos datos, siberia nada porque los rusos no liberan los datos). 

##librerias
library(dismo)

#################################################################
###########download ocurrences of Pinus genus as a whole#########
#################################################################

pines = gbif(genus="Pinus", #caracter. El nombre del género
    species="*", #Caraceter: nombre de la especie. Una "*" para descargar el género entero. Añade "*" al nombre de la especie para conseguir todas las variantes del nombre (por ejemplo, con y sin el nombre del autor) y los sub-taxa. 
    ext=NULL, #Extensión del objetvo para limtar la extensión geográfica de los datos. Una extension puede ser creada usando funcones como "drawExtent" and "extent"
    args=NULL, #Caracter. Argumentos adicionales para refinar la consulta Mira los parametro de consulta en http://www.gbif.org/developer/occurrence para más detalles. 
    geo=TRUE, #logical. si es TRUE, solos los records que estén georreferenciados (valores de longitud y latitd) serán descargados.
    sp=FALSE, #Logico. Si es TRUE entonces "geo" se pondrá automáticamente cmo TRUE y de se devolverá un objeto de la clase SpatialPointsDataFrame. Prefiero guardarlo como un data.frame normal. 
    removeZeros=TRUE, #logical. Si es TRUE, todos los records que tengan una latitud O longitud de cero serán eliminados si geo=TRUE, o se les dará el valor de NA si geo=FALSE. If False, solo los records que tienen TANTO la latitud como la longitud en cero seran eliminados o puestos como NA. De todas maneras, luego nosotros miramos si hay NAs en las columnas de longitud y latitud. 
    download=TRUE, #Lógico. Si es TRUE, los records serán descargados, si es FALSE simplemente se mostará el número de records existentes. 
    ntries=5, #Numero. Cuantas veces debería la función intentar descargar los datos si la base de datos lo deneiega (quizas porque el servidor de gbif está demasiado ocupado). 
    nrecs=300, #Número. Cuantas records se descargan en un única consulta (max is 300).
    start=1, #Núemro. Ocurrencia a partir del cual se empiezan a descargar los datos Te interesa empezar desde la primera ocurrencia. 
    end=Inf) #Número. Última ocurrencia que descargar. Nosotros queremos todas, no queremos que corte en un punto.   
    #NO FUNCIONA: En la ocurrencia 199800 acaba fallando, hay que descargar las ocurencias especie por especie.  
write.csv(pines, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/ocurrences/ocurrences/pines.csv")

#################################################################
###### Download the ocurrences for each species separately ######
#################################################################
##Descargamos las especies una a una con un loop, conforme vamos haciendo los modelos. 
#IMPORTANTE: BUSCA EN FLORA IBERICA Y ATALAS EUROPEO el nombre aceptado de cada espece cuando la busques. 
#Descargamos TOOODOS LOS DATOS de la especie, por eso usamos "*", queremos todas las variantes del nombre (con y sin autor, etc...). Si hay algo raro tanto a nivel taxoónicmo como de distribucion lo detectaremos en el depruado de ocurrencias con los mapas de Bianca. 
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

#remove jaliscana and tecunumanii. These species were included in the new phylogeny of Bianca, but during my stay we don't have neither phylogeny nor maps. These species were not included.
if("jaliscana" %in% epithet_species_list |  "tecunumanii" %in% epithet_species_list){
    epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("jaliscana"))]
    epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c( "tecunumanii"))]    
}
c("jaliscana", "tecunumanii") %in% epithet_species_list

#we add discolor because this species was used in the beginning
pos_discolor=28
epithet_species_list = c(epithet_species_list[1:(pos_discolor-1)], "discolor", epithet_species_list[pos_discolor:length(epithet_species_list)])
length(epithet_species_list) == 113
epithet_species_list[pos_discolor] == "discolor"

#load data frame with genus and especific_epithet of each species. 
list_species = cbind.data.frame(rep("Pinus", length(epithet_species_list)), epithet_species_list)
colnames(list_species) <- c("genus", "specific_epithet")
str(list_species)
head(list_species) #load data frame with genus and especific_epithet of each species. 

#Set problematic species: probably the problem will be in mugo or sylvestris
problem_species = list_species[list_species$specific_epithet == "mugo" | list_species$specific_epithet == "sylvestris",] #drop mugo and sylvestris because of the problems in download with the loop. I have discovered this after run the first loop, you can see the results of the gbif request in record_gbif_request.txt. This output was made with all species except pinus mugo and pinus sylvestris, which I will download separately.

#Set non-problematic species
list_species = list_species[!list_species$specific_epithet %in% problem_species$specific_epithet,]
nrow(list_species) == 113-nrow(problem_species) #the selected speceis should be equal to the difference between total speceis (113) and number of problem species
str(list_species)
problem_species$specific_epithet %in% list_species$specific_epithet #check that problematic species are not included in the list_species

#Download data of non-problematic species
require(dismo) #load the required package dismo
for (i in 1:nrow(list_species)){ 
       presences = #create a data.frame called presences
           gbif(genus="Pinus", #filled with a download from gbif of Pinus genus
               species=paste(list_species[i,]$genus, #select the species name from list_species. For the corresponding row (1 to 113), paste "genus"
                   paste(list_species[i,]$specific_epithet,"*", sep=""), #and the result of paste "specific_epithet" with "*"
                   sep=" "), #with a species
               ext=NULL, args=NULL, geo=TRUE, sp=FALSE, removeZeros=TRUE, download=TRUE, ntries=5, nrecs=300, start=1, end=Inf) 

       assign(paste(list_species[i,]$genus, list_species[i,]$specific_epithet, sep="_"), presences) #Change the name of presences by genus_specific_epithet. This create 113 data frames, one for each species with the name of the species. In case that we need to use the, before create the csv files. 
       write.csv(presences, paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species", paste(list_species[i,]$genus, paste(list_species[i,]$specific_epithet, "csv", sep="."), sep="_"), sep="/"), row.names=FALSE) #create a csv with presence. The path will be the result of paste all the path except the name and the result of paste genus with the result of paste specific_epithet with csv. This creates something like that: "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/ocurrences/ocurrences/Pinus_albicaulis.csv".
} #everything is ok


##Download data of problematic species
for (i in 1:nrow(problem_species)){ 
       presences = gbif(genus="Pinus", species=paste(problem_species[i,]$genus, paste(problem_species[i,]$specific_epithet,"*", sep=""), sep=" "), ext=NULL, args=NULL, geo=TRUE, sp=FALSE, removeZeros=TRUE, download=TRUE, ntries=5, nrecs=300, start=1, end=Inf) 
       assign(paste(problem_species[i,]$genus, problem_species[i,]$specific_epithet, sep="_"), presences) 
       write.csv(presences, paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species", paste(problem_species[i,]$genus, paste(problem_species[i,]$specific_epithet, "csv", sep="."), sep="_"), sep="/"), row.names=FALSE) 
} #we have the same problem with mugo

#Results of sylvestris
sylvestris = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_sylvestris.csv", header=T)
unique(sylvestris$speciesKey) #with sylvestris there is no problem. 

#Results of mugo
mugo = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_mugo_wrong.csv", header=T)
unique(mugo$speciesKey) #there is a lot of species. Here is the problem
nrow(mugo[mugo$speciesKey==5285385,]) == 21.882 #This is the number of ocurrences for pinus mugo according to the web page of gbif. The numbers are not equal.

#################################################################
########## Download ocurrences of sylvestris separately #########
#################################################################
pinus_sylvestris = gbif(genus="Pinus", species="Pinus sylvestris*", ext=NULL, args=NULL, geo=TRUE, sp=FALSE, removeZeros=TRUE, download=TRUE, ntries=5, nrecs=300, start=1, end=Inf) #download only sylvestris, There is no problem with this species
write.csv(pinus_sylvestris, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_sylvestris.csv", row.names=FALSE)

#################################################################
########## Download ocurrences of mugo separately ###############
#################################################################
#we are going to download directly from the gbif webpage. Only for Pinus mugo, and the only differences with our search in dismo is that we don't apply the removeZeros filter. 

#Load the file directly downloaded from gbif webpage. 
pinus_mugo_web_gbif = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_mugo_web_gbif.csv", header=TRUE, sep="\t")
str(pinus_mugo_web_gbif) #There are problem with the separation with excel, sep is tab and not commas. We cannot read it directly with excel, we have to read it with r (remember problems of strings in lon/lat). 

#Write the data of mugo in a csv sep with commas.
write.csv(pinus_mugo_web_gbif, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/ocurrences/ocurrences/tree_species/Pinus_mugo.csv")









