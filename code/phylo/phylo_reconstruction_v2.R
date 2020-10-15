###########################################
############ VARIABLE SELECTION ###########
###########################################

## Load ranked variables based on deviance, which is the % of variance explained of each variable for each species according to a glm and the deviance formula of Weisberg 1980; Guisan & Zimmermann 2000. For further information see "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/species_clustering.R".
d2_rank = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/clustering_results/d2_rank_species_variables.csv", header=TRUE)
str(d2_rank)

## create per variable : i) Sumatory of the rakings; ii) Number of species for whicht the variable is the top1,2 or 3.
ranks_compare = as.data.frame(matrix(NA, ncol=5, nrow=1))
colnames(ranks_compare) <- c("variable", "sum_rank", "number_ones", "number_twoes", "number_threes")
for(i in 2:ncol(d2_rank)){ #for each variable (first column is species names)
    
    #select the [i] variable
    selected_col = d2_rank[,i]
    
    #selcte the name of the [i] varaiable
    variable =  colnames(d2_rank)[i]
    
    #sum the ranks of the [i] variable across species
    sum_rank = sum(selected_col)

    #calculate the umber of species for whicht the variable is the top1,2 or 3.
    number_ones = length(which(selected_col == 1))
    number_twoes = length(which(selected_col == 2))
    number_threes = length(which(selected_col == 3)) 

    #bind all   
    ranks_compare = rbind.data.frame(ranks_compare, 
            cbind.data.frame(variable, sum_rank, number_ones, number_twoes, number_threes))
}

#delete row with NAs
ranks_compare = ranks_compare[-1,]

#check all variables are included
nrow(ranks_compare) == length(colnames(d2_rank)[which(!colnames(d2_rank) == "species")])

#set variable names in a vector and long variable names in another vector to add long vname variables to ranks_compare (variable will be the variable to merge)
variable = c(
    "bio1", 
    "bio2", 
    "bio3", 
    "bio4", 
    "bio5", 
    "bio6", 
    "bio7", 
    "bio8", 
    "bio9", 
    "bio10",
    "bio11",
    "bio12",
    "bio13",
    "bio14",
    "bio15",
    "bio16",
    "bio17",
    "bio18",
    "bio19",
    "clay", 
    "silt", 
    "sand", 
    "ph",   
    "cec",  
    "carbon",
    "depth")
long_var_names = c(
    "Annual Mean Temperature",
    "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
    "Isothermality (BIO2/BIO7) (* 100)",
    "Temperature Seasonality (standard deviation *100)",
    "Max Temperature of Warmest Month",
    "Min Temperature of Coldest Month",
    "Temperature Annual Range (BIO5-BIO6)",
    "Mean Temperature of Wettest Quarter",
    "Mean Temperature of Driest Quarter",
    "Mean Temperature of Warmest Quarter",
    "Mean Temperature of Coldest Quarter",
    "Total (annual) moisture",
    "Moisture of Wettest Month",
    "Moisture of Driest Month",
    "Moisture Seasonality (standard deviation)",
    "Moisture of Wettest Quarter",
    "Moisture of Driest Quarter",
    "Moisture of Warmest Quarter",
    "Moisture of Coldest Quarter",
    "Clay content (\\%)",
    "Silt content (\\%)",
    "Sand content (\\%)",
    "Ph (index * 10)",
    "Cation-exchange capacity (CEC; cmolc/kg)",
    "Organic carbon (g/kg)",
    "Absolute depth to bedrock (cm)")

#bind them
var_names = cbind.data.frame(variable, long_var_names)

#merge the long names with the ranks_compare data.frame
ranks_compare = merge(var_names, ranks_compare, by = "variable")

#order the rows in basis on sum_rank (increasing -> first most explicative variables)
ranks_compare = ranks_compare[order(ranks_compare$sum_rank, decreasing=FALSE),]
ranks_compare

#save results
write.table(ranks_compare, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/ranks_compare.csv", sep=",", row.names = FALSE, col.names = TRUE)

##TABLE 1
#version of sum rank for the paper
table_1 = ranks_compare[c(which(ranks_compare$variable %in% paste("bio", 1:11, sep="")), which(ranks_compare$variable %in% paste("bio", 12:19, sep="")), which(!ranks_compare$variable %in% paste("bio", 1:19, sep=""))), c(1:3)] #separate temperature, humidity and soil variables, but remains the order of ranks inside each group.
colnames(table_1)[1] <- "Bioclim abbreviations"
colnames(table_1)[2] <- "Variables"
colnames(table_1)[3] <- "Sum of ranks"

#reorder columns
table_1 = table_1[,c(2,1,3)]

#change to mayusculas bio abreviatures and remove soil name variables from columns abbreviature
new_abreviations = NULL
for(i in 1:nrow(table_1)){
    selected_row = table_1[i,]

    if(startsWith(as.vector(selected_row$"Bioclim abbreviations"), "bio")){
        new_abreviations = append(new_abreviations, paste("BIO", strsplit(as.vector(selected_row$"Bioclim abbreviations"), split="bio")[[1]][2], sep=""))
    } else { 
        new_abreviations = append(new_abreviations, NA)
    }
    print(paste(as.vector(new_abreviations[i]), "--", as.vector(selected_row$"Bioclim abbreviations"), sep=""))
}
table_1$"Bioclim abbreviations" <- new_abreviations
table_1

#save table as excel
write.table(table_1, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/table_1_paper.csv", sep=",", col.names=TRUE, row.names=FALSE)

#convert to a latex table
require(xtable)
print.xtable(xtable(table_1, align="lllc"), include.rownames=FALSE, NA.string="", floating = FALSE, sanitize.text.function=function(x) {x}, hline.after=c(-1, 0, 11, 19, 26)) #hline.after add a hline after the indicated row. 

#Copia la tabla en tables_latex.tex y corre este comando
system("cd /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering; pandoc -s tables_latex.tex -o tables_latex.doc")

## Selected variables
#-Selección de variables: He usado dos criterios: i) La sumatoria de la posición en el deviance ranking de una variable para todas las especies (más alto indica menos explicativo); ii) Número de especies para las que esa variable está en la posición 1,2 ó 3. Me gusta más el primero, porque tiene peso tanto las especies para las que explica mucho como para las que explica poco, por ejmplo: Una variable puede ser el top1 para muchas especies, pero para el resto no explicar nada, ese sería el caso de bio7, que es la variable más explicativa para más especies (25), pero luego su suma de los rankings es 160 mayor respecto de la primera variable (explica menos para muchas especies). Por tanto, la sumatoria de los rankings sería lo más idóneo para seleccionar una variable que se va a usar para reconstruir el estado de TODAS las especies.
    #-Temperatura -> bio4 (Temperature Seasonality). La sumatoria de los rankings es 593 frente a 740 de bio11 (Mean Temperature of Coldest Quarter), que es la segunda variable más alta de todas. bio4 es el top 1,2,3, para más especies y encima la suma del ranking es claramente menor. bio4 se usó para el cluster 1 y 2, mientras que la otras más cercanas no se usó para ninguno. Además son variables muy parecidas, así que me quedo con la primera según la sumatoria del rank. 
    #-Humedad -> bio17 (humedad del cuarto más seco). Sumatoria de 1109 frente a los 1161 de bio14 (humedad del mes más seco), que es la segunda variable de humedad más alta. En cuanto al número de especies para las que son el top, están muy igualadas. bio17 se usó para el cluster 2, la otra para ninguno. además son muy parecidas. Me quedo con bio17, que es la que tiene la sumatoria de nrakings más alta.
    #-Nota: la info sobre las variables usadas para cada cluster está en "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/variable_selection_inside_clusters/variable_selection_inside_clusters_v2.R".

#################################################################
############ EXTRACT P50 AND SE OF CLIMATIC VARIABLES ###########
#################################################################

#required packages
require(raster)

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

#function to calculate SE
se <- function(x) sd(x)/sqrt(length(x))

#load clay to mask variables
clay = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/clay.asc")

####################
##### bio4 #########
####################
bio4 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio4.asc")
res(bio4)#it is the sd of temperature (remember temepratures en wc 1.4 are multiplides by 10) multiplides by 100, so we add 1000. Therefore 9000 would be a sd of 9 centigrade degrees.

#Extract bio4 from species distribution
if(FALSE){
require(raster)
median_bio4 = NULL
sd_bio4 = NULL
se_bio4 = NULL
bio4_pines = stack() 
for (i in 1:length(epithet_species_list)){

    #select the species
    selected_epithet = epithet_species_list[i]

    #load the distribution buffer
    distri = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/", selected_epithet, "_distribution_buffer.asc", sep=""))

    #create a polygon from distributon
    polygon = rasterToPolygons(distri, fun=function(x){x==1}, n=16, dissolve=TRUE) #esta funcion de raster te transforma un raster completo o una parte del mismo en un poliogno. En nuestro caso solo queremos las celdas con valor=1, es decir, presencias. Por eso ponemos x==1. 

    #mask bio4 using the buffer to remove areas outside the extended distribution
    bio4_cropped = mask(bio4, polygon)

    #remove areas without soil
    bio4_cropped = mask(bio4_cropped, clay)    

    #save the raster in the empty stack
    bio4_pines = stack(bio4_pines, bio4_cropped) 
    names(bio4_pines)[i] <- paste("bio4_", selected_epithet, sep="")

    #extract humidity data and calculate the median
    median_bio4 = append(median_bio4, median(na.omit(getValues(bio4_cropped))))

    #extract sd
    sd_value = sd(na.omit(getValues(bio4_cropped)))

    #extract se
    se_value = sd_value/sqrt(length(na.omit(getValues(bio4_cropped))))

    #save both of them
    sd_bio4 = append(sd_bio4, sd_value)
    se_bio4 = append(se_bio4, se_value)
}

#check that that data from all species has been extracted
nlayers(bio4_pines) == 112 #112, on for each species
length(median_bio4) == 112
length(sd_bio4) == 112
length(se_bio4) == 112

#check that species names are ok
paste("bio4_", epithet_species_list, sep="") == names(bio4_pines)

#bind species names and medians
medians_bio4 = cbind.data.frame(epithet_species_list, median_bio4)
colnames(medians_bio4) <- c("species", "median_bio4")

#bind species and ranges
standard_dev_bio4 = cbind.data.frame(epithet_species_list, sd_bio4)
colnames(standard_dev_bio4) <- c("species", "sd_bio4")
standard_error_bio4 = cbind.data.frame(epithet_species_list, se_bio4)
colnames(standard_error_bio4) <- c("species", "se_bio4")

#save medians
write.table(medians_bio4, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_medians/medians_bio4.csv", sep=",", col.names = TRUE, row.names = FALSE)

#save sd
write.table(standard_dev_bio4, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/sd_bio4.txt", sep="\t", col.names = TRUE, row.names = FALSE)

#save se
write.table(standard_error_bio4, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/se_bio4.txt", sep="\t", col.names = TRUE, row.names = FALSE)

#plot
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_medians/plots_climate_ranges_bio4.pdf")
for(i in 1:length(epithet_species_list)){

    #select the species
    selected_epithet = epithet_species_list[i]

    #load distribution
    distri = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/", selected_epithet, "_distribution_buffer.asc", sep=""))

    #create a polygon from distributon
    polygon = rasterToPolygons(distri, fun=function(x){x==1}, n=16, dissolve=TRUE) #esta funcion de raster te transforma un raster completo o una parte del mismo en un poliogno. En nuestro caso solo queremos las celdas con valor=1, es decir, presencias. Por eso ponemos x==1. 

    #extract the raster of the [i] species
    raster_bio4 = bio4_pines[[which(names(bio4_pines) == paste("bio4_", selected_epithet, sep=""))]]    

    #crop the raster to improve the visualization
    raster_bio4 = crop(raster_bio4, polygon)

    #plot
    plot(raster_bio4, main=selected_epithet)
}
dev.off()

#save rasters
#writeRaster(bio4_pines, filename="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/phlyo/climatic_data_phylo/bio4_buffer.tif", options="INTERLEAVE=BAND", overwrite=TRUE)



####################
##### bio17 ########
####################
bio17 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio17.asc")
res(bio17)

#Extract bio17 from species distribution
require(raster)
median_bio17 = NULL
sd_bio17 = NULL
se_bio17 = NULL
bio17_pines = stack() 
for (i in 1:length(epithet_species_list)){

    #select the species
    selected_epithet = epithet_species_list[i]

    #load the distribution buffer
    distri = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/", selected_epithet, "_distribution_buffer.asc", sep=""))

    #create a polygon from distributon
    polygon = rasterToPolygons(distri, fun=function(x){x==1}, n=16, dissolve=TRUE) #esta funcion de raster te transforma un raster completo o una parte del mismo en un poliogno. En nuestro caso solo queremos las celdas con valor=1, es decir, presencias. Por eso ponemos x==1. 

    #mask bio4 using the buffer to remove areas outside the extended distribution
    bio17_cropped = mask(bio17, polygon)

    #remove areas without soil
    bio17_cropped = mask(bio17_cropped, clay)  

    #save the raster in the empty stack
    bio17_pines = stack(bio17_pines, bio17_cropped) 
    names(bio17_pines)[i] <- paste("bio17_", selected_epithet, sep="")

    #extract humidity data and calculate the median
    median_bio17 = append(median_bio17, median(na.omit(getValues(bio17_cropped))))

    #extract sd
    sd_value = sd(na.omit(getValues(bio17_cropped)))

    #extract se
    se_value = sd_value/sqrt(length(na.omit(getValues(bio17_cropped))))

    #save both of them
    sd_bio17 = append(sd_bio17, sd_value)
    se_bio17 = append(se_bio17, se_value)
}

#check that that data from all species has been extracted
nlayers(bio17_pines) == 112 #112, on for each species
length(median_bio17) == 112
length(se_bio17) == 112
length(sd_bio17) == 112

#check that species names are ok
paste("bio17_", epithet_species_list, sep="") == names(bio17_pines)

#bind species names and clim data 
medians_bio17 = cbind.data.frame(epithet_species_list, median_bio17)
colnames(medians_bio17) <- c("species", "median_bio17")

#bind species and ranges
standard_dev_bio17 = cbind.data.frame(epithet_species_list, sd_bio17)
colnames(standard_dev_bio17) <- c("species", "sd_bio17")
standard_error_bio17 = cbind.data.frame(epithet_species_list, se_bio17)
colnames(standard_error_bio17) <- c("species", "se_bio17")

#save medians
write.table(medians_bio17, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_medians/medians_bio17.csv", sep=",", col.names = TRUE, row.names = FALSE)

#save sd
write.table(standard_dev_bio17, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/sd_bio17.txt", sep="\t", col.names = TRUE, row.names = FALSE)

#save se
write.table(standard_error_bio17, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/se_bio17.txt", sep="\t", col.names = TRUE, row.names = FALSE)

#plot
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_medians/plots_climate_ranges_bio17.pdf")
for(i in 1:length(epithet_species_list)){

    #select the species
    selected_epithet = epithet_species_list[i]

    #load distribution
    distri = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/", selected_epithet, "_distribution_buffer.asc", sep=""))

    #create a polygon from distributon
    polygon = rasterToPolygons(distri, fun=function(x){x==1}, n=16, dissolve=TRUE) #esta funcion de raster te transforma un raster completo o una parte del mismo en un poliogno. En nuestro caso solo queremos las celdas con valor=1, es decir, presencias. Por eso ponemos x==1. 

    #extract the raster of the [i] species
    raster_bio17 = bio17_pines[[which(names(bio17_pines) == paste("bio17_", selected_epithet, sep=""))]]    

    #crop the raster to improve the visualization
    raster_bio17 = crop(raster_bio17, polygon)

    #plot
    plot(raster_bio17, main=selected_epithet)
}
dev.off()

#save rasters
#writeRaster(bio17_pines, filename="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/phlyo/climatic_data_phylo/bio17_buffer.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
}

#########################################
##### bind both variables of p50 ########
#########################################

#load bio4 p50 data
medians_bio4 = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_medians/medians_bio4.csv", sep=",", header=TRUE)

#load bio17 p50 data
medians_bio17 = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_medians/medians_bio17.csv", sep=",", header=TRUE)

#merge
climate_medians = merge(medians_bio4, medians_bio17, by="species")

#save
write.table(climate_medians, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_medians/climate_medians.csv", col.names = TRUE, row.names = FALSE, sep=",") 
write.table(climate_medians, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_medians/climate_medians.txt", col.names = TRUE, row.names = FALSE, sep="\t") 


#########################################
##### bind both variables of ranges ########
#########################################

#sd
sd_bio4 = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/sd_bio4.txt", sep="\t", header=TRUE)
sd_bio17 = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/sd_bio17.txt", sep="\t", header=TRUE)

#se
se_bio4 = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/se_bio4.txt", sep="\t", header=TRUE)
se_bio17 = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/se_bio17.txt", sep="\t", header=TRUE)
#merge
climate_sd = merge(sd_bio4, sd_bio17, by="species")
climate_se = merge(se_bio4, se_bio17, by="species")

#save
#sd
write.table(climate_sd, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/climate_sd.csv", col.names = TRUE, row.names = FALSE, sep=",") 
write.table(climate_sd, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/climate_sd.txt", col.names = TRUE, row.names = FALSE, sep="\t") 
#se
write.table(climate_se, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/climate_se.csv", col.names = TRUE, row.names = FALSE, sep=",") 
write.table(climate_se, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/climate_se.txt", col.names = TRUE, row.names = FALSE, sep="\t") 


#######################################################
########## RECONSTRUCTION OF ANCESTRAL STATE ##########
#######################################################

## load required packages
require(ape)
require(phytools)
require(geiger)
require(diversitree)

#load climate data
climate_medians = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_medians/climate_medians.csv", header=TRUE, sep=",") 
str(climate_medians)

#load climate data
climate_sd = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/climate_sd.csv", header=TRUE, sep=",") 
str(climate_sd)
climate_se = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/climate_se.csv", header=TRUE, sep=",")
str(climate_se)

## cargamos el arbol
tree<-read.nexus("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/phlyo/phylogeny/FBDl_MCC_commAnc.tre") 

## new species icnldued by bianca in tree that we have to drop, and also discolor
species_to_drop = tree$tip.label[which(!tree$tip.label %in% paste("Pinus_", climate_medians$species, sep=""))]

## prune the tree of speceis without seed mass data
tree_prunned = drop.tip(tree, species_to_drop)

## check
species_to_drop %in% tree_prunned$tip.label

##reorder rows of climate data in basis on tip labels 
climate_medians = climate_medians[match(tree_prunned$tip.label, paste("Pinus_", climate_medians$species,sep="")),]
climate_sd = climate_sd[match(tree_prunned$tip.label, paste("Pinus_", climate_sd$species,sep="")),]
climate_se = climate_se[match(tree_prunned$tip.label, paste("Pinus_", climate_se$species,sep="")),]

## save climatic variables in a vector with species names as names
bio4_vector = climate_medians$median_bio4
bio17_vector = climate_medians$median_bio17

## set names of these variables as species names
names(bio4_vector) <- paste("Pinus_", climate_medians$species, sep="")
names(bio17_vector) <- paste("Pinus_", climate_medians$species, sep="")

##intra variability as SE of all data across distribution
intra_var_bio4 = climate_se$se_bio4
names(intra_var_bio4) <- paste("Pinus_", climate_se$species, sep="")
intra_var_bio17 = climate_se$se_bio17
names(intra_var_bio17) <- paste("Pinus_", climate_se$species, sep="")

##check order
names(bio4_vector) == tree_prunned$tip.label
names(bio17_vector) == tree_prunned$tip.label
names(intra_var_bio4) == tree_prunned$tip.label
names(intra_var_bio17) == tree_prunned$tip.label

############################################
#### Correlation between bio17 and bio4 ####
############################################

#correlation without PIC
test_no_pic = cor.test(bio4_vector, bio17_vector, method="spearman")

#correlation with PIC
pic_bio4 <- pic(bio4_vector, tree_prunned)
pic.bio17 <- pic(bio17_vector, tree_prunned) 
test_pic = cor.test(pic.bio17, pic_bio4, method="spearman")

#plot both variables
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_medians/cors_bio4_bio17.pdf")
par(mfrow=c(2,2))

#NO pic
plot(bio17_vector~bio4_vector, climate_medians, xlab="Median BIO4", ylab="Median BIO17")
estimate_no_pic = bquote(italic(rho) == .(format(round(test_no_pic$estimate,2))))
text(x=14000, y=-520, labels = estimate_no_pic, cex=1)
p_no_pic = bquote(italic(p.value) == .(format(round(test_no_pic$p.value,4))))
text(x=14000, y=-600, labels = p_no_pic, cex=1)

#YES PIC
plot(pic.bio17~pic_bio4, xlab="PIC Median BIO4", ylab="PIC Median BIO17")
estimate_pic = bquote(italic(rho) == .(format(round(test_pic$estimate,2))))
text(x=-1500, y=-100, labels = estimate_pic, cex=1)
p_no_pic = bquote(italic(p.value) == .(format(round(test_pic$p.value,4))))
text(x=-1500, y=-120, labels = p_no_pic, cex=1)

dev.off()

####################################
#### extintio - speciaiton rate ####
###################################
b_d = birthdeath(tree_prunned) #0.242418
b_d

####################################
#### times nodes (height) ####
###################################
branching.times(tree_prunned)

############################
#### Señal Filogenética ####
############################

## bio4
#lambda
phylosig(tree_prunned, bio4_vector, method="lambda", test=TRUE) #phytools without intravariability: 0.59, P=0.00831993
phylosig(tree_prunned, bio4_vector, method="lambda", test=TRUE, se=intra_var_bio4, nsim=6000) #phytools with intravariability: 0.59, P=0.008320014
fitContinuous(phy = tree_prunned, dat = bio4_vector, model = "lambda", control = list(niter = 100, CI = 0.95), SE=0) #diversitree without intravariability: 0.59
fitContinuous(phy = tree_prunned, dat = bio4_vector, model = "lambda", control = list(niter = 100, CI = 0.95), SE=intra_var_bio4) #diversitree with intravariability: 0.44

#k
phylosig(tree_prunned, bio4_vector, method="K", nsim=6000, test=TRUE) #0.1578271

#plot under BM
obj = contMap(tree_prunned, bio4_vector)
plot(obj, type="fan")

## bio17
#lambda 
phylosig(tree_prunned, bio17_vector, method="lambda", test=TRUE) #phytools without intravariability: 0.87; P=4.558942e-13
phylosig(tree_prunned, bio17_vector, method="lambda", test=TRUE, se=intra_var_bio17, nsim=6000) #phytools with intravariability: 0.87; P=4.550917e-13
fitContinuous(phy = tree_prunned, dat = bio17_vector, model = "lambda", control = list(niter = 100, CI = 0.95), SE=0) #diversitree without intravariability: 0.87
fitContinuous(phy = tree_prunned, dat = bio17_vector, model = "lambda", control = list(niter = 100, CI = 0.95), SE=intra_var_bio17) #diversitree with intravariability: 0.87

#k
phylosig(tree_prunned, bio17_vector, method="K", nsim=6000, test=TRUE) #0.2537013

#plot under BM
obj = contMap(tree_prunned, bio17_vector)
plot(obj, type="fan")

#########################################
####### Comparison AIC of evo models ####
#########################################

## CREATE a FUNCTION with fitContinuous (GEIGER) to COMPARE MODELS white, BM and OU models considering the existence and estimation of intrapespecific variation (SE = NA) and considenrg no existence of intra variation (SE = 0). 
fitClim=function(trait=c("median_bio4","median_bio17")){

    # select the trait
    trait=match.arg(trait, c("median_bio4","median_bio17"))

    # define set of models to compare
    models=c("white", "BM", "OU", "lambda")
    summaries=c("white noise", "Brownian motion", "Ornstein-Uhlenbeck", "lambda")
    
    ###### ESTIMATING measurement error ######
    #Set the intraspecific variance (SE) as NA, with the purpose to be estimate. It is to say, it could exist intraspecific variance, but we don't know it. 
    
    #empty vectors to sabe lnl and aic
    aic.se=numeric(length(models)) 
    lnl.se=numeric(length(models))

    # extract values of the trait and add species names
    climate_variable = climate_medians[,trait] 
    names(climate_variable) = paste("Pinus_", climate_medians$species, sep="")

    #for each model
    for(m in 1:length(models)){

        #Print the name of the model
        cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep=""), models[m], " with SE *** \n", sep="")
        
        #extract SE of the variable form data within species
        intra_variation = climate_se[, which(colnames(climate_se) == paste("se_", strsplit(trait, split="_")[[1]][2], sep=""))]

        #set names of SE
        names(intra_variation) <- paste("Pinus_", climate_se$species, sep="")

        #fit the model
        tmp=fitContinuous(phy = tree_prunned, dat = climate_variable, SE=intra_variation, model=models[m], control = list(niter = 1000, CI = 0.95, method = c("subplex","L-BFGS-B")), bounds=list(SE=c(0,0.5)), ncores=2)
            #SE is for intraspecific variation. 
            #Beware: difficulty in finding the optimal solution is determined by an interaction between the nature and complexity of the likelihood space (which is data- and model-dependent) and the numerical optimizer used to explore the space. There is never a guarantee that the optimal solution is found, but using many random starting points (control$niter) and many optimization methods (control$method) will increase these odds.
            #Bounds for the relevant parameters of the fitted model may be given through the bounds argument. Bounds may be necessary (particularly under the OU model) if the likelihood surface is characterized by a long, flat ridge which can be exceedingly difficult for optimization methods. Several bounds can be given at a time (e.g., bounds=list(SE=c(0,0.1),alpha=c(0,1)) would constrain measurement error as well as the ’constraint’ parameter of the Ornstein-Uhlenbeck model). Default bounds under the different models are given below.
            #Models. We only use white, BM and OU, following Quintero & Wiens 2013. None of these models have assumptions that we violated (in constrast with drif for example -> Not valid for ultrametric trees) 
                #BM is the Brownian motion model (Felsenstein 1973), which assumes the correlation structure among trait values is proportional to the extent of shared ancestry for pairs of species. Default bounds on the rate parameter are sigsq=c(min=exp(-500),max=exp(100)). The same bounds are applied to all other models, which also estimate sigsq. 
                #OU is the Ornstein-Uhlenbeck model (Butler and King 2004), which fits a random walk with a central tendency with an attraction strength proportional to the parameter alpha. The OU model is called the hansen model in ouch, although the way the parameters are fit is slightly different here. Default bounds are alpha = c(min = exp(-500), max = exp(1)). 
                #white is a white-noise (non-phylogenetic) model, which assumes data come from a single normal distribution with no covariance structure among species. The variance parameter sigsq takes the same bounds defined under the BM model

        #print results
        print(tmp)

        #save aic and lnL        
        aic.se[m]=tmp$opt$aicc
        lnl.se[m]=tmp$opt$lnL
    }

    ###### ASSUMING no measurement error ######
    #Assuming that there is NO intraespecific variance (SE = 0).
    
    #empty vectors to sabe lnl and aic
    aic=numeric(length(models))
    lnl=numeric(length(models))
    
    #for each model
    for(m in 1:length(models)){

        #Print the name of the model
        cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep=""), models[m], " *** \n", sep="")
        
        #fit the model
        tmp=fitContinuous(phy = tree_prunned, dat = climate_variable, SE=0, model=models[m], control = list(niter = 1000, CI = 0.95, method = c("subplex","L-BFGS-B")), ncores=2) #SE = 0: NO intraespecific variance
        
        #print results
        print(tmp)

        #save aic and lnL
        aic[m]=tmp$opt$aicc
        lnl[m]=tmp$opt$lnL
    }

    ###### COMPARE AIC ######
    #set the names of aic and lnL vectors
    names(aic.se)<-names(lnl.se)<-names(aic)<-names(lnl)<-models

    #create a function to calculate differences of aic
    delta_aic<-function(x) x-x[which(x==min(x))] #This function calcualate the difference in AIC of all models respect the models with lower AIC (min(x))
    
    # delta AIC with no measurement error
    daic=delta_aic(aic)
    cat("\n\n\n\t\t\t\t*** MODEL COMPARISON: ",trait," *** \n",sep="")
    cat("\tdelta-AIC values for models assuming no measurement error \t\t\t\t zero indicates the best model\n\n")
    print(daic, digits=2)
    
    # delta AIC with measurement error
    daic.se=delta_aic(aic.se)
    cat("\n\n\n\n\t\t\t\t*** MODEL COMPARISON: ",trait," ***\n",sep="")
    cat("\t\t delta-AIC values for models estimating SE \t\t\t\t zero indicates the best model\n\n")
    print(daic.se, digits=2)
    cat("\n\n\n")

    #bind all results 
    res_aicc=rbind(aic, aic.se, daic, daic.se)

    #set names
    rownames(res_aicc)=c("AICc","AICc_SE","dAICc", "dAICc_SE")

    #return those results
    return(res_aicc)
} #function taken from page 19 of "https://cran.r-project.org/web/packages/geiger/geiger.pdf"


## bio4 
res_bio4=fitClim(trait = "median_bio4")
print(res_bio4) #Best model is OU with a difference of AICc of 14 with the second best model, which is lambda (with and without SE (intraespecific variation)). The last model is BM but close to white noise (11) and lambda (16). 
    #alpha = 0.075742

#OU with SE from my maps
ou_bio4_se_mine = fitContinuous(phy = tree_prunned, dat = bio4_vector, SE=intra_var_bio4, model="OU", control = list(niter = 100, CI = 0.95, method = c("subplex","L-BFGS-B")), bounds=list(SE=c(0,0.5)), ncores=2) #alpha =  0.075744. Similar alpha and AIC than the OU models fitted with SE=NA and SE=0


## bio17
res_bio17=fitClim(trait = "median_bio17")
print(res_bio17) #Best model is OU with a difference of AICc of 5 with the second best model, which is lambda (with and without SE (intraespecific variation)). The next model is BM with a difference of 16 (10 respect to lambda) The last model is white noise, very far away (56).
    #alpha = 0.036495
#alpha values are equal with and without SE of intraespecific variability.

#OU with SE from my maps
ou_bio17_se_mine = fitContinuous(phy = tree_prunned, dat = bio17_vector, SE=intra_var_bio17, model="OU", control = list(niter = 100, CI = 0.95, method = c("subplex","L-BFGS-B")), bounds=list(SE=c(0,0.5)), ncores=2) #alpha = 0.036466. Similar alpha and AIC than the OU models fitted with SE=NA and SE=0


##TABLE 2 with bio4 and bio17 AICc from differents models. 

#bind results of both variables
table_2 = rbind.data.frame(res_bio4, res_bio17)

#drop row names with AICc, dAICc...
row.names(table_2) <- NULL

#create a vector with parameters names. This two-step phase is needed becasue binding two dataste created repetead row.names with .1. Two rows cannot have the same name. 
row_names =c(
    "AICc",
    "AICc SE",
    "\\textdelta AICc",
    "\\textdelta AICc SE")

#create a table with names of parameters and results
table_2 = cbind.data.frame(
    rep(row_names, 2),
    table_2)

#create a column for the variable
table_2$variable <- NA

#give it values
table_2$variable[1] <- "BIO4"
table_2$variable[5] <- "BIO17"

#reorder columns
table_2 = table_2[,c(6,1,2,5,3,4)]

#set final names
colnames(table_2) <- c("Variables", "Parameter", "White noise", "Lambda", "Brownian Motion", "Ornstein-Uhlenbeck")

#write to excel
write.table(table_2, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/tables_figures/table_2.csv", sep=",", row.names=FALSE, col.names=TRUE)

#convert to a latex table
require(xtable)
print.xtable(xtable(table_2, align="lllcccc"), include.rownames=FALSE, NA.string="", floating = FALSE, sanitize.text.function=function(x) {x})

#Copia la tabla en tables_latex.tex y corre este comando
system("cd /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/phylo_reconstruction/tables_figures; pandoc -s table2_latex.tex -o table2_latex.doc")

################################################
############  PHYLOGENETIC MONTECARLO ##########
################################################

#Una de las propuestas de Cooper et al 2016 es hacer bootstrap paramétricos simulando datos bajo los diferentes modelos de evolución a comparar. Vamos a usar un "phylogenetic Montecarlo". 

#Package for phylogenetic montecarlo
require(pmc) #Con este paquete vamos a correr un "Phylogenetic Monte Carlo" para la selección entre dos modelos: ModeloA más simple y modeloB más complejo. Primero, los parámetros de ambos modelos son estimados usando los datos originales. Entonces, se simulan n datasets (1000 en nuestro caso como hacen en el paper del paquete) siguiendo la evolución que dicta cada modelo con los parámetros que hemos obtenido previamente. Así obtrendemos un rasgo que evoluciona en nuestra filogenia bajo un modelo BM con el mismo sigma que el de WP, y otro rasgo que evoluciona en nuestra filogenia con el mismo signa, theta y alpha que el de WP. Con cada rasgo se reestiman los paramatroes de los DOS modelos (BM y OU), es decir, el rasgo que evoluciona segun BM se usa para ajustar un BM y un OU, mientras que el rasgo que evoluciona bajo OU se usa para ajustar otro BM y otro OU. Se hace un likelihood ratio test entre ambos modelos para cada rasgo, es decir obtendríamos dos valores de LRT en cada simulación: El valor null, ó hipotesis nula, que sería la diferencia de likelihood entre modelo BM y OU asjutados con un rasgo simulado bajo BM; y el valor test ó nuestro test de interés, que sería la diferencia de likelihood entre modelo BM y OU ambos ajustados con un rasgo que sigue evolución OU. Cuanto mayor sea LRT, más apoyo para el segundo modelo, más complejo, OU en nuestro caso. Esperaríamos, que si el rasgo evouciona por OU y no hay un bias del árbol a favor del modelo más complejo, el LRT del rasgo que evoluciona bajo BM será más bajo que el del rasgo que evoluciona bajo OU, es decir, el modelo OU no es mejor para el caso del rasgo que evoluciona bajo BM. Al final tendremos una distribución de LRT bajo BM y OU que podremos comparar. Ojo al detalle que este procesi implica 4 ajusted por maxima verosimiltud (maximum likelihood), mientras que para un AIC solo usas dos. La figura 2 de "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS" explica muy bien esto. 

############################
#### Geiger - OU.1 BIO4 ####
############################
#OU con un optimo, en este caso no tenemos estado depedencia
require(geiger)

### extract WP data with species as row.names
dat = climate_medians[,which(colnames(climate_medians) %in% c("species", "median_bio4"))]
row.names(dat) <- paste("Pinus_", dat$species, sep="")
dat[,which(colnames(dat) == "species")] <- NULL
str(dat)

### bind tree and data into tmp
tmp = treedata(tree_prunned, dat)

### extract phylogeny
phy = tmp$phy

### extract WP data with species as row.names
datos = tmp$data

### run the phylo montecarlo
simulations_ou.1_bio4 = pmc(phy, datos, "BM", "OU", nboot=100, mc.cores=3)

### save simulations
save(simulations_ou.1_bio4, file="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/rdata/bio_4_geiger_BM_OU.1_nboot_100.rda")

### load it
load("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/rdata/bio_4_geiger_BM_OU.1_nboot_2000.rda")

### plot likelihood ratio test between BM (modelA) and OU.1 (modelB) models fitted with data simulated under BM and OU.1 respectively
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/pmc/figures/bio_4_geiger_BM_OU_1_nboot_2000.pdf")
par(mfrow=c(2,2))

## alpha values for data simulated under OU.1
#comparison BB: BM vs OU.1 model fitted with OU.1-simulated data and the value of alpha (selection strength)
lr_ou = simulations_ou.1_bio4$test

#calculate density distribution
density_lr_OU = density(lr_ou)

#plot
plot(density_lr_OU, xlim=c(-0.5,75), ylim=c(0,0.9), xlab="Likelihood ratio test", cex.lab = 1.3, main="")

#add color to the full area under the curve
x1 <- min(which(density_lr_OU$x >= 0))  
x2 <- max(which(density_lr_OU$x <  75))
with(density_lr_OU, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("black",alpha.f=0.65)))

#plot 95CI
#quantiles_ou = quantile(lr_ou, prob=c(0.025, 0.975))
#abline(v=c(quantiles_ou[1], quantiles_ou[2]), col="black", lty=2) #No los añado porque podría haber problemas para calcularlos así y no encuentro la función de pmc para ello: "Given the noisy nature of parameters estimated from phylogenies, we recommend that confidence interval should routinely be reported, and to facilitate this, have implemented this as pmc::confidenceIntervals.pow. Confidence intervals could also be estimated from the curvature of the likelihood surface, but these can be unreliable and problematic to compute." From "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS". 


## alpha values for data simulated under BM
#comparison AB: BM vs OU.1 model fitted with BM-simulated data and the value of alpha (selection strength)
lr_bm = simulations_ou.1_bio4$null

#calculate density distribution
density_lr_BM = density(lr_bm)

#add density plot to the previous plot
lines(density_lr_BM)

#add color to the full area under the curve
x1 <- min(which(density_lr_BM$x >= -0.8))  
x2 <- max(which(density_lr_BM$x <  13))
with(density_lr_BM, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("gray",alpha.f=0.65)))

#plot 95CI
#quantiles_BM = quantile(lr_bm, prob=c(0.025, 0.975))
#abline(v=c(quantiles_BM[1], quantiles_BM[2]), col="black", lty=2)

#add legen
#legend(53, 0.85, legend=c("BM", "OU.1"), fill=c("gray", adjustcolor("black",alpha.f=0.65)))

#add LRT value obtained from real data
abline(v=c(simulations_ou.1_bio4$lr), col="black", lty=5, lwd=2) #Este valor tiene que estar por encima del percentil 95 del LRT del modeloA (null) y por debajo del percentil 95 del LRT del modeloB (test). Mira la sección Methods-Model selection en "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS".
#add 95% LR test
#abline(v=c(quantile(lr_ou, probs=0.95)), col="black", lty=1)
#abline(v=c(quantile(lr_bm, probs=0.95)), col="black", lty=1) #La proporción de valores simulados mayore que LRT de los datos reales nos da una especie de P.value para el test de selección de modelos. La probabilidad de que la diferencia de likelihood entre modelos (LRT) observada se de bajo el modelo cero. Si cogemos el valor de LRT que es mayor que el 95% de los LRT bajo modeloA, si el LRT observado es mayor que ese threshold, entonces es poco probable que el valor obtenido se haya dado bajo el modeloA (más simple). También podemos calcular el poder del test, la probabilidad de acertar rechazando el modeloA porque los datos vienen del B. Para eso tenemos que usar la distribución de LRTs bajo modelo 1. Si como antes cogemos el valor de LRT mayor que el 95% de LRT simulados bajo modeloB, la cantidad de distribución que queda a la izquierda de ese threshold se aproxima a la probabiidad de rechaza el mdoeloA cuando los datos son producidos por el modeloB. 


### plot alpha between data simulated under BM (modelA) and OU.1 (modelB)

## alpha values for data simulated under OU.1
#comparison BB (BM vs OU.1 model fitted with OU.1-simulated data) and the value of alpha (selection strength)
alpha_ou = simulations_ou.1_bio4$par_dists$value[which(simulations_ou.1_bio4$par_dists$comparison == "BB" & simulations_ou.1_bio4$par_dists$parameter == "alpha")] 

#calculate density distribution
density_OU = density(alpha_ou)

#plot
plot(density_OU, xlim=c(0,0.39), ylim=c(0,140), , xlab=expression(alpha), cex.lab = 1.3, main="")

#add color to the full area under the curve
x1 <- min(which(density_OU$x >= 0))  
x2 <- max(which(density_OU$x <  0.2))
with(density_OU, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("black",alpha.f=0.65)))

#plot 95CI
#quantiles_ou = quantile(alpha_ou, prob=c(0.025, 0.975))
#abline(v=c(quantiles_ou[1], quantiles_ou[2]), col="black", lty=2) #No los añado porque podría haber problemas para calcularlos así y no encuentro la función de pmc para ello: "Given the noisy nature of parameters estimated from phylogenies, we recommend that confidence interval should routinely be reported, and to facilitate this, have implemented this as pmc::confidenceIntervals.pow. Confidence intervals could also be estimated from the curvature of the likelihood surface, but these can be unreliable and problematic to compute." From "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS". 

## alpha values for data simulated under BM
#comparison AB (BM vs OU.1 model fitted with BM-simulated data) and the value of alpha (selection strength)
alpha_BM = simulations_ou.1_bio4$par_dists$value[which(simulations_ou.1_bio4$par_dists$comparison == "AB" & simulations_ou.1_bio4$par_dists$parameter == "alpha")]

#calculate density distribution
density_BM = density(alpha_BM)

#add density plot to the previous plot
lines(density_BM)

#add color to the full area under the curve
x1 <- min(which(density_BM$x >= -0.2))  
x2 <- max(which(density_BM$x <  0.2))
with(density_BM, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("gray",alpha.f=0.65)))

#plot 95CI
#quantiles_BM = quantile(alpha_bm, prob=c(0.025, 0.975))
#abline(v=c(quantiles_BM[1], quantiles_BM[2]), col="black", lty=2)

#add legen
legend(0.27, 133, legend=c("BM", "OU"), fill=c("gray", adjustcolor("black",alpha.f=0.65)))

#add alpha value obtained from the real data
abline(v=c(simulations_ou.1_bio4$B$opt$alpha), col="black", lty=5, lwd=2) #Este valor tiene que estar por encima del percentil 95 del LRT del modeloA (null) y por debajo del percentil 95 del LRT del modeloA (test). Mira la sección Methods-Model selection en "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS".
#add 95% alpha
#abline(v=c(quantile(alpha_ou, probs=0.95)), col="black", lty=1)
#abline(v=c(quantile(alpha_BM, probs=0.95)), col="black", lty=1) #Quito estas lineas porque no se si tiene sentido aplicar el approach de poner el 95% de cada parametros en cada modelo, el paper de coop usan este approach con el likelihood ratio test nada más. 
dev.off()


#############################
#### Geiger - OU.1 BIO17 ####
#############################
#OU con un optimo, en este caso no tenemos estado depedencia
require(geiger)

### extract WP data with species as row.names
dat = climate_medians[,which(colnames(climate_medians) %in% c("species", "median_bio17"))]
row.names(dat) <- paste("Pinus_", dat$species, sep="")
dat[,which(colnames(dat) == "species")] <- NULL
str(dat)

### bind tree and data into tmp
tmp = treedata(tree_prunned, dat)

### extract phylogeny
phy = tmp$phy

### extract WP data with species as row.names
datos = tmp$data

### run the phylo montecarlo
simulations_ou.1_bio17 = pmc(phy, datos, "BM", "OU", nboot=100, mc.cores=3)

### save simulations
save(simulations_ou.1_bio17, file="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/rdata/bio_17_geiger_BM_OU.1_nboot_100.rda")

### load it
load("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/rdata/bio_17_geiger_BM_OU.1_nboot_2000.rda")


### plot likelihood ratio test between BM (modelA) and OU.1 (modelB) models fitted with data simulated under BM and OU.1 respectively
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/pmc/figures/bio_17_geiger_BM_OU_1_nboot_2000.pdf")
par(mfrow=c(2,2))

## alpha values for data simulated under OU.1
#comparison BB (BM vs OU.1 model fitted with OU.1-simulated data) and the value of alpha (selection strength)
lr_ou = simulations_ou.1_bio17$test

#calculate density distribution
density_lr_OU = density(lr_ou)

#plot
plot(density_lr_OU, xlim=c(-0.8,45), ylim=c(0,0.9), xlab="Likelihood ratio test", cex.lab = 1.3, main="")

#add color to the full area under the curve
x1 <- min(which(density_lr_OU$x >= -0.8))  
x2 <- max(which(density_lr_OU$x <  75))
with(density_lr_OU, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("black",alpha.f=0.65)))

#plot 95CI
#quantiles_ou = quantile(lr_ou, prob=c(0.025, 0.975))
#abline(v=c(quantiles_ou[1], quantiles_ou[2]), col="black", lty=2) #No los añado porque podría haber problemas para calcularlos así y no encuentro la función de pmc para ello: "Given the noisy nature of parameters estimated from phylogenies, we recommend that confidence interval should routinely be reported, and to facilitate this, have implemented this as pmc::confidenceIntervals.pow. Confidence intervals could also be estimated from the curvature of the likelihood surface, but these can be unreliable and problematic to compute." From "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS". 


## alpha values for data simulated under BM
#comparison AB (BM vs OU.1 model fitted with BM-simulated data) and the value of alpha (selection strength)
lr_bm = simulations_ou.1_bio17$null

#calculate density distribution
density_lr_BM = density(lr_bm)

#add density plot to the previous plot
lines(density_lr_BM)

#add color to the full area under the curve
x1 <- min(which(density_lr_BM$x >= -0.8))  
x2 <- max(which(density_lr_BM$x <  13))
with(density_lr_BM, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("gray",alpha.f=0.65)))

#plot 95CI
#quantiles_BM = quantile(lr_bm, prob=c(0.025, 0.975))
#abline(v=c(quantiles_BM[1], quantiles_BM[2]), col="black", lty=2)

#add legen
#legend(53, 0.85, legend=c("BM", "OU.1"), fill=c("gray", adjustcolor("black",alpha.f=0.65)))

#value obtained from real data
abline(v=c(simulations_ou.1_bio17$lr), col="black", lty=5, lwd=2) #Este valor tiene que estar por encima del percentil 95 del LRT del modeloA (null) y por debajo del percentil 95 del LRT del modeloB (test). Mira la sección Methods-Model selection en "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS".
#add 95% LR test
#abline(v=c(quantile(lr_ou, probs=0.95)), col="black", lty=1)
#abline(v=c(quantile(lr_bm, probs=0.95)), col="black", lty=1) #La proporción de valores simulados mayore que LRT de los datos reales nos da una especie de P.value para el test de selección de modelos. La probabilidad de que la diferencia de likelihood entre modelos (LRT) observada se de bajo el modelo cero. Si cogemos el valor de LRT que es mayor que el 95% de los LRT bajo modeloA, si el LRT observado es mayor que ese threshold, entonces es poco probable que el valor obtenido se haya dado bajo el modeloA (más simple). También podemos calcular el poder del test, la probabilidad de acertar rechazando el modeloA porque los datos vienen del B. Para eso tenemos que usar la distribución de LRTs bajo modelo 1. Si como antes cogemos el valor de LRT mayor que el 95% de LRT simulados bajo modeloB, la cantidad de distribución que queda a la izquierda de ese threshold se aproxima a la probabiidad de rechaza el mdoeloA cuando los datos son producidos por el modeloB. 


### plot alpha between data simulated under BM (modelA) and OU.1 (modelB)

## alpha values for data simulated under OU.1
#comparison BB (OU.1 model fitted with OU.1-simulated data) and the value of alpha (selection strength)
alpha_ou = simulations_ou.1_bio17$par_dists$value[which(simulations_ou.1_bio17$par_dists$comparison == "BB" & simulations_ou.1_bio17$par_dists$parameter == "alpha")] 

#calculate density distribution
density_OU = density(alpha_ou)

#plot
plot(density_OU, xlim=c(0,0.12), ylim=c(0,137), , xlab=expression(alpha), cex.lab = 1.3,, main="")

#add color to the full area under the curve
x1 <- min(which(density_OU$x >= 0))  
x2 <- max(which(density_OU$x <  0.2))
with(density_OU, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("black",alpha.f=0.65)))

#plot 95CI
#quantiles_ou = quantile(alpha_ou, prob=c(0.025, 0.975))
#abline(v=c(quantiles_ou[1], quantiles_ou[2]), col="black", lty=2) #No los añado porque podría haber problemas para calcularlos así y no encuentro la función de pmc para ello: "Given the noisy nature of parameters estimated from phylogenies, we recommend that confidence interval should routinely be reported, and to facilitate this, have implemented this as pmc::confidenceIntervals.pow. Confidence intervals could also be estimated from the curvature of the likelihood surface, but these can be unreliable and problematic to compute." From "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS". 

## alpha values for data simulated under BM
#comparison AB (BM vs OU.1 model fitted with BM-simulated data) and the value of alpha (selection strength)
alpha_BM = simulations_ou.1_bio17$par_dists$value[which(simulations_ou.1_bio17$par_dists$comparison == "AB" & simulations_ou.1_bio17$par_dists$parameter == "alpha")]

#calculate density distribution
density_BM = density(alpha_BM)

#add density plot to the previous plot
lines(density_BM)

#add color to the full area under the curve
x1 <- min(which(density_BM$x >= -0.2))  
x2 <- max(which(density_BM$x <  0.2))
with(density_BM, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("gray",alpha.f=0.65)))

#plot 95CI
#quantiles_BM = quantile(alpha_bm, prob=c(0.025, 0.975))
#abline(v=c(quantiles_BM[1], quantiles_BM[2]), col="black", lty=2)

#add legen
legend(0.085, 130, legend=c("BM", "OU"), fill=c("gray", adjustcolor("black",alpha.f=0.65)))

#add alpha vlaue form real data
abline(v=c(simulations_ou.1_bio17$B$opt$alpha), col="black", lty=5, lwd=2) #Este valor tiene que estar por encima del percentil 95 del LRT del modeloA (null) y por debajo del percentil 95 del LRT del modeloA (test). Mira la sección Methods-Model selection en "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS".
#add 95% alpha
#abline(v=c(quantile(alpha_ou, probs=0.95)), col="black", lty=1)
#abline(v=c(quantile(alpha_BM, probs=0.95)), col="black", lty=1) #Quito estas lineas porque no se si tiene sentido aplicar el approach de poner el 95% de cada parametros en cada modelo, el paper de coop usan este approach con el likelihood ratio test nada más. 
dev.off()

####### FINAL Figure for the paper ########
### plot likelihood ratio test between BM (modelA) and OU.1 (modelB) models fitted with data simulated under BM and OU.1 respectively
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/pmc/figures/bio_4_17_geiger_BM_OU_1_nboot_2000.pdf")
par(mfrow=c(2,2))

#######BIO4
### load it
load("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/rdata/bio_4_geiger_BM_OU.1_nboot_2000.rda")

## alpha values for data simulated under OU.1
#comparison BB (OU.1 model fitted with OU.1-simulated data) and the value of alpha (selection strength)
lr_ou = simulations_ou.1_bio4$test

#calculate density distribution
density_lr_OU = density(lr_ou)

#plot
plot(density_lr_OU, xlim=c(-0.5,75), ylim=c(0,0.9), xlab="Likelihood ratio test", cex.lab = 1.3, main="")

#add title
mtext(text=expression(bold("Temperature Seasonality  (BIO4)")), side=3, at=95, line=1.4)

#add color to the full area under the curve
x1 <- min(which(density_lr_OU$x >= 0))  
x2 <- max(which(density_lr_OU$x <  75))
with(density_lr_OU, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("black",alpha.f=0.65)))

#plot 95CI
#quantiles_ou = quantile(lr_ou, prob=c(0.025, 0.975))
#abline(v=c(quantiles_ou[1], quantiles_ou[2]), col="black", lty=2) #No los añado porque podría haber problemas para calcularlos así y no encuentro la función de pmc para ello: "Given the noisy nature of parameters estimated from phylogenies, we recommend that confidence interval should routinely be reported, and to facilitate this, have implemented this as pmc::confidenceIntervals.pow. Confidence intervals could also be estimated from the curvature of the likelihood surface, but these can be unreliable and problematic to compute." From "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS". 


## alpha values for data simulated under BM
#comparison AB (OU.1 model fitted with BM-simulated data) and the value of alpha (selection strength)
lr_bm = simulations_ou.1_bio4$null

#calculate density distribution
density_lr_BM = density(lr_bm)

#add density plot to the previous plot
lines(density_lr_BM)

#add color to the full area under the curve
x1 <- min(which(density_lr_BM$x >= -0.8))  
x2 <- max(which(density_lr_BM$x <  13))
with(density_lr_BM, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("gray",alpha.f=0.65)))

#plot 95CI
#quantiles_BM = quantile(lr_bm, prob=c(0.025, 0.975))
#abline(v=c(quantiles_BM[1], quantiles_BM[2]), col="black", lty=2)

#add legen
#legend(53, 0.85, legend=c("BM", "OU.1"), fill=c("gray", adjustcolor("black",alpha.f=0.65)))

#add 95% LR test
abline(v=c(simulations_ou.1_bio4$lr), col="black", lty=5, lwd=2) #Este valor tiene que estar por encima del percentil 95 del LRT del modeloA (null) y por debajo del percentil 95 del LRT del modeloB (test). Mira la sección Methods-Model selection en "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS".
#abline(v=c(quantile(lr_ou, probs=0.95)), col="black", lty=1)
#abline(v=c(quantile(lr_bm, probs=0.95)), col="black", lty=1) #La proporción de valores simulados mayore que LRT de los datos reales nos da una especie de P.value para el test de selección de modelos. La probabilidad de que la diferencia de likelihood entre modelos (LRT) observada se de bajo el modelo cero. Si cogemos el valor de LRT que es mayor que el 95% de los LRT bajo modeloA, si el LRT observado es mayor que ese threshold, entonces es poco probable que el valor obtenido se haya dado bajo el modeloA (más simple). También podemos calcular el poder del test, la probabilidad de acertar rechazando el modeloA porque los datos vienen del B. Para eso tenemos que usar la distribución de LRTs bajo modelo 1. Si como antes cogemos el valor de LRT mayor que el 95% de LRT simulados bajo modeloB, la cantidad de distribución que queda a la izquierda de ese threshold se aproxima a la probabiidad de rechaza el mdoeloA cuando los datos son producidos por el modeloB. 


### plot alpha between data simulated under BM (modelA) and OU.1 (modelB)

## alpha values for data simulated under OU.1
#comparison BB (OU.1 model fitted with OU.1-simulated data) and the value of alpha (selection strength)
alpha_ou = simulations_ou.1_bio4$par_dists$value[which(simulations_ou.1_bio4$par_dists$comparison == "BB" & simulations_ou.1_bio4$par_dists$parameter == "alpha")] 

#calculate density distribution
density_OU = density(alpha_ou)

#plot
plot(density_OU, xlim=c(0,0.39), ylim=c(0,140), , xlab=expression(alpha), cex.lab = 1.3, main="")

#add color to the full area under the curve
x1 <- min(which(density_OU$x >= 0))  
x2 <- max(which(density_OU$x <  0.2))
with(density_OU, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("black",alpha.f=0.65)))

#plot 95CI
#quantiles_ou = quantile(alpha_ou, prob=c(0.025, 0.975))
#abline(v=c(quantiles_ou[1], quantiles_ou[2]), col="black", lty=2) #No los añado porque podría haber problemas para calcularlos así y no encuentro la función de pmc para ello: "Given the noisy nature of parameters estimated from phylogenies, we recommend that confidence interval should routinely be reported, and to facilitate this, have implemented this as pmc::confidenceIntervals.pow. Confidence intervals could also be estimated from the curvature of the likelihood surface, but these can be unreliable and problematic to compute." From "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS". 

## alpha values for data simulated under BM
#comparison AB (OU.1 model fitted with BM-simulated data) and the value of alpha (selection strength)
alpha_BM = simulations_ou.1_bio4$par_dists$value[which(simulations_ou.1_bio4$par_dists$comparison == "AB" & simulations_ou.1_bio4$par_dists$parameter == "alpha")]

#calculate density distribution
density_BM = density(alpha_BM)

#add density plot to the previous plot
lines(density_BM)

#add color to the full area under the curve
x1 <- min(which(density_BM$x >= -0.2))  
x2 <- max(which(density_BM$x <  0.2))
with(density_BM, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("gray",alpha.f=0.65)))

#plot 95CI
#quantiles_BM = quantile(alpha_bm, prob=c(0.025, 0.975))
#abline(v=c(quantiles_BM[1], quantiles_BM[2]), col="black", lty=2)

#add legen
legend(0.27, 133, legend=c("BM", "OU"), fill=c("gray", adjustcolor("black",alpha.f=0.65)))

#add 95% LR test
abline(v=c(simulations_ou.1_bio4$B$opt$alpha), col="black", lty=5, lwd=2) #Este valor tiene que estar por encima del percentil 95 del LRT del modeloA (null) y por debajo del percentil 95 del LRT del modeloA (test). Mira la sección Methods-Model selection en "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS".
#abline(v=c(quantile(alpha_ou, probs=0.95)), col="black", lty=1)
#abline(v=c(quantile(alpha_BM, probs=0.95)), col="black", lty=1) #Quito estas lineas porque no se si tiene sentido aplicar el approach de poner el 95% de cada parametros en cada modelo, el paper de coop usan este approach con el likelihood ratio test nada más. 


####BIO17 

### load it
load("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/rdata/bio_17_geiger_BM_OU.1_nboot_2000.rda")


## alpha values for data simulated under OU.1
#comparison BB (OU.1 model fitted with OU.1-simulated data) and the value of alpha (selection strength)
lr_ou = simulations_ou.1_bio17$test

#calculate density distribution
density_lr_OU = density(lr_ou)

#plot
plot(density_lr_OU, xlim=c(-0.8,45), ylim=c(0,0.9), xlab="Likelihood ratio test", cex.lab = 1.3, main="")

#add title
mtext(text=expression(bold("Moisture of Driest Quarter  (BIO17)")), side=3, at=57, line=1.4)

#add color to the full area under the curve
x1 <- min(which(density_lr_OU$x >= -0.8))  
x2 <- max(which(density_lr_OU$x <  75))
with(density_lr_OU, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("black",alpha.f=0.65)))

#plot 95CI
#quantiles_ou = quantile(lr_ou, prob=c(0.025, 0.975))
#abline(v=c(quantiles_ou[1], quantiles_ou[2]), col="black", lty=2) #No los añado porque podría haber problemas para calcularlos así y no encuentro la función de pmc para ello: "Given the noisy nature of parameters estimated from phylogenies, we recommend that confidence interval should routinely be reported, and to facilitate this, have implemented this as pmc::confidenceIntervals.pow. Confidence intervals could also be estimated from the curvature of the likelihood surface, but these can be unreliable and problematic to compute." From "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS". 


## alpha values for data simulated under BM
#comparison AB (OU.1 model fitted with BM-simulated data) and the value of alpha (selection strength)
lr_bm = simulations_ou.1_bio17$null

#calculate density distribution
density_lr_BM = density(lr_bm)

#add density plot to the previous plot
lines(density_lr_BM)

#add color to the full area under the curve
x1 <- min(which(density_lr_BM$x >= -0.8))  
x2 <- max(which(density_lr_BM$x <  13))
with(density_lr_BM, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("gray",alpha.f=0.65)))

#plot 95CI
#quantiles_BM = quantile(lr_bm, prob=c(0.025, 0.975))
#abline(v=c(quantiles_BM[1], quantiles_BM[2]), col="black", lty=2)

#add legen
#legend(53, 0.85, legend=c("BM", "OU.1"), fill=c("gray", adjustcolor("black",alpha.f=0.65)))

#add 95% LR test
abline(v=c(simulations_ou.1_bio17$lr), col="black", lty=5, lwd=2) #Este valor tiene que estar por encima del percentil 95 del LRT del modeloA (null) y por debajo del percentil 95 del LRT del modeloB (test). Mira la sección Methods-Model selection en "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS".
#abline(v=c(quantile(lr_ou, probs=0.95)), col="black", lty=1)
#abline(v=c(quantile(lr_bm, probs=0.95)), col="black", lty=1) #La proporción de valores simulados mayore que LRT de los datos reales nos da una especie de P.value para el test de selección de modelos. La probabilidad de que la diferencia de likelihood entre modelos (LRT) observada se de bajo el modelo cero. Si cogemos el valor de LRT que es mayor que el 95% de los LRT bajo modeloA, si el LRT observado es mayor que ese threshold, entonces es poco probable que el valor obtenido se haya dado bajo el modeloA (más simple). También podemos calcular el poder del test, la probabilidad de acertar rechazando el modeloA porque los datos vienen del B. Para eso tenemos que usar la distribución de LRTs bajo modelo 1. Si como antes cogemos el valor de LRT mayor que el 95% de LRT simulados bajo modeloB, la cantidad de distribución que queda a la izquierda de ese threshold se aproxima a la probabiidad de rechaza el mdoeloA cuando los datos son producidos por el modeloB. 


### plot alpha between data simulated under BM (modelA) and OU.1 (modelB)

## alpha values for data simulated under OU.1
#comparison BB (OU.1 model fitted with OU.1-simulated data) and the value of alpha (selection strength)
alpha_ou = simulations_ou.1_bio17$par_dists$value[which(simulations_ou.1_bio17$par_dists$comparison == "BB" & simulations_ou.1_bio17$par_dists$parameter == "alpha")] 

#calculate density distribution
density_OU = density(alpha_ou)

#plot
plot(density_OU, xlim=c(0,0.12), ylim=c(0,137), , xlab=expression(alpha), cex.lab = 1.3,, main="")

#add color to the full area under the curve
x1 <- min(which(density_OU$x >= 0))  
x2 <- max(which(density_OU$x <  0.2))
with(density_OU, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("black",alpha.f=0.65)))

#plot 95CI
#quantiles_ou = quantile(alpha_ou, prob=c(0.025, 0.975))
#abline(v=c(quantiles_ou[1], quantiles_ou[2]), col="black", lty=2) #No los añado porque podría haber problemas para calcularlos así y no encuentro la función de pmc para ello: "Given the noisy nature of parameters estimated from phylogenies, we recommend that confidence interval should routinely be reported, and to facilitate this, have implemented this as pmc::confidenceIntervals.pow. Confidence intervals could also be estimated from the curvature of the likelihood surface, but these can be unreliable and problematic to compute." From "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS". 

## alpha values for data simulated under BM
#comparison AB (OU.1 model fitted with BM-simulated data) and the value of alpha (selection strength)
alpha_BM = simulations_ou.1_bio17$par_dists$value[which(simulations_ou.1_bio17$par_dists$comparison == "AB" & simulations_ou.1_bio17$par_dists$parameter == "alpha")]

#calculate density distribution
density_BM = density(alpha_BM)

#add density plot to the previous plot
lines(density_BM)

#add color to the full area under the curve
x1 <- min(which(density_BM$x >= -0.2))  
x2 <- max(which(density_BM$x <  0.2))
with(density_BM, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("gray",alpha.f=0.65)))

#plot 95CI
#quantiles_BM = quantile(alpha_bm, prob=c(0.025, 0.975))
#abline(v=c(quantiles_BM[1], quantiles_BM[2]), col="black", lty=2)

#add legen
#legend(0.085, 130, legend=c("BM", "OU"), fill=c("gray", adjustcolor("black",alpha.f=0.65)))

#add 95% LR test
abline(v=c(simulations_ou.1_bio17$B$opt$alpha), col="black", lty=5, lwd=2) #Este valor tiene que estar por encima del percentil 95 del LRT del modeloA (null) y por debajo del percentil 95 del LRT del modeloA (test). Mira la sección Methods-Model selection en "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS".
#abline(v=c(quantile(alpha_ou, probs=0.95)), col="black", lty=1)
#abline(v=c(quantile(alpha_BM, probs=0.95)), col="black", lty=1) #Quito estas lineas porque no se si tiene sentido aplicar el approach de poner el 95% de cada parametros en cada modelo, el paper de coop usan este approach con el likelihood ratio test nada más. 
dev.off()



################################
#### phylogenetic half-life ####
################################

#max height of the tree
height_tree = max(nodeHeights(tree_prunned)[,2])

#t[1/2] for BIO4 and BIO17
t_1_2_bio4 = log(2)/0.076
t_1_2_bio17 = log(2)/0.036

#difference between t[1/2] and max height
t_1_2_bio4/height_tree
t_1_2_bio17/height_tree #El timepo que tarda una especie que entra en nu nuevo nicho en llegar a la mitad de camino hacia su nuevo óptimo esperado. Si ese valor es muy pequeo en relación con la altura del árbol, entonces la evolución hacia el rago optimo es rápida y las correlaciones filogeneitcas residuales son debiles y hay poca influencia de los valores antiguos del rasgo. A nosotros nos sale un valor muy bajo respecto de la altura del árbol, pero hay que tener en cuenta que la estrcutura del arbul afecta al valor de alfa, ciertas formas de árbol pueden dar lugar a valores más altos de alfa independientemente del rasgo que estés mirando, y t_1_2 se calcula con alfa, así que tenemos esos sesgos. Cooper et al 2016 vieron que con alfa menor de 1, BM y OU son prácticamente indiferenciables. 


#####################################
####### Extract ancestal state ######
#####################################

####### OU and BM in Comare 4.6 #########
#extract species names
species_names = paste("Pinus_", climate_medians$species, sep="")

#write species names
write.table(species_names, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/compare4.6_data/compare_species_names.txt", sep="\t", col.names = FALSE, row.names =  FALSE, quote=FALSE) #We don't want quotes because these files will be copied to compare. quote: a logical value (‘TRUE’ or ‘FALSE’) or a numeric vector.  If ‘TRUE’, any character or factor columns will be surrounded by double quotes.  If a numeric vector, its elements are taken as the indices of columns to quote.  In both cases, row and column names are quoted if they are written.  If ‘FALSE’, nothing is quoted.

#bind climatic data to SE estimates (0 in our cases for all speices because we have not intraspecies data)
climate_medians$species == climate_se$species
compare_bio4 = paste(climate_medians$median_bio4, "<", climate_se$se_bio4, ">", sep="")
compare_bio17 = paste(climate_medians$median_bio17, "<", climate_se$se_bio17, ">", sep="")

#write variables
write.table(compare_bio4, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/compare4.6_data/compare_bio4.txt", sep="\t", col.names = FALSE, row.names =  FALSE, quote=FALSE) #YOU HAVE TO DROP " FROM THE FILE
write.table(compare_bio17, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/compare4.6_data/compare_bio17.txt", sep="\t", col.names = FALSE, row.names =  FALSE, quote=FALSE) #YOU HAVE TO DROP " FROM THE FILE

#save the tree prunned for compare 4.6
write.tree(tree_prunned, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/compare4.6_data/FBDl_MCC_commAnc_prunned.tree") #You have to copy hasta el ";", no etas espacio, sino no funciona


#En "/Volumes/GoogleDrive/My Drive/science/software/Comp46bExec" pinchas startForm.class para activar compare. Le das a main window. 

#Copias todos estos archivos creados con al código de justo arriba en Taxon names (names), Taxon means (rasgos) y Enter Phylogeny. Indicas que son 112 taxa, 1 rasgo y 1 filogenia. Hay que indicar que si queremos usar los SE dentro de especie ó asumir que la variación dentro de especie es desconocida. En este caso incluímos los valores de SE como variabilidad intraespecífica. Seleccionar PGLS-ancestros y exeecute.   

#Correr modelo: 
    #BM: Linear model which is the option by default. Not select specyfing alpha. Nothing else. 
    #OU:Luego Exponential model, specyfing Alpha, pon el valor de alfa de bio4  ó bio17 redondeados a dos decimales (0.076 y 0.036 respectivamente), 100 iteraciones y run (he comprobado el resultado con 1000 interaciones en ambas variablws y sale exactamente lo mismo; REVISADO JULIO 2019). Asú se corre un OU en comapre4.6. Esto se ha seguido de Guerrero et al.,... Wiens ., 2013 ("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3710863/")

#Le das a guardar en un archivo llamado "compare_res_bioX_XX.txt". De ese archivo copias la parte de "Trait #1: Ancestral state estimates" y la pegas en un excel, le das a pegar con el importador de datos (se hace en el boton que surge al pegar como el de mantener-quitar formato). Así te separará cada columna. Solo falta añadir a "Adj." el "SE" que queda en la siguiente columna (es SE adjusted) y guardar como .csv.

#load results of OU with SE intraespecífica (alpha = 0.076 for bio4 and 0.036 for bio17)
anc_ou_bio4 = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/compare_results/anc_bio4_ou.csv", sep=",", header=TRUE)
str(anc_ou_bio4)
anc_ou_bio17 = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/compare_results/anc_bio17_ou.csv", sep=",", header=TRUE)
str(anc_ou_bio17)

#change names of root by 112+1 (as notation of ape)
#bio4
first_node = length(tree_prunned$tip.label) + 1
levels(anc_ou_bio4$Node) = c(levels(anc_ou_bio4$Node), eval(first_node))
anc_ou_bio4$Node[which(anc_ou_bio4$Node == "Root")] <-  length(tree_prunned$tip.label) + 1
anc_ou_bio4$Node = droplevels(anc_ou_bio4$Node)
#bio17
first_node = length(tree_prunned$tip.label) + 1
levels(anc_ou_bio17$Node) = c(levels(anc_ou_bio17$Node), eval(first_node))
anc_ou_bio17$Node[which(anc_ou_bio17$Node == "Root")] <-  length(tree_prunned$tip.label) + 1
anc_ou_bio17$Node = droplevels(anc_ou_bio17$Node)

#change the name of Node to nodo1 for mergin with data.frame of node-species numbers
colnames(anc_ou_bio4)[which(colnames(anc_ou_bio4) == "Node")] <- "nodo1" 
colnames(anc_ou_bio17)[which(colnames(anc_ou_bio17) == "Node")] <- "nodo1" 

##### BM y OU con ace (ape) #######
## bio4
# BM
ml_recon_bio4 = ace(x=bio4_vector, phy=tree_prunned, type="continuous", method = "ML", CI = TRUE, model="BM", marginal=FALSE) 
    #If type = "continuous", the default model is Brownian motion where characters evolve randomly following a random walk. This model can be fitted by residual maximum likelihood (the default), maximum likelihood (Felsenstein 1973, Schluter et al. 1997), least squares (method = "pic", Felsenstein 1985), or generalized least squares (method = "GLS", Martins and Hansen 1997, Cunningham et al. 1998). In the last case, the specification of phy and model are actually ignored: it is instead given through a correlation structure with the option corStruct. 
    #In the setting method = "ML" and model = "BM" (this used to be the default until ape 3.0-7) the maximum likelihood estimation is done simultaneously on the ancestral values and the variance of the Brownian motion process; these estimates are then used to compute the confidence intervals in the standard way. The REML method first estimates the ancestral value at the root (aka, the phylogenetic mean), then the variance of the Brownian motion process is estimated by optimizing the residual log-likelihood. The ancestral values are finally inferred from the likelihood function giving these two parameters. If method = "pic" or "GLS", the confidence intervals are computed using the expected variances under the model, so they depend only on the tree. 
    #It could be shown that, with a continous character, REML results in unbiased estimates of the variance of the Brownian motion process while ML gives a downward bias. Therefore the former is recommanded.
    #By default, the likelihood of the different ancestral states of discrete characters are computed with a joint estimation procedure using a procedure similar to the one described in Pupko et al. (2000). If marginal = TRUE, a marginal estimation procedure is used (this was the only choice until ape 3.1- 1). With this method, the likelihood values at a given node are computed using only the information from the tips (and branches) descending from this node. With the joint estimation, all information is used for each node. The difference between these two methods is further explained in Felsenstein (2004, pp. 259-260) and in Yang (2006, pp. 121-126). The present implementation of the joint estimation uses a “two-pass” algorithm which is much faster than stochastic mapping while the estimates of both methods are very close.
reml_recon_bio4 = ace(x=bio4_vector, phy=tree_prunned, type="continuous", method = "REML", CI = TRUE, model="BM", marginal=FALSE)
pgls_recon_bio4 = ace(x=bio4_vector, phy=tree_prunned, type="continuous", method = "GLS", CI = TRUE, model="BM", marginal=FALSE, corStruct = corBrownian(1, phy = tree_prunned))
    #ML no estima el 95CI y además el manual de ape dice que ML puede dar estimar con bias en rasgos continuos. Ademá

# BM with reconstruct (APE)
ml_recontruct_bio4 = reconstruct(bio4_vector, tree_prunned, method = "ML", alpha = NULL, CI = TRUE)
reml_recontruct_bio4 = reconstruct(bio4_vector, tree_prunned, method = "REML", alpha = NULL, CI = TRUE)
pgl_recontruct_bio4 = reconstruct(bio4_vector, tree_prunned, method = "GLS", alpha = NULL, CI = TRUE)

#compare ancestral states of BM all models
ml_recon_bio4$ace - reml_recon_bio4$ace
ml_recon_bio4$ace - pgls_recon_bio4$ace
reml_recon_bio4$ace - pgls_recon_bio4$ace

ml_recontruct_bio4$ace - reml_recontruct_bio4$ace
ml_recontruct_bio4$ace - pgl_recontruct_bio4$ace
pgl_recontruct_bio4$ace - reml_recontruct_bio4$ace

ml_recontruct_bio4$ace - ml_recon_bio4$ace
reml_recontruct_bio4$ace - reml_recon_bio4$ace
pgl_recontruct_bio4$ace - pgls_recon_bio4$ace #REML y PGLS y ML (con reconstruct) iguales en todos los casos.

#OU 
pgls_ou_bio4 = reconstruct(bio4_vector, tree_prunned, method = "GLS_OU", alpha = 0.076, CI = TRUE)
pgls_ous_bio4 = reconstruct(bio4_vector, tree_prunned, method = "GLS_OUS", alpha = 0.076, CI = TRUE) #"GLS_OU" and "GLS_OUS" differs in the fact that "GLS_OUS" assume that the process starts from the optimum, while the root state has to be estimated for "GLS_OU", which may rise some issues (see reconstruct 219 Royer-Carenzi and Didier, 2016). Users may provide the attractive strength parameter alpha, for these two models. Users may provide the attractive strength parameter alpha, for these two models. "GLS_ABM", "GLS_OU" and "GLS_OUS" are all fitted by generalized least squares (Royer-Carenzi and Didier, 2016).
    #Note: GLS_OU may lead to aberrant reconstructions.

## bio17
# BM
ml_recon_bio17 = ace(x=bio17_vector, phy=tree_prunned, type="continuous", method = "ML", CI = TRUE, model="BM", marginal=FALSE) 
reml_recon_bio17 = ace(x=bio17_vector, phy=tree_prunned, type="continuous", method = "REML", CI = TRUE, model="BM", marginal=FALSE)
pgls_recon_bio17 = ace(x=bio17_vector, phy=tree_prunned, type="continuous", method = "GLS", CI = TRUE, model="BM", marginal=FALSE, corStruct = corBrownian(1, phy = tree_prunned))
    #ML no estima el 95CI y además el manual de ape dice que ML puede dar estimar con bias en rasgos continuos: "It could be shown that, with a continous character, REML results in unbiased estimates of the variance of the Brownian motion process while ML gives a downward bias. Therefore the former is recommanded". 

# BM with reconstruct (APE)
ml_recontruct_bio17 = reconstruct(bio17_vector, tree_prunned, method = "ML", alpha = NULL, CI = TRUE)
reml_recontruct_bio17 = reconstruct(bio17_vector, tree_prunned, method = "REML", alpha = NULL, CI = TRUE)
pgl_recontruct_bio17 = reconstruct(bio17_vector, tree_prunned, method = "GLS", alpha = NULL, CI = TRUE)

#compare ancestral states of BM all models
ml_recon_bio17$ace - reml_recon_bio17$ace
ml_recon_bio17$ace - pgls_recon_bio17$ace
reml_recon_bio17$ace - pgls_recon_bio17$ace

ml_recontruct_bio17$ace - reml_recontruct_bio17$ace
ml_recontruct_bio17$ace - pgl_recontruct_bio17$ace
pgl_recontruct_bio17$ace - reml_recontruct_bio17$ace

ml_recontruct_bio17$ace - ml_recon_bio17$ace
reml_recontruct_bio17$ace - reml_recon_bio17$ace
pgl_recontruct_bio17$ace - pgls_recon_bio17$ace #REML y PGLS y ML (con reconstruct) iguales en todos los casos. ML en ace no calcular el 95CI y además según el paquete mejor REML que ML.

#OU 
pgls_ou_bio17 = reconstruct(bio17_vector, tree_prunned, method = "GLS_OU", alpha = 0.036, CI = TRUE)
pgls_ous_bio17 = reconstruct(bio17_vector, tree_prunned, method = "GLS_OUS", alpha = 0.036, CI = TRUE) #Mas alto ape

# OU en ape vs OU en compare 4.6 (SE and Adj.SE) in bio17
abs(pgls_ous_bio17$CI95[,1] + (-pgls_ous_bio17$CI95[,2])) - anc_ou_bio17$SE
abs(pgls_ous_bio17$CI95[,1] + (-pgls_ous_bio17$CI95[,2])) - anc_ou_bio17$Adj.SE #range of pgls_ous are broader than SE and Adj.SE

#differences between SE and Adj.SE in bio17
median(anc_ou_bio17$SE)
median(anc_ou_bio17$Adj.SE)
wilcox.test(anc_ou_bio17$SE, anc_ou_bio17$Adj.SE)
plot(density(anc_ou_bio17$SE - anc_ou_bio17$Adj.SE)) 
quantile(anc_ou_bio17$SE - anc_ou_bio17$Adj.SE, probs=c(0.025, 0.975))
#CI95 overlap with zero, although significant the wilcoxon test

# OU en ape vs OU en compare 4.6 (SE and Adj.SE) in bio4
abs(pgls_ous_bio4$CI95[,1] + (-pgls_ous_bio4$CI95[,2])) - anc_ou_bio4$SE
abs(pgls_ous_bio4$CI95[,1] + (-pgls_ous_bio4$CI95[,2])) - anc_ou_bio4$Adj.SE #range of pgls_ous are broader than SE and Adj.SE

#differences between SE and Adj.SE in bio4
median(anc_ou_bio4$SE)
median(anc_ou_bio4$Adj.SE)
wilcox.test(anc_ou_bio4$SE, anc_ou_bio4$Adj.SE)
plot(density(anc_ou_bio4$SE - anc_ou_bio4$Adj.SE))
quantile(anc_ou_bio4$SE - anc_ou_bio4$Adj.SE, probs=c(0.025, 0.975)) #CI95 does not overlap with zero, and significant the wilcoxon test


#decision
#Adj.SE has lower values than SE, but as we are only to use BM, and we only want SE of OU to compare with SE of BM we are going to use the regular SE of compare, which probably is more similar to the SE of reconstruct According to comapre, SE is estimate with the Schluter et al. method, which may differ only slightly from Hansen method, the method used for reconstruct. Adj.SE correspond to the method of Rohlf (2001, Evolution 55:2143-2150), who proposed another way of estimating standard errors for ancestral states. We will use SE because of that.

#### final reconstructed variables ####
anc_ou_bio4
anc_ou_bio17
pgls_recon_bio4
pgls_recon_bio17 #For PGLS under BM we have to add corBrownian(1, phy = tree_prunned)


#calculate quantiles for the ancestral state of bio4 and the CI
quantile_ace_bio4 = quantile(pgls_recon_bio4$ace, probs=c(0.025, 0.5, 0.975))
quantile_ace_ci_bio4 = quantile(abs(pgls_recon_bio4$CI95[,1]-pgls_recon_bio4$CI95[,2]), probs=c(0.025, 0.5, 0.975))
quantile_ace_bio4
quantile_ace_ci_bio4

#calculate percentage of BIO4 that correspond to variability. For that we comapre the medians of ace BIO and medians of the difference CI BIO4
(quantile_ace_ci_bio4[2]*100)/quantile_ace_bio4[2]


#calculate quantiles for the ancestrlaa state of bio17 and the CI
quantile_ace_bio17 = quantile(pgls_recon_bio17$ace, probs=c(0.025, 0.5, 0.975))
quantile_ace_ci_bio17 = quantile(abs(pgls_recon_bio17$CI95[,1]-pgls_recon_bio17$CI95[,2]), probs=c(0.025, 0.5, 0.975))
quantile_ace_bio17
quantile_ace_ci_bio17

#calculate percentage of bio17 that correspond to variability. For that we comapre the medians of ace BIO and medians of the difference CI bio17
(quantile_ace_ci_bio17[2]*100)/abs(quantile_ace_bio17[2])

#### plots the final models

##check that the order of rows is correct
summary(names(pgls_recon_bio4$ace) == anc_ou_bio4$nodo1)
summary(names(pgls_recon_bio17$ace) == anc_ou_bio17$nodo1)

## modify function to add error bars of ancestral state (errorbar.contMap), because of a problem in the calculation of ii and jj. There are problems when the lower limit of CI95 and the lower limit of the range of current values of the trait are both negatives. The function calculates the difference, but the resulting number is negative, so I have added an abs()
#in "ancestral.states" you have to include an object with the state as ace and the CI95 as CI95. Typical ace object. 
errorbar_contMap_modified = function (obj, user=FALSE, anc.states=NULL, ...){
    if (hasArg(x)){
        x <- list(...)$x
    } else{
        x <- setNames(sapply(1:Ntip(obj$tree), function(x, obj) {
        ii <- which(obj$tree$edge[, 2] == x)
        ss <- names(obj$tree$maps[[ii]][length(obj$tree$maps[[ii]])])
        obj$lims[1] + as.numeric(ss)/(length(obj$cols) - 1) * 
            diff(obj$lims)
        }, obj = obj), obj$tree$tip.label)
    }
    if (hasArg(scale.by.ci)) {
        scale.by.ci <- list(...)$scale.by.ci
    } else {
        scale.by.ci <- TRUE
    }
    if (hasArg(lwd)){
        lwd <- list(...)$lwd
    } else {
        lwd <- 14
    }
    tree <- obj$tree
    if(user==FALSE){
        aa <- fastAnc(tree, x, CI = TRUE)
    } else {
        aa = anc.states
    }
    xlim <- range(aa$CI95)
    if (xlim[2] > obj$lims[2] || xlim[1] < obj$lims[1]) {
        cat(paste("  -----\n  The range of the contMap object, presently (", 
            round(obj$lims[1], 4), ",", round(obj$lims[2], 4), 
            "), should be equal to\n  or greater than the range of the CIs on ancestral states: (", 
            round(xlim[1], 4), ",", round(xlim[2], 4), ").\n", 
            sep = ""))
        cat(paste("  To ensure that your error bars are correctly plotted, please recompute your\n", 
            "  contMap object and increase lims.\n  -----\n", 
            sep = ""))
    }
    d <- diff(obj$lims)
    if (scale.by.ci) {
        v <- aa$CI95[, 2] - aa$CI95[, 1]
        v <- v/max(v)
    } else {
        v <- rep(0.5, tree$Nnode)
    }    
    n <- length(obj$cols) - 1
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    h <- max(nodeHeights(tree))
    for (i in 1:tree$Nnode) {
        ii <- round((abs(aa$CI95[i, 1] - obj$lims[1]))/d * n)
        jj <- round((abs(aa$CI95[i, 2] - obj$lims[1]))/d * (n + 1))
        cols <- obj$cols[ii:jj]
        add.color.bar(leg = 0.1 * h * v[i], cols = cols, prompt = FALSE, 
            x = lastPP$xx[i + Ntip(tree)] - 0.05 * h * v[i], 
            y = lastPP$yy[i + Ntip(tree)], title = "", subtitle = "", 
            lims = NULL, lwd = lwd)
    }
}

##create a tree with species names as "P. XXXXX" to save space in the plots
#before check that the order of the species names from climate_medians and tip.label of the tree is the same. 
summary(tree_prunned$tip.label == paste("Pinus_", climate_medians$species, sep=""))
#copy the tree
tree_prunned_to_plot = tree_prunned
#add reduced names
tree_prunned_to_plot$tip.label <- paste("P. ", climate_medians$species, sep="")

##BM bio4
#create a vector with the variable and species names as "P. XXXXX", similar to the tip.labels of tree_prunned_to_plot
summary(names(bio4_vector) == paste("Pinus_", climate_medians$species, sep="")) #before check that the order of the species names from climate_medians and bio4_vector is the same. 
bio4_vector_to_plot = bio4_vector #copy the vector of the trait
names(bio4_vector_to_plot) <- paste("P. ", climate_medians$species, sep="") #add reduced names

#open the pdf of the plot
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/plot_anc/plot_ancestral_states_bm_bio4.pdf", width = 24, height = 24)

#create the tree interpolated
obj_bm_bio4 = contMap(tree_prunned_to_plot, bio4_vector_to_plot, method="user", anc.states=pgls_recon_bio4$ace, plot=FALSE, lwd = 2, res=500) #aquí anc.states is a file with only ace, without 95CI.

#plot
plot(obj_bm_bio4, type="fan", fsize=c(2.5, 2.3)) #fsize tiene el primer valor del tamaño de las tip labels, el segundo es para la legenda.

#add the error bars
errorbar_contMap_modified(obj=obj_bm_bio4, user=TRUE, anc.states=pgls_recon_bio4, scale.by.ci=TRUE) #aquí anc.states is the complete file with ace and 95CI.
    #scale.by.ci=TRUE: that determines whether or not the length of the error bars will be scaled by the CI width
dev.off()


##OU bio4
#create a vector with the variable and species names as "P. XXXXX", similar to the tip.labels of tree_prunned_to_plot
summary(names(bio4_vector) == paste("Pinus_", climate_medians$species, sep="")) #before check that the order of the species names from climate_medians and bio4_vector is the same. 
bio4_vector_to_plot = bio4_vector #copy the vector of the trait
names(bio4_vector_to_plot) <- paste("P. ", climate_medians$species, sep="") #add reduced names

#create a list with ancestral state and CI95 of OU to plot
ace_bio4_ou = anc_ou_bio4$State #extract acenstral state
names(ace_bio4_ou) <- anc_ou_bio4$nodo1 #set names as node1
CI95_bio4_ou = matrix(NA, ncol=2, nrow=nrow(anc_ou_bio4)) #create a matrix with 95CI
CI95_bio4_ou[,1] <- anc_ou_bio4$State - anc_ou_bio4$SE #first column lower limit
CI95_bio4_ou[,2] <- anc_ou_bio4$State + anc_ou_bio4$SE #second column upper limit
row.names(CI95_bio4_ou) <- anc_ou_bio4$nodo1 #node1 as row.names
anc_ou_bio4_to_plot = list() #bind all
anc_ou_bio4_to_plot[["ace"]] <- ace_bio4_ou
anc_ou_bio4_to_plot[["CI95"]] <- CI95_bio4_ou

#create the tree interpolated
obj_ou_bio4 = contMap(tree_prunned_to_plot, bio4_vector_to_plot, method="user", anc.states=anc_ou_bio4_to_plot$ace, plot=FALSE, lwd = 2, res=500) #aquí anc.states is a file with only ace, without 95CI.

#open the pdf of the plot
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/plot_anc/plot_ancestral_states_ou_bio4.pdf", width = 24, height = 24)

#plot
plot(obj_ou_bio4, type="fan", fsize=c(2.5, 2.3)) #fsize tiene el primer valor del tamaño de las tip labels, el segundo es para la legenda.

#add the error bars
errorbar_contMap_modified(obj=obj_ou_bio4, user=TRUE, anc.states=anc_ou_bio4_to_plot, scale.by.ci=TRUE) #aquí anc.states is the complete file with ace and 95CI.
    #scale.by.ci=TRUE: that determines whether or not the length of the error bars will be scaled by the CI width
dev.off()

#plot for BM vs. OU comparison
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/plot_anc/anc_bio4_bm_ou.pdf", width = 24, height = 11)
par(mfrow=c(1,2))
plot(obj_bm_bio4, type="phylogram", fsize=c(1, 2.3), xlim=c(-10,155), ftype="off",offset = 2.2)
plot(obj_ou_bio4, type="phylogram", fsize=c(1, 2.3), xlim=c(-20,155), ftype="off", direction="leftwards") #fsize tiene el primer valor del tamaño de las tip labels, el segundo es para la legenda.
dev.off()


##BM bio17
#create a vector with the variable and species names as "P. XXXXX", similar to the tip.labels of tree_prunned_to_plot
summary(names(bio17_vector) == paste("Pinus_", climate_medians$species, sep="")) #before check that the order of the species names from climate_medians and bio17_vector is the same. 
bio17_vector_to_plot = bio17_vector #copy the vector of the trait
names(bio17_vector_to_plot) <- paste("P. ", climate_medians$species, sep="") #add reduced names

#open the pdf of the plot
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/plot_anc/plot_ancestral_states_bm_bio17.pdf", width = 24, height = 24)

#create the tree interpolated
obj_bm_bio17 = contMap(tree_prunned_to_plot, bio17_vector_to_plot, method="user", anc.states=pgls_recon_bio17$ace, plot=FALSE, lwd = 2, res=500) #aquí anc.states is a file with only ace, without 95CI.

#plot
plot(obj_bm_bio17, type="fan", fsize=c(2.5, 2.3)) #fsize tiene el primer valor del tamaño de las tip labels, el segundo es para la legenda.

#add the error bars
errorbar_contMap_modified(obj=obj_bm_bio17, user=TRUE, anc.states=pgls_recon_bio17, scale.by.ci=TRUE) #aquí anc.states is the complete file with ace and 95CI.
    #scale.by.ci=TRUE: that determines whether or not the length of the error bars will be scaled by the CI width
dev.off()


##OU bio17
#create a vector with the variable and species names as "P. XXXXX", similar to the tip.labels of tree_prunned_to_plot
summary(names(bio17_vector) == paste("Pinus_", climate_medians$species, sep="")) #before check that the order of the species names from climate_medians and bio17_vector is the same. 
bio17_vector_to_plot = bio17_vector #copy the vector of the trait
names(bio17_vector_to_plot) <- paste("P. ", climate_medians$species, sep="") #add reduced names

#create a list with ancestral state and CI95 of OU to plot
ace_bio17_ou = anc_ou_bio17$State #extract acenstral state
names(ace_bio17_ou) <- anc_ou_bio17$nodo1 #set names as node1
CI95_bio17_ou = matrix(NA, ncol=2, nrow=nrow(anc_ou_bio17)) #create a matrix with 95CI
CI95_bio17_ou[,1] <- anc_ou_bio17$State - anc_ou_bio17$SE #first column lower limit
CI95_bio17_ou[,2] <- anc_ou_bio17$State + anc_ou_bio17$SE #second column upper limit
row.names(CI95_bio17_ou) <- anc_ou_bio17$nodo1 #node1 as row.names
anc_ou_bio17_to_plot = list() #bind all
anc_ou_bio17_to_plot[["ace"]] <- ace_bio17_ou
anc_ou_bio17_to_plot[["CI95"]] <- CI95_bio17_ou

#create the tree interpolated
obj_ou_bio17 = contMap(tree_prunned_to_plot, bio17_vector_to_plot, method="user", anc.states=anc_ou_bio17_to_plot$ace, plot=FALSE, lwd = 2, res=500) #aquí anc.states is a file with only ace, without 95CI.

#open the pdf of the plot
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/plot_anc/plot_ancestral_states_ou_bio17.pdf", width = 24, height = 24)

#plot
plot(obj_ou_bio17, type="fan", fsize=c(2.5, 2.3)) #fsize tiene el primer valor del tamaño de las tip labels, el segundo es para la legenda.

#add the error bars
errorbar_contMap_modified(obj=obj_ou_bio17, user=TRUE, anc.states=anc_ou_bio17_to_plot, scale.by.ci=TRUE) #aquí anc.states is the complete file with ace and 95CI.
    #scale.by.ci=TRUE: that determines whether or not the length of the error bars will be scaled by the CI width
dev.off()

#plot for BM vs. OU comparison
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/plot_anc/anc_bio17_bm_ou.pdf", width = 24, height = 11)
par(mfrow=c(1,2))
plot(obj_bm_bio17, type="phylogram", fsize=c(1.45, 2.3), xlim=c(-10,155), ftype="off", offset = 2.2)
plot(obj_ou_bio17, type="phylogram", fsize=c(1.45, 2.3), xlim=c(-20,155), ftype="off", direction="leftwards") #fsize tiene el primer valor del tamaño de las tip labels, el segundo es para la legenda.
dev.off()


#Si recuerdas estoy estimando el estado ancestral bajo BM usando la función ace de APE. Pues bien, los intervalos de confianza de la estimas son extremadamente diferentes entre tipos de ajuste. Según PGLS los intervalos son muuuy estrechos, mientras que para REML son mucho más amplios (anchura del intervalo de 10 unidades en PGLS frente 300 en REML; con ML no se pueden llegar a calcular). Te adjunto dos plots con el estado ancestral de bio17 bajo REML y PGLS, las barritas en cada nodo indican el 95CI. He probado a calcular los estado ancestrales con la otra función de APE para esto, reconstruct, la cual se diferencia de ace en que "computations are not performed by numerical optimisation but through matrix calculus", pero se pueden usar los mismos tipos de ajuste (ML, REML y PGLS). En este caso salen los intervalos amplios para los tres casos. Por tanto, los intervalos muuuy estrechos solo salen en PGLS con ace.
#No he encontrado reconstrucciones con tan poca incertidumbre en ninguno de los ejemplos que tiene Revell por la red, esto unido al hecho de que no se replican las diferencias entre ajustes cuando se repite todo con "reconstruct", me hace tener muy poca seguridad para usar esos intervalos "estrechos". Ante esto hay cuatro opciones: 
    #Usar solo los modelos con intervalos estrechos. No me fío de esos intervalos, el extremo del intervalo podría estar ahí ó mucho más lejos.  
    #Usar los intervalos de confianza obtenidos con PGLS en reconstruct, ó con REML en cualquiera de las dos funciones (ace - reconstruct). Estos intervalos me dan algo más seguridad, pero al haber tanta incertidumbre (intervalos amplios) el valor actual de muchas especies va a caer dentro del intervalo y por tanto no se va a poder hacer corregir la idoneidad por la filogenia (si no tenemos certeza de que el valor actual se diferencie del valor del último nodo, la amplitud del rango filogenético es 0).
    #Usar todos los modelos (intervalos amplios y estrechos). Creo que esto no tiene mucho sentido, porque si para los casos con intervalos amplios no se puede corregir la idoneidad de hábitat por la filogenia, sería como estar usando solo los estrechos.
    #Calcular la media de los dos intervalos (REML vs PGLS). Si el extremo superior de ambos modelos es 500 y 300 respectivamente, el "consenso" sería 400. No me convence, porque no me da seguridad ningún intervalo, así que el consenso de ellos no debe ser muy fiable. Para eso, aunque ampliemos más el rango filo, prefiero la siguiente opción. 
    #NO usar los intervalos de confianza de las reconstrucciones. Es una putada no usar medida de incertidumbre pero creo que es la mejor opción dada la poca seguridad que tenemos sobre las 95CI. He comparados los valores de cada nodo (sin intervalos) y son prácticamente iguales entre todos los modelos y funciones (ace - reconstruct). Solo hay discrepancia entre ML de ace con el resto, pero esa reconstrucción era muy mala (ni se estima el 95CI) y de hecho en el manual de ape se desaconseja usarla. Si tomamos esta opción, el rango filogenético iría desde el valor del último nodo (sin incluirlo) hasta el valor actual (incluyéndolo). 
#Todo esto aplica para bio17, en el caso de bio4 hay menos diferencias entre PGLS y REML, pero estas siguen siendo grandes: Anchura del intervalo de 10 unidades frente a 150. Por tanto, no podemos tener certeza de cuales son los intervalos. 


#############################################
####### Comparsion of ancestral states ######
#############################################
### Comparison of OU in comapre with SE of intraespecific variation and BM in ace without intraespecific variation

#number of nodes
nrow(anc_ou_bio4)
nrow(anc_ou_bio17)
nrow(pgls_recon_bio4$CI95)
nrow(pgls_recon_bio17$CI95)

#check that the order of nodes is the same
summary(anc_ou_bio4$nodo1 == row.names(pgls_recon_bio4$CI95))
summary(anc_ou_bio17$nodo1 == row.names(pgls_recon_bio17$CI95))

###bio4
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/comparison_bm_ou/comparsion_bm_ou_bio4.pdf")
par(mfrow=c(2,2))

##pannel 1
#anc states  
plot(anc_ou_bio4$State ~ pgls_recon_bio4$ace, xlab="Ancestral states under BM", ylab="Ancestral states under OU", main="")
#cor between anc states
spearman<-cor.test(anc_ou_bio4$State, pgls_recon_bio4$ace, method="spearman")
mylabel = bquote(italic(rho) == .(format(spearman$estimate, digits = 3)))
text(x=9500, y=5800, labels = mylabel, cex=1)
if(!spearman$p.value == 0){ #if p.value is <2.2e-16, spearman save it as 0
    mylabel.p = bquote(italic(p.value) == .(format(spearman$p.value, digits = 3)))
} else {
    mylabel.p = bquote(italic(p.value) < .(format(2.2e-16)))
}
text(x=9500, y=5500, labels = mylabel.p, cex=1)
mylabel.t = bquote(italic(S) == .(format(spearman$statistic, digits = 3)))
text(x=9500, y=5200, labels = mylabel.t, cex=1)

##pannel 2
if(FALSE){
#densitplot of differences between ancestral states
    plot(density(abs(anc_ou_bio4$State - pgls_recon_bio4$ace)), main="Absolute difference ancestral states")
    p<-wilcox.test(anc_ou_bio4$State, pgls_recon_bio4$ace)$p.value
    mylabel.p = bquote(italic(p.value) == .(format(p, digits = 3)))
    text(x=3000, y=0.0005, labels = mylabel.p, cex=1)
}

##pannel 3
#SE of OU versus CI95 of BM
plot(anc_ou_bio4$SE~abs(pgls_recon_bio4$CI95[,1] - pgls_recon_bio4$CI95[,2]), xlab="Standard errors under BM", ylab="Standard errors under OU", main="")
#cors
spearman<-cor.test(anc_ou_bio4$SE, abs(pgls_recon_bio4$CI95[,1] - pgls_recon_bio4$CI95[,2]), method="spearman")
mylabel = bquote(italic(rho) == .(format(spearman$estimate, digits = 3)))
text(x=10, y=2700, labels = mylabel, cex=1)
if(!spearman$p.value == 0){ #if p.value is <2.2e-16, spearman save it as 0
    mylabel.p = bquote(italic(p.value) == .(format(spearman$p.value, digits = 3)))
} else {
    mylabel.p = bquote(italic(p.value) < .(format(2.2e-16)))
}
text(x=10, y=2400, labels = mylabel.p, cex=1)
mylabel.t = bquote(italic(S) == .(format(spearman$statistic, digits = 3)))
text(x=10, y=2100, labels = mylabel.t, cex=1)


##pannel 4
if(FALSE){
    #densitplot of differences between error of ancestral state estimates
    plot(density(abs(anc_ou_bio4$SE - abs(pgls_recon_bio4$CI95[,1] - pgls_recon_bio4$CI95[,2]))), main="Absolute difference standard errors")
    p<-wilcox.test(anc_ou_bio4$SE, abs(pgls_recon_bio4$CI95[,1] - pgls_recon_bio4$CI95[,2]))$p.value
    mylabel.p = bquote(italic(p.value) == .(format(p, digits = 3)))
    text(x=4800, y=0.0004, labels = mylabel.p, cex=1)
}
dev.off()

###bio17
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/comparison_bm_ou/comparsion_bm_ou_bio17.pdf")
par(mfrow=c(2,2))

##pannel 1
#anc states  
plot(anc_ou_bio17$State ~ pgls_recon_bio17$ace, xlab="Ancestral states under BM", ylab="Ancestral states under OU", main="")
#cor between anc states
spearman<-cor.test(anc_ou_bio17$State, pgls_recon_bio17$ace, method="spearman")
mylabel = bquote(italic(rho) == .(format(spearman$estimate, digits = 3)))
text(x=-600, y=-300, labels = mylabel, cex=1)
if(!spearman$p.value == 0){ #if p.value is <2.2e-16, spearman save it as 0
    mylabel.p = bquote(italic(p.value) == .(format(spearman$p.value, digits = 3)))
} else {
    mylabel.p = bquote(italic(p.value) < .(format(2.2e-16)))
}
text(x=-600, y=-340, labels = mylabel.p, cex=1)
mylabel.t = bquote(italic(S) == .(format(spearman$statistic, digits = 3)))
text(x=-600, y=-380, labels = mylabel.t, cex=1)

##pannel 2
if(FALSE){
    #densitplot of differences between ancestral states
    plot(density(abs(anc_ou_bio17$State - pgls_recon_bio17$ace)), main="Absolute difference ancestral states")
    p<-wilcox.test(anc_ou_bio17$State, pgls_recon_bio17$ace)$p.value
    mylabel.p = bquote(italic(p.value) == .(format(p, digits = 3)))
    text(x=100, y=0.02, labels = mylabel.p, cex=1)
} 

##pannel 3
#SE of OU versus CI95 of BM
plot(anc_ou_bio17$SE~abs(pgls_recon_bio17$CI95[,1] - pgls_recon_bio17$CI95[,2]), xlab="Standard errors under BM", ylab="Standard errors under OU", main="")
#cors
spearman<-cor.test(anc_ou_bio17$SE, abs(pgls_recon_bio17$CI95[,1] - pgls_recon_bio17$CI95[,2]), method="spearman")
mylabel = bquote(italic(rho) == .(format(spearman$estimate, digits = 3)))
text(x=11, y=100, labels = mylabel, cex=1)
if(!spearman$p.value == 0){ #if p.value is <2.2e-16, spearman save it as 0
    mylabel.p = bquote(italic(p.value) == .(format(spearman$p.value, digits = 3)))
} else {
    mylabel.p = bquote(italic(p.value) < .(format(2.2e-16)))
}
text(x=11, y=85, labels = mylabel.p, cex=1)
mylabel.t = bquote(italic(S) == .(format(spearman$statistic, digits = 3)))
text(x=11, y=70, labels = mylabel.t, cex=1)


##pannel 4
if(FALSE){
    #densitplot of differences between error of ancestral state estimates
    plot(density(abs(anc_ou_bio17$SE - abs(pgls_recon_bio17$CI95[,1] - pgls_recon_bio17$CI95[,2]))), main="Absolute difference standard errors")
    p<-wilcox.test(anc_ou_bio17$SE, abs(pgls_recon_bio17$CI95[,1] - pgls_recon_bio17$CI95[,2]))$p.value
    mylabel.p = bquote(italic(p.value) == .(format(p, digits = 3)))
    text(x=230, y=0.008, labels = mylabel.p, cex=1)
}
dev.off()

######Conclusions
#Bio17: Como era de esperar por su valor de lambda tan alto (0.87), a penas muestra diferencias entre modelos. Los estados ancestrales están muy correlacionados (rho=0.975)  y no son significativamente diferentes (p = 0.6; panel 1 y 2). Los errores también están muy correlacionados, pero si hay diferencias entre modelos, siendo los errores de OU más grandes (esto creo que se debe más al programa [compare vs ace] que a los modelos). 
#Bio4: En bio4 está la cosa un poco menos clara, lo cual es lógico viendo el valor de lambda no tan alto para esta variable (0.57). Sigue habiendo una correlación alta entre estados ancestrales (rho = 0.87), pero sí hay diferencias significativas entre modelos aunque con una p no muy baja (p = 0.03). Los errores de las estimaciones son significativamente muuucho más grandes en OU.
#Yo creo que hay una correlación clara entre ambos modelos, se parecen mucho, y es lógico teniendo en cuenta los valores de lambda altos (sobre todo en bio17) y los valores de alfa tan bajos (0.076 para bio4 y 0.036 para bio17). Sin embargo, me preocupan un poco las diferencias en valor absoluto de los errores estándar para bio4, pero como los errores no se van a usar, solo los estados ancestrales, tampoco veo mucho problema. Voy a seguir para adelante con BM. en caso de que hubiese que dar marcha atrás no va a ser mucho problema de trabajo, solo un poquito de tiempo para correr de nuevo los análisis.
#Lo DEJO TODO PREPARADO PARA USAR OU TAMBIÉN, PERO NO LO VOY A USAR EN LAS FIGURAS, DEJO ABIERTA LA OPCIÓN PARA MÁS ADELANTE. 


#####################################
####### calculate phylo ranges ######
#####################################

## Usamos PGLS aunque tengan interalos de confianza raros para BM, como no veíamos un simulitud entre los interalos de los diferentes ajustes para BM pasamos de usar interalos y usamos solo los estados ancestrales, los cuales son muuy parecidos entre ajusted bajo BM y además son muuy parecidos entre OU y BM. 

## Calculo los rangos para OU también pero no lo vamos a usar finalmente. 


## check that all files (from compare and ape) has the same order and names of nodes
summary(anc_ou_bio4$nodo1 == names(pgls_recon_bio4$ace))
summary(anc_ou_bio17$nodo1 == names(pgls_recon_bio17$ace))

#extract intial and final node of all branches
ramas = as.data.frame(tree_prunned$edge)
colnames(ramas) <- c("nodo1", "nodo2")

## create a variable final nodes of each species 
#select branch in which are implicate terminal nodes (speices)
final_ramas = ramas[which(ramas$nodo2 %in% 1:length(epithet_species_list)),]
#reoder tree labels in basis on nodo2 positions
species_ramas = tree_prunned$tip.label[final_ramas$nodo2]
#bind
final_ramas = cbind.data.frame(final_ramas, species_ramas)

## check that each node2 has the correct species name
summary(names(bio17_vector)[final_ramas$nodo2] == final_ramas$species_ramas) #all true

#plot the tree to check that each pair of species in the plot correspond with species with the same nodo1
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_recons/FBDl_big_names", width=12, height=16)
plot(tree_prunned)
dev.off()

## create a new data.frame with current values of climate variables to merge with ancestral data
climate_medians_merge = climate_medians
climate_medians_merge$species = paste("Pinus_", climate_medians_merge$species, sep="")
colnames(climate_medians_merge)[which(colnames(climate_medians_merge) == "species")] <- "species_ramas"
climate_medians_merge_bio4 = climate_medians_merge[,which(colnames(climate_medians_merge) %in% c("species_ramas", "median_bio4"))]
climate_medians_merge_bio17 = climate_medians_merge[,which(colnames(climate_medians_merge) %in% c("species_ramas", "median_bio17"))]


#### bind final ramas with ancestral reconstructions from OU (compare) ####
## merge species names and anc data
final_anc_ou_bio4 = merge(final_ramas, anc_ou_bio4, by="nodo1")
final_anc_ou_bio17 = merge(final_ramas, anc_ou_bio17, by="nodo1")

## check that each node has the correct ancestral value comparing with the raw data.frame of ancestrla reconstructions
summary(anc_ou_bio4[which(anc_ou_bio4$nodo1 %in% final_anc_ou_bio4$nodo1),]$State == final_anc_ou_bio4[which(!duplicated(final_anc_ou_bio4$nodo1)),]$State)
summary(anc_ou_bio17[which(anc_ou_bio17$nodo1 %in% final_anc_ou_bio17$nodo1),]$State == final_anc_ou_bio17[which(!duplicated(final_anc_ou_bio17$nodo1)),]$State)

## additional check to test if species that share the last node have the same ancestral state for bio4

#which species had duplicated nodes (i.e. share their last node with another species)
duplicated_nodo1 = which(duplicated(final_anc_ou_bio4$nodo1))

#bind the row number of species duplicated and their pairs (they are the previous row). Then sort the rows
species_pairs = sort(c(duplicated_nodo1-1, duplicated_nodo1))

#for each pair, sequence from the 1 to the total number of species with paris (selecting X and X+1)
test_bio4 = NULL
for(i in seq(1, length(species_pairs), 2)){

    #extract the state of the first species of the pair
    state_first_species = final_anc_ou_bio4[species_pairs[i],]$State

    #extract the state of the second species of the pair
    state_second_species = final_anc_ou_bio4[species_pairs[i+1],]$State

    #check that both share the nodo1
    print(final_anc_ou_bio4[species_pairs[i],]$nodo1 == final_anc_ou_bio4[species_pairs[i+1],]$nodo1)

    #test the existence of differences of ancestral state between them
    test_bio4 = append(test_bio4, state_first_species == state_second_species)
}
length(test_bio4) == length(duplicated_nodo1)
summary(test_bio4)

#for each pair, sequence from the 1 to the total number of species with paris (selecting X and X+1)
test_bio17 = NULL
for(i in seq(1, length(species_pairs), 2)){

    #extract the state of the first species of the pair
    state_first_species = final_anc_ou_bio4[species_pairs[i],]$State

    #extract the state of the second species of the pair
    state_second_species = final_anc_ou_bio4[species_pairs[i+1],]$State

    #check that both share the nodo1
    print(final_anc_ou_bio4[species_pairs[i],]$nodo1 == final_anc_ou_bio4[species_pairs[i+1],]$nodo1)

    #test the existence of differences of ancestral state between them
    test_bio17 = append(test_bio17, state_first_species == state_second_species)
}
length(test_bio17) == length(duplicated_nodo1)
summary(test_bio17)

## bind current climatic data
final_anc_ou_bio4 = merge(final_anc_ou_bio4, climate_medians_merge_bio4, by="species_ramas")
final_anc_ou_bio17 = merge(final_anc_ou_bio17, climate_medians_merge_bio17, by="species_ramas")
row.names(final_anc_ou_bio4) == final_anc_ou_bio4$nodo2
row.names(final_anc_ou_bio17) == final_anc_ou_bio17$nodo2


## check that current value is added correctly
summary(climate_medians_merge[match(climate_medians_merge$species_ramas, final_anc_ou_bio17$species_ramas),]$species_ramas == final_anc_ou_bio17$species_ramas)
summary(climate_medians_merge[match(climate_medians_merge$species_ramas, final_anc_ou_bio4$species_ramas),]$species_ramas == final_anc_ou_bio4$species_ramas)

## reorder columns and set new col names
#bio4
final_anc_ou_bio4 = final_anc_ou_bio4[,c(1,3,7,2,4,5,6)]
colnames(final_anc_ou_bio4) <- c("species", "node_species", "current_value", "node_antecesor", "ace", "SE", "Adj.SE")
str(final_anc_ou_bio4)
head(final_anc_ou_bio4, 10)
#bio17
final_anc_ou_bio17 = final_anc_ou_bio17[,c(1,3,7,2,4,5,6)]
colnames(final_anc_ou_bio17) <- c("species", "node_species", "current_value", "node_antecesor", "ace", "SE", "Adj.SE")
str(final_anc_ou_bio17)
head(final_anc_ou_bio17, 10)

#### bind final ramas with ancestral reconstructions from BM (ape) ####
##bio4
#Bind node number,95CI and ace,
final_anc_bm_bio4 = cbind.data.frame(row.names(pgls_recon_bio4$CI95), pgls_recon_bio4$ace, pgls_recon_bio4$CI95)
colnames(final_anc_bm_bio4) <- c("nodo1", "State", "low_bound", "upper_bound")
#merge anc data with species names
final_anc_bm_bio4 = merge(final_ramas, final_anc_bm_bio4, by="nodo1")
#reorder rows in basis on node numbers
final_anc_bm_bio4 = final_anc_bm_bio4[order(final_anc_bm_bio4$nodo1, decreasing=FALSE),]
#check that the anc data has been correctly extraced
summary(final_anc_bm_bio4[which(!duplicated(final_anc_bm_bio4$nodo1)),c(5,6)] == pgls_recon_bio4$CI95[which(row.names(pgls_recon_bio4$CI95) %in% final_ramas$nodo1),])
#add current values
final_anc_bm_bio4 = merge(final_anc_bm_bio4, climate_medians_merge_bio4, by="species_ramas")
#check that current value is added correctly
summary(climate_medians_merge[match(final_anc_bm_bio4$species_ramas, climate_medians_merge$species_ramas),]$species_ramas == final_anc_bm_bio4$species_ramas)
#reorder final data
final_anc_bm_bio4 = final_anc_bm_bio4[,c(1,3,7,2,4,5,6)]
colnames(final_anc_bm_bio4) = c("species", "node_species", "current_value", "node_antecesor", "ace", "ace_low_bound", "ace_upper_bound")
str(final_anc_bm_bio4)
head(final_anc_bm_bio4, 20)

##bio17
#Bind node number,95CI and ace,
final_anc_bm_bio17 = cbind.data.frame(row.names(pgls_recon_bio17$CI95), pgls_recon_bio17$ace, pgls_recon_bio17$CI95)
colnames(final_anc_bm_bio17) <- c("nodo1", "State", "low_bound", "upper_bound")
#merge anc data with species names
final_anc_bm_bio17 = merge(final_ramas, final_anc_bm_bio17, by="nodo1")
#reorder rows in basis on node numbers
final_anc_bm_bio17 = final_anc_bm_bio17[order(final_anc_bm_bio17$nodo1, decreasing=FALSE),]
#check that the anc data has been correctly extraced
summary(final_anc_bm_bio17[which(!duplicated(final_anc_bm_bio17$nodo1)),c(5,6)] == pgls_recon_bio17$CI95[which(row.names(pgls_recon_bio17$CI95) %in% final_ramas$nodo1),])
#add current values
final_anc_bm_bio17 = merge(final_anc_bm_bio17, climate_medians_merge_bio17, by="species_ramas")
#check that current value is added correctly
summary(climate_medians_merge[match(final_anc_bm_bio17$species_ramas, climate_medians_merge$species_ramas),]$species_ramas == final_anc_bm_bio17$species_ramas)
#reorder final data
final_anc_bm_bio17 = final_anc_bm_bio17[,c(1,3,7,2,4,5,6)]
colnames(final_anc_bm_bio17) = c("species", "node_species", "current_value", "node_antecesor", "ace", "ace_low_bound", "ace_upper_bound")
str(final_anc_bm_bio17)
head(final_anc_bm_bio17, 20)


## additional check to test if species that share the last node have the same ancestral state for bio4

#which species had duplicated nodes (i.e. share their last node with another species)
duplicated_nodo1 = which(duplicated(final_anc_bm_bio4$node_antecesor))

#for each pair
test_bio4_bm = NULL
for(i in 1:length(duplicated_nodo1)){

    #select the second species
    second_species = final_anc_bm_bio4[duplicated_nodo1[i],]

    #select the first species
    first_species = final_anc_bm_bio4[which(final_anc_bm_bio4$node_antecesor == second_species$node_antecesor)[1],]

    #check that the ace and CI intervals are the same
    test_bio4_bm = append(test_bio4_bm, ifelse(first_species$ace == second_species$ace & first_species$ace_low_bound == second_species$ace_low_bound & first_species$ace_upper_bound == second_species$ace_upper_bound, TRUE, FALSE))
}
length(test_bio4_bm) == length(duplicated_nodo1)
summary(test_bio4_bm)

## additional check to test if species that share the last node have the same ancestral state for bio17

#which species had duplicated nodes (i.e. share their last node with another species)
duplicated_nodo1 = which(duplicated(final_anc_bm_bio17$node_antecesor))

#bind the row number of species duplicated and their pairs (they are the previous row). Then sort the rows
species_pairs = sort(c(duplicated_nodo1-1, duplicated_nodo1))

#for each pair
test_bio17_bm = NULL
for(i in 1:length(duplicated_nodo1)){

    #select the second species
    second_species = final_anc_bm_bio17[duplicated_nodo1[i],]

    #select the first species
    first_species = final_anc_bm_bio17[which(final_anc_bm_bio17$node_antecesor == second_species$node_antecesor)[1],]

    #check that the ace and CI intervals are the same
    test_bio17_bm = append(test_bio17_bm, ifelse(first_species$ace == second_species$ace & first_species$ace_low_bound == second_species$ace_low_bound & first_species$ace_upper_bound == second_species$ace_upper_bound, TRUE, FALSE))
}
length(test_bio17_bm) == length(duplicated_nodo1)
summary(test_bio17_bm)


#### adittional check of correct selection of ancestral state per node #####
## bio4
#BM
raw_anc_bm_bio4 = pgls_recon_bio4$CI95[which(row.names(pgls_recon_bio4$CI95) %in% final_anc_bm_bio4$node_antecesor),] #select nodes of species from raw dataset
test_results_anc_bm_bio4 = NULL
test_results_current_bm_bio4 = NULL
for(i in 1:length(final_anc_bm_bio4$species)){

    #selected node
    selected_species = final_anc_bm_bio4$species[i]

    #select the row
    selected_row_final_data = final_anc_bm_bio4[which(final_anc_bm_bio4$species == selected_species),]

    #extract anc state from final data.set 
    anc_final_dataset = selected_row_final_data[,which(colnames(final_anc_bm_bio4) %in% c("ace_low_bound", "ace_upper_bound"))]

    #extract anc state from raw data.set     
    anc_raw_dataset = raw_anc_bm_bio4[which(row.names(raw_anc_bm_bio4) == selected_row_final_data$node_antecesor),]

    #test anc
    test_results_anc_bm_bio4 = append(test_results_anc_bm_bio4, anc_final_dataset == anc_raw_dataset)

    #extract current state from final data.set
    current_final_dataset = selected_row_final_data[, which(colnames(final_anc_bm_bio4) %in% c("species", "current_value"))]
    
    #extract current state from raw data.set
    current_raw_dataset = as.vector(bio4_vector[which(names(bio4_vector) == current_final_dataset$species)])

    #test current state
    test_results_current_bm_bio4 = append(test_results_current_bm_bio4, current_final_dataset$current_value == current_raw_dataset)
}
length(test_results_anc_bm_bio4) == 112*2 #two limits per species (95CI of anc test)
length(test_results_current_bm_bio4) == 112 #one current value per species


#OU
raw_anc_ou_bio4 = anc_ou_bio4[which(anc_ou_bio4$nodo1 %in% final_anc_ou_bio4$node_antecesor),] #select nodes of species from raw dataset
test_results_anc_ou_bio4 = NULL
test_results_current_ou_bio4 = NULL
for(i in 1:length(final_anc_ou_bio4$species)){

    #selected node
    selected_species = final_anc_ou_bio4$species[i]

    #select the row
    selected_row_final_data = final_anc_ou_bio4[which(final_anc_ou_bio4$species == selected_species),]

    #extract anc state from final data.set 
    anc_final_dataset = selected_row_final_data[,which(colnames(final_anc_ou_bio4) %in% c("ace", "SE", "Adj.SE"))]

    #extract anc state from raw data.set     
    anc_raw_dataset = raw_anc_ou_bio4[which(raw_anc_ou_bio4$nodo1 == selected_row_final_data$node_antecesor), which(colnames(raw_anc_ou_bio4) %in% c("State", "SE", "Adj.SE"))]

    #test anc
    test_results_anc_ou_bio4 = append(test_results_anc_ou_bio4, anc_final_dataset == anc_raw_dataset)

    #extract current state from final data.set
    current_final_dataset = selected_row_final_data[, which(colnames(final_anc_ou_bio4) %in% c("species", "current_value"))]
    
    #extract current state from raw data.set
    current_raw_dataset = as.vector(bio4_vector[which(names(bio4_vector) == current_final_dataset$species)])

    #test current state
    test_results_current_ou_bio4 = append(test_results_current_ou_bio4, current_final_dataset$current_value == current_raw_dataset)
}
length(test_results_anc_ou_bio4) == 112*3 #two SE varaibles (SE and Adj.SE) per species and one ace (anc test)
length(test_results_current_bm_bio4) == 112 #one current value per species


## bio17
#BM
raw_anc_bm_bio17 = pgls_recon_bio17$CI95[which(row.names(pgls_recon_bio17$CI95) %in% final_anc_bm_bio17$node_antecesor),] #select nodes of species from raw dataset
test_results_anc_bm_bio17 = NULL
test_results_current_bm_bio17 = NULL
for(i in 1:length(final_anc_bm_bio17$species)){

    #selected node
    selected_species = final_anc_bm_bio17$species[i]

    #select the row
    selected_row_final_data = final_anc_bm_bio17[which(final_anc_bm_bio17$species == selected_species),]

    #extract anc state from final data.set 
    anc_final_dataset = selected_row_final_data[,which(colnames(final_anc_bm_bio17) %in% c("ace_low_bound", "ace_upper_bound"))]

    #extract anc state from raw data.set     
    anc_raw_dataset = raw_anc_bm_bio17[which(row.names(raw_anc_bm_bio17) == selected_row_final_data$node_antecesor),]

    #test anc
    test_results_anc_bm_bio17 = append(test_results_anc_bm_bio17, anc_final_dataset == anc_raw_dataset)

    #extract current state from final data.set
    current_final_dataset = selected_row_final_data[, which(colnames(final_anc_bm_bio17) %in% c("species", "current_value"))]
    
    #extract current state from raw data.set
    current_raw_dataset = as.vector(bio17_vector[which(names(bio17_vector) == current_final_dataset$species)])

    #test current state
    test_results_current_bm_bio17 = append(test_results_current_bm_bio17, current_final_dataset$current_value == current_raw_dataset)
}
length(test_results_anc_bm_bio17) == 112*2 #two limits per species (95CI of anc test)
length(test_results_current_bm_bio17) == 112 #one current value per species


#OU
raw_anc_ou_bio17 = anc_ou_bio17[which(anc_ou_bio17$nodo1 %in% final_anc_ou_bio17$node_antecesor),] #select nodes of species from raw dataset
test_results_anc_ou_bio17 = NULL
test_results_current_ou_bio17 = NULL
for(i in 1:length(final_anc_ou_bio17$species)){

    #selected node
    selected_species = final_anc_ou_bio17$species[i]

    #select the row
    selected_row_final_data = final_anc_ou_bio17[which(final_anc_ou_bio17$species == selected_species),]

    #extract anc state from final data.set 
    anc_final_dataset = selected_row_final_data[,which(colnames(final_anc_ou_bio17) %in% c("ace", "SE", "Adj.SE"))]

    #extract anc state from raw data.set     
    anc_raw_dataset = raw_anc_ou_bio17[which(raw_anc_ou_bio17$nodo1 == selected_row_final_data$node_antecesor), which(colnames(raw_anc_ou_bio17) %in% c("State", "SE", "Adj.SE"))]

    #test anc
    test_results_anc_ou_bio17 = append(test_results_anc_ou_bio17, anc_final_dataset == anc_raw_dataset)

    #extract current state from final data.set
    current_final_dataset = selected_row_final_data[, which(colnames(final_anc_ou_bio17) %in% c("species", "current_value"))]
    
    #extract current state from raw data.set
    current_raw_dataset = as.vector(bio17_vector[which(names(bio17_vector) == current_final_dataset$species)])

    #test current state
    test_results_current_ou_bio17 = append(test_results_current_ou_bio17, current_final_dataset$current_value == current_raw_dataset)
}
length(test_results_anc_ou_bio17) == 112*3 #two SE varaibles (SE and Adj.SE) per species and one ace (anc test)
length(test_results_current_bm_bio17) == 112 #one current value per species

#summary
summary(test_results_anc_bm_bio4)
summary(test_results_current_bm_bio4)
summary(test_results_anc_ou_bio4)
summary(test_results_current_bm_bio4)
summary(test_results_anc_bm_bio17)
summary(test_results_current_bm_bio17)
summary(test_results_anc_ou_bio17)
summary(test_results_current_bm_bio17) #ALL test correct, GO ON!

#### Final dataset  ####
str(final_anc_ou_bio4)
head(final_anc_ou_bio4, 20)
str(final_anc_ou_bio17)
head(final_anc_ou_bio17, 20)
str(final_anc_bm_bio4)
head(final_anc_bm_bio4, 20)
str(final_anc_bm_bio17)
head(final_anc_bm_bio17, 20)

## save it
write.table(final_anc_ou_bio4, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_recons/final_anc_ou_bio4.csv", sep=",", col.names = TRUE, row.names = FALSE)
write.table(final_anc_ou_bio17, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_recons/final_anc_ou_bio17.csv", sep=",", col.names = TRUE, row.names = FALSE)
write.table(final_anc_bm_bio4, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_recons/final_anc_bm_bio4.csv", sep=",", col.names = TRUE, row.names = FALSE)
write.table(final_anc_bm_bio17, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_recons/final_anc_bm_bio17.csv", sep=",", col.names = TRUE, row.names = FALSE)

## NOTA: SI HUBIERA QUE COGER EL ESTADO ANCESTRAL EN UNA MISMA FECHA PARA TODAS LAS ESPECIES ANTES DEL ÚLTIMO NODO RESIVA LA FUNCIÓN CONTMAP DE PHYTOOLS, DE AHÍ PUEDES SACAR EL CÓDIGO PARA HACERLO. Mira en "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/phylo/contmap_phytools.R". Si no te aclaras mira este tutoria ("http://www.phytools.org/eqg/Exercise_5.2/"), la parte del traitgram. 

#save image
save.image("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/rdata/tests_phylo.RData")