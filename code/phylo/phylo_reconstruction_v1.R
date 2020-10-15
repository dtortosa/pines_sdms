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

## Selected variables
#-Selección de variables: He usado dos criterios: i) La sumatoria de la posición en el deviance ranking de una variable para todas las especies (más alto indica menos explicativo); ii) Número de especies para las que esa variable está en la posición 1,2 ó 3. Me gusta más el primero, porque tiene peso tanto las especies para las que explica mucho como para las que explica poco, por ejmplo: Una variable puede ser el top1 para muchas especies, pero para el resto no explicar nada, ese sería el caso de clay, que es la variable más explicativa para más especies (26), pero luego su suma de los rankings es el doble respecto de la primera variable (explica poco para muchas especies). Por tanto, la sumatoria de los rankings sería lo más idóneo para seleccionar una variable que se va a usar para reconstruir el estado de TODAS las especies.
    #-Temperatura -> bio4 (Temperature Seasonality). La sumatoria de los rankings es 635 frente a 724 de bio3 (Isothermality), que es la segunda variable más alta de todas. bio3 es el top1,2,3 para más especies (5), pero hay una diferencia considerable en la sumatoria. bio4 se usó para el cluster 3 y 4, mientras que la otras no se usó para ninguno. Además son variables muy parecidas, así que me quedo con la primera según la sumatoria del rank. 
    #-Humedad -> bio17 (humedad del cuarto más seco). Sumatoria de 1208 frente a los 1342 de bio14 (humedad del mes más seco), que es la segunda variable de humedad más alta. En cuanto al número de especies para las que son el top, están muy igualadas. bio17 se usó para el cluster 2 y 3, la otra para ninguno. además son muy parecidas. Me quedo con bio17, que es la que tiene la sumatoria de nrakings más alta.
    #-Nota: la info sobre las variables usadas para cada cluster está en "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/variable_selection_inside_clusters/variable_selection_inside_clusters.R".

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

#drop discolor (problem tazonomy, no diferetiaced from cembriodes)
if("discolor" %in% epithet_species_list){
    epithet_species_list = epithet_species_list[-which(epithet_species_list == "discolor")]
}

#check it
!"discolor" %in% epithet_species_list
length(epithet_species_list) == 112

#function to calculate SE
se <- function(x) sd(x)/sqrt(length(x))

#unzip the selected variables bio4 and bio17
unzip("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals.zip", files = c("finals/bio4.asc", "finals/bio17.asc"), list=FALSE, exdir = "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/phlyo/climatic_data_phylo", junkpaths = TRUE) #junkpaths = TRUE indicated that we don't want all directoy and subdirectories, only the files. 

####################
##### bio4 #########
####################
bio4 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/phlyo/climatic_data_phylo/bio4.asc")
res(bio4)

#Extract bio4 from species distribution
require(raster)
median_bio4 = NULL
sd_bio4 = NULL
se_bio4 = NULL
bio4_pines = stack() 
for (i in 1:length(epithet_species_list)){

    #select the species
    selected_epithet = epithet_species_list[i]

    #unzip species distribution with buffer
    buffer_path = unzip("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences.zip",  list=FALSE, exdir = "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/phlyo/buffers", files = paste("ocurrences/", selected_epithet, "_distribution_buffer.asc", sep=""), junkpaths = TRUE)

    #load it
    distri = raster(buffer_path)

    #create a polygon from distributon
    polygon = rasterToPolygons(distri, fun=function(x){x==1}, n=16) #esta funcion de raster te transforma un raster completo o una parte del mismo en un poliogno. En nuestro caso solo queremos las celdas con valor=1, es decir, presencias. Por eso ponemos x==1. 

    #crop bio4
    bio4_cropped = mask(bio4, polygon)

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

#save rasters
#writeRaster(bio4_pines, filename="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/phlyo/climatic_data_phylo/bio4_buffer.tif", options="INTERLEAVE=BAND", overwrite=TRUE)



####################
##### bio17 ########
####################
bio17 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/phlyo/climatic_data_phylo/bio17.asc")
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

    #unzip species distribution with buffer
    buffer_path = unzip("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences.zip",  list=FALSE, exdir = "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/phlyo/buffers", files = paste("ocurrences/", selected_epithet, "_distribution_buffer.asc", sep=""), junkpaths = TRUE)

    #load it
    distri = raster(buffer_path)

    #create a polygon from distributon
    polygon = rasterToPolygons(distri, fun=function(x){x==1}, n=16) #esta funcion de raster te transforma un raster completo o una parte del mismo en un poliogno. En nuestro caso solo queremos las celdas con valor=1, es decir, presencias. Por eso ponemos x==1. 

    #crop bio17
    bio17_cropped = mask(bio17, polygon)

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

#save rasters
#writeRaster(bio17_pines, filename="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/phlyo/climatic_data_phylo/bio17_buffer.tif", options="INTERLEAVE=BAND", overwrite=TRUE)


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
text(x=14000, y=-320, labels = estimate_no_pic, cex=1)
p_no_pic = bquote(italic(p.value) == .(format(round(test_no_pic$p.value,4))))
text(x=14000, y=-400, labels = p_no_pic, cex=1)

#YES PIC
plot(pic.bio17~pic_bio4, xlab="PIC Median BIO4", ylab="PIC Median BIO17")
estimate_pic = bquote(italic(rho) == .(format(round(test_pic$estimate,2))))
text(x=-1500, y=70, labels = estimate_pic, cex=1)
p_no_pic = bquote(italic(p.value) == .(format(round(test_pic$p.value,4))))
text(x=-1500, y=50, labels = p_no_pic, cex=1)

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
phylosig(tree_prunned, bio4_vector, method="lambda", test=TRUE) #phytools: 0.44
fitContinuous(phy = tree_prunned, dat = bio4_vector, model = "lambda", control = list(niter = 100, CI = 0.95), SE=0) #diversitree without intravariability: 0.44
fitContinuous(phy = tree_prunned, dat = bio4_vector, model = "lambda", control = list(niter = 100, CI = 0.95), SE=intra_var_bio4) #diversitree with intravariability: 0.44

#k
phylosig(tree_prunned, bio4_vector, method="K", nsim=200000, test=TRUE) #0.1556661

#plot under BM
obj = contMap(tree_prunned, bio4_vector)
plot(obj, type="fan")

## bio17
#lambda 
phylosig(tree_prunned, bio17_vector, method="lambda", test=TRUE) #phytools: 0.94
fitContinuous(phy = tree_prunned, dat = bio17_vector, model = "lambda", control = list(niter = 100, CI = 0.95), SE=0) #diversitree without intravariability: 0.94
fitContinuous(phy = tree_prunned, dat = bio17_vector, model = "lambda", control = list(niter = 100, CI = 0.95), SE=intra_var_bio17) #diversitree with intravariability: 0.94

#k
phylosig(tree_prunned, bio17_vector, method="K", nsim=10000, test=TRUE) #0.3076953

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
print(res_bio4) #Best model is OU with a difference of AICc of 15 with the second best model, which is lambda (with and without SE (intraespecific variation)). The last model is BM but close to white noise and lambda (11). 
    #alpha = 0.076650

## bio17
res_bio17=fitClim(trait = "median_bio17")
print(res_bio17) #Best model is OU with a difference of AICc of 10 with the second best model, which is BM (with and without SE (intraespecific variation)). The last model is white noise, very far away (63). Lambda is close to BM (more AICc but only 1 unit). 
    #alpha = 0.027591
#alpha values are equal with and without SE of intraespecific variability.


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
t_1_2_bio4 = log(2)/0.08
t_1_2_bio17 = log(2)/0.03

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
    #OU:Luego Exponential model, specyfing Alpha, pon el valor de alfa de bio4  ó bio17 redondeados a dos decimales (0.08 y 0.03 respectivamente), 100 iteraciones y run (he comprobado el resultado con 1000 interaciones en ambas variablws y sale exactamente lo mismo). Asú se corre un OU en comapre4.6. Esto se ha seguido de Guerrero et al.,... Wiens ., 2013 ("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3710863/")

#Le das a guardar en un archivo llamado "compare_res_bioX_XX.txt". De ese archivo copias la parte de "Trait #1: Ancestral state estimates" y la pegas en un excel, le das a pegar con el importador de datos (se hace en el boton que surge al pegar como el de mantener-quitar formato). Así te separará cada columna. Solo falta añadir a "Adj." el "SE" que queda en la siguiente columna (es SE adjusted) y guardar como .csv.

#load results of OU with SE intraespecífica (alpha = 0.07 for bio4 and 0.03 for bio17)
anc_ou_bio4 = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/compare_results/anc_bio4_ou.csv", sep=",", header=TRUE)
str(anc_ou_bio4)
anc_ou_bio17 = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/compare_results/anc_bio17_ou.csv", sep=",", header=TRUE)
str(anc_ou_bio17)

#change names of root by 112+1 (ass notation of ape)
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
pgls_ou_bio4 = reconstruct(bio4_vector, tree_prunned, method = "GLS_OU", alpha = 0.076650, CI = TRUE)
pgls_ous_bio4 = reconstruct(bio4_vector, tree_prunned, method = "GLS_OUS", alpha = 0.076650, CI = TRUE) #"GLS_OU" and "GLS_OUS" differs in the fact that "GLS_OUS" assume that the process starts from the optimum, while the root state has to be estimated for "GLS_OU", which may rise some issues (see reconstruct 219 Royer-Carenzi and Didier, 2016). Users may provide the attractive strength parameter alpha, for these two models. Users may provide the attractive strength parameter alpha, for these two models. "GLS_ABM", "GLS_OU" and "GLS_OUS" are all fitted by generalized least squares (Royer-Carenzi and Didier, 2016).
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
pgls_ou_bio17 = reconstruct(bio17_vector, tree_prunned, method = "GLS_OU", alpha = 0.027591, CI = TRUE)
pgls_ous_bio17 = reconstruct(bio17_vector, tree_prunned, method = "GLS_OUS", alpha = 0.027591, CI = TRUE) #Mas alto ape

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


#calculate quantiles for the ancestrlaa state of bio4 and the CI
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
obj_bm_bio4 = contMap(tree_prunned_to_plot, bio4_vector_to_plot, method="user", anc.states=pgls_recon_bio4$ace, plot=FALSE, lwd = 2) #aquí anc.states is a file with only ace, without 95CI.

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
obj_ou_bio4 = contMap(tree_prunned_to_plot, bio4_vector_to_plot, method="user", anc.states=anc_ou_bio4_to_plot$ace, plot=FALSE, lwd = 2) #aquí anc.states is a file with only ace, without 95CI.

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
obj_bm_bio17 = contMap(tree_prunned_to_plot, bio17_vector_to_plot, method="user", anc.states=pgls_recon_bio17$ace, plot=FALSE, lwd = 2, res=10) #aquí anc.states is a file with only ace, without 95CI.

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
obj_ou_bio17 = contMap(tree_prunned_to_plot, bio17_vector_to_plot, method="user", anc.states=anc_ou_bio17_to_plot$ace, plot=FALSE, lwd = 2) #aquí anc.states is a file with only ace, without 95CI.

#open the pdf of the plot
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/plot_anc/plot_ancestral_states_ou_bio17.pdf", width = 24, height = 24)

#plot
plot(obj_ou_bio17, type="fan", fsize=c(2.5, 2.3)) #fsize tiene el primer valor del tamaño de las tip labels, el segundo es para la legenda.

#add the error bars
errorbar_contMap_modified(obj=obj_ou_bio17, user=TRUE, anc.states=anc_ou_bio17_to_plot, scale.by.ci=TRUE) #aquí anc.states is the complete file with ace and 95CI.
    #scale.by.ci=TRUE: that determines whether or not the length of the error bars will be scaled by the CI width
dev.off()

#open the pdf of the plot
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/plot_anc/plot_ancestral_states_bm_REML_bio17.pdf", width = 24, height = 24)

#create the tree interpolated
obj_bm_bio17 = contMap(tree_prunned_to_plot, bio17_vector_to_plot, method="user", plot=FALSE, lwd = 2, anc.states=reml_recon_bio17$ace) #aquí anc.states is a file with only ace, without 95CI.

#plot
plot(obj_bm_bio17, type="fan", fsize=c(2.5, 2.3)) #fsize tiene el primer valor del tamaño de las tip labels, el segundo es para la legenda.

#add the error bars
errorbar_contMap_modified(obj=obj_bm_bio17, user=TRUE, anc.states=reml_recon_bio17, scale.by.ci=TRUE) #aquí anc.states is the complete file with ace and 95CI.
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
text(x=-50, y=-300, labels = mylabel, cex=1)
if(!spearman$p.value == 0){ #if p.value is <2.2e-16, spearman save it as 0
    mylabel.p = bquote(italic(p.value) == .(format(spearman$p.value, digits = 3)))
} else {
    mylabel.p = bquote(italic(p.value) < .(format(2.2e-16)))
}
text(x=-50, y=-340, labels = mylabel.p, cex=1)
mylabel.t = bquote(italic(S) == .(format(spearman$statistic, digits = 3)))
text(x=-50, y=-380, labels = mylabel.t, cex=1)

##pannel 2
if(FALSE){
    #densitplot of differences between ancestral states
    plot(density(abs(anc_ou_bio17$State - pgls_recon_bio17$ace)), main="Absolute difference ancestral states")
    p<-wilcox.test(anc_ou_bio17$State, pgls_recon_bio17$ace)$p.value
    mylabel.p = bquote(italic(p.value) == .(format(p, digits = 3)))
    text(x=3000, y=0.0005, labels = mylabel.p, cex=1)
} 

##pannel 3
#SE of OU versus CI95 of BM
plot(anc_ou_bio17$SE~abs(pgls_recon_bio17$CI95[,1] - pgls_recon_bio17$CI95[,2]), xlab="Standard errors under BM", ylab="Standard errors under OU", main="")
#cors
spearman<-cor.test(anc_ou_bio17$SE, abs(pgls_recon_bio17$CI95[,1] - pgls_recon_bio17$CI95[,2]), method="spearman")
mylabel = bquote(italic(rho) == .(format(spearman$estimate, digits = 3)))
text(x=11, y=80, labels = mylabel, cex=1)
if(!spearman$p.value == 0){ #if p.value is <2.2e-16, spearman save it as 0
    mylabel.p = bquote(italic(p.value) == .(format(spearman$p.value, digits = 3)))
} else {
    mylabel.p = bquote(italic(p.value) < .(format(2.2e-16)))
}
text(x=11, y=65, labels = mylabel.p, cex=1)
mylabel.t = bquote(italic(S) == .(format(spearman$statistic, digits = 3)))
text(x=11, y=50, labels = mylabel.t, cex=1)


##pannel 4
if(FALSE){
    #densitplot of differences between error of ancestral state estimates
    plot(density(abs(anc_ou_bio17$SE - abs(pgls_recon_bio17$CI95[,1] - pgls_recon_bio17$CI95[,2]))), main="Absolute difference standard errors")
    p<-wilcox.test(anc_ou_bio17$SE, abs(pgls_recon_bio17$CI95[,1] - pgls_recon_bio17$CI95[,2]))$p.value
    mylabel.p = bquote(italic(p.value) == .(format(p, digits = 3)))
    text(x=4800, y=0.0004, labels = mylabel.p, cex=1)
}
dev.off()

######Conclusions
#Bio4: Como era de esperar por su valor de lambda tan alto (0.94), bio17 apenas muestra diferencias entre modelos. Los estados ancestrales están muy correlacionados (rho=0.98)  y no son significativamente diferentes (p = 0.4; panel 1 y 2). Los errores también están muy correlacionados, pero si hay diferencias entre modelos, siendo los errores de OU más grandes (esto creo que se debe más al programa [compare vs ace] que a los modelos). 
#Bio17: En bio4 está la cosa un poco menos clara, lo cual es lógico viendo el valor de lambda no tan alto para esta variable (0.47). Sigue habiendo una correlación alta entre estados ancestrales (rho = 0.88), pero sí hay diferencias significativas entre modelos aunque con una p no muy baja (p = 0.03). Los errores de las estimaciones son significativamente muuucho más grandes en OU.
#Yo creo que hay una correlación clara entre ambos modelos, se parecen mucho, y es lógico teniendo en cuenta los valores de lambda altos (sobre todo en bio17) y los valores de alfa tan bajos (0.08 para bio4 y 0.03 para bio17). Sin embargo, me preocupan un poco las diferencias en valor absoluto de los errores estándar para bio4, que al final es el parámetro que usamos y nos está dando resultados muy diferentes. La verdad es que no lo tengo claro, pero voy a seguir para adelante con BM. en caso de que hubiese que dar marcha atrás no va a ser mucho problema de trabajo, solo un poquito de tiempo para correr de nuevo los análisis.
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

## bind current climatic data
final_anc_ou_bio4 = merge(final_anc_ou_bio4, climate_medians_merge_bio4, by="species_ramas")
final_anc_ou_bio17 = merge(final_anc_ou_bio17, climate_medians_merge_bio17, by="species_ramas")

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


#######################################################
################ EXAMPLE P.HALEPENSIS #################
#######################################################

###########################################################################
######### Extraction and preparation of climatic and species data #########
###########################################################################

#load packages
library(raster)
library(rgeos)

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

#load bio4 currently
bio4 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/phlyo/climatic_data_phylo/bio4.asc")
res(bio4)

#load bio17 currently
bio17 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/phlyo/climatic_data_phylo/bio17.asc")
res(bio17)

#plot bio4 and distribution polygons
#plot(crop(bio4, polygon_plot_buffer))
#plot(polygon_high_precision_points, lty=1, add=TRUE)
#plot(ocurrences_buffer_polygon, lty=3, add=TRUE)

#plot bio17 and distribution polygons
#plot(crop(bio17, polygon_plot_buffer))
#plot(polygon_high_precision_points, lty=1, add=TRUE)
#plot(ocurrences_buffer_polygon, lty=3, add=TRUE)

#list of continuos projections and binary projections for all scenarios and climatic models
climatic_scenarios = c("bc26", "bc45", "bc60", "bc85", "cc26", "cc45", "cc60", "cc85", "gs26", "gs45", "gs60", "gs85", "he26", "he45", "he60", "he85", "ip26", "ip45", "ip60", "ip85", "mg26", "mg45", "mg60", "mg85", "mr26", "mr45", "mr60", "mr85")

#load projections of variables for all scenarios
stack_bio17 = stack(paste("/Users/diegosalazar/phd_big_documents/pines_niche/climate_proj_phylo/bio17_", climatic_scenarios, ".asc", sep=""))
stack_bio4 = stack(paste("/Users/diegosalazar/phd_big_documents/pines_niche/climate_proj_phylo/bio4_", climatic_scenarios, ".asc", sep=""))

## load ancestral reconstruction
final_anc_ou_bio4 = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_recons/final_anc_ou_bio4.csv", sep=",", header=TRUE)
final_anc_ou_bio17 = read.table( "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_recons/final_anc_ou_bio17.csv", sep=",", header=TRUE)
final_anc_bm_bio4 = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_recons/final_anc_bm_bio4.csv", sep=",", header=TRUE)
final_anc_bm_bio17 = read.table( "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_recons/final_anc_bm_bio17.csv", sep=",", header=TRUE)

#take a look
str(final_anc_ou_bio4)
str(final_anc_ou_bio17)
str(final_anc_bm_bio4)
str(final_anc_bm_bio17)

#list of evolution models: ONLY BM
#list_models_bio4 = c("final_anc_ou_bio4", "final_anc_bm_bio4")
#list_models_bio17 = c("final_anc_ou_bio17", "final_anc_bm_bio17")
list_models_bio4 = c("final_anc_bm_bio4")
list_models_bio17 = c("final_anc_bm_bio17")

#set the species
species = "halepensis"

#load buffer of pseudo absences (broader) to crop climatic variables
PA_buffer = raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/pa_buffers/", species, "_PA_buffer.asc", sep="")) #usamos PA buffer que es más amplio para la ver el rango filo por si salen sitios interesantes lejos para migración asistida, pero los plots en general se hacen con un area más pequeña. 

#crop climatic variables with PA buffer
stack_bio17_cropped = crop(stack_bio17, PA_buffer)
stack_bio4_cropped = crop(stack_bio4, PA_buffer)

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
polygon_plot_buffer = gBuffer(ocurrences_buffer_polygon, byid=FALSE, id=NULL, width=15, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0) 

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

#crop future sutiabiity aroun the whole PA buffer not for plotting but for analsis
ensamble_suitability = crop(projected_suit, PA_buffer)

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

#create a spatial point data frame
coors_high_precision_points = high_precision_ocurrences[,which(colnames(high_precision_ocurrences) %in% c("longitude", "latitude"))]
coordinates(coors_high_precision_points) = c("longitude", "latitude")

#extrac id cells from that spatial point data frame
id.cell <- extract(ocurrences_buffer_raster, coors_high_precision_points, cellnumbers=TRUE)[,1]

#create a polygon with those id cell, which will include all high precision points
raster_high_precision_points = raster(extent(ocurrences_buffer_raster), resolution = res(ocurrences_buffer_raster)) #empty raster with the same extent and res than ocurrences_buffer_raster
raster_high_precision_points[id.cell] <- ocurrences_buffer_raster[id.cell] #fill the raster with the values from the cells ocurrences_buffer_raster with presence of the species
polygon_high_precision_points = rasterToPolygons(raster_high_precision_points, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to polygon

#loop for comparing projected climate and phylo range
phylo_rasters_bio17 = stack()
phylo_rasters_bio17_proportion = stack()
for(s in 1:length(climatic_scenarios)){

    #select the climatic scenario
    selected_scenario = climatic_scenarios[s]

    #create a empty raster with the same extent and resolution of the raster layer of the [s] IPCC scenario
    raster_subsetted = raster(extent(stack_bio17_cropped[[s]]), resolution = res(stack_bio17_cropped[[s]]))

    #add to the empty raster those cells of the [s] IPCC raster in which the habitat suitability is higher than 25 and lower than 75
    raster_subsetted[which(getValues(ensamble_suitability) > 25 & getValues(ensamble_suitability) < 75)] <- stack_bio17_cropped[[s]][which(getValues(ensamble_suitability) > 25 & getValues(ensamble_suitability) < 75)] #suitability map has the same resolution and extent than IPCC maps because the suitability map was obtained from these raster of ICPP scenarios. 

    #for each evolution model 
    for(m in 1:length(list_models_bio17)){

        #select the [m] evolution model 
        selected_model = list_models_bio17[m]

        #extract data of [m] model
        model = get(selected_model)

        #select the row of the corresponding species
        model = model[which(model$species == paste("Pinus_", species, sep="")),]

        #extract all cell values from the raster with climatic data of the [s] scenario only in those areas with uncertainty (raster_subsetted)
        cell_values = getValues(raster_subsetted)

        #extract ID of those cells without NA
        cells_withot_NA = which(!is.na(cell_values))

        #extract, from all cells withput NA, those whose value is inside the phylogenetic range (including the current value but not including the ancestral). For that we used is.between function, created by me. 
        cell_inside_phylo_range = which(is.between(cell_value = na.omit(getValues(raster_subsetted)), ancestral = model$ace, current = model$current_value))

        #from ID of cells without NA, select the ID of those whose vale is inside of the phylo range
        final_cells = cells_withot_NA[cell_inside_phylo_range]

        #create a empty raster with the same extent and resolution than the [s] IPCC raster
        final_raster = raster(extent(stack_bio17_cropped[[s]]), resolution = res(stack_bio17_cropped[[s]]))
        final_raster_proportion = raster(extent(stack_bio17_cropped[[s]]), resolution = res(stack_bio17_cropped[[s]]))

        #fill the raster with zeros
        final_raster[] <- 0
        final_raster_proportion[] <- 0

        #add to these final cells a value of suitability without and with proportion
        final_raster[final_cells] <- 1
        final_raster_proportion[final_cells] <- phylo_proportion(x=raster_subsetted[final_cells], ancestral_value=model$ace, current_value=model$current_value)

        #add the name of the raster
        names(final_raster) <- paste(selected_scenario, "_", strsplit(selected_model, split="_")[[1]][3], "_", strsplit(selected_model, split="_")[[1]][4], sep="")
        names(final_raster_proportion) <- paste(selected_scenario, "_", strsplit(selected_model, split="_")[[1]][3], "_", strsplit(selected_model, split="_")[[1]][4], sep="")        

        #save the raster into a stack
        phylo_rasters_bio17 = stack(phylo_rasters_bio17, final_raster)
        phylo_rasters_bio17_proportion = stack(phylo_rasters_bio17_proportion, final_raster_proportion)
    }    
}

#check that all scenarios has been included
nlayers(phylo_rasters_bio17) == length(climatic_scenarios)
nlayers(phylo_rasters_bio17_proportion) == length(climatic_scenarios)

#sum suitability across IPCC scenarios for bio17
sum_phylo_bio17 = calc(phylo_rasters_bio17, function(x) (sum(x)))

#loop for comparing projected climate and phylo range
phylo_rasters_bio4 = stack()
phylo_rasters_bio4_proportion = stack()
for(s in 1:length(climatic_scenarios)){

    #select the climatic scenario
    selected_scenario = climatic_scenarios[s]

    #create a empty raster with the same extent and resolution of the raster layer of the [s] IPCC scenario
    raster_subsetted = raster(extent(stack_bio4_cropped[[s]]), resolution = res(stack_bio4_cropped[[s]]))

    #add to the empty raster those cells of the [s] IPCC raster in which the habitat suitability is higher than 25 and lower than 75
    raster_subsetted[which(getValues(ensamble_suitability) > 25 & getValues(ensamble_suitability) < 75)] <- stack_bio4_cropped[[s]][which(getValues(ensamble_suitability) > 25 & getValues(ensamble_suitability) < 75)] #suitability map has the same resolution and extent than IPCC maps because the suitability map was obtained from these raster of ICPP scenarios. 

    #for each evolution model 
    for(m in 1:length(list_models_bio4)){

        #select the [m] evolution model 
        selected_model = list_models_bio4[m]

        #extract data of [m] model
        model = get(selected_model)

        #select the row of the corresponding species
        model = model[which(model$species == paste("Pinus_", species, sep="")),]

        #extract all cell values from the raster with climatic data of the [s] scenario only in those areas with uncertainty (raster_subsetted)
        cell_values = getValues(raster_subsetted)

        #extract ID of those cells without NA
        cells_withot_NA = which(!is.na(cell_values))

        #extract, from all cells withput NA, those whose value is inside the phylogenetic range (including the current value but not including the ancestral). For that we used is.between function, created by me. 
        cell_inside_phylo_range = which(is.between(cell_value = na.omit(getValues(raster_subsetted)), ancestral = model$ace, current = model$current_value))

        #from ID of cells without NA, select the ID of those whose vale is inside of the phylo range
        final_cells = cells_withot_NA[cell_inside_phylo_range]

        #create a empty raster with the same extent and resolution than the [s] IPCC raster
        final_raster = raster(extent(stack_bio4_cropped[[s]]), resolution = res(stack_bio4_cropped[[s]]))
        final_raster_proportion = raster(extent(stack_bio4_cropped[[s]]), resolution = res(stack_bio4_cropped[[s]]))

        #fill the raster with zeros
        final_raster[] <- 0
        final_raster_proportion[] <- 0

        #add to these final cells a value of suitability without and with proportion
        final_raster[final_cells] <- 1
        final_raster_proportion[final_cells] <- phylo_proportion(x=raster_subsetted[final_cells], ancestral_value=model$ace, current_value=model$current_value)

        #add the name of the raster
        names(final_raster) <- paste(selected_scenario, "_", strsplit(selected_model, split="_")[[1]][3], "_", strsplit(selected_model, split="_")[[1]][4], sep="")
        names(final_raster_proportion) <- paste(selected_scenario, "_", strsplit(selected_model, split="_")[[1]][3], "_", strsplit(selected_model, split="_")[[1]][4], sep="")        

        #save the raster into a stack
        phylo_rasters_bio4 = stack(phylo_rasters_bio4, final_raster)
        phylo_rasters_bio4_proportion = stack(phylo_rasters_bio4_proportion, final_raster_proportion)
    }    
}

#check that all scenarios has been included
nlayers(phylo_rasters_bio4) == length(climatic_scenarios)
nlayers(phylo_rasters_bio4_proportion) == length(climatic_scenarios)

#sum suitability across IPCC scenarios for bio4
sum_phylo_bio4 = calc(phylo_rasters_bio4, function(x) (sum(x)))

#calculate intersection between sum both suitaiblity maps (bio17, bio4)
intersection_ensamble_phylo = sum_phylo_bio17 * sum_phylo_bio4 #this will be used to exclude areas not shared between variables from the final ensamble. 

#bind bio4 and bio17 rasters without proportions
phylo_rasters = stack(phylo_rasters_bio4, phylo_rasters_bio17)
nlayers(phylo_rasters) == length(climatic_scenarios)*2

#bind bio4 and bio17 rasters with proportions
phylo_rasters_proportion = stack(phylo_rasters_bio4_proportion, phylo_rasters_bio17_proportion)
nlayers(phylo_rasters_proportion) == length(climatic_scenarios)*2

#calculate the proportion of cells that fall inside the differents phylo ranges across IPCC scenarios and global circulation models without and with proportions
ensamble_phylo = calc(phylo_rasters, function(x) (sum(x)*0.5)/nlayers(phylo_rasters)) #En tanto por 0.5, para luego sumarle 0.5 y que todos los valores estén por encima de 0.5. Así ganamos contraste y los mapas son comparables entre especies. 
ensamble_phylo_proportion = calc(phylo_rasters_proportion, function(x) (sum(x)*0.5)/nlayers(phylo_rasters_proportion)) #En tanto por 0.5, para luego sumarle 0.5 y que todos los valores estén por encima de 0.5. Así ganamos contraste y los mapas son comparables entre especies. 

#sumamos 0.5 para tener valores de idoneidad de 0.5 a 1 (mayor contraste)
ensamble_phylo = calc(ensamble_phylo, function(x) (x+0.5))
ensamble_phylo_proportion = calc(ensamble_phylo_proportion, function(x) (x+0.5))

#ponemos como cero los casos con el valor mínimo (sería el 0 en el raster inicial antes de sumar), para que así queden transparentes
second_ensamble_phylo = ensamble_phylo
second_ensamble_phylo[which(getValues(second_ensamble_phylo) == min(getValues(second_ensamble_phylo)))] <- 0
second_ensamble_phylo_proportion = ensamble_phylo_proportion
second_ensamble_phylo_proportion[which(getValues(second_ensamble_phylo_proportion) == min(getValues(second_ensamble_phylo_proportion)))] <- 0

#check that raster withput proportians has a equal or higher suitability than the raster with proportions
median(getValues(second_ensamble_phylo)[which(!getValues(second_ensamble_phylo)==0)]) >= median(getValues(second_ensamble_phylo_proportion)[which(!getValues(second_ensamble_phylo_proportion)==0)])

#select only those areas suitable for at least one scenario for each variables in the suitability without proportions
#bio17 and bio4 suitaiblity don't overlap in areas with zero in intersection_ensamble_phylo
second_ensamble_phylo[intersection_ensamble_phylo == 0] <- 0
second_ensamble_phylo_proportion[intersection_ensamble_phylo == 0] <- 0

#save second ensables
writeRaster(second_ensamble_phylo, paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/phylo_ensamble/without_proportions/", species, "_phylo_ensamble_without_proportions.asc", sep=""), overwrite=TRUE)
writeRaster(second_ensamble_phylo_proportion, paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/phylo_ensamble/with_proportions/", species, "_phylo_ensamble_with_proportions.asc", sep=""), overwrite=TRUE)

#creamos un raster con todas las áreas que caen en algún ranog filo para algún escenario para el raster sin proporciones. PARA plotear
final_ensamble_phylo = raster(extent(second_ensamble_phylo), resolution = res(second_ensamble_phylo))
final_ensamble_phylo[]<-0
final_ensamble_phylo[which(getValues(second_ensamble_phylo) > 0)] <- 1

#creamos un raster con todas las áreas que caen en algún ranog filo para algún escenario para el raster sin proporciones. PARA plotear. 
final_ensamble_phylo_proportion = raster(extent(second_ensamble_phylo_proportion), resolution = res(second_ensamble_phylo_proportion))
final_ensamble_phylo_proportion[]<-0
final_ensamble_phylo_proportion[which(getValues(second_ensamble_phylo_proportion) > 0)] <- 1

###plot final con todos los paneles sin usar proporciones
#ploteamos seas áreas sonbre la idoneidad de hábitat, pero con un nivel de transparencia dependiente del número casos en los que la celda ha caído dentro del rango filo (incertidumbre; alpha=second_ensamble_phylo). Por tanto, aquellas zonas que han caído dentro del rango para pocos scenarios ó solo para una de las varuables se ven poco (podría verse con solo una variable dentro del rango sin entrease en el rango bajo muuuchos escenarios).
pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/final_figures/suple/without_proportions/", species, "_without_proportions.pdf", sep=""), width=12, height = 12)
par(oma=c(0,0,2.7,2))
par(mfcol=c(2,2))

###Pannel 1
plot(current_suit,  main="", axis.args=list(cex.axis=1.5))
mtext(text="Predicted habitat suitability currently", side=3,line=2, outer=FALSE, cex=1.4, font = 2)
#add acuatic bodies
plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
#plot low precision points
if(nrow(low_precision_ocurrences)>0){
    points(low_precision_ocurrences$longitude, low_precision_ocurrences$latitude, cex=0.5, col="white", bg=NA, lwd=0.6, type="p", pch=21)
}
#plot high precision points
if(nrow(high_precision_ocurrences)>0){
    points(high_precision_ocurrences$longitude, high_precision_ocurrences$latitude, cex=0.5, col="black", bg=NA, lwd=0.2, type="p", pch=21)
} 
#plot the legend    
legend("top", legend=c("High precision points", "Low precision points"), fill=c("black", "white"), cex=1.25, bty="o", horiz=TRUE, bg="white")

###Pannel 2
plot(current_suit,  main="", axis.args=list(cex.axis=1.5))
mtext(text="Predicted habitat suitability currently", side=3,line=2, outer=FALSE, cex=1.4, font = 2)
#add acuatic bodies
plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
#add lines around high precision points
plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=2)
#add polygon of cirtifield distribution + buffer
plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=2)
#plot the legend    
legend("top", legend=c("High precision points", "Low precision points"), lty=c(1,3), cex=1.25, bty="o", lwd=2, horiz=TRUE, bg="white")

###Pannel 3
plot(projected_suit, main="", , axis.args=list(cex.axis=1.5))
mtext(text="Predicted habitat suitability in 2070", side=3,line=2, outer=FALSE, cex=1.4, font = 2)
#add acuatic bodies
plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
#add lines around high precision points
plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=2)
#add polygon of cirtifield distribution + buffer
plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=2)
#plot the legend    
legend("top", legend=c("High precision points", "Low precision points"), lty=c(1,3), cex=1.25, bty="o", lwd=2, horiz=TRUE, bg="white")

###Pannel 4
plot(projected_suit, main="", axis.args=list(cex.axis=1.5))
mtext(text="Predicted habitat suitability in 2070", side=3, line=3, outer=FALSE, cex=1.4, font = 2)
mtext(text="+", side=3,line=2.2, outer=FALSE, cex=1, font = 2)
mtext(text="Phylo-corrected suitability", side=3,line=1, outer=FALSE, cex=1.4, font = 2)
#add acuatic bodies
plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
#add lines around high precision points
plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=2)
#add polygon of cirtifield distribution + buffer
plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=2)
#add info phlo
plot(final_ensamble_phylo, col="#00ffe9", alpha=second_ensamble_phylo, add=TRUE, legend=FALSE) 
    #col="#0000FF": color azul oscuro; ALTERNATIVA el negro ("#0000FF"). 
    #alpha=second_ensamble_phylo: transparencia dependiente de la incertidumbre (número casos en los que la celda ha caído dentro del rango filo)
    #add=TRUE: Add to the previous plot
    #legend=FALSE: Not add legend, in that way the legen of the first plot will not be hidden
#plot the legend of precision    
legend("top", legend=c("High precision points", "Low precision points"), lty=c(1,3), cex=1.25, bty="o", lwd=2, horiz=TRUE, bg="white")
#add legend of acnestral
legend("bottom", legend="Phylo-corrected suitability", fill="#00ffe9", bty="o", cex=1.25, bg="white")

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
plot(current_suit,  main="", axis.args=list(cex.axis=1.5))
mtext(text="Predicted habitat suitability currently", side=3,line=2, outer=FALSE, cex=1.4, font = 2)
#add acuatic bodies
plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
#plot low precision points
if(nrow(low_precision_ocurrences)>0){
    points(low_precision_ocurrences$longitude, low_precision_ocurrences$latitude, cex=0.5, col="white", bg=NA, lwd=0.6, type="p", pch=21)
}
#plot high precision points
if(nrow(high_precision_ocurrences)>0){
    points(high_precision_ocurrences$longitude, high_precision_ocurrences$latitude, cex=0.5, col="black", bg=NA, lwd=0.2, type="p", pch=21)
} 
#plot the legend    
legend("top", legend=c("High precision points", "Low precision points"), fill=c("black", "white"), cex=1.25, bty="o", horiz=TRUE, bg="white")

###Pannel 2
plot(current_suit,  main="", axis.args=list(cex.axis=1.5))
mtext(text="Predicted habitat suitability currently", side=3,line=2, outer=FALSE, cex=1.4, font = 2)
#add acuatic bodies
plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
#add lines around high precision points
plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=2)
#add polygon of cirtifield distribution + buffer
plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=2)
#plot the legend    
legend("top", legend=c("High precision points", "Low precision points"), lty=c(1,3), cex=1.25, bty="o", lwd=2, horiz=TRUE, bg="white")

###Pannel 3
plot(projected_suit, main="", , axis.args=list(cex.axis=1.5))
mtext(text="Predicted habitat suitability in 2070", side=3,line=2, outer=FALSE, cex=1.4, font = 2)
#add acuatic bodies
plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
#add lines around high precision points
plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=2)
#add polygon of cirtifield distribution + buffer
plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=2)
#plot the legend    
legend("top", legend=c("High precision points", "Low precision points"), lty=c(1,3), cex=1.25, bty="o", lwd=2, horiz=TRUE, bg="white")

###Pannel 4
plot(projected_suit, main="", axis.args=list(cex.axis=1.5))
mtext(text="Predicted habitat suitability in 2070", side=3, line=3, outer=FALSE, cex=1.4, font = 2)
mtext(text="+", side=3,line=2.2, outer=FALSE, cex=1, font = 2)
mtext(text="Phylo-corrected suitability", side=3,line=1, outer=FALSE, cex=1.4, font = 2)
#add acuatic bodies
plot(aquatic_bodies, col = "steelblue2", add=TRUE, legend=FALSE)
#add lines around high precision points
plot(polygon_high_precision_points, lty=1, add=TRUE, lwd=2)
#add polygon of cirtifield distribution + buffer
plot(ocurrences_buffer_polygon, lty=3, add=TRUE, lwd=2)
#add info phlo
plot(final_ensamble_phylo_proportion, col="#00ffe9", alpha=second_ensamble_phylo_proportion, add=TRUE, legend=FALSE) 
    #col="#0000FF": color azul oscuro; ALTERNATIVA el negro ("#0000FF"). 
    #alpha=second_ensamble_phylo: transparencia dependiente de la incertidumbre (número casos en los que la celda ha caído dentro del rango filo)
    #add=TRUE: Add to the previous plot
    #legend=FALSE: Not add legend, in that way the legen of the first plot will not be hidden
#plot the legend of precision    
legend("top", legend=c("High precision points", "Low precision points"), lty=c(1,3), cex=1.25, bty="o", lwd=2, horiz=TRUE, bg="white")
#add legend of acnestral
legend("bottom", legend="Phylo-corrected suitability", fill="#00ffe9", bty="o", cex=1.25, bg="white")

#### main title
title("") #this is used because mtext only does not work (see http://stackoverflow.com/questions/12895783/r-language-mtext-not-working-with-image-plot-array)
mtext(bquote(italic('Pinus') ~italic(.(species))), outer = TRUE, cex = 2.5, font=2) #bquote is used to convert to italic the specific epithet (see http://stackoverflow.com/questions/27266398/using-italics-in-the-title-on-an-object-from-a-dataframe)
dev.off()

#save image
save.image("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/phylo_reconstruction/rdata/tests_phylo.RData")