#!/usr/bin/env Rscript

#we use the shebang "#!/usr/bin/env Rscript" to facilitate portability between machines. "env" offers more portability because it searches for the first occurrence of the R interpreter in the environment where we are running. It should thus look for R version installed in the container
    #https://www.r-bloggers.com/2019/11/r-scripts-as-command-line-tools/
    #https://www.baeldung.com/linux/bash-shebang-lines

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

#pass command line arguments
require(optparse)
    #we are going to use the pythonic way thanks to optparse, but you can also use base for this by using "commandArgs" 
        #https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/
option_list = list(
    make_option(opt_str=c("-s", "--species"), type="character", default="halepensis,radiata", help="Species to be analyzed. Strings separated by commas WITHOUT spaces between them [default=%default]", metavar="character"),
    make_option(opt_str=c("-b", "--batch"), type="character", default="batch_1", help="Name of the batch [default=%default]", metavar="character"))
    #opt_str: string indicating the short and long flags for the argument
    #type: string indicating if integer, logical, character...
    #default: default value
    #help: A character string describing the option to be used by 'print_help' in generating a usage message.
        #run the script with the flag --help
    #metavar: A character string that stands in for the option argument when printing help text
        #this changes the type of the argument when printing the help
    #dest: indicates the name of the field of the output list storing the value of the argument. By default, it uses the name of the long flag
        #for example, "batch" will store the batch number
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
    #this takes a "OptionParser" instance as input and generates as output a list with the values of the arguments
    #there are options to deal with positional arguments, etc..
#get a vector with the arguments
species_to_analyze=strsplit(opt$species, split=",")[[1]]
batch_number=opt$batch

#packages
require(raster)

#set the seed for reproducibility
set.seed(4357430)
    #we have a random selection of occurrences within the 50x50km cells

#create some folders
system("mkdir -p ./results/global_test_phylo_current/exsitu_occurrences")
system("mkdir -p ./results/global_test_phylo_current/env_predictors")
system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_no_phylo/", sep=""))
system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_phylo/", sep=""))

#load species names
list_species = read.table("./code/presences/species.txt", sep="\t", header=TRUE)

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
buffer_albicaulis = raster(paste("./results/ocurrences/albicaulis_distribution_buffer", ".asc", sep=""))
#prepare fact value for aggregate
fact_value=res(buffer_albicaulis)[1]/res(environment_var)[1]
    #in raster::aggregate, fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
    #we calculate it dividing the target (coarser) resolution by the original (finer) resolution, so we can know by how much we have to increase the number of cells in order to reach the target (coarser) resolution
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
if(!file.exists("./results/global_test_phylo_current/exsitu_occurrences/pinus_occurrences_fordryad_exoticonly.csv")){
    system(paste(" \\
        unzip \\
            -o \\
                ./datos/phylo/method_validation/doi_10_5061_dryad_1hr1n52__v20181213.zip \\
                pinus_occurrences_fordryad_exoticonly.csv \\
            -d ./results/global_test_phylo_current/exsitu_occurrences/", sep="")) 
}
naturalized_occurrences=read.csv("./results/global_test_phylo_current/exsitu_occurrences/pinus_occurrences_fordryad_exoticonly.csv", header=TRUE)
summary(naturalized_occurrences)
#check
if(nrow(naturalized_occurrences)!=597){
    stop("ERROR! FALSE! WE DO NOT HAVE ALL NATURALIZED OCCURRENCES WE SHOULD HAVE ACCORDING TO THE PERRET'S MANUSCRIPT")
}
if(sum(unique(naturalized_occurrences$species) %in% epithet_species_list)!= length(unique(naturalized_occurrences$species))){
    stop("ERROR! FALSE! WE DO NOT HAVE THE SAME SPECIES NAMES IN PERRET DATA")
}
if(
    class(naturalized_occurrences$species)!="character" | 
    class(naturalized_occurrences$longitude)!="numeric" |
    class(naturalized_occurrences$latitude)!="numeric"){
    stop("ERROR! FALSE! WE HAVE A PROBLEM WITH THE DATA TYPES IN THE OCCURRENCE DATA")
}

#plot minimum number of naturalized occurrences vs number of species retained
#calculate the number of occurrences per species
#unique_species=unique(naturalized_occurrences$species)[1]
n_naturalized=data.frame(unique_species=NA, n_nat_occurrences=NA)
for(unique_species in unique(naturalized_occurrences$species)){

    #calculate the number of occurrences of the selected species
    n_nat_occurrences=nrow(naturalized_occurrences[which(naturalized_occurrences$species==unique_species),])

    #save with the species name
    n_naturalized=rbind.data.frame(n_naturalized, cbind.data.frame(unique_species, n_nat_occurrences))
}
#remove first row with NA
n_naturalized=n_naturalized[which(apply(is.na(n_naturalized), 1, sum)!=ncol(n_naturalized)),]
    #remove row for which the sum of TRUEs of is.na() is equal to the number of columns, i.e., row with NA for all columns
#update the number of rows
row.names(n_naturalized)=1:nrow(n_naturalized)
#reorder base on the number of occurrences
n_naturalized=n_naturalized[order(n_naturalized$n_nat_occurrences, decreasing=FALSE),]
print(n_naturalized)

#plot the threshold of number of occurrences against the remaining number of species
threshold=1:15
    #number of thresholds to test
plot(1, type="n", xlab="", ylab="", xlim=c(0, length(threshold)), ylim=c(0, nrow(n_naturalized)), xaxt="n", yaxt="n")
    #open empty plot
        #https://stackoverflow.com/a/23409042
axis(1, at=seq(0, length(threshold), by=1), las=2)
axis(2, at=seq(0, nrow(n_naturalized), by=1), las=2)
    #add x and y ticks
threshold_results=data.frame(threshold_value=NA, species_remaining=NA)
    #open DF to save results
#threshold_value=threshold[1]
for(threshold_value in threshold){

    #count number of species remaining after applying the threshold
    species_remaining=nrow(n_naturalized[which(n_naturalized$n_nat_occurrences>=threshold_value),])

    #plot
    points(x=threshold_value, y=species_remaining)

    #save results
    threshold_results=rbind.data.frame(threshold_results, cbind.data.frame(threshold_value,species_remaining))
}
dev.off()
    #close the plot
threshold_results=threshold_results[which(apply(is.na(threshold_results), 1, sum)!=ncol(threshold_results)),]
    #remove first row with NAs
row.names(threshold_results)=1:nrow(threshold_results)
    #update row names
print(threshold_results)
    #Great decrease from 1 to 5, but then from 5 to 7 only 2 species lost, while from 8 to 10 only one specie more lost. We will select 10.
    #We have to bear in mind that the final number of occurrences can be lower because we do resampling and also remove occurrences inside the PA buffer, even if they are outside of the natural distribution to avoid using any area considered during modeling.
        #remember we should use new data not seen by the model. An occurrence point within the PA buffer could fall on a pseudo-absence, i.e., a spatial point already considered by the model. It do not think we should consider that are for an independent evaluation, better to go to other continents.
        #also, doing this we avoid any potential risk of considering a real natural population that is not included in the natural range for any reason. We are being here more stringent, selecting occurrences naturalized for sure.
    #Therefore, we could also lose more species. Using a threshold of five, we would likely lose muricata, rigida, torreyana and clausa, having a total 15 species.

#set the threshold of FINAL number of naturalized occurrences we have
nat_threshold=5
    #the boyce index is specially well suited to dataset where absences are not reliable, so I understand we can use it even if we have very few presences.
        #if only 5 presence exits and from the whole globe it fall exactly in a place with high suitability, it is telling you something. it is not so strong than having 50 positive hits, but it is something. Maybe, other occurrences would not fall within high suitability but this one did among many many other pixels that are not suitable
        #The same applies for the opposite. If you only have 50 presences and it falls in low-suitability bin, for now you cannot be confidence about your model.
    #I still want to avoid those cases with just 1 or 2 occurrences, because the number is extremely low and these cases give me strange errors in the code. 5 seems to be a good balance.
    #Once the analyses are done, we will check the final number of occurrences, and take that in consideration when discussing the results.
    #Note that the aspect directly mentioned about boyce that can be problematic, i.e., the existence of suitability bins with low number of pixels, IS CONSIDERED. If a bin has less than 60 pixels (30 is the default threshold for the boyce function), we will stop the analyses (see below).

#NOTE: In general, I have used extract to obtain any information about the points: elevation, environmental data and number of cell. Extract only consider that a point falls inside a cell if its center is inside that cell. It is important consider this. 




###################################################
##### SELECT OCCURRENCES OUTSIDE DISTRIBUTION #####
###################################################

#packages
require(raster)
require(sf)

#species="radiata"
exsitu_occurrences=function(species){

    #start
    print(paste("STARTING exsitu_occurrences FOR ", species), sep="")

    ##momentarily set the flag for the next step as FALSE
    next_step=FALSE
        #it will change to TRUE if the whole script runs

    ##stop if pumila, just in case, because this species has part of its distribution in a corner and we should check it carefully
    if(species=="pumila"){
        stop("ERROR! FALSE! YOU ARE TRYING TO ANALYZE PUMILA, AND WE HAVE NOT CHECKED THIS SCRIPT FOR THAT AS THIS SPECIES HAS PART OF ITS DISTRIBUTION AS AN INDEPENDENT PATCH IN THE MAP (THE LITTLE PEAK OF SIBERIA IN THE LEFT CORNER")
    }

    ##open folders
    system(paste("mkdir -p ./results/global_test_phylo_current/exsitu_occurrences/", species, sep=""))

    ##load the distribution buffer
    #read the raster
    distribution_buffer = raster(paste("./results/ocurrences/", species, "_distribution_buffer.asc", sep=""))

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
    if((nrow(occurrences)-nrow(environment_presences))>2){
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
    n_points_in_out_pa_buffer=cbind.data.frame(nrow(points_and_values_outside_pa_buffer), length(points_inside_pa_buffer), length(points_outside_pa_buffer))
    names(n_points_in_out_pa_buffer)=c("total_points", "points_inside", "points_outside")
    if(n_points_in_out_pa_buffer$points_inside > (n_points_in_out_pa_buffer$total_points*0.20)){
        print(paste("WARNING! WE HAVE A PROBLEM, MORE THAN 20% OF OCCURRENCES OF PERRET FALL WITHIN PA BUFFER FOR SPECIES ", species, sep=""))
        print((n_points_in_out_pa_buffer$points_inside/n_points_in_out_pa_buffer$total_points)*100)
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
        if(length(subset_2)>0 | length(subset_3)>0){
            checks[[2]]=FALSE
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
            #i guess the previous commented lines were used in the original script for taking a look to some examples of cells with high precision points: less than 3 or more than 3 in the cell, respectively.
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
                
                #stop this because we should not have cases with precision weight==1
                stop(paste("ERROR! FALSE! WE HAVE A PROBLEM, OCURRENCES WITH PRECISION WEIGHT==1 FOR SPECIES ", species, sep=""))

                #open png
                png(paste("./results/global_test_phylo_current/exsitu_occurrences/", species, "/", species, "_check_altitudinal_sampling.png", sep=""), width=2200, height=1600, res=300)

                #plot elevation
                if(species=="radiata"){
                    #if radiata, plot the canary islands where it has been introduced and where we can see selection of high precision points across the Teide
                    canariensis_buffer=raster(paste("./results/ocurrences/canariensis_distribution_buffer", ".asc", sep=""))
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

    #add the number of occurrences (outside PA buffer) after resampling in the counting file
    if(exists("ultimate_ocurrences")){
        n_points_in_out_pa_buffer$points_outside_after_resampling=nrow(ultimate_ocurrences)
    } else {
        n_points_in_out_pa_buffer$points_outside_after_resampling=0
    }

    #change the status of the flag for the next step to TRUE
    next_step=TRUE

    #end
    print(paste("ENDING exsitu_occurrences FOR ", species), sep="")

    #return the flag for the next step and the number of points before resampling
    return(list(next_step, n_points_in_out_pa_buffer))
}

#run it for one species
#exsitu_occurrences("radiata")




##############################################
##### PREDICT AND EVALUATE WITHOUT PHYLO #####
##############################################

#packages
require(randomForest)
require(gam)
require(modEvA)
    #for boyce index. I am using this instead of ecospat.boyce because the latter gives me errors when loading the predicted raster
require(terra)
    #to convert RasterLayers to SpatRasters, which is the required input format of the boyce index function

#species="radiata"
predict_eval_no_phylo = function(species, status_previous_step){

    #starting
    print(paste("STARTING predict_eval_no_phylo FOR ", species), sep="")

    #open folder
    system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/", sep=""))

    #momentarily set the flag for the next step as NO
    next_step=FALSE
        #it will change to TRUE if the whole script is run and everything is OK

    #check if we have output from the previous step
    presence_file_exist=file.exists(paste("./results/global_test_phylo_current/exsitu_occurrences/", species, "/", species, "_final_presences.tsv", sep=""))

    #check if the number of presences is higher than the threshold of naturalized occurrences
    #do it only if we have a presence file indeed
    if(presence_file_exist){
        check_n_presences=nrow(read.table(paste("./results/global_test_phylo_current/exsitu_occurrences/", species, "/", species, "_final_presences.tsv", sep=""), sep="\t", header=TRUE))>=nat_threshold
    } else{
        #if no file is present, just set the check as FALSE
        check_n_presences=FALSE
    }

    #do stuff if we have the presence file and more than 10 naturalized occurrences. Also we need green light from the previous step
    if(presence_file_exist & check_n_presences & status_previous_step){

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
   
        #save the environmental stack
        system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/env_stack", sep=""))
        writeRaster(variables_stack_masked, filename=paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/env_stack/", species, "_env_stack", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        print(system(paste("
            cd ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/; 
            tar -zcvf ", species, "_env_stack.tar.gz ./env_stack; 
            rm -rf ./env_stack", sep=""), intern=TRUE))
            #we force the output to be sent to the output file of the species using intern=TRUE that generates a vector with the string output and the force printing with print
            #print should not be necessary, but this is only the way I have found to show the system's output while parallelizing


        ###predict suitability and evaluate
        ##load models
        #open folder to save the decompressed files
        system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/decompressed_models", sep=""))

        #decompress into the folder of the species
        print(system(paste("unzip -o ./results/models/models_", species, ".zip -d ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/decompressed_models/", sep=""), intern=TRUE))
            #the models have been obtained from Rafa PRO: "/Users/dsalazar/nicho_pinus/results/final_analyses/models". This is the path indicated in "/home/dftortosa/diego_docs/science/phd/nicho_pinus/code/models/fit_predict_current_suit/fit_predict_current_suit_v5.R"

        #load the RDA files of the models
        load(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/decompressed_models/", species, "_glm_model.rda", sep=""))
        load(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/decompressed_models/", species, "_gam_model.rda", sep=""))
        load(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/decompressed_models/", species, "_rf_model.rda", sep=""))


        ##extract thresholds
        #These thresholds were previously calculated using the evaluation data of each partition hence selecting the threshold maximizing TSS for each data partition.
        #decompress thresholds
        system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/thresholds/", sep=""))
        print(system(paste("unzip -o ./results/threshold/threshold_", species, ".zip -d ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/thresholds/", sep=""), intern=TRUE))
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
        if(length(glm_resample)!=12 | length(gam_resample)!=12 | length(rf_resample)!=12 | length(glm_threshold)!=12 | length(gam_threshold)!=12){
            stop("ERROR! FALSE! WE DO NOT HAVE 12 PARITIONS OF THE DATA, I.E., WE DO NOT HAVE 12 MODELS OR 12 THRESHOLDS FOR EACH ALGORITHM TYPE")
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

        #set the check about min sample size in the bins of suitability as TRUE. It will be changed to FALSE if we have suitability bins with less than 30 pixels (see below)
        check_min_bin_sample_size=TRUE

        #run a loop using the models trained on the 12 data partitions
        #k=1
        for(k in 1:length(glm_resample)){

            ##predict obtaining a continuous probability of presence for each cell outside of the PA_buffer, which will be the input for the Boyce index
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
                #calculate also predictions using the binary results of RF so we can use it later when ensembling predictions. We get binary predictions using response for this model because we used classification trees (see above).


            ##calculate the boyce index
            #run the function removing or not duplicates (see below)
            jpeg(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/plots/", species, "_part_", k, "_boyce_index_plot.jpeg", sep=""), width=960, height=960, pointsize=24)
            par(mfcol=c(3,2))
            glm_boyce=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(glm_predict[[k]]), n.bins=NA, bin.width="default", res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=FALSE, na.rm=TRUE, plot=TRUE, main=paste("GLM", sep=""))
            mtext(paste(species, " - part ", k, sep=""), side=3, line=-1.9, outer=TRUE, cex=1, font=2)
            gam_boyce=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(gam_predict[[k]]), n.bins=NA, bin.width="default", res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=FALSE, na.rm=TRUE, plot=TRUE, main=paste("GAM", sep=""))
            rf_boyce=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(rf_predict[[k]]), n.bins=NA, bin.width="default", res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=FALSE, na.rm=TRUE, plot=TRUE, main=paste("RF", sep=""))
            glm_boyce_no_dup=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(glm_predict[[k]]), n.bins=NA, bin.width="default", res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=TRUE, na.rm=TRUE, plot=TRUE, main=paste("GLM - no duplicates", sep=""))
            gam_boyce_no_dup=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(gam_predict[[k]]), n.bins=NA, bin.width="default", res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=TRUE, na.rm=TRUE, plot=TRUE, main=paste("GAM - no duplicates", sep=""))
            rf_boyce_no_dup=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(rf_predict[[k]]), n.bins=NA, bin.width="default", res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=TRUE, na.rm=TRUE, plot=TRUE, main=paste("RF - no duplicates", sep=""))
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
                            #0.5/0.5=1 is the P/E ratio
                        #now imagine that we remove the duplicates, so only 1 of the two presences remains
                            #1/3=0.33 is the proportion of presences in the bin
                            #4/8=0.5 is the proportion of cells/pixels in the bin
                            #0.33/0.5=0.66 is the P/E ratio
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
                        #The information I got from copilot assumes that the duplicated points are indeed the same, just duplicates and non-independent, but this is not the case. We have resampled the occurrences to have only 3 occurrences per 50x50km cell, and we have checked that we do not have more than 3 occurrences in any of these cells. Therefore, we have independent points.
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
                            #In their opinion, ten classes (0.1) is a good optimum. So we are going to use the default, i.e., 1/10 or 10 classes.
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
                        #these bins with the same value will decrease the correlation between P/E ratio and suitability if they all have the same value! Remember that each new window has a higher average of suitability (e.g., we move from 0.70-0.80 to 0.71-81 of suitability values), so it should have a higher proportion of presences respect to the expected ratio strengthening the correlation. If we do not have the expected increase in presences, the correlation will be weaken.
                        #If we remove these bins that are contiguous and have the same P/E values:
                            #1) we will lose resolution, in the sense we will not know if the model is able or not to distinguish between different levels of suitability within this range with repeated values. If you remove 0.70-0.8 and 0.71-0.81, you do not know whether the model discriminate well between these classes at that fine resolution. In other words, pooling data reduces the resolution of the model’s predictions. The resolution refers to the model’s ability to distinguish between many different classes (levels) of suitability. Because of this, ecospat docs says that setting this as TRUE "lower the assessment of the evaluation of the model resolution".
                            #2) We will focuses more on the discriminative aspect of the predictions. We are just evaluating whether more broad suitability classes (levels) differ in their proportion of presences. In other words, we focus the evaluation on the ability of the model to discriminate areas of presence from areas of absence or lower presence probability, taking less care of discriminating between more classes (0.70-0.80 vs 0.71-81 vs 0.72-0.82...). 
                                #And this is what we want, just check how well the model differentiate between areas with low and high suitability, predicting naturalized occurrences in these high-suitable areas.
                        #Let's talk about the case of GAM for radiata in partition 1.
                            #Without removing duplicated classes, we get Boyce=0.0002, while removing the duplicated classes gives Boyce=0.5
                            #We have all presences in just two suitability levels, around 0 and around 1. This is in line with the fact that most of the suitability values are around 0 and 1 with the fact that the threshold maximizing TSS in that partition is 0.78, so most of the suitable areas are only those veery suitable, there is no a gradient of suitability values. Despite being continuous, it is very similar to binary.
                            #The first suitability bin (0-0.1) has a proportion of presences of presences respect to the total of presences of 0.64, while the proportion of pixels is 0.89. In other words, the proportion of presences is lower than the random expectation. The proportion of presences is a 29% less than the random expectation (P/E ratio=0.71).
                            #The last bin (0.891-0.991) has 0.015 proportion of presences while the proportion of pixels is 0.016. Again, the proportion of presences is lower than the random expectation BUT we are not much closer. Now the proportion of presences is lower only in a 5% (P/E ratio=0.002).
                            #Despite being the proportion of presences lower than the random expectation in both bins, there is a clear increase in the proportion of presences in the bin of high suitability.
                            #if you consider the ability of the model to discriminate between multiple suitability levels, then the model is doing pretty bad as you have maaaany suitability levels with the same P/E ratio, so there is no a continuous increase in the proportion of presences as predicted suitability increases. In line with this, Boyce=0.002.
                                #As explained in Harzel et al 2006, this information can be relevant for managment. If you plot the curves across partitions, you can get a median and confidence interval, so you can see in which levels of suitability you model works consistenly well across partitions and consider that when making management decisions. But we are not interested in this right now, see next line.
                            #if you remove all the suitability levels in the middle that have the same P/E ratio and hence, you only consider the ability of the model to discriminate between low and high suitability (0-0.1 vs 0.89-0.991), then the model is doing well. In, areas with a high suitability, the proportion of presences tend to higher with respect to the random expectation under an uniform distribution of presences. The proportion of presences is still lower than the proportion of pixels in the last bin, but that presence proportion gets much closer to the random expection in the high suitability bin. Consequently, we get Boyce=0.5.
                                #Indeed, if we check the plot with the presences on top of the predictions of gam for partition 1, we can see presences in areas with 0 suitbility, but also some of them in areas with high suitability.
                                #We are going to consider the median of Boyce across partitions, anyways, so if, like in radiata, some partitions are very high but in generla we have this situaton of extreme suitabilities, others will have very low boyce index, and the median will be reduded.
                        #The default of ecospat is TRUE for this argument, and this was the option used by Nick, but it is not the default for modEvA.
                        #Note that the duplicated classes remain in "$bins" of the Boyce results, you have to remove them if you want to replicate the Byce value we used here.
                    #na.rm: 
                        #Logical value indicating if missing values should be removed from computations. The default is TRUE.
                        #TRUE, because we want to avoid cells without suitability values. There is nothing to do there.
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
  

            ##check ptsrast2obspred, which is internally used by boyce, and the number of points within 10x10km cells for each algorithm
            #it should be the same in all models, but we are checking it just in case
            #algorithm="glm"
            for(algorithm in c("glm", "gam", "rf")){

                #select the raster with the predictions of the selected model
                selected_predict_raster=eval(parse(text=paste(algorithm, "_predict[[k]]", sep="")))

                #run the function to calculate DF with presence/absence and predicted probability for each model
                obspred_dup=ptsrast2obspred(rst=terra::rast(selected_predict_raster), pts=presences[,c("longitude", "latitude")], rm.dup=FALSE, na.rm=TRUE, verbosity=2)
                obspred_no_dup=ptsrast2obspred(rst=terra::rast(selected_predict_raster), pts=presences[,c("longitude", "latitude")], rm.dup=TRUE, na.rm=TRUE, verbosity=2)

                #the number of presences/absences after removing dups (second point in the same cell) should be the same than the total number of cell in the raster of predictions as we have now one value per cell, no more.
                if(nrow(obspred_no_dup) != length(which(!is.na(getValues(selected_predict_raster))))){
                    stop(paste("ERROR! FALSE! WE ARE NOT CORRECTLY CALCULATING THE PRESENCES/ABSENCES FOR EVALUATION IN SPECIES ", species, " AND MODEL ", algorithm, sep=""))
                }
                    #the additional presences in the same cell are included as additional rows in the data.frame
                    #these will be also considered when counting the number of presences in each suitability bin and when counting the total number of presences, but this is ok, because these are legit and independent observations following our definition (see above).

                #check if we have too many points within the same 10x10 cells
                if((nrow(obspred_dup)-nrow(obspred_no_dup))>(nrow(presences)*0.35)){

                    #print the percentage of occurrences in the same 10x10km cells
                    print(((nrow(obspred_dup)-nrow(obspred_no_dup))/nrow(presences))*100)
                    
                    #stop
                    stop(paste("ERROR! FALSE! WE HAVE A PROBLEM WITH THE NUMBER OF PRESENCES WITIN THE SAME 10x10km CELL. MORE THAN 35% OF THE PRESENCES ARE PRESENT IN THE SAME CELL WITH OTHER PRESENCES FOR SPECIES ", species, " AND MODEL ", algorithm, sep=""))
                }
                    #the difference in points between obspred_dup and obspred_no_dup is the number of presences being in 10x10km cells with other presences, as the latter was created by removing these presences, in contrast with the former.
            }


            ##check we do not have a low number of pixels in any bin
            #In the docs of "modEvA::Boyce", they say the following: "In bins with overly small sample sizes, the comparison between median prediction and random expectation may not be meaningful, although these bins will equally contribute to the overall Boyce index. When there are bins with less than 30 values, a warning is emitted and their points are plotted in red, but mind that 30 is a largely arbitrary number. See the $bins$bin.N section of the console output, and use the 'bin.width' argument to enlarge the bins."
            #For species with a small natural range, it would be possible that no pixel has high suitability values, making a low sample size for high suitability bins. But this should not be very problematic in our case, because we are predicting across the whole world, increasing the probability to find enough pixels with low and high suitability. For example, clausa has a small range but its minimum pixel number per bin is 70. We are going to check that just in case
            min_bin_n_glm=min(glm_boyce[[1]]$bin.N) 
            min_bin_n_gam=min(gam_boyce[[1]]$bin.N) 
            min_bin_n_rf=min(rf_boyce[[1]]$bin.N) 
                #calculate the min number of points across bins in each model. You can have different sample size between models because these are bins based on suitability, and the suitability is differently predicted by each model
            #we are going stop the analyses if the minimum number of pixels per bin is 60, which is the double of the default value set by "modEvA::Boyce" to print a warning (i.e., 30)
            if((min_bin_n_glm<60) | (min_bin_n_gam<60) | (min_bin_n_rf<60)){
                print(paste("ERROR! FALSE! WE HAVE AT LEAST 1 NON-PHYLO SUITABILITY BIN WITH LESS THAN 60 DATA.POINTS FOR SPECIES ", species, ". THIS CAN MAKE THAT THE COMPARISON BETWEEN MEDIAN PREDICTION AND RANDOM EXPECTATION MAY NOT BE MEANINGFUL", sep=""))
                check_min_bin_sample_size=FALSE
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
            #using this approach makes all numeric variables to be converted to character. Look at the next major step for approach to change if needed. Not required right now, as we are directly saving the table to TSV, and the numbers are saved as num and not characters


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
        n_layers_binary_predictions=nlayers(binary_predictions)

        #calculate the percentage of models for which a pixel is suitable
        ensamble_predictions_bin = calc(binary_predictions, function(x) (sum(x)*100)/n_layers_binary_predictions)
            #ensamble_predictions_bin2 = calc(binary_predictions, function(x) (sum(x)*100)/nlayers(binary_predictions))
            #identical(getValues(ensamble_predictions_bin), getValues(ensamble_predictions_bin2))


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
        print(system(paste("
            cd ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/; 
            tar -zcvf ", species, "_prediction_rasters.tar.gz ./prediction_rasters; 
            rm -rf ./prediction_rasters", sep=""), intern=TRUE))
            #https://unix.stackexchange.com/a/93158
    
        #if we have passed the check about min sample size in suitability bins
        if(check_min_bin_sample_size){

            #green light for the next step
            next_step=TRUE
        }
    } else {
        #if we cannot do stuff, remove the folder of the species
        system(paste("rm -rf ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/", sep=""))
    }

    #ending
    print(paste("ENDING predict_eval_no_phylo FOR ", species), sep="")

    #return flag for the next step
    return(next_step)
}

#run it for one species
#predict_eval_no_phylo("radiata")




###########################################
##### PREDICT AND EVALUATE WITH PHYLO #####
###########################################

#packages
require(modEvA)
require(terra)
require(raster)
require(gtools)

#species="radiata"
predict_eval_phylo = function(species, status_previous_step){

    #staring
    print(paste("STARTING predict_eval_phylo FOR ", species), sep="")

    #open folder
    system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_phylo/", species, "/", sep=""))

    #check the outputs of the previous step exists
    pred_raster_check=file.exists(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/", species, "_prediction_rasters.tar.gz", sep=""))
    boyce_table_check=file.exists(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/", species, "_boyce_table.tsv.gz", sep=""))

    #if the previous outputs exists and we have green light from the previous step, do operations
    if(pred_raster_check & boyce_table_check & status_previous_step){
        
        ##perform the phylogenetic correction
        #load bio4 and bio17
        bio4_raw = raster("./datos/finals/bio4.asc")
        bio17_raw = raster("./datos/finals/bio17.asc")

        #load ancestral reconstruction
        final_anc_ou_bio4 = read.table("./results/phylo_reconstruction/final_recons/final_anc_ou_bio4.csv", sep=",", header=TRUE) # nolint
        final_anc_ou_bio17 = read.table( "./results/phylo_reconstruction/final_recons/final_anc_ou_bio17.csv", sep=",", header=TRUE) # nolint
        final_anc_bm_bio4 = read.table("./results/phylo_reconstruction/final_recons/final_anc_bm_bio4.csv", sep=",", header=TRUE) # nolint
        final_anc_bm_bio17 = read.table( "./results/phylo_reconstruction/final_recons/final_anc_bm_bio17.csv", sep=",", header=TRUE) # nolint
            #we set nolint for these objects along with the list of models and bio4-bio17 rasters to avoid the warning indicating that the variables are not being used. These objects are used but with get() and eval(parse()), so the linter does not detect them. The ou models will be not used, only BM.

        #list of evolution models: ONLY BM
        #list_models_bio4 = c("final_anc_ou_bio4", "final_anc_bm_bio4")
        #list_models_bio17 = c("final_anc_ou_bio17", "final_anc_bm_bio17")
        list_models_bio4 = c("final_anc_bm_bio4") # nolint
        list_models_bio17 = c("final_anc_bm_bio17") # nolint

        #remove the areas inside the PA buffer
        #we have selected naturalized occurrences outside the PA buffer, see occurrences script for further details
        raster_pa_buffer = raster(paste("./results/pseudo_absences/", species, "_PA_buffer.asc", sep=""))

        #convert to polygon 
        polygon_raster_pa_buffer = rasterToPolygons(raster_pa_buffer, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon

        #IMPORTANT STEP: do an inverse mask of the environmental variable using the PA buffer as mask. This will give us all areas outside of the PA buffer.
        bio4=mask(bio4_raw, polygon_raster_pa_buffer, inverse=TRUE) # nolint
        bio17=mask(bio17_raw, polygon_raster_pa_buffer, inverse=TRUE) # nolint
            #we do not need the area inside the PA buffer, as we have removed all occurrences occurring there. We do not want to use regions considered during the modeling.
            #also, the boyce index (evaluation metric for presence-only data) uses the whole raster of predictions, considering as absences everything that does not have a presence. It would consider as absences the PA buffer if we do not remove this area because we have removed the occurrences inside that buffer in these exsitu analyses.

        #load ensamble of binary suitability
        print(system(paste("tar --directory ./results/global_test_phylo_current/predict_eval_phylo/", species, " -zxvf ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/", species, "_prediction_rasters.tar.gz ./prediction_rasters/", species, "_ensemble.grd ./prediction_rasters/", species, "_ensemble.gri", sep=""), intern=TRUE))
            #--directory has to go first
        ensamble_suitability = raster(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters/", species, "_ensemble.grd", sep=""))

        #open stacks to save rasters
        phylo_rasters_subset = stack()
        phylo_rasters_proportion_subset = stack()
        phylo_rasters_no_subset = stack()
        phylo_rasters_proportion_no_subset = stack()

        #vector with env variable names
        environmental_variables=c("bio4", "bio17")

        #check what cells fall within the phylogenetic range for bio4 and bio17 and save the result as a stack (with and without proportion)
        #env_var="bio4"
        for(env_var in environmental_variables){

            #get the selected variable
            selected_env_var=get(env_var)

            #create a empty raster with the same extent and resolution of the raster layer of the [s] IPCC scenario
            raster_subsetted = raster(extent(selected_env_var), resolution=res(selected_env_var))
                #we can just use the environmental variable because this comes from "finals", which is the folder of the variables used to the create the input stack used for global predictions (see above).

            #add to the empty raster those cells of the [s] IPCC raster in which the habitat suitability is higher than 25 and lower than 75
            raster_subsetted[which(getValues(ensamble_suitability) > 25 & getValues(ensamble_suitability) < 75)] <- selected_env_var[which(getValues(ensamble_suitability) > 25 & getValues(ensamble_suitability) < 75)]
                #it is ok to select only regions 25-75 considering all models (ensemble) and then apply the phylo-correction in these regions within the prediction of a single model (e.g., glm). We need regions with uncertainty and you can detect uncertainty when combining different partitions and models.
                #In other words, we have areas of disagreement between partitions, so we use the phylogenetic correction to solve the conflict.
                #Regions within the phylogenetic correction are considered suitable independently of the algorithm. This is useful because we can see what happens with RF and GLM when applying the phylogenetic correction separately. Maybe the correction is not useful for some algorithms, but it could be for others. Maybe, for GLM there is an area that is not suitable but the correction makes it suitable, improving a lot boyce, while that area is already suitable in RF, or RF has very low correlation in other areas so this change does not have impact.
                #As a final step after running all the analyses, we are going to model non-phylo boyce MINUS phylo boyce as a function of algorithm type, partition and species, so we will be able to see the impact of each of these factors on performance of the models. For example, an interaction of species X algorithm would mean that the phylogenetic correction has an impact for some combinations of species and algorithms but not for others. Maybe it is more difficult to see the impact if only look for main effects, without separating between algorithm types.
 
            #create also a raster with all values of the variable, so we can apply the correction across the whole world
            raster_no_subsetted = selected_env_var

            #get the list of models
            list_selected_models=get(paste("list_models_", env_var, sep=""))

            #for each evolution model
            #m=1
            for(m in 1:length(list_selected_models)){ # nolint

                #select the [m] evolution model 
                selected_model = list_selected_models[m]

                #extract data of [m] model
                model = get(selected_model)

                #select the row of the corresponding species
                model = model[which(model$species == paste("Pinus_", species, sep="")),]

                #extract all cell values from the raster with climatic data of the selected varaible only in those areas with uncertainty (raster_subsetted)
                cell_values_subset = getValues(raster_subsetted)
                #extract ID of those cells without NA
                cells_withot_NA_subset = which(!is.na(cell_values_subset))
                #check that the index of the cells without NA is correct and the same used as input for is.between()
                if(sum(getValues(raster_subsetted)[cells_withot_NA_subset]!=na.omit(getValues(raster_subsetted)))!=0){
                    stop("ERROR! FALSE! WE HAVE A PROBLEM SELECTING CELLS INSIDE THE PHYLO RANGE")
                }

                #do the same but considering all cells independenlty of suitability
                cell_values_no_subset = getValues(raster_no_subsetted)
                #extract ID of those cells without NA
                cells_withot_NA_no_subset = which(!is.na(cell_values_no_subset))
                #check that the index of the cells without NA is correct and the same used as input for is.between()
                if(sum(getValues(raster_no_subsetted)[cells_withot_NA_no_subset]!=na.omit(getValues(raster_no_subsetted)))!=0){
                    stop("ERROR! FALSE! WE HAVE A PROBLEM SELECTING CELLS INSIDE THE PHYLO RANGE")
                }

                #extract, from all cells withput NA, those whose value is inside the phylogenetic range (including the current value but not including the ancestral). For that we used is.between function, created by me.
                #set function to check if cell values is between ancestral and current value
                #cell_value=1000; ancestral=900; current=1100
                #cell_value=889; ancestral=900; current=1100
                #cell_value=1101; ancestral=900; current=1100
                #cell_value=900; ancestral=900; current=1100
                #cell_value=1100; ancestral=900; current=1100
                is.between <- function(cell_value, ancestral, current) {
                    
                    #select the cell_value
                    selected_value = cell_value

                    #test if the [v] cell value is between ancestral an current values
                    test = (selected_value - ancestral)  *  (current - selected_value) #code taken from "https://stat.ethz.ch/pipermail/r-help/2008-August/170749.html". The order is irrelevant, ancestral can be higher or lower than current value. Idem for the sign of numbers, it works with only negative, only positive and negative-positive numbers.  

                    #if test is not zero 
                    if(!test == 0){

                        #test if test is lower or higher than zero to know is the [v] cell value is between current and ancestral values. then save
                        result = test > 0
                    } else { #if not, then [v] cell value is equal to the current or ancestral value, but we only want TRUE if the value is equal to the current value. 

                        #If the [v] cell value is equal to the current value
                        if(current == selected_value){

                            #result is TRUE
                            result = TRUE
                        } else { #if not

                            #if the [v] cell value is equal to ancestral value
                            if(ancestral == selected_value){

                                #result is FALSE
                                result = FALSE
                            } else {

                                #result is NA, problem
                                result = NA
                            }
                        }
                    }

                    #return results
                    return(result)
                } 
                    #Is very important to add ancestral first, and second current value, becuase TRUE will be returned if the cell_values is equal to "current" (second argument), but FALSE if it is equal to ancestral (first argument)
                cell_inside_phylo_range_subset=which(apply(X=array(na.omit(getValues(raster_subsetted))), MARGIN=1, FUN=is.between, ancestral=model$ace, current=model$current_value))
                cell_inside_phylo_range_no_subset=which(apply(X=array(na.omit(getValues(raster_no_subsetted))), MARGIN=1, FUN=is.between, ancestral=model$ace, current=model$current_value))
                    #apply the function over a vector with the value of selected cells
                        #you have to convert to array the vector in order to use it as input for apply. I have checked this does not change the values.
                        #a vector converted to array has only 1 dimension, so margin has to be 1.
                        #you could parallelize here if needed

                #from ID of cells without NA, select the ID of those whose vale is inside of the phylo range
                final_cells_subset = cells_withot_NA_subset[cell_inside_phylo_range_subset]
                final_cells_no_subset = cells_withot_NA_no_subset[cell_inside_phylo_range_no_subset]
                    #we can use the list of cells within the phylogenetic range directly on "cells_without NA" because we checked what cells where inside the range using na.omit(raster_subsetted) and na.omit(raster_no_subsetted), so the ID of cells is the same in both cases.

                #create a empty raster with the same extent and resolution than the selected env variable
                final_raster_subset = raster(extent(selected_env_var), resolution=res(selected_env_var))
                final_raster_proportion_subset = raster(extent(selected_env_var), resolution=res(selected_env_var))
                final_raster_no_subset = raster(extent(selected_env_var), resolution=res(selected_env_var))
                final_raster_proportion_no_subset = raster(extent(selected_env_var), resolution=res(selected_env_var))

                #fill the raster with zeros
                final_raster_subset[] <- 0
                final_raster_proportion_subset[] <- 0
                final_raster_no_subset[] <- 0
                final_raster_proportion_no_subset[] <- 0

                #add to these final cells a value of suitability without proportion
                final_raster_subset[final_cells_subset] <- 1
                final_raster_no_subset[final_cells_no_subset] <- 1

                #add value of suitability with proportion
                #set function to covert the suitibalityi phylo correct to a proportion from 0 to 1
                #x=1000; ancestral_value=900; current_value=1100
                #x=1100; ancestral_value=900; current_value=1100
                #x=900; ancestral_value=900; current_value=1100
                phylo_proportion = function(x, ancestral_value, current_value){

                    #calculate the maximum distance to the ancestal value (i.e the current value)
                    range_length = abs(current_value - ancestral_value)

                    #calculate between the cell value and the ancestral value
                    distance_to_ancestral = abs(x - ancestral_value)

                    #if range_length is the 1, distance_to_ancestral will be x; so x = (distance_to_ancestral*1)/range_length 
                    proportion=distance_to_ancestral/range_length
                    return(proportion)
                } 
                    #Like in the latter function, the order is key. The proportion will have 1 as value is close to the second argument (current value).
                final_raster_proportion_subset[final_cells_subset] = apply(X=array(raster_subsetted[final_cells_subset]), MARGIN=1, FUN=phylo_proportion, ancestral_value=model$ace, current_value=model$current_value)
                final_raster_proportion_no_subset[final_cells_no_subset] = apply(X=array(raster_no_subsetted[final_cells_no_subset]), MARGIN=1, FUN=phylo_proportion, ancestral_value=model$ace, current_value=model$current_value)
                    #here we do not need to apply na.omit because we are already have the ID of cells inside the phylo range considering the whole raster (NA and non-NA)
                        #for is.between(), we tested whether each cell was inside the range considering only non-NA cells
                        #then, we extract the ID of these cells using a vector with the global ID (i.e., considering the whole raster) of non-NA cells.
                        #this number now can be used to select cells in the whole raster
                        #using it as input here will let us to select the cells inside the range in the whole raster, and to convert them to proportion.

                #add the name of the raster
                layer_name=paste(strsplit(selected_model, split="_")[[1]][3], "_", strsplit(selected_model, split="_")[[1]][4], sep="")
                names(final_raster_subset) <- paste(layer_name, "_subset", sep="")
                names(final_raster_proportion_subset) <- paste(layer_name, "_proportion_subset", sep="")      
                names(final_raster_no_subset) <- paste(layer_name, "_no_subset", sep="")
                names(final_raster_proportion_no_subset) <- paste(layer_name, "_proportion_no_subset", sep="")

                #save the raster into a stack
                phylo_rasters_subset = stack(phylo_rasters_subset, final_raster_subset)
                phylo_rasters_proportion_subset = stack(phylo_rasters_proportion_subset, final_raster_proportion_subset)
                phylo_rasters_no_subset = stack(phylo_rasters_no_subset, final_raster_no_subset)
                phylo_rasters_proportion_no_subset = stack(phylo_rasters_proportion_no_subset, final_raster_proportion_no_subset)
            }
        }

        #check we have correct number of layers, i.e., 2 because we have 2 environmental variables.
        if(
            (nlayers(phylo_rasters_subset)!=2) | #nolint
            (nlayers(phylo_rasters_proportion_subset)!=2) | #nolint
            (nlayers(phylo_rasters_no_subset)!=2) | #nolint
            (nlayers(phylo_rasters_proportion_no_subset)!=2)){
            stop(paste("ERROR! FALSE! WE HAVE NOT ANALYZED ALL VARIABLES-MODELS FOR ", species, sep=""))
        }
            #the warning here is caused because we are using "|" insted of "||"
                #c(TRUE, FALSE, TRUE) | c(FALSE, TRUE, FALSE) gives c(TRUE, TRUE, TRUE), because it makes a element-wise comparison.
                #c(TRUE, FALSE, TRUE) || c(FALSE, TRUE, FALSE) gives TRUE, because it only compares the first element of each vector, and returns TRUE if at least one of them is TRUE.
                #in general, we should use the first approach because "if" requires 1 boolean value. We are going to leave it as it is because we are sure that we are only doing one comparison, as nlayers returns only one value.

        #check we have the correct names in each phylo stack
        if(
            (FALSE %in% grepl(pattern="_subset", x=names(phylo_rasters_subset), ignore.case=FALSE, fixed=TRUE)) | #nolint
            (FALSE %in% grepl(pattern="_proportion_subset", x=names(phylo_rasters_proportion_subset), ignore.case=FALSE, fixed=TRUE)) | #nolint
            (FALSE %in% grepl(pattern="_no_subset", x=names(phylo_rasters_no_subset), ignore.case=FALSE, fixed=TRUE)) | #nolint
            (FALSE %in% grepl(pattern="_proportion_no_subset", x=names(phylo_rasters_proportion_no_subset), ignore.case=FALSE, fixed=TRUE))){
            stop(paste("ERROR! FALSE! WE HAVE NOT ANALYZED ALL VARIABLES-MODELS FOR ", species, sep=""))
        }
            #grepl
                #pattern: string containing regular expression or just a string if "fixed=TRUE"
                #fixed: if TRUE, the pattern is a string to matched, not regular expression.
                    #we want to look for a string, not regular expression
                #x: a character vector where we look for matches
                #ignore.case: if "TRUE", case is ignored. Default is FALSE.
            #same warning as previous lines, the script ensures we only have 1 boolean value, so we are going to leave it as it is.

        #combine the phylogenetic correction of all env variables into one single ensemble (excluding and not excluding the areas inside the phylogenetic range of one variable but outside of the other phylo-ranges)
        #open stack to save the ensembles (i.e., considering all env variables) but across different phylo-approaches (e.g., with and without proportion...)
        ensamble_phylo_no_inter=list()
        ensamble_phylo_inter=list()
        #create a vector with the names of the input phylo stacks
        stack_name=c("phylo_rasters_subset", "phylo_rasters_proportion_subset", "phylo_rasters_no_subset", "phylo_rasters_proportion_no_subset")
        #for each phylo-stack
        #k=1
        for(k in 1:length(stack_name)){ #nolint

            #extract the name of the selected stack
            selected_stack_name=stack_name[k]

            #get the actuals stack
            selected_stack=get(selected_stack_name)
            
            #check the name of the layers matches the name of the stack we are working with
            if(!identical(names(selected_stack), paste("bm_", environmental_variables, "_", strsplit(selected_stack_name, split="phylo_rasters_")[[1]][2], sep=""))){
                stop("ERROR! FALSE! WE HAVE A PROBLEM DOING THE STACK OF PHYLO CORRECTIONS, WE ARE NOT SELECTING THE CORRECT PHYLO APPROACH")
            }

            #count the number of layers in the stack
            n_layers_stack=nlayers(selected_stack)
            if(n_layers_stack!=2){
                stop("ERROR! FALSE! WE HAVE MORE ENVIRONMENTAL PHYLO RASTERS THAT WE SHOULD")
            }

            #calculate the average phylo-suitability per cell, i.e., sum the result of the phylo-correction for a cell across environmental variables (phylo correction of all variables) and then divided by the total of layers
            ensamble_phylo = calc(selected_stack, function(x) sum(x)/n_layers_stack)
                #it is much faster if we previously calculate the number of layers. I guess that it calculates the number of layers for each cell otherwise, but this is not necessary, because the number of layers is constant within the stack. A stack has rasters with the same ext, all cells have the same number of layers.
                #ensamble_phylo2 = calc(selected_stack, function(x) (sum(x))/nlayers(selected_stack))
                #identical(getValues(ensamble_phylo), getValues(ensamble_phylo2))
                #TRUE
                #Note about abrupt changes in phylo ensemble:
                    #I have detected in strobus abrupt changes in the ensemble of phylo-corrected suitability proportion so subset.
                    #This is NOT caused by any kind of abrupt change in bio17 (humidity of the driest quarter), but because bio4, which temperature seasonality, a variables without any abrupt change whatsoever. Therefore, this is not a problem of the input climatic data.
                    #The reason behind the abrupt change is in the approach of calculating phylo-proportion only in regions within the phylogenetic range. The first calculation of the phylogenetic range (is.between) is binomial, i.e., the cell is whithin the phylogenetic yes or no. Therefore you are going to have cell in and out of the range together without any transition. The propotion within the range will be continuous, but once you reach the first boundary, the phylo-proportion suitability will go to zero. This is what we see at the longitude around Moscow.
                    #This should not be a problem for model evaluation
                        #TSS is calculated in the continuous predictions of GLM and GAM looking for a threshold that maximize TSS. Therefore, we are binarizying and then evaluating, looking for the threshold best separating presences and absences.
                        #We, indeed, have abrupt changes in the ensembles of continuous predictions without phylo like in strobus. In central europe, you can see abrupt reductions of suitability from green (above 0.5) to gray (0).
                        #In the case of boyce, we could have some occurrences that whether they fall at one side or other of the boundary (inside or outside of the phylogenetic correction), will end up in high or low-suitability bins, making abrupt changes in the values.
                        #This could explain, as we will see in the results, the great increase in suitability for phylo non subset in strobus, but this is not an artifact.
                        #As done for evaluating the original models, we have predictions that tend to be abrupt, and can capture or not occurrences. If we just extend these suitable regions without capturing any presences, boyce index will be low, so we are not artifically inflating this evaluation metric. And we are seeing positive impact in several models, not just one case.... We could have a problem if the positive impact comes only for one model/species, because maybe it was just chances and then amplified by the abrupt change, but this is happening for other models/species.
                        #And again, we have also abrupt changes in our continuous predictions without phylo, and boyce is doing ok. We are not using boyce in a binary raster, just a continuous raster with a little bit less of variation in some regions.

            #add the name of the layer indicating that we are NOT excluding areas that do not fall within all phylogenetic ranges and save
            names(ensamble_phylo) = paste(selected_stack_name, "_no_inter", sep="")
            ensamble_phylo_no_inter[[k]] = ensamble_phylo

            #calculate the intersection between the phylogenetic correction of all variables
            #get a stack with the rasters of each variable separated
            bio4_raster=selected_stack[[which(grepl("bio4", names(selected_stack), fixed=TRUE))]]
            bio17_raster=selected_stack[[which(grepl("bio17", names(selected_stack), fixed=TRUE))]]

            #stop if we do NOT have 1 layer in each stack, as we should have only 1 model
            if(nlayers(bio4_raster)!=1 | nlayers(bio17_raster)!=1){ #nolint
                stop("ERROR! FALSE! WE HAVE A PROBLEM MAKING THE INTERSTCION OF BIO17 AND BIO4, WE HAVE MORE THAN 1 LAYER PER VARIABLE, BUT THIS SCRIPT ONLY CAN DEAL WITH 1 PER VARIABLE. ERROR. YOU HAVE TO SUM ALL RASTERS OF EACH ENV BIO STACK SO YOU HAVE ONE VALUE PER CELL AND KNOW WHAT CELLS ARE ZERO FOR ALL CLIMATIC PROJECTIONS FOR THAT ENV VARIABLE. In that way, when we multiply bio4*bio17, we get zero only in those cells where all bio4 and bio17 rasters have 0, leaving only cells with value for at least one phylo rasters in bio4 AND bio17")
            }

            #multiply both stacks to get the intersection
            #IMPORTANT: if having several rasters of the same env variable (e.g., future climate projections), then you should sum all rasters of the env bio stack (e.g., bio4)
            intersection=bio4_raster*bio17_raster
                #a cell being zero in any of the raster will be zero in the rest as zeros propagates

            #those cases with zero in the intersection, i.e., areas outside of the phylo range in at least one phylogenetic range
            ensamble_phylo[intersection == 0] <- 0

            #rename the ensemble phylo now indicating that we have selected only regions within the intersection of all env variables, then save
            names(ensamble_phylo) = paste(selected_stack_name, "_inter", sep="")
            ensamble_phylo_inter[[k]] = ensamble_phylo
        }

        #check we have the correct names in the new rasters, having one for each type of phylo-approach
        if(
            !identical(sapply(ensamble_phylo_no_inter, {function(x){names(x)}}), paste(stack_name, "_no_inter", sep="")) | #nolint
            !identical(sapply(ensamble_phylo_inter, {function(x){names(x)}}), paste(stack_name, "_inter", sep=""))
            ){
            stop("ERROR! FALSE! WE HAVE A PROBLEM GENERATING THE STACK OF PHYLO CORRECTIONS")
        }

        #check that non-proportion rasters are indeed not proportion, i.e., they only have two values for intersection (0 and 1) or three values for non-intersection (0, 0.5 and 1). This is because:
            #with intersection, regions with zero for at least one raster of bio4 or bio17 are removed. Given that we only have one raster for bio17 and one for bio4, if a cell falls withing the phylo-range of bio4, it will also fall within the range of bio17. In other words, we do bio4*bio17, so cells with zero in any of the two rasters are going to be zero. Therefore, the result will be 1 or 0, fall in both rasters o in no one.
            #without intersection, we can have a cell falling within one range but not in other. In that case, the phylo-suitability will be the average of 1 and 0, i.e., 0.5. Therefore, we can have 0, 1 and 0.5
        #raster_layer=ensamble_phylo_no_inter[[1]]
        #raster_layer=ensamble_phylo_inter[[1]]
        fun_check_data_range=function(raster_layer){
            not_proportion=!grepl("proportion", names(raster_layer), fixed=TRUE)
            not_inter=grepl("no_inter", names(raster_layer), fixed=TRUE)
            if(not_proportion){
                if(not_inter){
                    if(FALSE %in% (unique(getValues(raster_layer)) %in% c(0, 0.5, 1))){
                        stop("ERROR! FALSE! WE HAVE A PROBLEM WHEN CALCULATING THE PHYLOGENETIC CORRECTION ACROSS ENV VARIABLES")
                    }
                } else {
                    if(FALSE %in% (unique(getValues(raster_layer)) %in% c(0, 1))){
                        stop("ERROR! FALSE! WE HAVE A PROBLEM WHEN CALCULATING THE PHYLOGENETIC CORRECTION ACROSS ENV VARIABLES")
                    }   
                }
            }
        }
        sapply(ensamble_phylo_no_inter, fun_check_data_range)
        sapply(ensamble_phylo_inter, fun_check_data_range)

        #create stacks
        ensamble_phylo_no_inter_stack=stack(ensamble_phylo_no_inter)
        ensamble_phylo_inter_stack=stack(ensamble_phylo_inter)

        #open folder to save the stacks
        system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_phylo/", species, "/phylo_stacks/" , sep=""))

        #save them
        writeRaster(ensamble_phylo_no_inter_stack, filename=paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/phylo_stacks/", species, "_phylo_stacks_no_inter", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
            #ensamble_phylo_no_inter_stack=stack(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/phylo_stacks/", species, "_phylo_stacks_no_inter", sep=""))
        writeRaster(ensamble_phylo_inter_stack, filename=paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/phylo_stacks/", species, "_phylo_stacks_inter", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
            #ensamble_phylo_inter_stack=stack(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/phylo_stacks/", species, "_phylo_stacks_inter", sep=""))


        ##evaluate predictions with phylo-correction
        #decompress the predictions of glm, gam and rf
        print(system(paste(" \\
            tar \\
                --verbose \\
                --wildcards \\
                --directory ./results/global_test_phylo_current/predict_eval_phylo/", species, " \\
                --file ./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/", species, "_prediction_rasters.tar.gz \\
                --gunzip \\
                --extract \\
                ./prediction_rasters/", species, "_predictions_glm* \\
                ./prediction_rasters/", species, "_predictions_gam* \\
                ./prediction_rasters/", species, "_predictions_rf* \\
                ./prediction_rasters/", species, "_predictions_bin_glm* \\
                ./prediction_rasters/", species, "_predictions_bin_gam* \\
                ./prediction_rasters/", species, "_predictions_bin_rf*", sep=""), intern=TRUE))
            #--verbose: Verbosely list files processed
            #--wildcards: use wildcards
                #we are going to decompress the grd and gri files for each model (e.g., glm*...)
            #--directory: select where to decompress
            #--file: Use archive file
            #--gunzip: filter through gunzip
            #--extract: Extract files from an archive. Arguments are optional. When given, they specify names of the archive members to be extracted.
                #https://www.man7.org/linux/man-pages/man1/tar.1.html

        #load presences
        presences=read.table(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/", species, "_env_predictors.tsv.gz", sep=""), sep="\t", header=TRUE)

        #make folders to save results
        system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/plots/", sep=""))
        system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/boyce_partitions", sep=""))
        system(paste("mkdir -p ./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters_phylo/", sep=""))

        #load the thresholds used to binarize GLM and GAM
        load(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/thresholds/", species, "_glm_threshold.rda", sep=""))
        load(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/thresholds/", species, "_gam_threshold.rda", sep=""))

        #load Boyce non-phylo so we can now where the last bin ends and we can put cells within the phylo-range within that bin (see below)
        load(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/boyce_partitions/", species, "_glm_boyce.rda", sep=""))
        load(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/boyce_partitions/", species, "_gam_boyce.rda", sep=""))
        load(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/boyce_partitions/", species, "_rf_boyce.rda", sep=""))

        #create a vector with the names of all phylo-approaches
        phylo_models=c(
            paste(stack_name, "_inter", sep=""),
            paste(stack_name, "_no_inter", sep="")
            )
        if(length(phylo_models)!=8){
            stop("ERROR! FALSE! WE ARE NOT ANALYZING THE CORRECT NUMBER OF PHYLO-APPROACHES")
        }

        #set a sequence with the number of partitions
        data_partitions=1:12

        #open DF to save boyce results. First column is partition:
        final_boyce_results=data.frame(
            partition=c(
                seq(1,12,1), 
                "percentile_2.5", 
                "percentile_50", 
                "percentile_97.5") )

        #empty lists for raster predictions. We only need one list for all phylo models, because we will only save the rasters for the phylo approach used in the manuscript
        glm_phylo_list=list()
        gam_phylo_list=list()
        rf_phylo_list=list()
        glm_phylo_bin_list=list()
        gam_phylo_bin_list=list()
        rf_phylo_bin_list=list()

        #for each phylo-approach
        #l=1
        for(l in 1:length(phylo_models)){ # nolint

            #select the [l] phylo-approach
            selected_phylo_model=phylo_models[l]
            
            #extract the corresponding raster layer
            if(grepl("no_inter", selected_phylo_model, fixed=TRUE)){
                selected_ensemble_phylo=ensamble_phylo_no_inter_stack[[which(names(ensamble_phylo_no_inter_stack)==selected_phylo_model)]]
            } else {
                selected_ensemble_phylo=ensamble_phylo_inter_stack[[which(names(ensamble_phylo_inter_stack)==selected_phylo_model)]]  
            }

            #empty list to save the results of boyce for each model
            glm_phylo_boyce_list=list()
            gam_phylo_boyce_list=list()
            rf_phylo_boyce_list=list()

            #for each data partition
            #k=1
            for(k in data_partitions){
                
                #we have the same order of partitions here as we created predictions_ iterating over the list of partitions from 1 to 12
                glm_predict = stack(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters/", species, "_predictions_glm.grd", sep=""), bands=k)[[1]]
                gam_predict = stack(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters/", species, "_predictions_gam.grd", sep=""), bands=k)[[1]]
                rf_predict = stack(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters/", species, "_predictions_rf.grd", sep=""), bands=k)[[1]]
                    #we use the number of layer to extract it, getting a stack of 1 layer. Then we have to extract it to get the raster layer instead of a stack

                #get a value of suitability between the end of the penultimate suitability bin and the last bin, but as close as possible to the latter
                #model_name="glm"
                suit_value_phylo=function(model_name){
                    
                    #select the boyce eval for the corresponding model and partition
                    selected_boyce=eval(parse(text=paste(model_name, "_eval[[k]]$bins", sep="")))
                        #careful if you use the data stored in bins, as this include all continuous bins with the same P/E ratio that we have NOT considered, as we used rm.dup.class=TRUE

                    #calculate the number of rows, i.e., bins
                    n_bins=nrow(selected_boyce)
                    
                    #get suitability at the last two bins
                    end_last_bin=selected_boyce[n_bins,]$bin.max
                    end_penultimate_bin=selected_boyce[n_bins-1,]$bin.max
                    
                    #calculate the difference between the two ends
                    distance_last_bins=end_last_bin-end_penultimate_bin
                    
                    #add to the penultimate end the 99% of the difference so we get very close to the end of the last bin but without touching it
                    suit_value_for_phylo=end_penultimate_bin+(distance_last_bins*0.99)
                    return(suit_value_for_phylo)
                }
                    #we are going to select the a suitability value as high as possible without surpassing the end of the last bind of each Boyce evaluation. This will be the suitability value of cells withing the phylogenetic range.
                    #In this way, we ensure that all cells within the phylo range are indeed considered in Boyce (the last bin usually ends at 0.99), while having a high value of suitability (around 0.99). So, if for example a cell has suitability of 1 according to the SDMs, it will continue to be highly suitable. If the cell is non-suitable, then the phylo correction is going to increase a lot its suitability, and Boyce is going to be affected by that having a last bin with more bins (and likely more occurrences if the phylo-correction is working properly for a given species).
                #apply the function for all models
                suit_value_phylo_across_models=sapply(c("glm", "gam", "rf"), suit_value_phylo)
                    #we are going to use the same function for the three models, so we are going to use sapply to apply the function to each model. We are going to use the result to select the value of suitability for cells within the phylogenetic range.

                #phylo-conversion of the models
                glm_phylo=glm_predict
                glm_phylo[selected_ensemble_phylo>=0.1]=suit_value_phylo_across_models["glm"]
                gam_phylo=gam_predict
                gam_phylo[selected_ensemble_phylo>=0.1]=suit_value_phylo_across_models["gam"]
                rf_phylo=rf_predict
                rf_phylo[selected_ensemble_phylo>=0.1]=suit_value_phylo_across_models["rf"]
                    #we can use selected_ensemble_phylo as condition because the phylo rasters were calculated using the res and extent of the ensemble of binary predictions of the three models and, of course, these binary predictions come from the continuos predictions we are using here. Therefore, there is a match between the cells of these rasters
                    #IMPORTANT: if you change this step, the way to combine phylo and models, then maybe you have to change the way you calculate the binary RF phylo below because you just directly take the binary RF prediction and convert to 1 those cells with phylo-suit above 0.1
                    #we are not completely binarizying as we will still have intermediate values, which are required to calculate a correlation between suitability and number of presence points
                    #this is also the closest thing to what we do in the analyses of the MS
                        #remember that we cannot just calculate the average because we have zero for most of the regions in the phylo, so the average would be very low in almost every place. 
                        #also, phylo-suitability is not the same than the SDMs suitability! Here, a value above 0 means that the climate is inside the phylogenetic range! yes, maybe far away from the current status, but still after the first split with the closest sister species.

                #Important notes about the approach to combine SDMs and phylo
                    #we are changing the value of suitability of cells within the phylogenetic range to ~1 (0.99), moving cells from intermediate suitability bins to the last bin. Therefore, we are moving cells from bins of intermediate suitability to bins with max suitability. If these cells have presence points, boyce index will increase, but if they don't, then the boyce index will decrease. In that latter case, we are putting noise in the correlation between suitability and the probability of presence because now more cells with high suitability do not have points.
                    #For an example of what is happening, we can explain the cases of P. strobus for proportion_no_subset_inter, which shows the greatest impact of the phylogenetic correction. We will explain GAM, but GLM and RF also show a lot of change with phylo.
                        #GAM without phylo shows a peak of proportion presences and P/Eratio around 0.31 of suitabilty. From there, the proportion of presences and P/E ratio decreases steadily, generating a negative correlation. Lowergam_phylo probability to find presences as suitability increases.
                        #GAM with phylo, occurrences that where from the bin 63 (median suitability 0.60) were all moved to to high suitability (yes, the phylo-correction has correctly selected all cells with ocurrences in intermediate suitable bins! that is really cool). This makes the proportion of presences for these bins zero and P/E ratio is also zero (you are dividing by zero). Now, we have instead of a steadly decrease of P/E ratio until 1.79, you stop at a P/E ratio of 4.47. Then you have a bunch of P/E values of zero, from which only 1 is counted, because we are only interested in the ability of the model to discriminate low vs high suitability. If you have a bunch of low P/E values around 0.6 suitability, you just need to take one to know the model is not predicting well there. Then, you have a point of very high proportion of presences and P/E ratio at the last bin where occurrences has been sent, but only if you set the suitability for cells within phylo as 0.99, because the last bin ends at 0.99000.... instead of one. This is only one point in the correlation, one cell, because of this, sending the cells to 1 or 0.99 does not make a lot of change. Again, this is really cool, because it confirms the correction have been collecting occurrences from bins of intermediate suitability and then accumulate them in the last bin, increasing the proportion of presences more than the proportion of pixels in that bin, leading to a peak in P/E ratio. In other words, although the number of pixel in that bin increases (because we are increasing the number of high suitability cells), the number of presences increases even more in relation to the total of presences. Putting the occurrences into the last bin makes the correlation increase from 0.39 to 0.4.
                        #Note that in general, the P/E ratio decreases (y axis), and that makes sense because we are removing presences from intermediate suitability to bins to the last bin. That is ok, becuase the important thing here is the pattern, the correlation between P/E ratio and suitability.
                    #Strobus RF
                        #Same story, there is a steadily decrease of P/E ratio from bin 61 to the end, generating a negative correlation between P/E and suitability.
                        #the correction collects occurrences from bin 61, leading to zero presence probability until bin 89. These are all have the same value and are contiguous, so only one value is considered.
                        #From bin 89, we still have occurrences, but the proportion is smaller compared to the original model.
                        #if we set suitability of cells within the phylo range to 0.99, the last bin show a peak of proportion presences and P/E ratio, increasing the correlation from 0.35 to 0.39.
                        #Again, we have a decrease of the range of P/E ratio, but this is not a problem because the pattern between P/E ratio and suitability improves.
                    #sylvestris GLM
                        #The base model is good (Boyce=0.771), but phylo improves a bit (Boyce=0.825) because it reduces the P/E ratio of bins around 0.2, making the correlation between P/E and suitability more positive.
                        #In other words, it is correctly reducing the proportion of presences respect to the proportion of pixels in low-suitabilty bins where, indeed, we expect a lower P/E ratio is the model works well.
                        #Putting occurrences in the last bin or not does not change much, in both cases we are around 0.82.
                    #sylvestris RF
                        #The base model is good (Boyce=0.52), but phylo improves a bit (Boyce=0.55) because it reduces the peak of P/E ratio of bins around 0.4, making the correlation between P/E and suitability more positive.
                        #Again, it is reducing the proportion of presences respect to the proportion of pixels in low-suitabilty bins where, indeed, we expect a lower P/E ratio is the model works well.
                        #Putting occurrences in the last bin or not does not change much, in both cases we are around 0.55-0.56.
                    #This is totally correct
                        #we are taking cells from intermediate suitability bins that were making the model fail. These presences were distributed in a way that made the P/E ratio to decrease as suitability of the bin increases, i.e., the model was not able to correctly increase the probability of presence as suitability increases.
                        #what the PA correction has done is to take ALL cells with occurrences in these intermediate bins and move them to the last, highest suitability bin. SO IT HAS BEEN CORRECTLY SELECTING CELLS WITH OCCURRENCES and moving them to high suitability, THIS IS GREAT! Because of this we have now a bettere correlation, we have made a cleaner correlation.
                        #Also note that we are applying the same phylogenetic-map to the predictions of the three models, so the new boyce values are correlated, but we can still have differences between models of the same species. A region that is included in the phylogenetic correction could be missed by GLM while at the same time it could be considered as suitable by RF. In the first case the phylo correction would improve Boyce, but not in the second case, as the the cell remain with the same suitability and hence it remains in the same bin.
                        #Also note that we do not have a problem with other extreme, i.e., suitability around 0. We are not moving cells to low-suitability giving them exactly zero, thus the value of low-suitability cells will be 0.00000.... and that is the value for start of the first bin. So it is very likely that very low-suitable cells are included in the first bin. We are not biasing the results by adding cells exact 0, leaving them outside of the first bin, while including high suitable cells in the last bin. This is not the case.

                #add the names of the rasters
                names(glm_phylo)=paste("glm_", selected_phylo_model, "_part_", k, sep="")
                names(gam_phylo)=paste("gam_", selected_phylo_model, "_part_", k, sep="")
                names(rf_phylo)=paste("rf_", selected_phylo_model, "_part_", k, sep="")

                #check whether the suitability value given to cells within the phylo-range is still within the last bin of boyce after generating the new phylo-raster
                #define function to do the check
                #model_to_check="gam"
                check_last_bin_boyce=function(model_to_check){
                    
                    #get the phylo and non-phylo rasters for the selected model
                    selected_non_phylo_raster=eval(parse(text=paste(model_to_check, "_predict", sep="")))
                    selected_phylo_raster=eval(parse(text=paste(model_to_check, "_phylo", sep="")))

                    #if the min or the max of suitability in the new phylo raster are not same compared to the original prediction raster, then the range of pred is changed and hence the bins in boyce will be different, so the suitability value given to cells within the phylo range could be again outside of the last bin, and hence it will not be considered. The min/max can change if all the cells with the lowest/highest suitaibility value are within the phylo range. 
                        #For example, if the highest suitability is 1 being in three cells. These three cells are within the phylogenetic range, so we give them the new and lower suitability value so they can enter in the last bin of boyce (remember, 1 is outside of the last bin in Boyce). Now these cells have a value lower than 1 and, given they were the ones with the highest suitability, now the maximum suitability of the raster has decreased. They can still be the highest value or other cell (in or outside the phylo range) could take the first place. It does not matter what cells take the first place. When we calculate Boyce in this new data, we will have difference bins because the range of suitability is different now, that is the thing. The max is smaller, so the max of the last bin will be smaller, so the suitabiltiy value previously given to the phylo-cells can be now higher than the last bin of Boyce.
                        #The opposite case is less likely because we would need to have ALL cels with the minimum suitability (likely 0) within the phylogenetic range, but we are also accounting for that anyways. If the min value of suitability changes, the range of suitability changes, so we recalculate the bins using Boyce
                        #In both cases, if we calculate Boyce using the new range, we will get the final bins, and we can use them to select the suitability for phylo-cells to enter the last bin.
                    if(
                        max(getValues(selected_non_phylo_raster), na.rm=TRUE) != max(getValues(selected_phylo_raster), na.rm=TRUE) |  # nolint: vector_logic_linter.
                        min(getValues(selected_non_phylo_raster), na.rm=TRUE) != min(getValues(selected_phylo_raster), na.rm=TRUE)
                    ){

                        #run Boyce 
                        new_last_bin_boyce=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(selected_phylo_raster), n.bins=NA, bin.width="default", res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=FALSE, na.rm=TRUE, plot=FALSE)

                        #calculate the number of rows, i.e., bins
                        n_bins=nrow(new_last_bin_boyce$bins)
                        
                        #get suitability at the last two bins
                        end_last_bin=new_last_bin_boyce$bins[n_bins,]$bin.max
                        end_penultimate_bin=new_last_bin_boyce$bins[n_bins-1,]$bin.max
                        
                        #calculate the difference between the two ends
                        distance_last_bins=end_last_bin-end_penultimate_bin
                        
                        #add to the penultimate end the 99% of the difference so we get very close to the end of the last bin but without touching it
                        new_suit_value_for_phylo=end_penultimate_bin+(distance_last_bins*0.99)

                        #give the new suitability value to those cells within the phylo-range
                        selected_phylo_raster[selected_ensemble_phylo>=0.1]=new_suit_value_for_phylo

                        #save the new phylo raster overwritting the previous phylo_raster in the global environment
                        assign(paste(model_to_check, "_phylo", sep=""), selected_phylo_raster, envir=.GlobalEnv)
                    }

                    #calculate again boyce to check whether the phylo-cells are within the last bin
                    final_last_bins=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(selected_phylo_raster), n.bins=NA, bin.width="default", res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=FALSE, na.rm=TRUE, plot=FALSE)

                    #get the suitability of phylo-cells
                    final_suitability_phylo=unique(selected_phylo_raster[selected_ensemble_phylo>=0.1])
                    if(length(final_suitability_phylo)!=1){
                        stop(paste("ERROR! FALSE! WE HAVE A PROBLEM SELECTING THE NEW SUITAIBLITY FOR PHYLO-CELLS, ALL THEY SHOULD HAVE THE SAME SUITABILITY, BUT THIS IS NOT CASE: ", species, "-", selected_phylo_model, "-", k, "-", model_to_check, sep=""))
                    }

                    #update the final_suitability value for phylo-cells
                    if(final_suitability_phylo!=suit_value_phylo_across_models[model_to_check]){
                        suit_value_phylo_across_models[model_to_check]=final_suitability_phylo
                    }
                    if(
                        final_suitability_phylo!=suit_value_phylo_across_models[model_to_check] | # nolint: vector_logic_linter.
                        suit_value_phylo_across_models[model_to_check]<0.95
                    ){
                        stop(paste("ERROR! FALSE! WE HAVE A PROBLEM SELECTING THE NEW SUITAIBLITY FOR PHYLO-CELLS, NOT CORRECTLY CHANGED OR TOO LOW, I.E., BELOW 0.95: ", species, "-", selected_phylo_model, "-", k, "-", model_to_check, sep=""))
                    }

                    #now do the check about whether the new suitability value for phylo-cells is indeed in the last bin of boyce
                    if(
                        !(final_suitability_phylo > final_last_bins$bins[99,]$bin.max &  # nolint: vector_logic_linter.
                        final_suitability_phylo < final_last_bins$bins[100,]$bin.max)
                    ){
                        print(paste("WARNING! WE HAVE BEEN UNABLE TO PUT THE CELLS OF THE PHYLOGENETIC RANGE IN THE LAST BIN FOR: ", species, "-", selected_phylo_model, "-", k, "-", model_to_check, ". The problem is that the cells inside the phylo-range are the one with the highest value of suitability. When we give them a new, lower value to enter within the last bin of Boyce, we reduce the max value of suitability for the whole raster. The script has tried to calculate again the bins considering the range of suitability values (new max), but the problem persists, likely because the phylo-cells are still the ones with the highest value, so the second time we gave them a new suitability value, we reduced again the max suitability of the raster. We could continue doing this until a cells outside the phylo-range has a higher suitability but this has the risk of reducing too much the suitability of phylo-cells. From what I have seen, the greatest impact of the correction comes from the removel of occurrences from low and intermediate suitability bins. Adding these occurrences to the last, high suitability bin usutually increases Boyce a bit, but it is much less important that the previous removal of occurrences. Considering this, and the fact that doing this process recursively could put again phylo-cells in intermediate or low bins, lead me to stop here.", sep=""))
                    }
                }
                #apply function across models
                sapply(c("glm", "gam", "rf"), check_last_bin_boyce)

                #calculate the binary predictions for ensemble ONLY FOR THE PHYLO MODEL WE ARE GOING TO USE IN THE PAPER. WE WILL HAVE BOYCE FOR THE REST OF MODELS, BUT NOT ENSEMBLES SO WE SAVE SPACE AND TIME
                if(selected_phylo_model=="phylo_rasters_proportion_subset_inter"){
                    
                    #add cells within the phylo range to the binary projections, i.e., value of 1 instead of 0
                    #glm
                    glm_predict_bin = stack(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters/", species, "_predictions_bin_glm.grd", sep=""), bands=k)[[1]]
                    glm_phylo_bin = glm_predict_bin
                    glm_phylo_bin[selected_ensemble_phylo>=0.1]=1
                        #we can use selected_ensemble_phylo as condition because the phylo rasters were calculated using the res and extent of the ensemble of binary predictions of the three models, thus there is a match between the cells of these rasters

                    #gam
                    gam_predict_bin = stack(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters/", species, "_predictions_bin_gam.grd", sep=""), bands=k)[[1]]
                    gam_phylo_bin = gam_predict_bin
                    gam_phylo_bin[selected_ensemble_phylo>=0.1]=1

                    #rf
                    rf_predict_bin = stack(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters/", species, "_predictions_bin_rf.grd", sep=""), bands=k)[[1]]
                    rf_phylo_bin = rf_predict_bin
                    rf_phylo_bin[selected_ensemble_phylo>=0.1]=1
                    
                    #IMPORTANT HERE: 
                        #We are directly using the binary prediction of RF and this is ok.
                        #think what we do with glm.
                            #take the continuous prediction, convert to high suitability those cells with phylo-suitability > 0.1.
                            #then apply the threshold to binarize
                            #if the suitability value given to cells within the phylo-range is equal or above the threshold
                                #these cells are going to be considered 1 
                                #therefore, it is the same than if we binarize glm and then convert to 1 those cells with phylo-suitability>0.1
                            #if the suitability value given to cells within the phylo-range is below the threshold, then we could have differences
                                #if that is the case, we would have to use the suitability value given to phylo-cells as threshold. All cells above that value would be considered 1, even those outside of the phylogenetic range and below the original threshold, if they still are above the suitability value given to cells within the phylo-range.
                            #our current approach, i.e., applying phylo after binarization, is cleaner. Doing that, outside the phylo-range, only cells above the original threshold are considered as suitable, and then those witin phylo are considered as suitable even if their suibility value is below the threshold. This is what we want
                                #Note that the suitability value given to cells within the phylo-range will be very high anyways, as we selected a value very close to the end of the last suitability bin in Boyce calculations, i.e., very close to 0.99...

                    #check whether the threshold used to create binary predictions is above the value of suitability given to cells within the phylo range
                        #if it is, just print a warning, but not problem because we are applying the phylo-correction after binarizing, including the phylo-suitable areas in that way
                        #if it is not, we can check whether applying before or after gives the same result
                    #model_to_check="glm"
                    check_threshold_vs_phylo=function(model_to_check){
                        
                        #select the threshold of the model
                        selected_threshold=eval(parse(text=paste(model_to_check, "_threshold[[k]][2,2]", sep="")))
                        
                        #if the threshold is above the suitability of cells within the phylo-range     
                        if(selected_threshold > suit_value_phylo_across_models[model_to_check]){

                            #print warning
                            print(paste("WARNING! The threshold for ", model_to_check, " is above the suitability given to the cells inside the phylo-range for ", species, " in ", selected_phylo_model, sep=""))

                            #print the threshold and phylo-suitability values
                            print(paste("The threshold value is ", selected_threshold, " while the suitability value given to cells inside the phylo-range is ", suit_value_phylo_across_models[model_to_check], sep=""))
                        } else { 

                            #means that the suitability value given to cells within the phylo-range is equal or higher than the threshold so, in that case, applying the phylo correction before or after binarization should give the same result

                            #get the original phylo bin raster
                            phylo_raster_bin_check_1=eval(parse(text=paste(model_to_check, "_phylo_bin", sep="")))

                            #get the original phylo continuous raster
                            phylo_raster_continuous_check=eval(parse(text=paste(model_to_check, "_phylo", sep="")))

                            #binarize the continuous phylo raster
                            phylo_raster_bin_check_2 = phylo_raster_continuous_check
                            phylo_raster_bin_check_2[phylo_raster_bin_check_2>=selected_threshold,] <- 1
                            phylo_raster_bin_check_2[phylo_raster_bin_check_2<selected_threshold,] <- 0

                            #check both approaches are identical
                            if(!identical(getValues(phylo_raster_bin_check_1), getValues(phylo_raster_bin_check_2))){
                                stop("ERROR! FALSE! WE HAVE A PROBLEM DURING THE BINARIZATION OF PHYLO-MODELS")
                            }
                        }
                    }
                    check_threshold_vs_phylo("glm")
                    check_threshold_vs_phylo("gam")

                    #add names to the rasters
                    names(glm_phylo_bin)=paste("glm_", selected_phylo_model, "_bin_part_", k, sep="")
                    names(gam_phylo_bin)=paste("gam_", selected_phylo_model, "_bin_part_", k, sep="")
                    names(rf_phylo_bin)=paste("rf_", selected_phylo_model, "_bin_part_", k, sep="")

                    #save the rasters the corresponding lists
                    glm_phylo_list[[k]]=glm_phylo
                    gam_phylo_list[[k]]=gam_phylo
                    rf_phylo_list[[k]]=rf_phylo
                    glm_phylo_bin_list[[k]]=glm_phylo_bin
                    gam_phylo_bin_list[[k]]=gam_phylo_bin
                    rf_phylo_bin_list[[k]]=rf_phylo_bin
                }

                #calculate the boyce index (see previous major step for further details)
                jpeg(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/plots/", species, "_part_", k, "_", selected_phylo_model, "_boyce_index_plot.jpeg", sep=""), width=960, height=960, pointsize=24)
                par(mfcol=c(2,2))
                glm_phylo_boyce=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(glm_phylo), n.bins=NA, bin.width="default", res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=FALSE, na.rm=TRUE, plot=TRUE, main=paste("GLM", sep=""))
                mtext(paste(species, " - part ", k, " - ", selected_phylo_model, sep=""), side=3, line=-1.6, outer=TRUE, cex=1, font=2)
                gam_phylo_boyce=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(gam_phylo), n.bins=NA, bin.width="default", res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=FALSE, na.rm=TRUE, plot=TRUE, main=paste("GAM", sep=""))
                rf_phylo_boyce=modEvA::Boyce(obs=presences[,c("longitude", "latitude")], pred=terra::rast(rf_phylo), n.bins=NA, bin.width="default", res=100, method="spearman", rm.dup.classes=TRUE, rm.dup.points=FALSE, na.rm=TRUE, plot=TRUE, main=paste("RF", sep=""))
                dev.off()

                #add the partition to the boyce results
                glm_phylo_boyce$partition=k
                gam_phylo_boyce$partition=k
                rf_phylo_boyce$partition=k

                #save the boyce results
                glm_phylo_boyce_list[[k]]=glm_phylo_boyce
                gam_phylo_boyce_list[[k]]=gam_phylo_boyce
                rf_phylo_boyce_list[[k]]=rf_phylo_boyce


                ##check we do not have a low number of pixels in any bin
                #In the docs of "modEvA::Boyce", they say the following: "In bins with overly small sample sizes, the comparison between median prediction and random expectation may not be meaningful, although these bins will equally contribute to the overall Boyce index. When there are bins with less than 30 values, a warning is emitted and their points are plotted in red, but mind that 30 is a largely arbitrary number. See the $bins$bin.N section of the console output, and use the 'bin.width' argument to enlarge the bins."
                #For species with a small natural range, it would be possible that no pixel has high suitability values, making a low sample size for high suitability bins. But this should not be very problematic in our case, because we are predicting across the whole world, increasing the probability to find enough pixels with low and high suitability. For example, clausa has a small range but its minimum pixel number per bin is 70. We are going to check that just in case
                min_bin_n_glm_phylo=min(glm_phylo_boyce[[1]]$bin.N) 
                min_bin_n_gam_phylo=min(gam_phylo_boyce[[1]]$bin.N) 
                min_bin_n_rf_phylo=min(rf_phylo_boyce[[1]]$bin.N) 
                    #calculate the min number of points across bins in each model. You can have different sample size between models because these are bins based on suitability, and the suitability is differently predicted by each model
                #we are going stop the analyses if the minimum number of pixels per bin is 60, which is the double of the default value set by "modEvA::Boyce" to print a warning (i.e., 30)
                if((min_bin_n_glm_phylo<60) | (min_bin_n_gam_phylo<60) | (min_bin_n_rf_phylo<60)){
                    print(paste("ERROR! FALSE! WE HAVE AT LEAST 1 PHYLO SUITABILITY BIN WITH LESS THAN 60 DATA.POINTS FOR PHYLO MODEL ", selected_phylo_model, " AND SPECIES ", species, ". THIS CAN MAKE THAT THE COMPARISON BETWEEN MEDIAN PREDICTION AND RANDOM EXPECTATION MAY NOT BE MEANINGFUL", sep=""))
                }
            }
            
            #check
            if(
                (length(glm_phylo_boyce_list)!=length(data_partitions)) |
                (length(gam_phylo_boyce_list)!=length(data_partitions)) |
                (length(rf_phylo_boyce_list)!=length(data_partitions))){
                stop("ERROR! FALSE! WE DO NOT HAVE BOYCE INDEX FOR ALL PARTITIONS")
            }

            #save evaluations
            save(glm_phylo_boyce_list, file=paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/boyce_partitions/", species, "_", selected_phylo_model, "_glm_boyce.rda", sep=""))
            save(gam_phylo_boyce_list, file=paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/boyce_partitions/", species, "_", selected_phylo_model, "_gam_boyce.rda", sep=""))
            save(rf_phylo_boyce_list, file=paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/boyce_partitions/", species, "_", selected_phylo_model, "_rf_boyce.rda", sep=""))

            #check we have the same partitions in all lists
            #x=glm_phylo_boyce_list[[1]]
            glm_partitions=sapply(glm_phylo_boyce_list, function(x){return(x[[3]])})
            gam_partitions=sapply(gam_phylo_boyce_list, function(x){return(x[[3]])})
            rf_partitions=sapply(rf_phylo_boyce_list, function(x){return(x[[3]])})
            if((!identical(glm_partitions, gam_partitions)) | (!identical(glm_partitions, rf_partitions))){ # nolint: vector_logic_linter.
                stop("ERROR! FALSE! WE HAVE A PROBLEM WITH THE CALCULATION OF BOYCE ACROSS PARTITIONS")
            }

            #calculate a table with boyce index across models and partitions
            #add the boyce index of each combination as a row
            #x=glm_phylo_boyce_list[[1]]
            boyce_table=cbind.data.frame(
                cbind.data.frame(glm_partitions),
                cbind.data.frame(sapply(glm_phylo_boyce_list, function(x){return(x[[2]])})), 
                cbind.data.frame(sapply(gam_phylo_boyce_list, function(x){return(x[[2]])})),
                cbind.data.frame(sapply(rf_phylo_boyce_list, function(x){return(x[[2]])})))
            names(boyce_table)=c("partition", paste(selected_phylo_model, c("_glm_", "_gam_", "_rf_"), "boyce", sep=""))
                #for each partition, get the second element which is the boyce index, obtaining a vector with all of them. 
                #Then, convert to a column of DF. 
                #repeate for all models and combine the different columns in a DF along with the number of the partition.

            #add the confidence interval as a new row
            #calculate the different percentiles per column
            low_interval_boyce=apply(X=boyce_table[,which(colnames(boyce_table) != "partition")], MARGIN=2, FUN=quantile, probs=0.025, na.rm=TRUE)
            median_boyce=apply(X=boyce_table[,which(colnames(boyce_table) != "partition")], MARGIN=2, FUN=quantile, probs=0.5, na.rm=TRUE)
            high_interval_boyce=apply(X=boyce_table[,which(colnames(boyce_table) != "partition")], MARGIN=2, FUN=quantile, probs=0.975, na.rm=TRUE)
                #calculate the median of each column (excluding the first one with the number of the partition)
                #margin lets you select if you want to apply the function across columns (2) or rows (1)
                #before calculating the quantiles of each, remove NAs (na.rm=TRUE).
                    #We can have NA for boyce if, for example, the phylo correction captures ALL occurrences in low and intermediate suitability bins, and put them above the last bin of Boyce, because the suitability value given to them is still above that bin.
                    #Remember that we try to give phylo-cells a suitability within the last bin of Boyce, but sometimes is not possible because the highest suitability value is still within these phylo-cells, so any change we do to them will change the max value of suitability of the raster and hence the range used for Boyce to calculate the bins. Doing that recursively wcould reduce a lot the suitability of the phylo cells, somethign that would be misleading as the greatest impact of the correction comes from the removal of occurrences from intermediate suitabiltiy bins. See above for code delaing with this.

            #add them as new rows. Do it in two steps to avoid converting numeric into string, as we have to add a string ("percentile_XX") and numbers (the actual percentile value)
            boyce_table[nrow(boyce_table) + 1, 1]="percentile_2.5"
            boyce_table[nrow(boyce_table), 2:ncol(boyce_table)]=low_interval_boyce
            boyce_table[nrow(boyce_table) + 1, 1]="percentile_50"
            boyce_table[nrow(boyce_table), 2:ncol(boyce_table)]=median_boyce
            boyce_table[nrow(boyce_table) + 1, 1]="percentile_97.5"
            boyce_table[nrow(boyce_table), 2:ncol(boyce_table)]=high_interval_boyce
                #https://stackoverflow.com/a/44150746

            #check
            if(
                (ncol(boyce_table)!=4) |
                (nrow(boyce_table)!=length(data_partitions)+3)){
                stop("ERROR! FALSE! WE HAVE A PROBLEM WITH THE BOYCE TABLE")
            }

            #merge the boyce table with the empty data.frame to save results across phylo-approaches
            final_boyce_results=merge(final_boyce_results, boyce_table, by="partition")
        }

        #check we have generated all boyce plots
        n_boyce_plots=as.numeric(system(paste("ls ./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/plots/ | awk 'END{print NR}'", sep=""), intern=TRUE))
        if(n_boyce_plots!=length(phylo_models)*length(data_partitions)){
            stop("ERROR! FALSE! WE HAVE NOT ALL THE BOYCE PLOTS WE SHOULD HAVE")
        }

        #check we have the correct number of rasters in all lists. It should be only 12, because we are saving rasters only for the first phylo model, so 12 partitions, 12 rasters
        if(
            length(glm_phylo_list)!=length(data_partitions) |
            length(gam_phylo_list)!=length(data_partitions) |
            length(rf_phylo_list)!=length(data_partitions) |
            length(glm_phylo_bin_list)!=length(data_partitions) |
            length(gam_phylo_bin_list)!=length(data_partitions) |
            length(rf_phylo_bin_list)!=length(data_partitions)){
            stop(paste("ERROR! FALSE! WE HAVE NOT ANALYZED ALL PARTTITIONS FOR ", species, sep=""))
        }

        #check that we have only selected the first phylo-approach (i.e., the one used in the MS) and we have analyzed all 12 partitions
        if(
            !identical(sapply(glm_phylo_list, {function(x) names(x)}), paste("glm_phylo_rasters_proportion_subset_inter_part_", c(1:12), sep="")) |
            !identical(sapply(gam_phylo_list, {function(x) names(x)}), paste("gam_phylo_rasters_proportion_subset_inter_part_", c(1:12), sep="")) |
            !identical(sapply(rf_phylo_list, {function(x) names(x)}), paste("rf_phylo_rasters_proportion_subset_inter_part_", c(1:12), sep="")) |
            !identical(sapply(glm_phylo_bin_list, {function(x) names(x)}), paste("glm_phylo_rasters_proportion_subset_inter_bin_part_", c(1:12), sep="")) |
            !identical(sapply(gam_phylo_bin_list, {function(x) names(x)}), paste("gam_phylo_rasters_proportion_subset_inter_bin_part_", c(1:12), sep="")) |
            !identical(sapply(rf_phylo_bin_list, {function(x) names(x)}), paste("rf_phylo_rasters_proportion_subset_inter_bin_part_", c(1:12), sep=""))
        ){
            stop("ERROR! FALSE! THIS SCRIPT IS PREPARED TO DEAL WITH ONLY 1 PHYLO MODEL, SPECIFICALLY 'phylo_rasters_proportion_subset_inter'")
        }

        #make stacks
        predictions_glm=stack(glm_phylo_list[grep("glm_phylo_rasters_proportion_subset_inter_part_", sapply(glm_phylo_list, {function(x) names(x)}), fixed=TRUE)])
        predictions_gam=stack(gam_phylo_list[grep("gam_phylo_rasters_proportion_subset_inter_part_", sapply(gam_phylo_list, {function(x) names(x)}), fixed=TRUE)])
        predictions_rf =stack(rf_phylo_list[grep("rf_phylo_rasters_proportion_subset_inter_part_", sapply(rf_phylo_list, {function(x) names(x)}), fixed=TRUE)])
        predictions_glm_bin=stack(glm_phylo_bin_list[grep("glm_phylo_rasters_proportion_subset_inter_bin_part_", sapply(glm_phylo_bin_list, {function(x) names(x)}), fixed=TRUE)])
        predictions_gam_bin=stack(gam_phylo_bin_list[grep("gam_phylo_rasters_proportion_subset_inter_bin_part_", sapply(gam_phylo_bin_list, {function(x) names(x)}), fixed=TRUE)])
        predictions_rf_bin=stack(rf_phylo_bin_list[grep("rf_phylo_rasters_proportion_subset_inter_bin_part_", sapply(rf_phylo_bin_list, {function(x) names(x)}), fixed=TRUE)])

        #check
        if(
            !identical(names(predictions_glm), paste("glm_phylo_rasters_proportion_subset_inter_part_", c(1:12), sep="")) |
            !identical(names(predictions_gam), paste("gam_phylo_rasters_proportion_subset_inter_part_", c(1:12), sep="")) |
            !identical(names(predictions_rf), paste("rf_phylo_rasters_proportion_subset_inter_part_", c(1:12), sep="")) |
            !identical(names(predictions_glm_bin), paste("glm_phylo_rasters_proportion_subset_inter_bin_part_", c(1:12), sep="")) |
            !identical(names(predictions_gam_bin), paste("gam_phylo_rasters_proportion_subset_inter_bin_part_", c(1:12), sep="")) |
            !identical(names(predictions_rf_bin), paste("rf_phylo_rasters_proportion_subset_inter_bin_part_", c(1:12), sep=""))
        ){
            stop("ERROR! FALSE! THIS SCRIPT IS PREPARED TO DEAL WITH ONLY 1 PHYLO MODEL, SPECIFICALLY 'phylo_rasters_proportion_subset_inter'")
        }

        #prepare ensemble
        #stack all binary predictions 
        binary_predictions = stack(predictions_glm_bin, predictions_gam_bin, predictions_rf_bin)
        n_layers_binary_predictions=nlayers(binary_predictions)

        #calculate the percentage of models for which a pixel is suitable
        ensamble_predictions_bin = calc(binary_predictions, function(x) (sum(x)*100)/n_layers_binary_predictions)
            #ensamble_predictions_bin2 = calc(binary_predictions, function(x) (sum(x)*100)/nlayers(binary_predictions))
            #identical(getValues(ensamble_predictions_bin), getValues(ensamble_predictions_bin2))

        #plot ensembles with and without phylo
        pdf(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/", species, "_ensembles_plot.pdf", sep=""), width=12, height=12)
        par(mfcol=c(2,1), mar=c(3, 4, 3, 2) + 0.1)
        plot(ensamble_suitability, main="Without phylo")
        points(presences[,c("longitude", "latitude")], pch=20, col="red", cex=0.5)
        plot(ensamble_predictions_bin, main="With phylo")
        points(presences[,c("longitude", "latitude")], pch=20, col="red", cex=0.5)
        dev.off()

        #save stacks
        writeRaster(predictions_glm, filename=paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters_phylo/", species, "_phylo_glm", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        writeRaster(predictions_gam, filename=paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters_phylo/", species, "_phylo_gam", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        writeRaster(predictions_rf, filename=paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters_phylo/", species, "_phylo_rf", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        writeRaster(predictions_glm_bin, filename=paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters_phylo/", species, "_phylo_glm_bin", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        writeRaster(predictions_gam_bin, filename=paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters_phylo/", species, "_phylo_gam_bin", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        writeRaster(predictions_rf_bin, filename=paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters_phylo/", species, "_phylo_rf_bin", sep=""), options="COMPRESS=LZW", overwrite=TRUE)
        writeRaster(ensamble_predictions_bin, filename=paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/prediction_rasters_phylo/", species, "_phylo_ensemble", sep=""), options="COMPRESS=LZW", overwrite=TRUE)

        #prepare final boyce table
        #sort by numerical order in partition to avoid problems when comparing with non-phylo results
        final_boyce_results=final_boyce_results[mixedorder(final_boyce_results$partition),]
        rownames(final_boyce_results)=1:nrow(final_boyce_results)
            #this applies natural sorting
            #mixedsort(final_boyce_results$partition)
                #https://stackoverflow.com/a/2778060
        #check we have the correct number of columns and rows
        if(ncol(final_boyce_results)!=length(phylo_models)*3+1 | nrow(final_boyce_results)!=length(data_partitions)+3){
            stop("ERROR! FALSE! WE HAVE NOT ANALYZED ALL PHYLO MODELS AND PARTTITIONS")
        }


        ###PROBLEM here with strobus
        #check problem with correlation proporion non-proportion in strobus


        #check for each algorithm and phylo-approach that proportion and non-proportion is the same
        #model_type="glm"
        for(model_type in c("glm", "gam", "rf")){
            
            #phylo_option="subset_inter"
            for(phylo_option in c("subset_inter", "no_subset_inter", "subset_no_inter", "no_subset_no_inter")){

                #extract the proportion and non-proportion columns for the selected algorithm and phylo-approach
                non_proportion_column=final_boyce_results[
                    which(!final_boyce_results$partition %in% c("percentile_2.5", "percentile_50", "percentile_97.5")),
                    which(colnames(final_boyce_results)==paste("phylo_rasters_", phylo_option, "_", model_type, "_boyce", sep=""))]
                proportion_column=final_boyce_results[
                    which(!final_boyce_results$partition %in% c("percentile_2.5", "percentile_50", "percentile_97.5")),
                    which(colnames(final_boyce_results)==paste("phylo_rasters_proportion_", phylo_option, "_", model_type, "_boyce", sep=""))]
                
                #if all observations are NA
                if(sum(is.na(non_proportion_column))==length(non_proportion_column) | sum(is.na(proportion_column))==length(proportion_column)){
                    
                    #do not calculate the correlation
                    correlation=list()
                    correlation$estimate=NA
                } else {
                    correlation=cor.test(non_proportion_column, proportion_column, method="spearman")
                }

                #check
                if(!is.na(correlation$estimate)){
                    if(correlation$estimate<0.90){
                        print(paste("WARNING! PROPORTION AND NON-PROPORTION ARE NOT SIMILAR AS EXPECTED. CORRELATION IS ", correlation$estimate, " FOR SPECIES ", species, " model ", model_type, " phylo option ", phylo_option, sep=""))
                    }
                } else { #we also print the warning if the correlation is NA
                    print(paste("WARNING! PROPORTION AND NON-PROPORTION ARE NOT SIMILAR AS EXPECTED. CORRELATION IS ", correlation$estimate, " FOR SPECIES ", species, " model ", model_type, " phylo option ", phylo_option, sep="")) 
                }
            }
        }

        #make a reduced version of the table with only the first two phylo-approaches
        final_boyce_results_reduced=final_boyce_results[,which(
            colnames(final_boyce_results)=="partition" | 
            grepl("phylo_rasters_subset_inter_", colnames(final_boyce_results), fixed=TRUE) | 
            grepl("phylo_rasters_proportion_subset_inter_", colnames(final_boyce_results), fixed=TRUE))]
        #reorder columns by algorithm
        final_boyce_results_reduced=final_boyce_results_reduced[,c(
            which(colnames(final_boyce_results_reduced)=="partition"),
            which(grepl("glm", colnames(final_boyce_results_reduced), fixed=TRUE)),
            which(grepl("gam", colnames(final_boyce_results_reduced), fixed=TRUE)),
            which(grepl("rf", colnames(final_boyce_results_reduced), fixed=TRUE)))]
        #make a version with phylo-approaches as rows and partitions as columns
        final_boyce_results_transponse=as.data.frame(t(final_boyce_results[,which(colnames(final_boyce_results)!="partition")]))
        colnames(final_boyce_results_transponse)=final_boyce_results$partition
        final_boyce_results_transponse$model=row.names(final_boyce_results_transponse)
        row.names(final_boyce_results_transponse)=NULL
        final_boyce_results_transponse = final_boyce_results_transponse[, c(
            which(colnames(final_boyce_results_transponse)=="model"),
            which(colnames(final_boyce_results_transponse)!="model"))]

        #save median boyce
        write.table(final_boyce_results_transponse, gzfile(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/", species, "_boyce_table.tsv.gz", sep="")), sep="\t", col.names=TRUE, row.names=FALSE)
        write.table(final_boyce_results_reduced, gzfile(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/", species, "_boyce_table_reduced.tsv.gz", sep="")), sep="\t", col.names=TRUE, row.names=FALSE)

        #calculate the difference with non-phylo
        #load non-phylo boyce
        non_phylo_boyce=read.table(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/", species, "_boyce_table.tsv.gz", sep=""), sep="\t", header=TRUE)
        #merge phylo and non-phylo
        merged_boyce=merge(non_phylo_boyce, final_boyce_results, by="partition", sort=FALSE)
            #avoid sort because we have already sorted the rows using the partition and follwoing natural sorting, 1,2,3....
        #remove percentile rows
        merged_boyce=merged_boyce[which(!merged_boyce$partition %in% c("percentile_2.5", "percentile_50", "percentile_97.5")),]
        if(nrow(merged_boyce)!=12){
            stop("ERROR! FALSE! WE HAVE A PROBLE WHEN MERGING THE BOYCE TABLES WITH PHYLO AND NON-PHYLO")
        }
        #calculate the difference phylo vs non-phylo for each phylo-approach and algorithm
        #selected_phylo_model=phylo_models[1]
        for(selected_phylo_model in phylo_models){

            #algorithm="glm"
            for(algorithm in c("glm", "gam", "rf")){

                #extract boyce values without phylo
                boyce_without_phylo=merged_boyce[,which(colnames(merged_boyce)==paste(algorithm, "_eval", sep=""))]

                #extract boyce values with phylo
                boyce_with_phylo=merged_boyce[,which(colnames(merged_boyce)==paste(selected_phylo_model, "_", algorithm, "_boyce", sep=""))]

                #calculate the difference
                phylo_vs_non_phylo=boyce_without_phylo-boyce_with_phylo
                
                #save the difference
                merged_boyce[, paste(selected_phylo_model, "_", algorithm, "_diff", sep="")] = phylo_vs_non_phylo
            }
        }

        #check we have correctly calculated the difference between phylo and non-phylo
        #algorithm="glm"
        for(algorithm in c("glm", "gam", "rf")){

            #extract the boyce index for non-phylo
            non_phylo_selected=non_phylo_boyce[
                which((!non_phylo_boyce$partition %in% c("percentile_2.5", "percentile_50", "percentile_97.5")) & (non_phylo_boyce$partition %in% final_boyce_results$partition)),
                which(colnames(non_phylo_boyce)==paste(algorithm, "_eval", sep=""))]
                #avoid percentiles and select those partitions present in phylo boyce

            #extract the boyce index for phylo
            phylo_selected=final_boyce_results[
                which(!final_boyce_results$partition %in% c("percentile_2.5", "percentile_50", "percentile_97.5")),
                which(grepl(paste(algorithm, sep=""), colnames(final_boyce_results), fixed=TRUE))]
                #select all phylo approaches for the selected algorithm

            #extract the difference non-phylo vs phylo previously calculated
            diff_selected=merged_boyce[
                which(!merged_boyce$partition %in% c("percentile_2.5", "percentile_50", "percentile_97.5")),
                which(grepl(paste(algorithm, "_diff", sep=""), colnames(merged_boyce), fixed=TRUE))]
                #select all phylo approaches for the selected algorithm

            #calculate again the difference and check
            if(FALSE %in% c(diff_selected==non_phylo_selected-phylo_selected)){
                stop("ERROR! FALSE! WE HAVE NOT PROPERLY CALCULATED THE DIFFERENCE BETWEEN PHYLO AND NON-PHYLO")
            }
        }

        #calculate percentiles
        low_interval_phylo_diff=apply(X=merged_boyce[,which(colnames(merged_boyce) != "partition")], MARGIN=2, FUN=quantile, probs=0.025, na.rm=TRUE)
        median_phylo_diff=apply(X=merged_boyce[,which(colnames(merged_boyce) != "partition")], MARGIN=2, FUN=quantile, probs=0.5, na.rm=TRUE)
        high_interval_phylo_diff=apply(X=merged_boyce[,which(colnames(merged_boyce) != "partition")], MARGIN=2, FUN=quantile, probs=0.975, na.rm=TRUE)
            #calculate the median of each column (excluding the first one with the number of the partition)
            #margin lets you select if you want to apply the function across columns (2) or rows (1)
            #we can have NAs for some phylo approaches (see above), so remove NAs to avoid problems

        #add them as new rows. Do it in two steps to avoid converting numeric into string
        merged_boyce[nrow(merged_boyce) + 1, 1]="percentile_2.5"
        merged_boyce[nrow(merged_boyce), 2:ncol(merged_boyce)]=low_interval_phylo_diff
        merged_boyce[nrow(merged_boyce) + 1, 1]="percentile_50"
        merged_boyce[nrow(merged_boyce), 2:ncol(merged_boyce)]=median_phylo_diff
        merged_boyce[nrow(merged_boyce) + 1, 1]="percentile_97.5"
        merged_boyce[nrow(merged_boyce), 2:ncol(merged_boyce)]=high_interval_phylo_diff
            #https://stackoverflow.com/a/44150746

        #put the phylo-models as rows and partitions as columns
        merged_boyce_transponse=as.data.frame(t(merged_boyce[,which(colnames(merged_boyce)!="partition")]))
        colnames(merged_boyce_transponse)=final_boyce_results$partition
        merged_boyce_transponse$model=row.names(merged_boyce_transponse)
        row.names(merged_boyce_transponse)=NULL
        merged_boyce_transponse=merged_boyce_transponse[, c(
            which(colnames(merged_boyce_transponse)=="model"),
            which(colnames(merged_boyce_transponse)!="model"))]

        #save
        write.table(merged_boyce_transponse, gzfile(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/", species, "_boyce_table_non_phylo_vs_phylo.tsv.gz", sep="")), sep="\t", col.names=TRUE, row.names=FALSE)

        #compress and remove folders
        print(system(paste(" \\
            cd ./results/global_test_phylo_current/predict_eval_phylo/", species, "/; 
            tar \\
                --verbose \\
                --create ./phylo_stacks \\
                --file ./", species, "_phylo_stacks.tar.gz \\
                --gzip; \\
            rm -rf ./phylo_stacks", sep=""), intern=TRUE))
            #verbose: list files being compressed
            #create: create an archive from the files in XXX
            #file: store the output as a file named as XXX
            #use gzip to compres
                #https://unix.stackexchange.com/a/93158
        print(system(paste(" \\
            cd ./results/global_test_phylo_current/predict_eval_phylo/", species, "/; 
            tar \\
                --verbose \\
                --create ./prediction_rasters_phylo \\
                --file ./", species, "_prediction_rasters_phylo.tar.gz \\
                --gzip; \\
            rm -rf ./prediction_rasters_phylo", sep=""), intern=TRUE))
        system(paste(" \\
            cd ./results/global_test_phylo_current/predict_eval_phylo/", species, "/; 
            rm -rf ./prediction_rasters", sep=""))
    } else {
        #if we cannot do stuff, remove the folder of the species
        system(paste("rm -rf ./results/global_test_phylo_current/predict_eval_phylo/", species, "/", sep=""))
    }

    #ending
    print(paste("ENDING predict_eval_phylo FOR ", species), sep="")
}

#run it for one species
#predict_eval_phylo("radiata")
  



##########################
##### RUN EVERYTHING #####
##########################

#make folder to save outputs of species
system("mkdir -p ./results/global_test_phylo_current/species_output_files")

##prepare master function
#better this way, so if a species is too slow it will not get the other species stuck in the first step, as the other ones will continue running
master_processor=function(species){

    #send output to a specific file for the species
    output_file=file(paste("./results/global_test_phylo_current/species_output_files/", species, ".txt", sep=""), open = "wt")
    sink(output_file, type="output")
    sink(output_file, type="message")
        #https://stackoverflow.com/a/75991645

    #start
    print(paste("STARTING ", species), sep="")

    #check species
    if(!species %in% unique(naturalized_occurrences$species)){
        stop(paste("ERROR! FALSE! WE HAVE USED A WRONG SPECIES NAME: ", species, sep=""))
    }
        #if we put here before creating the output file, the stop message is sent to the output file of other species that are running right now

    #occurrences preparation
    output_first_step=exsitu_occurrences(species)
    status_step_1=output_first_step[[1]]
    n_points_before_resampling=output_first_step[[2]]

    #predict and evaluate without phylo
    status_step_2=predict_eval_no_phylo(species, status_step_1)

    #predict and evaluate with phylo
    predict_eval_phylo(species, status_step_2)

    #finish
    print(paste("ENDING ", species), sep="")

    #return the number of points before resampling
    return(n_points_before_resampling)

    #stop writing to the file
    sink(type="message")
    sink(type="output")
    close(output_file)
        #https://stackoverflow.com/a/75991645
}


##parallelize the function
#load packages
require(foreach)
require(doParallel)

#set up cluster
clust <- makeCluster(length(species_to_analyze), outfile="")
    #You can usually figure out why the worker died by using the makeCluster "outfile" option so that the error message generated by the worker isn't thrown away. I usually recommend using outfile=""
        #https://stackoverflow.com/a/24352032
registerDoParallel(clust)

#packages
packages_parallel=c("raster", "sf", "randomForest", "gam", "modEvA", "gtools")

#run the function in parallel
n_points_before_resampling = foreach(i=species_to_analyze, .packages=packages_parallel) %dopar% {
    master_processor(species=i)
}

#stop the cluster 
stopCluster(clust)


##create table with number of occurences before resampling
#convert to DF the list with the number of data.points
n_points_before_resampling_df=do.call(rbind.data.frame, n_points_before_resampling)
    #do.call pass elements of your_list as arguments to rbind. It's equivalent of rbind(your_list[[1]], your_list[[2]], your_list[[3]], ....., your_list[[length of your_list]])
        #https://stackoverflow.com/a/4227483

#add column with species
n_points_before_resampling_df$species=species_to_analyze

#put species column first
n_points_before_resampling_df=n_points_before_resampling_df[,c(
    which(colnames(n_points_before_resampling_df)=="species"), 
    which(colnames(n_points_before_resampling_df)!="species"))]
print(n_points_before_resampling_df)

#check
if(FALSE %in% (n_points_before_resampling_df$total_points==n_points_before_resampling_df$points_outside+n_points_before_resampling_df$points_inside)){
    stop("ERROR! FALSE! WE HAVE A PROBLEM CALCULATING THE NUMBER OF OCCURENCES INSIDE AND OUTSIDE THE PA BUFFER")
}

#save
system("mkdir -p ./results/global_test_phylo_current/n_points_before_resampling/")
write.table(n_points_before_resampling_df, paste("./results/global_test_phylo_current/n_points_before_resampling/n_points_before_resampling_", batch_number, ".tsv", sep=""), sep="\t", row.names=FALSE)




##################
##### FINISH #####
##################

#finish the script
print("## FINISH ##")

#singularity exec 01_global_test_phylo_ubuntu_20_04_v1.sif ./code/phylo/global_test_phylo_current_v1.R --species="halepensis,radiata" --batch="batch_1" 2>&1 ./code/phylo/global_test_phylo_current_v1_batch_1.Rout




######################
##### NEXT STEPS #####
######################

#took old version of _check_altitudinal_sampling.png in results/occurrences folders for radiata

#check the thing about glm non binary!
    #mira pagina 284 Nivk's book, he seems to do as we




#check in github fast the differences between the version finished before debugging and the last version
##SEND AGAIN NECESARRY FILES




#script to check general output and outputs per species
    #table with n_points_before_resampling_df
        #CHECK SPECIES WITH LOW OCCURRENCES (e.g., clausa)
            #check pattern, these cases have higher or lower boyce? to see if they are biasing our results, our glmm, see below.
        #check how many naturalized presences inside PA buffer across species in "n_points_before_resampling". Use this to answer comment 8 of review about buffer size too big.


#after running everything, quick check how other phylo-approaches work in the species with more impact, i.e., radiata, sylvestris, strobus...
    #I have already done it with the frist run and subset_no_inter is a bit better, but probably not enough to change the whole phylo approach of the main text
    #If the other phylo-approaches work much better, you have several options
        #split perret dataset (maybe 75-25), compare the different phylo approaches in the 75% dataset, select the best one and then evaluate in the 25% dataset
            #if you just select the best option performing in the whole Perret dataset, then you are not testing independently, and there is risk of overfitting.
        #combine all phylo approaches into a single ensemble and use it as correction
        #remember, you cannot just chcedk what is the best phylo approach in current distirbutions because everyhing close to current will be 1 in the correction, so the phylo correction is not going to have a great impact there.
