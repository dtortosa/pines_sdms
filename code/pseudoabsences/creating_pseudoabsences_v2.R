#code for create pseudoabsences for all species using a function parallelizable

####define work directory
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus")

#Libraries
require(raster)
require(rgeos) #for creating buffers
require(foreach) #for repeat a process several times 
require(doParallel) #for parallel

#set the seed to have reproducible results
set.seed(36743)

###list ocurrences
list_species = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/species.txt", sep="\t", header=T)
list_species

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

##load variables selected for making the pseudo-absences
bio6 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio6.asc") #mininum tempreature of the coldest month (Bio6) 
bio18 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio18.asc") #the moisture index of the warmest quarter (bio18)
clay = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/clay.asc") #the percentage of clay. We will use this variable, for drop the areas without soil data for the environmental categories raster. In this way, PAs will not be created in areas without soil data. 
bio15 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio15.asc")#initial variable candidate that was removed

#Note: On the contrary to presences, in this case there is not points outside of area of environmental variables (areas without data) because these points are created with environmental sampling of the spaces using these variables with the same resolution that resolution used in the analysis. All pseudoabsences have values of environmental variables. 

####create environmental categoires around all the world#######
require(ecospat)

#mask climatic variables with soil for avoid sea areas
bio6 = mask(bio6, clay)
bio8 = mask(bio18, clay)

#calculate the minimun values of bio18 for avoiding negative values for the log. 
values_bio18 = getValues(bio18) #extract the values of bio18
global_sum_to_bio18 = abs(floor(min(na.omit(values_bio18)))) #select min without na and inf

#convert inf and -inf in NA for avoiding problems in log
bio18[is.infinite(bio18)] <- NA

##create 9 caterogires for our two environmental variables
##2 or 3 variables, only, temperature and moisture variables.. ad seasonality. 
global_variables_categories = ecospat.rcls.grd(bio6,9)+ecospat.rcls.grd(log(bio18+global_sum_to_bio18),9)*10 #ecospat.rcls.grd take a raster and divide it in several categories. We do it for all variables and sum them. We apply the log to bio18 for the non-normal distribution of the data, most of the values are included in some areas, and the rest of values are in very litle areas, which would have a the same number point than the largers. More info in the note_book. 
#hist(getValues(bio6))
#hist(getValues(bio18))
#hist(getValues(log(bio18+global_sum_to_bio18)))
#plot(getValues(bio6), getValues(bio18))
global_variables_categories
global_strata = table(getValues(global_variables_categories)) #the values of the raster are several categories. 
length(global_strata)

##plot
##create an adequate palette of colors
global_cspan<-maxValue(global_variables_categories)-minValue(global_variables_categories) #difference between the lowest and the highest category for calcualte the spectrum of rainbow colors
global_yb<-rainbow(100)[round(runif(global_cspan,0.5,100.5))] #take 80 (number of cspan) numbers from uniform distribution, without decimals (round) and use them for select colors from rainbow. In this way, we can see as different the categories.
plot(global_variables_categories, col=global_yb, asp=1)

###################################################
#####TEst different environmental categoires#######
###################################################
###create diffent environmental categories for calculate the optinmal numver of variables used. DON`T RUN. 

#bio6 - bio18
variables_categories_0 = ecospat.rcls.grd(bio6,9)+ecospat.rcls.grd(bio18,9)*10
variables_categories_0 = mask(variables_categories_0, clay)
cspan_0<-maxValue(variables_categories_0)-minValue(variables_categories_0)
yb_0<-rainbow(100)[round(runif(cspan_0,.5,100.5))]

#bio6 - bio18 - bio15
variables_categories_1 = ecospat.rcls.grd(bio6,9)+ecospat.rcls.grd(bio18,9)+ecospat.rcls.grd(bio15,9)*10
variables_categories_1 = mask(variables_categories_1, clay)
cspan_1<-maxValue(variables_categories_1)-minValue(variables_categories_1)
yb_1<-rainbow(100)[round(runif(cspan_1,.5,100.5))]

#bio6 - bio18 - clay
variables_categories_2 = ecospat.rcls.grd(bio6,9)+ecospat.rcls.grd(bio18,9)+ecospat.rcls.grd(clay,9)*10
cspan_2<-maxValue(variables_categories_2)-minValue(variables_categories_2)
yb_2<-rainbow(100)[round(runif(cspan_1,.5,100.5))]

#bio6 - bio18 - bio15 - clay
variables_categories_3 = ecospat.rcls.grd(bio6,9)+ecospat.rcls.grd(bio18,9)+ecospat.rcls.grd(bio15,9)+ecospat.rcls.grd(clay,9)*10
cspan_3<-maxValue(variables_categories_3)-minValue(variables_categories_3)
yb_3<-rainbow(100)[round(runif(cspan_3,.5,100.5))]

#plot everything
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences/environm_categor_with_different_variables.png", width=1300, height=1300)
par(mfcol=c(2,2))
plot(variables_categories_0,col=yb_0,main="bio6 - bio18 ", asp=1)
plot(variables_categories_1,col=yb_1,main="bio6 - bio18 - bio15", asp=1)
plot(variables_categories_2,col=yb_2,main="bio6 - bio18 - clay", asp=1)
plot(variables_categories_3,col=yb_3,main="bio6 - bio18 - bio15 - clay", asp=1)
dev.off()#the best option is bio6 and bio18

save.image("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/rdata/environm_categor_with_different_variables.RData")


############################
#######Function PAs#########
############################
pseudo_absence = function(i){
    
    #print the specific_epithet
    print(i)

    #load ocurrences
    data = read.csv(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences", paste(i, "final.presences.csv", sep="_"), sep="/"), header=TRUE)

    #load distribution
    distribution = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p", paste(i, "01.img", sep="_"), sep="_"))


    ################################
    ######polygon_PA_buffer#########
    ################################

    #read raster distribution buffer
    raster_buffer = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences", paste(i, "distribution_buffer.asc", sep="_"), sep="/"))

    #create a polygon from that raster
    polygon_buffer = rasterToPolygons(raster_buffer, fun=function(x){x==1}, dissolve=TRUE)
    
    ###Create the pseudoabsence buffer for crop the environmental variables. The difference between this buffer and the distribution_buffer will be the area with PAs.
    if(i == "pumila"){#if the species is pumila we load the PA buffer previously created with buffer in both sides of the map
        #load the polygon with polygons in both sides of the map previously created
        raster_PA_buffer = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences/pumila_PA_buffer.asc")

        #create the new polygon halepensis distribution without sea areas
        polygon_PA_buffer = rasterToPolygons(raster_PA_buffer, fun=function(x){x==1}, dissolve=TRUE)

        #crop the distribution raster with the polygon PA buffer
        distribution_crop =  crop(distribution, polygon_PA_buffer)         
    } else {#if not, we create the buffer here
        #create a polygon buffer around the distribution buffer
        polygon_PA_buffer = gBuffer(polygon_buffer, byid=FALSE, id=NULL, width=22.5, quadsegs=5, capStyle="ROUND", joinStyle="ROUND", mitreLimit=1.0)

        #crop the distribution raster with the polygon PA buffer
        distribution_crop =  crop(distribution, polygon_PA_buffer) 

        #create a raster from the PA buffer polygon 
        raster_PA_buffer = raster() 
        extent(raster_PA_buffer) = extent(distribution_crop) 
        res(raster_PA_buffer) = res(distribution) 
        raster_PA_buffer  = rasterize(polygon_PA_buffer, raster_PA_buffer) 
        #plot(crop(bio6, raster_PA_buffer), col="gray80")
        #plot(raster_PA_buffer, add=T)
        #plot(polygon_PA_buffer, add=T)
        #plot(polygon_buffer, add=T)

        #drop the sea areas
        raster_PA_buffer = distribution_crop*raster_PA_buffer 
        raster_PA_buffer[!is.na(raster_PA_buffer)] <- 1 
        #plot(crop(bio6, raster_PA_buffer), col="gray80")
        #plot(raster_PA_buffer, add=T, col="yellow")
        #plot(raster_buffer, add=T, col="blue")
        #legend(10, 57, legend=c("Pseudo-Absences buffer", "Presences buffer"), fill=c("yellow", "blue"))
        #We mantenain some water bodies inside the continents because we can't multiply our soil variables with our PA buffer, they have different resolution. Because of this, we will use part of vaule of environmental variables over water bodies for variables selection, I think it is not very important.
           
        #create the new polygon halepensis distribution without sea areas
        polygon_PA_buffer = rasterToPolygons(raster_PA_buffer, fun=function(x){x==1}, dissolve=TRUE)

        #write the PA buffer without water bodies.  
        writeRaster(raster_PA_buffer, paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences", paste(i, "PA_buffer.asc", sep="_"), sep="/"), format="ascii", overwrite=TRUE)              
    }
    
    #create a raster from polygon_PA_buffer with the same extent and resolution than climatic variables for masking them. Crop let a lot of space outside the buffer, mask remove all areas outside the buffer. Mask only can be applied between rasters with the same resolution and extent
    raster_PA_buffer_mask = raster() 
    extent(raster_PA_buffer_mask) = extent(bio6) 
    res(raster_PA_buffer_mask) = res(bio6) 
    raster_PA_buffer_mask  = rasterize(polygon_PA_buffer, raster_PA_buffer_mask) 

    #crop  environmental variables with the buffer of PA
    bio6_crop =  mask(bio6, raster_PA_buffer_mask)
    bio18_crop =  mask(bio18, raster_PA_buffer_mask) 
    clay_crop = mask(clay, raster_PA_buffer_mask) #we use the raster of the PA buffer to crop the rest of variables because y faster than using a polygon. 

    #mask the climatic raster for dropping areas without soil data
    bio6_crop = mask(bio6_crop, clay_crop)
    bio18_crop = mask(bio18_crop, clay_crop)

 
    ################################
    ######Create categories#########
    ################################
    
    #calculate the minimun values of bio18 for avoiding negative values for the log. 
    sum_to_bio18 = NULL
    if (bio18_crop@data@min==0){ #if the minimun vale is 0
        sum_to_bio18 = 1 #we will sum 1 to bio18
    } else { #if not
        if (bio18_crop@data@min<0){ #if the minimun value of bio18 is negative
            sum_to_bio18 = abs(floor(bio18_crop@data@min)) #round the number to the lower limite (because of this floor; eg -6.7, we take -7) and calculate the absolute value. 
        } else {
            sum_to_bio18 = 0
        }
    } #we will sum this value tu bio18_crop when we make the log of it. 

    ##create 9 caterogires for our two environmental variables
    ##2 or 3 variables, only, temperature and moisture variables.. ad seasonality. 
    variables_categories = ecospat.rcls.grd(bio6_crop,9)+ecospat.rcls.grd(log(bio18_crop+sum_to_bio18),9)*10 #ecospat.rcls.grd take a raster and divide it in several categories. We do it for all variables and sum them. We apply the log to bio18 for the non-normal distribution of the data, most of the values are included in some areas, and the rest of values are in very litle areas, which would have a the same number point than the largers. More info in the note_book and the beginning of this script. 
    #the number of clases (9) and the multplication was established by men in switezerlan, probably under the supervision of nick, becuase the help of ecospat.rcls.grd (made by nick and achilleas) uses that.
    variables_categories #the result map is a raster with a value for each strata. So all cells of a given strata have the same value. Therefore, the number more frequent across the raster will correspond to a bigger strata.
    strata = table(getValues(variables_categories)) #the values of the raster are several categories. 
    length(strata)

    ##create an adequate palette of colors
    cspan<-maxValue(variables_categories)-minValue(variables_categories) #difference between the lowest and the highest category for calcualte the spectrum of rainbow colors
    yb<-rainbow(100)[round(runif(cspan,.5,100.5))] #take 80 (number of cspan) numbers from uniform distribution, without decimals (round) and use them for select colors from rainbow. In this way, we can see as different the categories.

    #plot result
    png(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences", paste(i, "environmental_categories.png", sep="_"),sep="/"), width=1300, height=1300)
    par(mfcol=c(2,1))
    hist(variables_categories, breaks=100, col=heat.colors(cspan), main="Histogram values")
    plot(variables_categories,col=yb,main="Stratified map", asp=1)
    dev.off()

    ################################
    ######Make sampling#############
    ################################

    #calculate number PAs
    if(10*nrow(data)/length(strata)<30){ #if for each class there is less of 30 points with a number of PAs of (10*ocurrences): 

        number_PAs = length(strata)*30 #with this number of PAs we have 30 PAs for each strata. 

    } else { 
        number_PAs = 10*nrow(data)
    }

    #make the stratifief sampling in the enviromental categories: 
    envstrat_ecu = ecospat.recstrat_regl(in_grid = variables_categories, sample_no = number_PAs)
    envstrat_prp = ecospat.recstrat_prop(in_grid = variables_categories, sample_no = number_PAs) 

    #compare both types of methods
    png(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences", paste(i, "pseudo_absences_types_samplings_without_resampling.png", sep="_"),sep="/"), width=2000, height=2000)
    par(mfcol=c(2,2))
    plot(variables_categories, col=yb, main="Equal Sampling")
    points(envstrat_ecu$x,envstrat_ecu$y,pch=16,cex=.2,col="black")
    plot(variables_categories, col=yb, main="Proportional Sampling")
    points(envstrat_prp$x,envstrat_prp$y,pch=16,cex=.2,col="black")
    barplot(table(envstrat_ecu$class),col="slategray4", main="Equal point allocation")
    barplot(table(envstrat_prp$class),col="firebrick", main="Proportional point allocation")
    dev.off() #we will select the proportional, just in case there is categories very very bign and other very very small. In that case, with the regular sampling we will introduce the same number of points in the small and in the big, but this number it is not enought for cover all conditions inside the big cathegory. 
        #Nick says: 
            #No, each stratum has the same environmental range… this imposed by the stratification. But not all strata are geographically as abundant…
        #Lo que dice nick tiene sentido, un strato se supone que es climáticamente homogeneous. Lo de muestrear mas lo más grandes aun así tiene sentido, porque uno pequeño lo vas a cubir con menos ausencias. Esto no cambia nada los análisis y el texto Nick ha modificado y queda bien. Los geograficamente más amplios se muestrean más.

    ##select points inside PA buffer
    PA_buffer_values_of_points = extract(raster_PA_buffer, envstrat_prp[c("x", "y")]) #extract the values of each PAs in the raster PA buffer, the possibilities are two: 1 or NA. 1 inside PA buffer, Na outside PA buffer. We use the raster instead of the buffer because it is faster. 
    #head(PA_buffer_values_of_points)
    #length(PA_buffer_values_of_points)
    #nrow(envstrat_prp) 
    points_and_values_PA_buffer = cbind(envstrat_prp[c("x", "y")], PA_buffer_values_of_points) #bind these values with the corresponding PA points
    #str(points_and_values_PA_buffer)
    points_inside_PA_buffer = which(!is.na(points_and_values_PA_buffer$PA_buffer_values_of_points)) #obtain rows in which values of the PA buffer raster is 1, thus points inside the buffer
    points_outside_PA_buffer = which(is.na(points_and_values_PA_buffer$PA_buffer_values_of_points)) #obtain rows in which values of the PA buffer raster is NA, thus points outside the buffer

    ##Select points outside distribution buffer
    distribution_buffer_values_of_points = extract(raster_buffer, envstrat_prp[c("x", "y")]) #extract the values of each PAs in raster of buffer distribution, the possibilities are two: 1 or NA. 1 inside PA buffer, Na outside PA buffer. We use the raster instead of the buffer because it is faster. 
    #head(distribution_buffer_values_of_points)
    #length(distribution_buffer_values_of_points)
    #nrow(envstrat_prp) #equal length than points created in envstrat_prp, some of the points fall outside distribution buffer
    points_and_values_distribution_buffer = cbind(envstrat_prp[c("x", "y")], distribution_buffer_values_of_points) #bind these values with the corresponding PA points
    points_inside_distribution_buffer = which(!is.na(points_and_values_distribution_buffer$distribution_buffer_values_of_points)) #obtain rows in which values of the PA buffer raster is 1, thus points outside the buffer
    points_outside_distribution_buffer = which(is.na(points_and_values_distribution_buffer$distribution_buffer_values_of_points)) #obtain rows in which values of the PA buffer raster is NA, thus points inside the buffer

    #plot
    png(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences", paste(i, "pseudo_absences_inOut_buffer_without_resampling.png", sep="_"),sep="/"), width=1300, height=1300)
    plot(clay_crop, col="gray80")
    plot(polygon_PA_buffer, col="green", add=T)
    plot(polygon_buffer, col="red", add=T)
    points(envstrat_prp[c(points_outside_PA_buffer),]$x, envstrat_prp[c(points_outside_PA_buffer),]$y, pch=16,cex=.1,col="yellow") #points outside the PA buffer
    points(envstrat_prp[points_inside_distribution_buffer,]$x, envstrat_prp[points_inside_distribution_buffer,]$y,pch=16,cex=.1,col="blue") #points inside distribution buffer
    points(envstrat_prp[-c(points_outside_PA_buffer, points_inside_distribution_buffer),]$x,envstrat_prp[-c(points_outside_PA_buffer, points_inside_distribution_buffer),]$y, pch=16,cex=.1,col="black") #these points are the interesting for me. 
    legend("topright", legend=c(expression(bold("PA buffer")), expression(bold("Distribution buffer")), expression(bold("PAs outside PA buffer")), expression(bold("PAs inside distribution buffer")), expression(bold("PAs inside PA buffer and outside distribution buffer"))), fill=c("green", "red", "yellow", "blue", "black"), cex=1)
    dev.off()

    #create PA data frame with the same variables of presence
    PA = as.data.frame(envstrat_prp[-c(points_outside_PA_buffer, points_inside_distribution_buffer),]) #create data.frame from the coordinates of created points inside of the PA buffer (-points_outside_PA_buffer) and outside of distributon buffer (-points_inside_distribution_buffer).
    PA$precision_weight = rep(x=sum(data$precision_weight)/number_PAs, #repeat the sum of all weights of presences divided by the number of presences
        times=nrow(PA)) #as many times as pseudo-absences exists. This will be the weigh of each pseudo-absence. 
    PA$presence = rep(x=0, times=nrow(PA)) #repeat 0 (absence) as many times as pseudo-abscences exists. 
    names(PA)[1] <- "longitude"
    names(PA)[2] <- "latitude"
    #str(PA)
    #head(PA)

    #####PROBLEMA: porngo el precision weight en base al 10*nº ocurrencias porque eso es lo que conseguiremos al final, después del segundo muestro. PERO OJO: ten en cuenta que los datos se parten en dos para entrenamiento y evaluación. Por lo que si queremos que sea igualm para RF y GLM-GAM, tenemos que dividir el weight de las ocurrencias finales de cada set de entrenamiento entre el número de PAs en ese SET. Eso se ha solucionado en la version 5 de modelado de las condiciones actuales. En esa versión, además hemos indicado que todas las presencias de baja precision tendrán un weight de 1 si NO hay presencias de alta precisión. No tiene sentido darle menos peso cuando es lo unico que hay, y además así es como se hacíamos de forma manual para random forest con el sampleo.

    ######################################################################################################
    #######Reduce the number of points to selected number of PAs (number_PAs created before)##############
    ######################################################################################################

    ###Explanation about this step (extract from a mail to Rafi): I have explained the problem bad. Of course, we can obtain a proportional sampling in relation to the size of the strata from ecospat. The problem is that we cannot stablish the exact number of PAs, we can set 3000 PAs in the argument "sample_no", but If ecospat.recstrat_prop consider this number as too low (according to number of strata), it increases the number of PAs. I have seen increments of two times the number. Therefore, we don't have a fine control of the PAs numbers. For reducing number, you proposed to me make a stratified sampling on the PAs, which have been created with a proportional sampling previously with ecospat. I have sampled with more probability those PAs included in strata with higher number of PAs (larger strata according to proportional sampling of ecospat). It is to say, larger strata are more sampled than smaller strata. In this way, we can select the exact number of PAs that we want without lose proportional sampling made by ecospat. I have checked the result of this for several species and it is what we want: A reduction in the number of PAs, but keeping the proportional sampling of the strata (more points in larger strata). You can see the result of this in the plots that I sent you previously: Per each species, "without resampling" show the PAs created with ecospat without any change, whilst "with resampling" show the same PAs but after make a stratified sampling with the purpose of reducing number of PAs. Last, I have made a change, I have introduced a conditional in the selection of number of PAs. In the case of species with a low number of ocurrences, 10*(presences) is too low. You can see this in the attached plot of P. cubensis ("cubensis_PA_low.png"), the number of PAs is exactly 10*(presences), i think that the space is not sampled completely. Because of this, I have increased the number of PAs in these cases, check the result in ("cubensis_PA_high.png"). The condition minimun 30 points per class ((10*presences)/(number classes) higher than 30) is a way to determine if the number of PAs is enough in relation to the number of strata, of course the propotional sampling of ecospat will cause that small classes have less than 30 points, but I will be sure that the large classes are well sampled. 

    #the number of PAs for proportional sampling is higher than number_PAs, very higher. This is caused by ecospat.recstrat_prop function, that does not create the exact number of PAs that you indicate. 
    nrow(PA) == number_PAs

    ####Reduction the number of PAs of the proportional samplings###
    #we will make a sampling of the total of PAs of the proportional sampling, but only on the strata with a number of points higher than the mean of the points between all the classes

    #calculate the number of PAs in each strata (environmental categorie)
    table_classes = table(PA$class)
    classes = data.frame(class=names(table_classes), prob = as.vector(table_classes)) #give names to columns. We will use the number of PAs in each strata as probability of the sampling. We want that the points of larger stratas (categories with more points thanks to the proportional sampling of the space) have higher probability to be taken. In this way, the areas larger (higher number of points) will be more sampled, which was the objective to use the proportional sampling, but now we will have exactly 10*nrow(data) pseudoabsences. 
    summary(classes$class==names(table_classes)) #the names of classes between files are the same
    summary(classes$prob==table_classes) #the probs are the same in the corresponding categories

    #create a column with probablities as the number of points of the correspondin category. 
    PA$probs = NA #create a column empty in the data frame of PAs
    envstrat_prp_fixed = data.frame() #create a data frame empty

    for (k in unique(PA$class)){ #for each class (strata)
        subset = PA[PA$class==k,] #select all the PAs (rows) of this strata. 
        subset$probs = classes[classes$class==k,]$prob #the probs of these rows will be the number of PAs in the corresponding strata taken from classes data.frame
        envstrat_prp_fixed = rbind(envstrat_prp_fixed, subset) #bind these rows with the rows of the new data.frame.
    }

    test=NULL #test if the loops works 
    for (k in unique(PA$class)){ #for each class
        test = append(test, envstrat_prp_fixed[envstrat_prp_fixed$class==k,]$probs == classes[classes$class==k,]$prob) #the probs of the class [i] in envstrat_prp_fixed are the same than the probs of this class in classes data.frame. 
    }
    summary(test)

    #Create vector with the name of the clases with more and less PAs than the mean of the all classes.
    mean_points_per_class = as.vector(quantile(table(PA$class), probs=0.25)) #first quartile of the all classes
    small_classes = classes[classes$prob<=mean_points_per_class,]$class #classes wiht  a number of points lower than the mean of the classes
    big_classes = classes[classes$prob>mean_points_per_class,]$class #classes wiht  a number of points higher than the mean of the classes
    summary(small_classes %in% big_classes) #there are not coincidences. Must be FALSE. 
    length(small_classes) + length(big_classes) == nrow(classes) #the total of classes are the same than the sum of small and big classes

    #Separate envstrat_prp_fixed_big_classes into two big and small strata
    envstrat_prp_fixed_big_classes = envstrat_prp_fixed[envstrat_prp_fixed$class %in% big_classes,]
    envstrat_prp_fixed_small_classes = envstrat_prp_fixed[envstrat_prp_fixed$class %in% small_classes,]
    summary(envstrat_prp_fixed_big_classes$class %in% big_classes)
    summary(envstrat_prp_fixed_small_classes$class %in% small_classes) #the separation works

    #make the sampling of the large strata 
    sampling = envstrat_prp_fixed_big_classes[sample(1:nrow(envstrat_prp_fixed_big_classes), size=number_PAs-nrow(envstrat_prp_fixed_small_classes), prob=envstrat_prp_fixed_big_classes$prob),] #sample from 1 to the number of rows of large strata, only (number_PAs, frm before) less the number of points in small areas, which will be taken without sampling. The probability of take any points depends of the number of PAs in its strata (envstrat_prp_fixed_big_classes$prob).

    #bind the select points of large strata, with all the points of small strata
    final = rbind(sampling, envstrat_prp_fixed_small_classes)
    nrow(final) == number_PAs #the number of PAs is 10 times the number of presences
    
    #final test
    final_test = NULL
    for (k in names(table_classes)){
        final_test = append(final_test, unique(final[final$class==k,]$probs) == as.vector(table_classes[which(names(table_classes)==k)]))       
    }
    summary(final_test)

    #save the new PAs
    PA_final = final[,-which(names(final)=="probs")]
    nrow(PA_final)

    #plot
    png(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences", paste(i, "pseudo_absences_on_classes.png", sep="_"),sep="/"), width=2000, height=2000)
    par(mfcol=c(1,2))
    plot(variables_categories, col=yb, main="Without resampling")
    points(PA$longitude,PA$latitude,pch=16,cex=.2,col="black")
    plot(variables_categories, col=yb, main="With resampling")
    points(PA_final$longitude,PA_final$latitude,pch=16,cex=.2,col="black")
    dev.off()

    ##merge with presences
    head(PA_final)
    head(data)
    complete.presences = rbind(data, PA_final[,c("longitude", "latitude", "precision_weight", "presence")])
    #head(complete.presences)
    #str(complete.presences)

    #plot all 
    png(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences", paste(i, "PA_presences.png", sep="_"),sep="/"), width=2000, height=2000)
    plot(bio6_crop, col="gray80")
    points(complete.presences[complete.presences$presence==1,]$longitude, complete.presences[complete.presences$presence==1,]$latitude, cex=.2, col="red", pch=16)
    points(complete.presences[complete.presences$presence==0,]$longitude, complete.presences[complete.presences$presence==0,]$latitude, cex=.2, col="blue", pch=16)
    legend("topright", legend=c(expression(bold("Presences")), expression(bold("Pseudo-absences"))), cex=1, pch=16, col=c("red", "blue"))
    dev.off()

    write.csv(x=complete.presences, file=paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences", paste(i, "complete.presences.csv", sep="_"),sep="/"), row.names=FALSE)

    save.image(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/rdata", paste(i, "preparation_pseudo_Absences.RData", sep="_"),sep="/"))

    #delete all documents
    rm("big_classes",                                             
       "bio18_crop" ,                                                 
       "bio6_crop" ,                       "classes",                         
                                           "clay_crop" ,                      
       "cspan"  ,                          "data"   ,                         
       "distribution" ,                    "distribution_crop" ,              
       "envstrat_ecu" ,                    "envstrat_prp" ,                   
       "envstrat_prp"  ,             "envstrat_prp_fixed" ,             
       "envstrat_prp_fixed_big_classes",   "envstrat_prp_fixed_small_classes",
       "final",                            "final_test"  ,                    
                                             
       "mean_points_per_class",            "polygon_buffer" ,                 
       "polygon_PA_buffer",                "raster_buffer"  ,                 
       "raster_PA_buffer",                 "sampling"  ,                      
       "small_classes",                    "subset" ,                         
       "sum_to_bio18" ,                    "table_classes" ,                  
       "test",                             "variables_categories",            
       "yb")
    gc() 
}


#########################
#######Paralellize#######
#########################
#create a vector with all species
species = epithet_species_list

# set up cluster
clust <- makeCluster(2) #Wu use fork because we don't want that the files that are loaded in each loop (rasters of bio18, bio6 and clay), will be copied only one time, in the master cluster, and not in ALL cluster. This let save memorie, but increase processning time if the files are used a lot by the code. In our case, these raster are used only 1 time in each loop. 
registerDoParallel(clust)

# run for all species
foreach(i = species, .packages = c("raster", "ecospat", "rgeos")) %dopar% { 
    pseudo_absence(i = i)
} #the "stringr" package is used for "str_split_fixed" function

#stop the cluster 
stopCluster(clust)
