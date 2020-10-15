####Code for determine if coordinatePrecision and coordinateuncertaintyinmeters are useful to determine data quality.

#########################################
##########Load species data##############
#########################################

########Tree species#########
#species by species using the bianca's species of the tree and gbif function of dismo
list_ocurrences = list.files(path="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/final_species", pattern=".csv", full.names=TRUE)
length(list_ocurrences) #Load the path of each ocurrence file

precision_table_1 = NULL
precision_table_2 = NULL
precision_table_3 = NULL
precision_table_4 = NULL
for (i in 1:length(list_ocurrences)){
    table = read.csv(list_ocurrences[i], header=TRUE)
    if (length(table$coordinatePrecision) > 0 & length(table$coordinateUncertaintyInMeters) > 0 ){ #we include coordinateuncertaintyinmeters and coordinateprecision, just in case these variable have these names in some species. 
        precision_table_1=rbind(precision_table_1,table[,c("lon", "lat","coordinatePrecision", "coordinateUncertaintyInMeters")])
    } else {
        if (length(table$coordinatePrecision) > 0 ){
            precision_table_2=rbind(precision_table_2,table[,c("lon", "lat","coordinatePrecision")])
        } else {
            if (length(table$coordinateUncertaintyInMeters) > 0 ){
                precision_table_3=rbind(precision_table_3,table[,c("lon", "lat","coordinateUncertaintyInMeters")])
            } else {
                precision_table_4 = rbind(precision_table_4,table[,c("lon", "lat")])
            }
        }
    }
} #Note: for all species the notation is in this way coordinatePrecision, nor coordinateprecision. I have revised it. 

#Test: Each categorie carry out with its characteristics?
length(precision_table_1$coordinatePrecision) > 0 & length(precision_table_1$coordinateUncertaintyInMeters) > 0 #must have coordinatePrecision and coordinateUncertaintyInMeters
length(precision_table_2$coordinatePrecision) > 0 & length(precision_table_2$coordinateUncertaintyInMeters) == 0 #must heve only coordinatePrecision
length(precision_table_3$coordinatePrecision) == 0 & length(precision_table_3$coordinateUncertaintyInMeters) > 0 #must have only coordinateUncertaintyInMeters
length(precision_table_4$coordinatePrecision) == 0 & length(precision_table_4$coordinateUncertaintyInMeters) == 0 #must not have any of that variables. 

#create empty columns for that variables in the data frame without them
precision_table_2$coordinateUncertaintyInMeters = rep(NA, nrow(precision_table_2))
precision_table_3$coordinatePrecision = rep(NA, nrow(precision_table_3))
precision_table_4$coordinateUncertaintyInMeters = rep(NA, nrow(precision_table_4))
precision_table_4$coordinatePrecision = rep(NA, nrow(precision_table_4))

#see structure
str(precision_table_1)
str(precision_table_2)
str(precision_table_3)
str(precision_table_4)

#bind all data frame
precision_table = rbind(precision_table_1, precision_table_2, precision_table_3, precision_table_4)
str(precision_table)

#check that all rows has been included 
nrow(precision_table) == nrow(precision_table_1) + nrow(precision_table_2) + nrow(precision_table_3) + nrow(precision_table_4)

#see number of data coordinatePrecision and coordinateUncertaintyInMeters
nrow(precision_table[!is.na(precision_table$coordinatePrecision),]) #105644
nrow(precision_table[!is.na(precision_table$coordinateUncertaintyInMeters),]) #210626

####################################################################
########Calculate number of decimals of coordinates#################
####################################################################
require(DescTools)
precision_table$lon[1]
Ndec(as.character(precision_table$lon[1])) #we have to convert into character the integers, because these intergers are rounded, the new decimals that appears after covnert to character are not invented. For example, the first ocurrence, whici is the first ocurrence of the first species (P. abicaulis) has a longitude of -118.55102, however you can see that precision_table$lon[1] show -118.551 because makes a round. 

#create the variables of number of decimals of coordinates
precision_table$n.decimals_lon = Ndec(as.character(precision_table$lon))
precision_table$n.decimals_lat = Ndec(as.character(precision_table$lat))
str(precision_table)
head(precision_table)


#test of n.decimals variables creation
Ndec(as.character(precision_table$lon[1])) == precision_table$n.decimals_lon[1]
Ndec(as.character(precision_table$lon[2])) == precision_table$n.decimals_lon[2]
Ndec(as.character(precision_table$lon[3])) == precision_table$n.decimals_lon[3]
Ndec(as.character(precision_table$lon[4])) == precision_table$n.decimals_lon[4]
Ndec(as.character(precision_table$lon[5])) == precision_table$n.decimals_lon[5]
Ndec(as.character(precision_table$lon[456])) == precision_table$n.decimals_lon[456]
Ndec(as.character(precision_table$lon[842])) == precision_table$n.decimals_lon[842]
Ndec(as.character(precision_table$lon[236])) == precision_table$n.decimals_lon[236]


############################################
########coordinatePrecision#################
############################################
precision_table_precision = precision_table[!is.na(precision_table$coordinatePrecision),]

##2 or more decimals in both coordinates
nrow(precision_table_precision[precision_table_precision$n.decimals_lon>=2 & precision_table_precision$n.decimals_lat>=2,]) #number of ocurrences with coordinatePrecision that exhibit two or more decimals in longitude AND latitude
nrow(precision_table_precision) - nrow(precision_table_precision[precision_table_precision$n.decimals_lon>=2 & precision_table_precision$n.decimals_lat>=2,]) #ocurrences that don't exhibit two or more decimals in latitude AND longitude

##Different number of decimals in coordinates
nrow(precision_table_precision[precision_table_precision$n.decimals_lon == precision_table_precision$n.decimals_lat,]) #number of ocurrences with the same number of decimals in lon and lat
nrow(precision_table_precision[!precision_table_precision$n.decimals_lon == precision_table_precision$n.decimals_lat,])#number of ocurrences with different number of decimals in lon and lat

##2 or more decimanl in latitude or longitude
nrow(precision_table_precision[precision_table_precision$n.decimals_lon>=2 | precision_table_precision$n.decimals_lat>=2,]) #number of ocurrences with coordinatePrecision that exhibit two or more decimals in longitude OR latitude (some ocurrences have 
nrow(precision_table_precision) - nrow(precision_table_precision[precision_table_precision$n.decimals_lon>=2 | precision_table_precision$n.decimals_lat>=2,]) #ocurrences that don't exhibit two or more decimals in latitude OR longitude 

#loop for calculating the number of ocurrences with 2 or more decimals for each value of coordinatePrecision.
cPrecision = NULL
high_precision_2_coords = NULL
low_precision_2_coords = NULL
high_precision_1_coords = NULL
low_precision_1_coords = NULL

for (i in unique(precision_table_precision$coordinatePrecision)){ 
        cPrecision = append(cPrecision, i)
        high_precision_2_coords = append(high_precision_2_coords, nrow(precision_table_precision[precision_table_precision$coordinatePrecision == i & precision_table_precision$n.decimals_lon>=2 & precision_table_precision$n.decimals_lat>=2,]))
        high_precision_1_coords = append(high_precision_1_coords, nrow(precision_table_precision[precision_table_precision$coordinatePrecision == i & precision_table_precision$n.decimals_lon>=2 | precision_table_precision$coordinatePrecision == i & precision_table_precision$n.decimals_lat>=2,]))
        low_precision_2_coords = append(low_precision_2_coords, nrow(precision_table_precision[precision_table_precision$coordinatePrecision == i,]) - nrow(precision_table_precision[precision_table_precision$coordinatePrecision == i & precision_table_precision$n.decimals_lon>=2 & precision_table_precision$n.decimals_lat>=2,]))
        low_precision_1_coords = append(low_precision_1_coords, nrow(precision_table_precision[precision_table_precision$coordinatePrecision == i,]) - nrow(precision_table_precision[precision_table_precision$coordinatePrecision == i & precision_table_precision$n.decimals_lon>=2 | precision_table_precision$coordinatePrecision == i & precision_table_precision$n.decimals_lat>=2,]))
}

coordinatePrecision = as.data.frame(cbind(cPrecision, high_precision_2_coords, low_precision_2_coords, high_precision_1_coords, low_precision_1_coords))
str(coordinatePrecision)

coordinatePrecision = coordinatePrecision[order(coordinatePrecision$cPrecision),]

write.csv(coordinatePrecision, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinatePrecision/coordinatePrecision.csv", row.names=FALSE) #the cPrecison with more points with more than 2 decimals in the coordinates are in this order: 25, 1, 10 and 5. 

############################################
########coordinateUncertaintyInMeters#######
############################################

##Select ocurrences with values of coordinateUncertaintyInMeters from all ocurrences
precision_table_uncert = precision_table[!is.na(precision_table$coordinateUncertaintyInMeters),]

##2 or more decimals in both coordinates
nrow(precision_table_uncert[precision_table_uncert$n.decimals_lon>=2 & precision_table_uncert$n.decimals_lat>=2,]) #number of ocurrences with coordinateUncertaintyInMeters that exhibit two or more decimals in longitude AND latitude
nrow(precision_table_uncert) - nrow(precision_table_uncert[precision_table_uncert$n.decimals_lon>=2 & precision_table_uncert$n.decimals_lat>=2,]) #ocurrences that don't exhibit two or more decimals in latitude AND longitude

##Different number of decimals in coordinates
nrow(precision_table_uncert[precision_table_uncert$n.decimals_lon == precision_table_uncert$n.decimals_lat,]) #number of ocurrences with the same number of decimals in lon and lat
nrow(precision_table_uncert[!precision_table_uncert$n.decimals_lon == precision_table_uncert$n.decimals_lat,])#number of ocurrences with different number of decimals in lon and lat

##2 or more decimanl in latitude or longitude
nrow(precision_table_uncert[precision_table_uncert$n.decimals_lon>=2 | precision_table_uncert$n.decimals_lat>=2,]) #number of ocurrences with coordinateUncertaintyInMeters that exhibit two or more decimals in longitude OR latitude (some ocurrences have 
nrow(precision_table_uncert) - nrow(precision_table_uncert[precision_table_uncert$n.decimals_lon>=2 | precision_table_uncert$n.decimals_lat>=2,]) #ocurrences that don't exhibit two or more decimals in latitude OR longitude 


##number of coordinate decimals below and above of 5000. 
#above of 5000
a = nrow(precision_table_uncert[precision_table_uncert$coordinateUncertaintyInMeters > 5000 & precision_table_uncert$n.decimals_lon>=2 & precision_table_uncert$n.decimals_lat>=2,]) #2 or more decimals in BOTH coordinates for ocurrences with Uncertainty higher than 5000
b = nrow(precision_table_uncert[precision_table_uncert$coordinateUncertaintyInMeters > 5000,]) - a #ocurrences that don't have more than 2 decimals in both coordinates. they could hace more than two decimal in longitude OR latitude (not in both), or in any of them. With Uncertainty higher than 5000 

c = nrow(precision_table_uncert[precision_table_uncert$coordinateUncertaintyInMeters > 5000 & precision_table_uncert$n.decimals_lon>=2 | precision_table_uncert$coordinateUncertaintyInMeters > 5000 & precision_table_uncert$n.decimals_lat>=2,]) #2 or more decimals in longitude OR latitude for ocurrences with Uncertainty higher than 5000
d = nrow(precision_table_uncert[precision_table_uncert$coordinateUncertaintyInMeters > 5000,]) - c #less than 2 decimals in both coordinates. With Uncertainty higher than 5000 

#below of 5000
e = nrow(precision_table_uncert[precision_table_uncert$coordinateUncertaintyInMeters <= 5000 & precision_table_uncert$n.decimals_lon>=2 & precision_table_uncert$n.decimals_lat >= 2,]) #2 or more decimals in longitude AND latitude for ocurrences with Uncertainty lower than 5000
f = nrow(precision_table_uncert[precision_table_uncert$coordinateUncertaintyInMeters <= 5000,]) - e #ocurrences that don't have more than 2 decimals in both coordinates. they could hace more than two decimal in longitude OR latitude (not in both), or in any of them. With Uncertainty lower than 5000 

g = nrow(precision_table_uncert[(precision_table_uncert$coordinateUncertaintyInMeters <= 5000 & precision_table_uncert$n.decimals_lon>=2 | precision_table_uncert$coordinateUncertaintyInMeters <= 5000 & precision_table_uncert$n.decimals_lat >= 2),]) #2 or more decimals in longitude OR latitude for ocurrences with Uncertainty lower than 5000
h = nrow(precision_table_uncert[precision_table_uncert$coordinateUncertaintyInMeters <= 5000,]) - g#less than 2 decimals in both coordinates. With Uncertainty higher than 5000 with Uncertainty lower than 5000

coordinateUncertaintyInMeters = rbind(c(a,b,c,d),c(e,f,g,h))
rownames(coordinateUncertaintyInMeters) <- c(">5000", "<=5000")
colnames(coordinateUncertaintyInMeters) <- c(">=2 lon&lat", "<2 lon&lat;<2 lon|lat", ">=2 lon|lat", "<2 lon&lat")
coordinateUncertaintyInMeters 
    #rows: More than 5000 of coordinateUncertaintyInMeters; Lower than 5000 of coordinateUncertaintyInMeters. 
    #columns: More than 2 decimals in both coordinates; Less than 2 decimals in both coordinates and less than two decimals in lon OR lat; more than 2 decimals in one of the coordinates; less than 2 decimals in both coordinates. The most interesingt columns are the first and the fourth. 
    #In general there is a lot of points with more than 2 decimals in both categories (more and low 5000), but there are more in the group of less 5000. 

write.csv(coordinateUncertaintyInMeters, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinateUncertaintyInMeters/coordinateUncertaintyInMeters.csv", row.names=TRUE)

###En coordinatePreccision parece que de 1 a 25 hay muchas ocurrencias con más de 2 decimales, mientras que en coordinateUncertaintyInMeters hay bastantes más ocurrencias con dos decimales entre aquellas que tienen valores menores de 5000 metros. Sin embargo esto no es definitivo. 

#############################################################################
########Export to kml points with coordinatePrecision #######################
#############################################################################
library("sp")
library("rgdal")

#stablish seed 
#because we will use sample function and we want obtain the same result always 
set.seed(3000) 

###points with coordinatePrecision
cPrecision_Kml = precision_table[!is.na(precision_table$coordinatePrecision),][,c("lon", "lat", "coordinatePrecision")] #take the data frame without NA for coordinatePrecision
str(cPrecision_Kml)

#select points with the coordinatePrecision values with the highest number of points witrh more than two decimals in longitude and latitude (25,10,5,1)
cPrecision_Kml_25 = cPrecision_Kml[cPrecision_Kml$coordinatePrecision==25, c("lon", "lat")]
cPrecision_Kml_10 = cPrecision_Kml[cPrecision_Kml$coordinatePrecision==10, c("lon", "lat")]
cPrecision_Kml_5= cPrecision_Kml[cPrecision_Kml$coordinatePrecision==5, c("lon", "lat")]
cPrecision_Kml_1 = cPrecision_Kml[cPrecision_Kml$coordinatePrecision==1, c("lon", "lat")]

#plot the points
require(raster)
bio1 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio1.asc")
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinatePrecision/all_points_plotted.png", width=800, height=800, pointsize=30)
plot(crop(bio1,c(-30,40,25,75)), col="gray80")
points(cPrecision_Kml_25$lon, cPrecision_Kml_25$lat, col="black", cex=0.01)
points(cPrecision_Kml_10$lon, cPrecision_Kml_10$lat, col="red", cex=0.01)
points(cPrecision_Kml_5$lon, cPrecision_Kml_5$lat, col="green", cex=0.01)
points(cPrecision_Kml_1$lon, cPrecision_Kml_1$lat, col="yellow", cex=0.01)
legend("topright", legend=c(25,10,5,1), fill=c("black", "red", "green", "yellow"), cex=0.5)
dev.off()

#select 100 random points
cPrecision_Kml_25_sample = cPrecision_Kml_25[sample(x=1:nrow(cPrecision_Kml_25), size=100),]
cPrecision_Kml_10_sample =cPrecision_Kml_10[sample(x=1:nrow(cPrecision_Kml_10), size=100),]
cPrecision_Kml_5_sample =cPrecision_Kml_5[sample(x=1:nrow(cPrecision_Kml_5), size=100),]
cPrecision_Kml_1_sample =cPrecision_Kml_1[sample(x=1:nrow(cPrecision_Kml_1), size=100),]

#Convert to SpatialPoints
coordinates(cPrecision_Kml_25_sample) <- c("lon", "lat") 
coordinates(cPrecision_Kml_10_sample) <- c("lon", "lat") 
coordinates(cPrecision_Kml_5_sample) <- c("lon", "lat") 
coordinates(cPrecision_Kml_1_sample) <- c("lon", "lat") #indicte coordinates

#create the .kml que the KML function of Raster
KML(cPrecision_Kml_25_sample, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinatePrecision/coordinatePrecision_25.kml", overwrite=TRUE)
KML(cPrecision_Kml_10_sample, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinatePrecision/coordinatePrecision_10.kml", overwrite=TRUE)
KML(cPrecision_Kml_5_sample, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinatePrecision/coordinatePrecision_5.kml", overwrite=TRUE)
KML(cPrecision_Kml_1_sample, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinatePrecision/coordinatePrecision_1.kml", overwrite=TRUE) #In general, all cP exhibit a good precision. Their points are very close to forests. 

##Repeat all the process but for ALL values of coordinatePrecision. 
for (i in unique(cPrecision_Kml$coordinatePrecision)){
    table = cPrecision_Kml[cPrecision_Kml$coordinatePrecision==i, c("lon", "lat")]
    if (nrow(table)>100){
        table_sample = table[sample(x=1:nrow(table), size=100),]
    } else {
        table_sample = table
    }
    coordinates(table_sample) <- c("lon", "lat") 
    KML(table_sample, paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinatePrecision/coordinatePrecision", paste(i, "kml", sep="."),sep="_"), overwrite=TRUE)
} #In general, all cP exhibit a good precision. Their points are very close to forests. 


#######################################################################################
########Export to kml points with coordinateUncertaintyInMeters #######################
#######################################################################################
library("sp")
library("rgdal")

#stablish seed 
#because we will use sample function and we want obtain the same result always 
set.seed(3000) 

###points with coordinatePrecision
cUncertaintyInMeters_Kml = precision_table[!is.na(precision_table$coordinateUncertaintyInMeters),][,c("lon", "lat", "coordinateUncertaintyInMeters")] #take the data frame without NA for coordinatePrecision
str(cUncertaintyInMeters_Kml)

#select points with the coordinateUncertaintyInMeters higher and lower than 5000
cUncertaintyInMeters_less_5000_Kml = cUncertaintyInMeters_Kml[cUncertaintyInMeters_Kml$coordinateUncertaintyInMeters<=5000, c("lon", "lat")]
cUncertaintyInMeters_plus_5000_Kml = cUncertaintyInMeters_Kml[cUncertaintyInMeters_Kml$coordinateUncertaintyInMeters>5000, c("lon", "lat")]

#plot the points
require(raster)
bio1 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/climate/bioclim/wc2-5/bio1.asc")
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinateUncertaintyInMeters/all_points_plotted.png", width=1400, height=800, pointsize=30)
par(mfcol=c(1,2))
plot(bio1, col="gray80", main="coordinateUncertaintyInMeters less 5000")
points(cUncertaintyInMeters_less_5000_Kml$lon, cUncertaintyInMeters_less_5000_Kml$lat, col="black", cex=0.1)
plot(bio1, col="gray80", main="coordinateUncertaintyInMeters plus 5000")
points(cUncertaintyInMeters_plus_5000_Kml$lon, cUncertaintyInMeters_plus_5000_Kml$lat, col="black", cex=0.1)
dev.off()

#select 100 random points
cUncertaintyInMeters_less_5000_Kml_sample = cUncertaintyInMeters_less_5000_Kml[sample(x=1:nrow(cUncertaintyInMeters_less_5000_Kml), size=100),]
cUncertaintyInMeters_plus_5000_Kml_sample =cUncertaintyInMeters_plus_5000_Kml[sample(x=1:nrow(cUncertaintyInMeters_plus_5000_Kml), size=100),]

#Convert to SpatialPoints
coordinates(cUncertaintyInMeters_less_5000_Kml_sample) <- c("lon", "lat") 
coordinates(cUncertaintyInMeters_plus_5000_Kml_sample) <- c("lon", "lat") 

#create the .kml que the KML function of Raster
KML(cUncertaintyInMeters_less_5000_Kml_sample, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinateUncertaintyInMeters/cUncertaintyInMeters_less_5000_Kml_sample.kml", overwrite=TRUE)
KML(cUncertaintyInMeters_plus_5000_Kml_sample, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinateUncertaintyInMeters/cUncertaintyInMeters_plus_5000_Kml_sample.kml", overwrite=TRUE)

##### With pinus albicaulis #####
pinus_albicaulis = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_albicaulis.csv", header=TRUE)

#no coordinatePrecision
table(pinus_albicaulis$coordinatePrecision)

pinus_albicaulis_no_NA = pinus_albicaulis[!is.na(pinus_albicaulis$coordinateUncertaintyInMeters), ]
pinus_albicaulis_less_5000 = pinus_albicaulis_no_NA[pinus_albicaulis_no_NA$coordinateUncertaintyInMeters<=5000, c("lon", "lat")]
pinus_albicaulis_plus_5000 = pinus_albicaulis_no_NA[pinus_albicaulis_no_NA$coordinateUncertaintyInMeters>5000, c("lon", "lat")]

#plot points
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinateUncertaintyInMeters/pinus_albicaulis.pdf")
par(mfcol=c(1,2))
plot(crop(bio1, c(-150, -50, 20, 60)), col="gray80", main="Pinus_albicaulis less 5000")
points(pinus_albicaulis_less_5000$lon, pinus_albicaulis_less_5000$lat, col="black", cex=0.1)
plot(crop(bio1, c(-150, -50, 20, 60)), col="gray80", main="Pinus_albicaulis plus 5000")
points(pinus_albicaulis_plus_5000$lon, pinus_albicaulis_plus_5000$lat, col="black", cex=0.1)
dev.off()

#convert to spatial point object
coordinates(pinus_albicaulis_less_5000) <- c("lon", "lat") 
coordinates(pinus_albicaulis_plus_5000) <- c("lon", "lat") 

#save as kml
KML(pinus_albicaulis_less_5000, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinateUncertaintyInMeters/pinus_albicaulis_less_5000.kml", overwrite=TRUE)
KML(pinus_albicaulis_plus_5000, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinateUncertaintyInMeters/pinus_albicaulis_plus_5000.kml", overwrite=TRUE) #No se ven muchas diferencias, ambos tipos de puntos están cerca de bosque.

#### With pinus sylvestris ####
pinus_sylvestris = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/raw_ocurrences/ocurrences/tree_species/Pinus_sylvestris.csv", header=TRUE)

#there is data of cP
table(pinus_sylvestris$coordinatePrecision)

#also of coordinateUncertaintyInMeters
table(pinus_sylvestris$coordinateUncertaintyInMeters)

#summary
summary(pinus_sylvestris$coordinatePrecision)
summary(pinus_sylvestris$coordinateUncertaintyInMeters)
summary(pinus_sylvestris$lon)
summary(pinus_sylvestris$lat)

#drop rows with Na in lon, lat or cU
pinus_sylvestris_clean = pinus_sylvestris[!(is.na(pinus_sylvestris$lon) | is.na(pinus_sylvestris$lat) | is.na(pinus_sylvestris$coordinateUncertaintyInMeters)),]
summary(pinus_sylvestris_clean$coordinateUncertaintyInMeters)
summary(pinus_sylvestris_clean$lon)
summary(pinus_sylvestris_clean$lat) #no NA in any case

#select points with cU lower and higher than 5000
pinus_sylvestris_less_5000 = pinus_sylvestris_clean[pinus_sylvestris_clean$coordinateUncertaintyInMeters<=5000, c("lon", "lat")]
pinus_sylvestris_plus_5000 = pinus_sylvestris_clean[pinus_sylvestris_clean$coordinateUncertaintyInMeters>5000, c("lon", "lat")]

#plot points
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinateUncertaintyInMeters/pinus_sylvestris.pdf")
par(mfcol=c(1,2))
plot(crop(bio1, c(-10, 60, 35, 75)), col="gray80", main="Pinus_sylvestris less 5000")
points(pinus_sylvestris_less_5000$lon, pinus_sylvestris_less_5000$lat, col="black", cex=0.1)
plot(crop(bio1, c(-10, 60, 35, 75)), col="gray80", main="Pinus_sylvestris plus 5000")
points(pinus_sylvestris_plus_5000$lon, pinus_sylvestris_plus_5000$lat, col="black", cex=0.1)
dev.off()

#convert to spatial point object
coordinates(pinus_sylvestris_less_5000) <- c("lon", "lat") 
coordinates(pinus_sylvestris_plus_5000) <- c("lon", "lat") 

#save as kml
KML(pinus_sylvestris_less_5000, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinateUncertaintyInMeters/pinus_sylvestris_less_5000.kml", overwrite=TRUE)
KML(pinus_sylvestris_plus_5000, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/ocurrences_precision/coordinateUncertaintyInMeters/pinus_sylvestris_plus_5000.kml", overwrite=TRUE) ###points with lower than 5000 have forest closer than 5000 meters. 

#save work space
save.image("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/rdata/test_precision_points.RData")
