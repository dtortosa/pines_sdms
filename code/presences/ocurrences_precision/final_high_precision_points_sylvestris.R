####check that FINAL high precision points of sylvestris do not came from a atlas.
#I am concern about of the high precision points of syvlestris ein Sweden and norway, there are a lot and i cannot see precisely if there are a equally separated (data from a flora, like those of germany that has been considered low precision points) or not.  

#required packages
require(sp)
require(rgdal)
require(raster)

#selected species
species = "sylvestris"

#load complete occurrence data
ocurrence_data = read.csv(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ocurrences/", species, "_complete.presences.csv", sep=""), header=TRUE)

#select high precision points
high_precision_points = ocurrence_data[which(ocurrence_data$precision_weight==1), which(colnames(ocurrence_data) %in% c("longitude", "latitude"))]
str(high_precision_points)

#conver to spatialPoint
coordinates(high_precision_points) <- c("longitude", "latitude")
str(high_precision_points)

#plot
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/final_high_precision_points_sweden_sylvestris.pdf")
plot(crop(high_precision_points, c(0,25,55,70))) #it does not seem a flora, there is not a regular pattern, a similar separation between points.
dev.off() 

#convert to KML
KML(high_precision_points, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/ocurrences_precision/final_high_precision_points_sylvestris.kml", overwrite=TRUE) #Open with google earth. The high precision points inside the Iberian peninsula fall inside forest, indeed points close to Granada do not fall inside the Zubia pines (likely Pinaster), instead fall closer to Botanical garde of the Cortijuela. Points in Sweden fall inside or very close to forest, so great. 