#RUN ONE TIME: This code takes the new sylvestris maps (fixed by Bianca) and drop sea areas from them. For that uses the continental limits of other distribution map (albicaulis)

require(raster)

#load distributions
sylvestris = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/datos_brutos/Maps/p_sylvestris_michael.img")
albicaulis = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_albicaulis_01.img") #use the limite of othe map species

albicaulis[albicaulis==1] <- 0 #eliminate the distribution of albicaulis
plot(albicaulis)


sylvestris[sylvestris==0] <- NA #eliminat the areas outside the distribution
plot(sylvestris)

new_sylvestris = merge(sylvestris, albicaulis) #merge both raster
plot(new_sylvestris)

writeRaster(new_sylvestris, "/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_sylvestris_01.img", overwrite=TRUE) 