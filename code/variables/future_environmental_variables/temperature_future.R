#######Code for calculating tavg (average temperature) for all climatic scenarios of the future (2070) and crop solar dation raster to tavg extension


######################################################
#####Creating tavg for all climatics scenarios########
######################################################

#library
require(raster)

#load the names of all climatic scenarios
climatic_scenarios = c("bc26", "bc45", "bc60", "bc85", "cc26", "cc45", "cc60", "cc85", "gs26", "gs45", "gs60", "gs85", "he26", "he45", "he60", "he85", "ip26", "ip45", "ip60", "ip85", "mr26", "mr45", "mr60", "mr85", "mg26", "mg45", "mg60", "mg85")
    
#run a loop for each scenario and for each month 
for (k in climatic_scenarios){

    dir.create(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/tavg/", k, "tavg70", sep=""))

    for (i in 1:12){
        tmax = raster(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/tmax/", k, "tx70/", k, "tx70", i, ".tif", sep=""))
        tmin = raster(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/tmin/", k, "tn70/", k, "tn70", i, ".tif", sep=""))
        tavg = overlay(tmin, tmax, fun  = function(x) sum(x)/2)

        writeRaster(tavg, paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/tavg/", k, "tavg70/", k, "tavg70", i, ".tif", sep=""), overwrite=TRUE)
    }

}


######################################################
#####Crop solar radiation to tavg ####################
######################################################
#Data of solar radiation have antarctica, but we are not interested in this area, so we are going to drop it. 

#list all raster of solar radiation
list_solar = list.files("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/solar_rad", pattern=".tif", full.name=TRUE)

#stack all of them
solar_radiation = stack(list_solar)

#load one raster ot tavg (randomly selected)
one_raster_temperature = raster("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/tavg/bc26tavg70/bc26tavg701.tif")

#crop all solar radiation raster with tavg
solar_radiation_crop = crop(solar_radiation, one_raster_temperature)

#write the stack
writeRaster(solar_radiation_crop, filename="/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/solar_rad/srad_crop", overwrite=TRUE, bylayer=TRUE, suffix='numbers', format="ascii")

#save.image
save.image("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/rdata/calulation_temperature.RData")
