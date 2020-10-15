#Code for making calculations of variables: calculate mean and resample to get the same resolution. 

#######################################
###########Create variables############
#######################################
require(raster)

###list variables
list_bioclim_moisture = list.files("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/climate/bioclim_moisture", pattern=".asc", full.names=TRUE)
list_clay = list.files(path="/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/soil/physical/clay_content", pattern=".tif", full.names=TRUE)
list_silt = list.files(path="/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/soil/physical/silt_content", pattern=".tif", full.names=TRUE)
list_sand = list.files(path="/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/soil/physical/sand_content", pattern=".tif", full.names=TRUE)
list_ph = list.files(path="/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/soil/chemical/ph", pattern=".tif", full.names=TRUE)
list_cec = list.files(path="/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/soil/chemical/cec", pattern=".tif", full.names=TRUE)
list_carbon = list.files(path="/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/soil/chemical/carbon", pattern=".tif", full.names=TRUE)

###load variables
require(raster)
bioclim_moisture = stack(list_bioclim_moisture)
clay_layers = stack(list_clay)
silt_layers = stack(list_silt)
sand_layers = stack(list_sand)
ph_layers = stack(list_ph)
cec_layers = stack(list_cec)
carbon_layers = stack(list_carbon)

#calculate mean of soil variable for all depths 
clay = mean(clay_layers)
silt = mean(silt_layers) 
sand = mean(sand_layers) 
ph = mean(ph_layers) 
cec = mean(cec_layers) 
carbon = mean(carbon_layers) 
depth = raster("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/soil/site_characteristics/depth/BDTICM_M_1km_ll.tif")
soil = brick(clay, silt, sand, ph, cec, carbon, depth) #include all soil variables in a brick (faster computation)
names(soil) = c("clay", "silt", "sand", "ph", "cec", "carbon", "depth")
soil_mask = soil[["carbon"]] * soil[["cec"]] * soil[["clay"]] * soil[["depth"]] * soil[["ph"]] * soil[["sand"]] * soil[["silt"]] #mask to reproduce all NAs 
soil=mask(soil, soil_mask) #apply the maks to all soil variable to reproduce all NAs in all of them. 
soil.resampled = resample(x=soil, bioclim_moisture[[1]], method="bilinear") #Resample soil variable to the resolution of bioclim variables. bilinear method is better than nearest neigbourg for continous variables

##stack all variables in a unique stack
variables = stack(bioclim_moisture, soil.resampled)
writeRaster(variables, filename="/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/finals/", suffix=names(variables), bylayer=TRUE, format="ascii", overwrite=TRUE)
