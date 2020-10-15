#####Create bioclim variable from tmin, tmax and moisture index (instead of precipitation)

####define work directory
setwd("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus")

#Libraries
require(raster)

##Load all world raster of moisture index for each month 
moisture = stack("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/climate/moisture/wc2.1_5m_mind/mind.grd")

##Load tmin and tmax data at 5 min 
#List tmin raster for each month
list.tmin = list.files(path="/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/climate/worldclim/tmin_5m_bil", pattern=".bil", full.names=TRUE)

#List tmax raster for each month
list.tmax = list.files(path="/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/climate/worldclim/tmax_5m_bil", pattern=".bil", full.names=TRUE)

#create a stack with all raster for each varaible
tmin = stack(list.tmin)
tmax = stack(list.tmax) #create the stacks

#crop moisture variables becasue there is no data of tmin and tmax in antarctica
moisture = crop(moisture, tmin)

#test if all stacks have the same resolution, extent, ncol and nrow
xres(moisture)
xres(tmin) 
xres(tmax)
extent(moisture)
extent(tmin)
extent(tmax) 
ncol(moisture)
ncol(tmin)
ncol(tmax) 
nrow(moisture)
nrow(tmin)
nrow(tmax) 
ncell(moisture)
ncell(tmin)
ncell(tmax) #same resolution, extent, ncol, nrow and ncell

#plot all stacks
plot(moisture)
plot(tmin)
plot(tmax)

#create the biolcim variables of moisture with "biovars"
require(dismo)
bio_moisture = biovars(prec = moisture, #vector, matrix, or RasterStack/Brick of precipitation data. In our case is moisture index of each month for the whole globe
    tmin = tmin, #vector, matrix, or RasterStack/Brick of minimum temperature data
    tmax = tmax) #vector, matrix, or RasterStack/Brick of maximum temperature data

#problem bio15: 
plot(bio_moisture[["bio15"]]) #extrange results.
#bio15 is the Precipitation Seasonality (Coefficient of Variation), in our case moisture seasonality. Coefficient of Variation is sd/mean (https://en.wikipedia.org/wiki/Coefficient_of_variation), and it makes sense use it for precipitation, because this variable can not exhibit negative values (min value always is zero). However, our variable is moisture index, that can exhibit negative values, which indicate that evapotanspiration is higher than precipitation. In this case, the coefficiente of variation would be very negative. We will use tha same way than bioclim use for temperaturature seasonality, use only standart deviation. 
bio_moisture[["bio15"]] = calc(moisture, sd) #use calc for calculate the sd of the raster of the stack, and change bio15 form biovars by our variable. 
plot(bio_moisture[["bio15"]])

#save de rasters
writeRaster(x = bio_moisture, filename="/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/climate/bioclim_moisture/", suffix=names(bio_moisture), bylayer=TRUE, format="ascii", overwrite=TRUE)

save.image("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/rdata/bioclim_variables_with_moisture.RData")


