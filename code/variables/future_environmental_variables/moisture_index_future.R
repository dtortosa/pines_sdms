######## Code for creating moisture index of each scenario and month in the future (2070)

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



###########################################################
#####Create etpt and mind  in ONE step ####################
###########################################################

#library
require(foreach) #for repeat a process several times 
require(doParallel) #for parallel

#create variables of months and scenarios
months <- formatC(1:12) #12 months
scenarios <- c("bc26", "bc45", "bc60", "bc85", "cc26", "cc45", "cc60", "cc85", "gs26", "gs45", "gs60", "gs85", "he26", "he45", "he60", "he85", "ip26", "ip45", "ip60", "ip85", "mr26", "mr45", "mr60", "mr85", "mg26", "mg45", "mg60", "mg85") #28 scenarios
scenarios_months <- expand.grid(scenarios,months) #create all possible combinatios in a data frmae with tow columns
scenarios_months <- with(scenarios_months, interaction(Var1, Var2, sep="_")) #combine the two columns in a unique factor

#load the function to calculate evapotranspiration vectorized 
overfun <- Vectorize(function(x, y, z) {
    ifelse(z < 0, 0, (x / (x + 15)) * ((0.0239001 * y) + 50.0) * (4.0 / 30.0))
}) #selected from Rafi (see "/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/climate/moisture/wc2.1_5m_etpt/etpt.r") 

#load a function to calculate the number of days of a month (variable according to leap years or not). Number of days will be used in the calculation of moisture index
days.in.month <- function(month, year = NULL){
    month = as.integer(month)
    if (is.null(year))
        year = as.numeric(format(Sys.Date(), '%Y'))
        dt = as.Date(paste(year, month, '01', sep = '-'))
        dates = seq(dt, by = 'month', length = 2)
        as.numeric(difftime(dates[2], dates[1], units = 'days'))
} # from: http://stackoverflow.com/questions/6243088/find-out-the-number-of-days-of-a-month-in-r

#load the function to create all: evapotranspiration and moisture index. 
fefun <- function(index){

    #split the name "scenario"_"month" into these two element
    splitted_name <- str_split_fixed(index, "_", 2)

    #load tavg for the corresponding scenario and month
    ta <- raster(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/tavg/", splitted_name[1], "tavg70/", splitted_name[1], "tavg70", splitted_name[2], ".tif", sep="")) #it has been calculated previously as the mean of tmin and tmax
    
    #load solar radiation for the correspoinding month
    sr <- raster(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/solar_rad/srad_crop", splitted_name[2], ".asc", sep="")) #we assume that the radiation will be the same in 2070
    
    #calculate etpt with overfun function
    etpt <- overlay(ta, sr, ta, fun = overfun)

    #load precipitation 
    prec <- raster(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/pp/", splitted_name[1], "pr70/", splitted_name[1], "pr70", splitted_name[2], ".tif", sep=""))
    
    #calculate the number of the days of the corresponding month
    days <- days.in.month(splitted_name[2], year = 2070)
    
    #calculate moisture index
    mind <- prec - etpt / 10 * days

    ##export
    writeRaster(mind, filename = paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/mind/", splitted_name[1], "mind70", splitted_name[2], ".tif", sep=""), overwrite=TRUE)

    writeRaster(etpt, paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/etpt/", splitted_name[1], "etpt70", splitted_name[2], ".tif", sep=""), overwrite=TRUE)

    #delete all data for saving memory
    rm(ta, sr, etpt, prec, mind)
    gc()    

}

# set up cluster
clust <- makeCluster(4)
registerDoParallel(clust)

# run
foreach(i = scenarios_months, .packages = c("raster", "stringr")) %dopar% { 
    fefun(index = i)
} #the "stringr" package is used for "str_split_fixed" function

#stop the cluster 
stopCluster(clust)

################################################################
#####Create bioclim variables with moisture ####################
################################################################

#library
require(foreach) #for repeat a process several times 
require(doParallel) #for parallel

#create variables of months and scenarios
scenarios <- c("bc26", "bc45", "bc60", "bc85", "cc26", "cc45", "cc60", "cc85", "gs26", "gs45", "gs60", "gs85", "he26", "he45", "he60", "he85", "ip26", "ip45", "ip60", "ip85", "mr26", "mr45", "mr60", "mr85", "mg26", "mg45", "mg60", "mg85") #28 scenarios

biofun = function(index){

    #list of tmin raster for the corresponding scenario
    tnlist = list.files(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/tmin/", index, "tn70", sep=""), pattern=".tif", full.names=TRUE)

    #stack them
    tmin=stack(tnlist)

    #order to same order of months 
    tmin = tmin[[c(paste(index, "tn701", sep=""), paste(index, "tn702", sep=""), paste(index, "tn703", sep=""), paste(index, "tn704", sep=""), paste(index, "tn705", sep=""), paste(index, "tn706", sep=""), paste(index, "tn707", sep=""), paste(index, "tn708", sep=""), paste(index, "tn709", sep=""), paste(index, "tn7010", sep=""), paste(index, "tn7011", sep=""), paste(index, "tn7012", sep=""))]]

    #list of tmax raster for the corresponding scenario
    txlist = list.files(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/tmax/", index, "tx70", sep=""), pattern=".tif", full.names=TRUE)

    #stack them
    tmax=stack(txlist)

    #order to same order of months 
    tmax = tmax[[c(paste(index, "tx701", sep=""), paste(index, "tx702", sep=""), paste(index, "tx703", sep=""), paste(index, "tx704", sep=""), paste(index, "tx705", sep=""), paste(index, "tx706", sep=""), paste(index, "tx707", sep=""), paste(index, "tx708", sep=""), paste(index, "tx709", sep=""), paste(index, "tx7010", sep=""), paste(index, "tx7011", sep=""), paste(index, "tx7012", sep=""))]]

    #list of mind raster for the corresponding scenario
    mindlist = list.files("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/mind/", pattern=paste(index, "mind70", sep=""), full.names=TRUE)

    #stack them
    mind=stack(mindlist)

    #order to same order of months 
    mind = mind[[c(paste(index, "mind701", sep=""), paste(index, "mind702", sep=""), paste(index, "mind703", sep=""), paste(index, "mind704", sep=""), paste(index, "mind705", sep=""), paste(index, "mind706", sep=""), paste(index, "mind707", sep=""), paste(index, "mind708", sep=""), paste(index, "mind709", sep=""), paste(index, "mind7010", sep=""), paste(index, "mind7011", sep=""), paste(index, "mind7012", sep=""))]]

    #create bioclim variables with moisture 
    bio_moisture = biovars(prec = mind, tmin = tmin, tmax = tmax) 

    #problem bio15: 
    bio_moisture[["bio15"]] = calc(mind, sd) #for more info see ("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/code/variables_presences/bioclim_variables_with_moisture.R")

    #export bioclim variables 
    writeRaster(bio_moisture, filename=paste(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/bioclim_moisture/bio", 1:nlayers(bio_moisture), sep=""), index, sep="_"), format="ascii", bylayer=TRUE, overwrite=TRUE)

    #delete all data for saving memory
    rm(tmin, tmax, bio_moisture)
    gc()

}

# set up cluster
clust <- makeCluster(4)
registerDoParallel(clust)

# run
foreach(i = scenarios, .packages = c("raster", "dismo")) %dopar% { 
    biofun(index = i)
} #the "dismo" package is used for "biovars" function, which creates bioclim variables

#stop the cluster 
stopCluster(clust)

#Comprobation that all scenarios have 19 bioclim variables
test = data.frame()
for (i in scenarios){
    list = list.files("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/bioclim_moisture", pattern=i, full.names=TRUE)
    stack_variables = stack(list)
    bio_names = c(paste("bio1", i, sep="_"), paste("bio10", i, sep="_"), paste("bio11", i, sep="_"), paste("bio12", i, sep="_"), paste("bio13", i, sep="_"), paste("bio14", i, sep="_"), paste("bio15", i, sep="_"), paste("bio16", i, sep="_"), paste("bio17", i, sep="_"), paste("bio18", i, sep="_"), paste("bio19", i, sep="_"), paste("bio2", i, sep="_"), paste("bio3", i, sep="_"), paste("bio4", i, sep="_"), paste("bio5", i, sep="_"), paste("bio6", i, sep="_"), paste("bio7", i, sep="_"), paste("bio8", i, sep="_"), paste("bio9", i, sep="_"))

    test = rbind(test, cbind(i, names(stack_variables) == bio_names))
}
names(test) <- c("scenario", "test")
test #All have to be TRUE and this is the case 

