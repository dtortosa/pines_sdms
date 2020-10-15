######## Code for checking that mindez of w2 is similar to mindex of wc1.4
#I have only one thing for which I'm little worried, and I would like to comment to you: For calculating current bioclimatic variables in the niche paper, I used temperature variables from wc 1.4 and a moisture index calculated with wc2. I calculated the moisture index for future CC scenarios, but for current climate I used the index you calculated, which I think was calculates with wc2. I did not realize that when I was doing the analyses. As I mention before, in all cases units seems to be ok (check done with chelsa), but I don't know if other issues could raise from calculating bioclimatic variables with data from different versions.

#The future bioclimatic variables are calculated with wc1.4, but the current variables were calculated with temperature from wc1.4 and m.index from wc2.0. The SDMs were fitted with this mixed data wc1.4/wc.2. Therefore, suitability projections would be affected by any potential problem. As I mentioned before, the units are ok, so any potential problem related to different units between WC versions is unlikely. In addition, I forget an important detail, in any case we will always have to use part of wc2.0 because solar radiation is not included in wc.1.4. Therefore, the moisture index under future conditions is not full wc1.4 neither, because solar radiation used for the moisture index came from wc2.0. I guess I did not move to wc2.0 because that version doesn't include future projections yet. 


######################################################
#####Creating tavg for all climatics scenarios########
######################################################

#library
require(raster)

#unzip the monthly tmax and tmin from wc.14 (5 min resolution)
system("unzip /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmax_5m_bil.zip -d /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmax_5m_bil/")
system("unzip /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmin_5m_bil.zip -d /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmin_5m_bil/")

#create a directory to save average temperature
system("rm -rf /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tavg_5m_bil; mkdir /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tavg_5m_bil")

#loop to calculate average temperature for each month 
for (i in 1:12){

    #load tmax and tmin of the [i] month
    tmax = raster(paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmax_5m_bil/tmax", i, ".bil", sep=""))
    tmin = raster(paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmin_5m_bil/tmin", i, ".bil", sep=""))

    #calculate average temperature of the [i] month
    tavg = overlay(tmin, tmax, fun  = function(x) sum(x)/2)

    #save it
    writeRaster(tavg, paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tavg_5m_bil/tavg", i, ".tif", sep=""), overwrite=TRUE)
}

#remove folders with tmin and tmax. We only leave the zip files
system("rm -rf /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmax_5m_bil")
system("rm -rf /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmin_5m_bil")

#compare avg temperature calculated by me and avg downloaded from worldclim
system("unzip /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tavg_5m_bil_direct_worldclim.zip -d /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tavg_5m_bil_direct_worldclim/")

#compare tavg calculated by me with tavg downloaded directly from worldclim
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/mindex_wc1.4_wc_2_check/differences_tavg_calculated_by_me_&_wc.pdf")
par(mfcol=c(2,2))
for(i in 1:12){

    #load my tavg form month 1
    own_tavg = raster(paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tavg_5m_bil/tavg", i, ".tif", sep=""))
    
    #load tavg directly downloaded from worldclim 
    direct_downloaded_tavg = raster(paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tavg_5m_bil_direct_worldclim/tmean", i, ".bil", sep="")) 

    #plot both and the difference
    plot(own_tavg, main=paste("month ", i, " tavg calculated by me", sep=""))
    plot(direct_downloaded_tavg, main=paste("month ", i, " tavg downloaded from worldclim", sep=""))
    plot(own_tavg - direct_downloaded_tavg, main=paste("month ", i, " differences between both tavgs", sep=""))
    frame()
}
dev.off() #everything is ok, very similar. Differences from 0.04 to -0.1 degrees

#remove the unzipped folder with tavg directly downloaded from worldclim
system("rm -rf /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tavg_5m_bil_direct_worldclim/")


#########################################
#####Crop solar radiation to tavg #######
#########################################
#Data of solar radiation have antarctica, but we are not interested in this area, so we are going to drop it. 

#unzip solar radiation data from wc2 (5 min resolution)
system("unzip /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/wc2.0_5m_srad.zip -d /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/wc2.0_5m_srad/")

#list all raster of solar radiation from wc2
list_solar = list.files("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/wc2.0_5m_srad", pattern=".tif", full.name=TRUE)
length(list_solar) == 12

#stack all of them
solar_radiation = stack(list_solar)
nlayers(solar_radiation) == 12
names(solar_radiation)

#load one raster ot tavg (randomly selected)
one_raster_temperature = raster("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tavg_5m_bil/tavg1.tif")

#crop all solar radiation raster with tavg
solar_radiation_crop = crop(solar_radiation, one_raster_temperature)

#create a new directory
system("rm -rf /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/wc2.0_5m_srad_cropped/; mkdir /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/wc2.0_5m_srad_cropped")

#write the stack
writeRaster(solar_radiation_crop, filename="/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/wc2.0_5m_srad_cropped/wc2.0_5m_srad_cropped", overwrite=TRUE, bylayer=TRUE, suffix='numbers', format="GTiff")

#remove the first directory of solar radiation
system("rm -rf /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/wc2.0_5m_srad")

###########################################################
#####Create etpt and mind in ONE step #####################
###########################################################

#library
require(foreach) #for repeat a process several times 
require(doParallel) #for parallel

#unzip pp data from wc1.4 (5 min res)
system("unzip /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/prec_5m_bil.zip -d /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/prec_5m_bil/")

#create a directory for saving etpt and mindex data
system("mkdir /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/mindex_5m_bil")
system("mkdir /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/etpt_5m_bil")

#create variables of months and scenarios
month <- formatC(1:12) #12 months

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
fefun <- function(month){

    #load tavg for the corresponding scenario and month
    ta <- raster(paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tavg_5m_bil/tavg", month, ".tif", sep="")) #it has been calculated previously as the mean of tmin and tmax
    
    #load solar radiation for the correspoinding month
    sr <- raster(paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/wc2.0_5m_srad_cropped/wc2.0_5m_srad_cropped_", month, ".tif", sep="")) #we assume that the radiation will be the same in 2070
    
    #calculate etpt with overfun function
    etpt <- overlay(ta, sr, ta, fun = overfun)

    #load precipitation 
    prec <- raster(paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/prec_5m_bil/prec", month, ".bil", sep=""))
    
    #calculate the number of the days of the corresponding month
    days <- days.in.month(month, year = 1990)
    
    #calculate moisture index
    mind <- prec - etpt / 10 * days

    ##export
    writeRaster(mind, filename = paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/mindex_5m_bil/mindex_5m_bil", month, ".tif", sep=""), overwrite=TRUE)
    writeRaster(etpt, filename = paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/etpt_5m_bil/etpt_5m_bil", month, ".tif", sep=""), overwrite=TRUE)

    #delete all data for saving memory
    rm(ta, sr, etpt, prec, mind)
    gc()    
}

# set up cluster
clust <- makeCluster(2)
registerDoParallel(clust)

# run
foreach(i = month, .packages = c("raster")) %dopar% { 
    fefun(month = i)
} #the "stringr" package is used for "str_split_fixed" function

#stop the cluster 
stopCluster(clust)

#remove directory of precipitation, we only leave the zip file
system("rm -rf /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/prec_5m_bil")

#########################################
####comparison mindex wc1.4 and wc2

###load mindex calculated by Rafi (wc2). Extracted from "/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/climate.zip" with "extract_rafi_etpt_script_from_climate_zip.sh"
mindex_wc2 = stack("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/wc2.1_5m_mind/mind.grd")
nlayers(mindex_wc2) == 12
names(mindex_wc2)

###load mindex wc1.4
list_rasters_mindex_wc14 = list.files("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/mindex_5m_bil", full.names = TRUE)
mindex_wc14 = stack(list_rasters_mindex_wc14)
nlayers(mindex_wc14) == 12
names(mindex_wc14) #month names are not in order

#remove antarctica from wc2
mindex_wc2_crop = crop(mindex_wc2, mindex_wc14)

###for each month plot the moisture index from different sources
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/mindex_wc1.4_wc_2_check/comparsion_wc1.4_wc2_chelsa.pdf")
for(i in 1:12){

    #select the month [i]
    selected_month = month.name[i]

    #print the selected month
    print(selected_month)

    #select the raster of mindex of wc 1.4
    selected_mindex_wc1.4 = raster(paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/mindex_5m_bil/mindex_5m_bil", i, ".tif", sep=""))

    #select the raster of mindex of wc 2
    selected_mindex_wc2 = mindex_wc2_crop[[i]]

    #select the raster of mindex of chelsa
    if(i %in% 1:9){
        selected_mindex_chelsa = raster(paste("/Volumes/Maxtor/diego/science_big_documents/seed_mass_aridity/really_final_variables/mind_chelsa/unzipped_10x10km/10x10km_CHELSA_mind_0", i, "_V1.2_land.tif", sep=""))
    }else{
        selected_mindex_chelsa = raster(paste("/Volumes/Maxtor/diego/science_big_documents/seed_mass_aridity/really_final_variables/mind_chelsa/unzipped_10x10km/10x10km_CHELSA_mind_", i, "_V1.2_land.tif", sep=""))      
    }

    #Get values of cells without NAs for WC versions
    values_mindex_wc1.4 = getValues(selected_mindex_wc1.4)
    values_mindex_wc2 = getValues(selected_mindex_wc2)

    #extract min and max values of cells in moisture index of wc1.4 and wc2
    max_values_mindex_wc1.4 = max(na.omit(values_mindex_wc1.4))
    max_values_mindex_wc2 = max(na.omit(values_mindex_wc2))
    min_values_mindex_wc1.4 = min(na.omit(values_mindex_wc1.4))
    min_values_mindex_wc2 = min(na.omit(values_mindex_wc2))

    #print the differences in min and max values between wc versions
    max_vals_differs = abs(max_values_mindex_wc1.4 - max_values_mindex_wc2)
    min_vals_differs = abs(min_values_mindex_wc1.4 - min_values_mindex_wc2)
    print(paste("Difference in max values ", round(max_vals_differs, 2), sep=""))
    print(paste("Difference in min values ", round(min_vals_differs, 2), sep=""))

    #indicate if we need to change any of the version for the min and max values
    max_min_index = NULL
    version_index = NULL
    if(max_vals_differs > 400){#if the differences between max values are higher than 400

        #the maximum value is the target
        max_min_index = "max"

        #if the maximum value of wc1.4 is higher
        if(max_values_mindex_wc1.4 > max_values_mindex_wc2){

            #wc 1.4 is the target
            version_index = "mindex_wc1.4"
        } else {

            #if not, then the target is wc 2
            version_index = "mindex_wc2"
        }
    } else{

        #if not, and the min values differ in 400 units
        if(min_vals_differs > 400){

            #the target is the minimum value 
            max_min_index = "min"

            #if the minimun value of wc1.4 is higher
            if(min_values_mindex_wc1.4 > min_values_mindex_wc2){

                #wc 1.4 is the target
                version_index = "mindex_wc1.4"
            } else {
                
                #if not, then the target is wc 2                
                version_index = "mindex_wc2"
            }
        } else {

            #if min nor max values showed differences
            version_index = "none"
        }          
    }

    #version_index is not none, hence min or max values differ
    if(version_index != "none"){

        #copy the target raster (the version that we have to change)
        raster_to_rerange = eval(parse(text=paste("selected_", version_index, sep="")))

        #get the values of the target raster
        values_raster_to_rerange = getValues(raster_to_rerange)

        #if the max value is only included in ONE cell
        if(length(which(values_raster_to_rerange == max(na.omit(values_raster_to_rerange)))) == 1){

            #remove the value of that cell
            raster_to_rerange[which(values_raster_to_rerange == max(na.omit(values_raster_to_rerange)))] <- NA
        }
    }

    #plot the maps
    #if we changed mindex_wc1.4
    if(version_index == "mindex_wc1.4"){

        #plot mindex_wc1.4 as the reranged map and mindex_wc2 normal 
        par(mfcol=c(3,1), mar=c(1, 4, 4, 2) + 0.1)
        plot(raster_to_rerange, main=paste("WorldClim 1.4 ", selected_month, sep=""))   
        plot(selected_mindex_wc2, main=paste("WorldClim 2 ", selected_month, sep=""))
    } else {

        #if not and we changed selected_mindex_wc2
        if(version_index == "mindex_wc2"){

            #plot mindex_wc2 as the reranged map and mindex_wc1.4 normal 
            par(mfcol=c(3,1), mar=c(1, 4, 4, 2) + 0.1)
            plot(selected_mindex_wc1.4, main=paste("WorldClim 1.4 ", selected_month, sep=""))   
            plot(raster_to_rerange, main=paste("WorldClim 2 ", selected_month, sep=""))
        } else {

            #if no version was reranged, plot them without changes
            par(mfcol=c(3,1), mar=c(1, 4, 4, 2) + 0.1)            
            plot(selected_mindex_wc1.4, main=paste("WorldClim 1.4 ", selected_month, sep=""))              
            plot(selected_mindex_wc2, main=paste("WorldClim 2 ", selected_month, sep=""))            
        }   
    }

    #plot the chelsa version
    plot(selected_mindex_chelsa, main=paste("Chelsa ", selected_month, sep=""))
}
dev.off()


###for each month calculate and plot differences in the moisture index from different sources
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/mindex_wc1.4_wc_2_check/difference_wc1.4_wc2_chelsa.pdf")
for(i in 1:12){

    #select the month [i]
    selected_month = month.name[i]

    #print the selected month
    print(selected_month)

    #select the raster of mindex of wc 1.4
    selected_mindex_wc1.4 = raster(paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/mindex_5m_bil/mindex_5m_bil", i, ".tif", sep=""))

    #select the raster of mindex of wc 2
    selected_mindex_wc2 = mindex_wc2_crop[[i]]

    #select the raster of mindex of chelsa
    if(i %in% 1:9){
        selected_mindex_chelsa = raster(paste("/Volumes/Maxtor/diego/science_big_documents/seed_mass_aridity/really_final_variables/mind_chelsa/unzipped_10x10km/10x10km_CHELSA_mind_0", i, "_V1.2_land.tif", sep=""))
    }else{
        selected_mindex_chelsa = raster(paste("/Volumes/Maxtor/diego/science_big_documents/seed_mass_aridity/really_final_variables/mind_chelsa/unzipped_10x10km/10x10km_CHELSA_mind_", i, "_V1.2_land.tif", sep=""))      
    }

    #calculate rasters of differences
    w14_wc2 = selected_mindex_wc1.4 - selected_mindex_wc2
    chelsa_14 = selected_mindex_chelsa - selected_mindex_wc1.4    
    chelsa_w2 = selected_mindex_chelsa - selected_mindex_wc2    

    #Get values of cells without NAs for WC versions
    values_mindex_w14_wc2 = getValues(w14_wc2)
    values_mindex_chelsa_14 = getValues(chelsa_14)    
    values_mindex_chelsa_w2 = getValues(chelsa_w2)

    #extract min and max values of cells in moisture index of wc1.4 and wc2
    max_values_values_mindex_w14_wc2 = max(na.omit(values_mindex_w14_wc2))
    min_values_values_mindex_w14_wc2 = min(na.omit(values_mindex_w14_wc2))
    max_values_values_mindex_chelsa_14 = max(na.omit(values_mindex_chelsa_14))     
    min_values_values_mindex_chelsa_14 = min(na.omit(values_mindex_chelsa_14))    
    max_values_values_mindex_chelsa_w2 = max(na.omit(values_mindex_chelsa_w2))
    min_values_values_mindex_chelsa_w2 = min(na.omit(values_mindex_chelsa_w2))
    
    #print the differences in min and max values between wc versions
    range_w14_wc2 = abs(max_values_values_mindex_w14_wc2 - min_values_values_mindex_w14_wc2)
    range_chelsa_14 = abs(max_values_values_mindex_chelsa_14 - min_values_values_mindex_chelsa_14)
    range_chelsa_w2 = abs(max_values_values_mindex_chelsa_w2 - min_values_values_mindex_chelsa_w2)    
    print(paste("Range for differences between w14_wc2 ", round(range_w14_wc2, 2), sep=""))
    print(paste("Range for differences between chelsa_14 ", round(range_chelsa_14, 2), sep=""))
    print(paste("Range for differences between chelsa_w2 ", round(range_chelsa_w2, 2), sep=""))

    #indicate what version and the extreme we have to change
    version_index = NULL
    max_min_index = NULL      

    #if the differences between max values are higher than 400    
    if(range_w14_wc2 > 1000){

        #range of wc 1.4 - wc2 is the target
        version_index = append(version_index, "w14_wc2")

        #if the maximum value of wc1.4 is higher
        if(abs(max_values_values_mindex_w14_wc2) > abs(min_values_values_mindex_w14_wc2)){

            #max is the target
            max_min_index = append(max_min_index, "max")
        } else {

            #in is the target
            max_min_index = append(max_min_index, "min")
        }        
    }
    
    #if the differences between max values are higher than 400    
    if(range_chelsa_14 > 1000){

        #range of chelsa - wc 1.4 is the target
        version_index = append(version_index, "chelsa_14")

        #if the maximum value of wc1.4 is higher
        if(abs(max_values_values_mindex_chelsa_14) > abs(min_values_values_mindex_chelsa_14)){

            #max is the target
            max_min_index = append(max_min_index, "max")
        } else {

            #in is the target
            max_min_index = append(max_min_index, "min")
        }         
    }

    #if the differences between max values are higher than 400    
    if(range_chelsa_w2 > 1000){

        #range of chelsa - wc2 is the target
        version_index = append(version_index, "chelsa_w2")

        #if the maximum value of wc1.4 is higher
        if(abs(max_values_values_mindex_chelsa_w2) > abs(min_values_values_mindex_chelsa_w2)){

            #max is the target
            max_min_index = append(max_min_index, "max")
        } else {

            #in is the target
            max_min_index = append(max_min_index, "min")
        }           
    } 
  
    #if we have rasters to be reranged
    if(length(version_index) != 0){
        #for each version that has to be reranged
        for(j in 1:length(version_index)){

            #select the version [i]
            selected_version = version_index[j]

            #selected range extreme (min or max)
            selected_extreme = max_min_index[j]

            #copy the target raster (the version that we have to change)
            raster_to_rerange = eval(parse(text=paste(selected_version, sep="")))

            #get the values of the target raster
            values_raster_to_rerange = getValues(raster_to_rerange)

            #if the extreme value is only included in ONE cell
            if(selected_extreme == "max"){
                if(length(which(values_raster_to_rerange == max(na.omit(values_raster_to_rerange)))) == 1){

                    #remove the value of that cell
                    raster_to_rerange[which(values_raster_to_rerange == max(na.omit(values_raster_to_rerange)))] <- NA
                }
            } else{
                if(selected_extreme == "min"){
                    if(length(which(values_raster_to_rerange == min(na.omit(values_raster_to_rerange)))) == 1){

                        #remove the value of that cell
                        raster_to_rerange[which(values_raster_to_rerange == min(na.omit(values_raster_to_rerange)))] <- NA
                    }                
                }
            } 

            #assign the new reranged raster to the correspoding version
            assign(paste(selected_version, sep=""), raster_to_rerange)
        }
    }    

    #plot the maps
    par(mfcol=c(3,1), mar=c(1, 4, 4, 2) + 0.1)            
    plot(w14_wc2, main=paste("WorldClim 1.4 - WorldClim 2 ", selected_month, sep="")) 
    plot(chelsa_14, main=paste("Chelsa - WorldClim 1.4 ", selected_month, sep=""))
    plot(chelsa_w2, main=paste("Chelsa - WorldClim 2 ", selected_month, sep=""))    
}
dev.off()

#########################################
######comparison between bioclimatic variables

###load mindex wc1.4
list_rasters_mindex_wc14 = list.files("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/mindex_5m_bil", full.names = TRUE)
#reoder in numeric order (1,2,3...)
library(gtools)
list_rasters_mindex_wc14 = mixedsort(sort(list_rasters_mindex_wc14))
#stack
mindex_wc14 = stack(list_rasters_mindex_wc14)
nlayers(mindex_wc14) == 12
names(mindex_wc14) #month names are not in order

###load temperature
#unzip the monthly tmax and tmin from wc.14 (5 min resolution)
system("unzip /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmax_5m_bil.zip -d /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmax_5m_bil/")
system("unzip /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmin_5m_bil.zip -d /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmin_5m_bil/")

#load tmax
list_rasters_tmax_wc14 = list.files("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmax_5m_bil", full.names = TRUE, pattern=".bil")
#reoder in numeric order (1,2,3...)
library(gtools)
list_rasters_tmax_wc14 = mixedsort(sort(list_rasters_tmax_wc14))
#stack
tmax_wc14 = stack(list_rasters_tmax_wc14)
nlayers(tmax_wc14) == 12
names(tmax_wc14) #month names are not in order

#load tmin
list_rasters_tmin_wc14 = list.files("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmin_5m_bil", full.names = TRUE, pattern=".bil")
#reoder in numeric order (1,2,3...)
library(gtools)
list_rasters_tmin_wc14 = mixedsort(sort(list_rasters_tmin_wc14))
#stack
tmin_wc14 = stack(list_rasters_tmin_wc14)
nlayers(tmin_wc14) == 12
names(tmin_wc14) #month names are not in order

#create the biolcim variables of moisture with "biovars"
require(dismo)
bio_moisture = biovars(prec = mindex_wc14, #vector, matrix, or RasterStack/Brick of precipitation data. In our case is moisture index of each month for the whole globe
    tmin = tmin_wc14, #vector, matrix, or RasterStack/Brick of minimum temperature data
    tmax = tmax_wc14) #vector, matrix, or RasterStack/Brick of maximum temperature data

#problem bio15: 
plot(bio_moisture[["bio15"]]) #extrange results.
#bio15 is the Precipitation Seasonality (Coefficient of Variation), in our case moisture seasonality. Coefficient of Variation is sd/mean (https://en.wikipedia.org/wiki/Coefficient_of_variation), and it makes sense use it for precipitation, because this variable can not exhibit negative values (min value always is zero). However, our variable is moisture index, that can exhibit negative values, which indicate that evapotanspiration is higher than precipitation. In this case, the coefficiente of variation would be very negative. We will use tha same way than bioclim use for temperaturature seasonality, use only standart deviation. 
bio_moisture[["bio15"]] = calc(mindex_wc14, sd) #use calc for calculate the sd of the raster of the stack, and change bio15 form biovars by our variable. We don't multiply by 100, as Hiijmans make in Temperature seasonality, this is a method to use less memory. Rafi told that it is nor neccesary. 
plot(bio_moisture[["bio15"]])

#create a directory to save biovariables
system("rm -rf /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/biovars_5m_bil; mkdir /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/biovars_5m_bil")

#save rasters
writeRaster(x = bio_moisture, filename="/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/biovars_5m_bil/biovars_5m_bil", suffix=names(bio_moisture), bylayer=TRUE, format="ascii", overwrite=TRUE)

#remove folders with tmin and tmax. We only leave the zip files
system("rm -rf /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmax_5m_bil")
system("rm -rf /Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/tmin_5m_bil")

#required for friendly paletes
require(RColorBrewer)

###for each month calculate and plot differences in the biovariables from different sources
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/mindex_wc1.4_wc_2_check/difference_biovars_wc1.4_finals.pdf")
for(i in c(1:7, 10:17)){

    #select the month [i]
    selected_biovar = paste("bio", i)

    #print the selected month
    print(selected_biovar)

    #select the raster of biovars of wc 1.4
    selected_biovar_wc1.4 = crop(raster(paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/biovars_5m_bil/biovars_5m_bil_bio", i, ".asc", sep="")), c(-180, 180, -5, 90))

    #select the raster of biovars of finals_v0 (bioclim with errors)
    selected_bioclim_finals = crop(raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals_v0/bio", i, ".asc", sep="")), c(-180, 180, -5, 90))

    #select the raster of biovars of chelsa
    selected_bioclim_chelsa = crop(raster(paste("/Volumes/Maxtor/diego/science_big_documents/seed_mass_aridity/really_final_variables/bioclim_variables_chelsa/unzipped_10x10km/10x10km_bioclim_variables_bio", i, ".tif", sep="")), c(-180, 180, -5, 90))

    #calculate rasters of differences
    w14_finals = selected_biovar_wc1.4 - selected_bioclim_finals
    chelsa_14 = selected_bioclim_chelsa - selected_biovar_wc1.4    
    chelsa_finals = selected_bioclim_chelsa - selected_bioclim_finals    

    #Get values of cells without NAs for WC versions
    values_biovars_w14_finals = getValues(w14_finals)
    values_biovars_chelsa_14 = getValues(chelsa_14)    
    values_biovars_chelsa_finals = getValues(chelsa_finals)

    #extract min and max values of cells in moisture index of wc1.4 and finals
    max_values_values_biovars_w14_finals = max(na.omit(values_biovars_w14_finals))
    min_values_values_biovars_w14_finals = min(na.omit(values_biovars_w14_finals))
    max_values_values_biovars_chelsa_14 = max(na.omit(values_biovars_chelsa_14))     
    min_values_values_biovars_chelsa_14 = min(na.omit(values_biovars_chelsa_14))    
    max_values_values_biovars_chelsa_finals = max(na.omit(values_biovars_chelsa_finals))
    min_values_values_biovars_chelsa_finals = min(na.omit(values_biovars_chelsa_finals))
    
    #print the differences in min and max values between wc versions
    range_w14_finals = abs(max_values_values_biovars_w14_finals - min_values_values_biovars_w14_finals)
    range_chelsa_14 = abs(max_values_values_biovars_chelsa_14 - min_values_values_biovars_chelsa_14)
    range_chelsa_finals = abs(max_values_values_biovars_chelsa_finals - min_values_values_biovars_chelsa_finals)    
    print(paste("Range for differences between w14_finals ", round(range_w14_finals, 2), sep=""))
    print(paste("Range for differences between chelsa_14 ", round(range_chelsa_14, 2), sep=""))
    print(paste("Range for differences between chelsa_finals ", round(range_chelsa_finals, 2), sep=""))

    #indicate what version and the extreme we have to change
    version_index = NULL
    max_min_index = NULL      

    #if the differences between max values are higher than 400    
    if(range_w14_finals > 1000){

        #range of wc 1.4 - finals is the target
        version_index = append(version_index, "w14_finals")

        #if the maximum value of wc1.4 is higher
        if(abs(max_values_values_biovars_w14_finals) > abs(min_values_values_biovars_w14_finals)){

            #max is the target
            max_min_index = append(max_min_index, "max")
        } else {

            #in is the target
            max_min_index = append(max_min_index, "min")
        }        
    }
    
    #if the differences between max values are higher than 400    
    if(range_chelsa_14 > 1000){

        #range of chelsa - wc 1.4 is the target
        version_index = append(version_index, "chelsa_14")

        #if the maximum value of wc1.4 is higher
        if(abs(max_values_values_biovars_chelsa_14) > abs(min_values_values_biovars_chelsa_14)){

            #max is the target
            max_min_index = append(max_min_index, "max")
        } else {

            #in is the target
            max_min_index = append(max_min_index, "min")
        }         
    }

    #if the differences between max values are higher than 400    
    if(range_chelsa_finals > 1000){

        #range of chelsa - finals is the target
        version_index = append(version_index, "chelsa_finals")

        #if the maximum value of wc1.4 is higher
        if(abs(max_values_values_biovars_chelsa_finals) > abs(min_values_values_biovars_chelsa_finals)){

            #max is the target
            max_min_index = append(max_min_index, "max")
        } else {

            #in is the target
            max_min_index = append(max_min_index, "min")
        }           
    } 
  
    #if we have rasters to be reranged
    if(length(version_index) != 0){
        #for each version that has to be reranged
        for(j in 1:length(version_index)){

            #select the version [i]
            selected_version = version_index[j]

            #selected range extreme (min or max)
            selected_extreme = max_min_index[j]

            #copy the target raster (the version that we have to change)
            raster_to_rerange = eval(parse(text=paste(selected_version, sep="")))

            #get the values of the target raster
            values_raster_to_rerange = getValues(raster_to_rerange)

            #if the extreme value is only included in ONE cell
            if(selected_extreme == "max"){
                if(length(which(values_raster_to_rerange == max(na.omit(values_raster_to_rerange)))) == 1){

                    #remove the value of that cell
                    raster_to_rerange[which(values_raster_to_rerange == max(na.omit(values_raster_to_rerange)))] <- NA
                }
            } else{
                if(selected_extreme == "min"){
                    if(length(which(values_raster_to_rerange == min(na.omit(values_raster_to_rerange)))) == 1){

                        #remove the value of that cell
                        raster_to_rerange[which(values_raster_to_rerange == min(na.omit(values_raster_to_rerange)))] <- NA
                    }                
                }
            } 

            #assign the new reranged raster to the correspoding version
            assign(paste(selected_version, sep=""), raster_to_rerange)
        }
    }    

    #plot the maps
    par(mfcol=c(3,1), mar=c(1, 4, 4, 2) + 0.1)            
    plot(w14_finals, main=paste("WorldClim 1.4 - WorldClim 2 (used in modelling): ", selected_biovar, sep=""), col=brewer.pal(n=9, name="PuOr")) 
    plot(chelsa_14, main=paste("Chelsa - WorldClim 1.4: ", selected_biovar, sep=""), col=brewer.pal(n=9, name="PuOr"))
    plot(chelsa_finals, main=paste("Chelsa - WorldClim 2 (used in modelling): ", selected_biovar, sep=""), col=brewer.pal(n=9, name="PuOr"))    
}
dev.off()


#####comparison of the bioclimatic variables created in this code and the final used in the analyses, which were obviously copied from these and changed the names #######

#for each month plot
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/mindex_wc1.4_wc_2_check/comparison_new_bio_variables_and_final_used.pdf")
par(mfcol=c(2,2))
for(i in 1:19){

    #load the initial raster of [i] month
    initial_raster = raster(paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/biovars_5m_bil/biovars_5m_bil_bio", i, ".asc", sep=""))
    
    #load the final raster of [i] month (final version without errors)
    new_raster = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio", i, ".asc", sep=""))
    
    #plot
    plot(initial_raster, main=paste("bio ", i, sep=""))
    plot(new_raster, main=paste("bio ", i, sep=""))
    plot(initial_raster - new_raster, main=paste("difference between maps bio ", i, sep=""))
    frame()
}
dev.off()

#plot bioclimatic variables into one single pdf
pdf("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/biovars_5m_bil/bioclim_variables_plotted.pdf")
for(i in 1:19){

    #load the raster
    bio_clim_raster = raster(paste("/Volumes/Maxtor/diego/science_big_documents/pines_niche/climate/mindex_wc1.4_wc_2_check/biovars_5m_bil/biovars_5m_bil_bio", i, ".asc", sep=""))

    #plot raster
    plot(bio_clim_raster, main=paste("bio ", i, " 10x10km", sep=""))
}
dev.off()