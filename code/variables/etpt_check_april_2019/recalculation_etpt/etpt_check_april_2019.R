#library
require(raster)        

#load a function to calculate the number of days of a month (variable according to leap years or not). Number of days will be used in the calculation of moisture index
days.in.month <- function(month, year = NULL){
    month = as.integer(month)
    if (is.null(year))
        year = as.numeric(format(Sys.Date(), '%Y'))
        dt = as.Date(paste(year, month, '01', sep = '-'))
        dates = seq(dt, by = 'month', length = 2)
        as.numeric(difftime(dates[2], dates[1], units = 'days'))
} # from: http://stackoverflow.com/questions/6243088/find-out-the-number-of-days-of-a-month-in-r

#load the function to calculate evapotranspiration vectorized 
overfun <- Vectorize(function(x, y, z) {
    ifelse(z < 0, 0, (x / (x + 15)) * ((0.0239001 * y) + 50.0) * (4.0 / 30.0))
}) #selected from Rafi (see "/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/climate/moisture/wc2.1_5m_etpt/etpt.r") 

#load rasters of temperature for each month
tavg_list = list()#create an empty list to save rasters of tav, sr and pp. This list will be transformed into a stack then
sr_list = list()
pp_list = list()
etpt_list = list()
mind_list = list()

#open pdf to save figure womparing two different transformations of etpt
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/etpt_check_april_2019/recalculation_etpt/etpt_checks_current_conditions.pdf")
for(i in 1:12){

    ### load min and max temperature and calcualte average temperature
    #load temperature variables of the [i] month
    tmax = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/etpt_check_april_2019/recalculation_etpt/data/tmax_10m_bil/tmax", i, ".bil", sep=""))
    tmin = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/etpt_check_april_2019/recalculation_etpt/data/tmin_10m_bil/tmin", i, ".bil", sep=""))

    #calculate average temperature of the [i] month
    tavg = overlay(tmin, tmax, fun  = function(x) sum(x)/2)

    #set the name of the [i] layer
    names(tavg) <- month.name[i]

    #save average temperature in a stack
    tavg_list[[i]] = tavg


    ###load solar radiation
    #open solar radiation of the [i] month
    if(i %in% 1:9){
        sr = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/etpt_check_april_2019/recalculation_etpt/data/wc2.0_10m_srad/wc2.0_10m_srad_0", i, ".tif", sep=""))
    } else {
        sr = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/etpt_check_april_2019/recalculation_etpt/data/wc2.0_10m_srad/wc2.0_10m_srad_", i, ".tif", sep=""))
    }   

    #crop solar radiation with ta raster
    sr = crop(sr, tavg)

    #set the name of the [i] layer
    names(sr) <- month.name[i] 

    #save average temperature in a stack
    sr_list[[i]] = sr


    ###calculate etpt 
    #calculate etpt for ONE day of january
    etpt <- overlay(tavg, sr, tavg, fun = overfun)

    #set the name of the [i] layer
    names(etpt) <- month.name[i] 

    #save average temperature in a stack
    etpt_list[[i]] = etpt


    ###load precipitation
    #load pp of the [i] month
    pp = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/etpt_check_april_2019/recalculation_etpt/data/prec_10m_bil/prec", i, ".bil", sep=""))

    #set the name of the [i] layer
    names(pp) <- month.name[i]

    #save average temperature in a stack
    pp_list[[i]] = pp 


    ###calcualte moisture index
    #calculate the number of the days of the corresponding month
    days <- days.in.month(i, year = 2070)

    #calculate moisture index
    mind <- pp - etpt / 10 * days

    #set the name of the [i] layer
    names(mind) <- month.name[i] 

    #save average temperature in a stack
    mind_list[[i]] = mind

    #check that units of etpt are correct. Compare pp with (etpt / 10 * days) and (etpt * days)
    par(mfcol=c(2,2), mar=c(5, 2, 4, 6) + 0.1)
    plot(pp, main=paste(month.name[i], " Precipitation", sep=""))
    plot(etpt / 10 * days, , main=paste(month.name[i], " etpt [etpt / 10 * days]", sep=""))
    plot(pp, main=paste(month.name[i], " Precipitation", sep=""))
    plot(etpt * days, main=paste(month.name[i], " etpt [etpt / days]", sep=""))        
}
dev.off()

#check that each layer have the month name in order
for(i in 1:12){
    print(names(tavg_list[[i]]) == month.name[i])     
    print(names(sr_list[[i]]) == month.name[i])     
    print(names(pp_list[[i]]) == month.name[i])    
    print(names(etpt_list[[i]]) == month.name[i])     
    print(names(mind_list[[i]]) == month.name[i])
}

#convert the lists as stacks
tavg_stack = stack(tavg_list)
sr_stack = stack(sr_list)
pp_stack = stack(pp_list)
etpt_stack = stack(etpt_list)
mind_stack = stack(mind_list)

### calculate the sum of mindex across months (would be equivalent to bio12)
sum_mindex = calc(mind_stack, function(x) sum(x))
sum_mindex
plot(sum_mindex)
    #we compared this with BIO12 from chelsa ("GoogleDrive/My Drive/science/phd/nicho_pinus/datos/comparison_units_niche_seed_mass_papers")

### calculate the standard deviation of moisture index, which is equivalent to bio15. In that way we check that the calculation of this biovariable is correct. 
#bio15 is the Precipitation Seasonality (Coefficient of Variation), in our case moisture seasonality. Coefficient of Variation is sd/mean (https://en.wikipedia.org/wiki/Coefficient_of_variation), and it makes sense use it for precipitation, because this variable can not exhibit negative values (min value always is zero). However, our variable is moisture index, that can exhibit negative values, which indicate that evapotanspiration is higher than precipitation. In this case, the coefficiente of variation would be very negative. We will use tha same way than bioclim use for temperaturature seasonality, use only standart deviation
coef_var_mindex = calc(mind_stack, function(x) sd(x))
coef_var_mindex
plot(coef_var_mindex)#everything ok.
    #this is we are doing: sd(c(-1500, -1000, -200, 500, 700, 1000, 1200)) = 1032.796
    #with variation coefficente would be: sd(c(-1500, -1000, -200, 500, 700, 1000, 1200))/mean(c(-1500, -1000, -200, 500, 700, 1000, 1200)) = 10.32796

#plot it
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/etpt_check_april_2019/recalculation_etpt/sum_moisture_wc1_4_current_across_months.pdf")
par(mfcol=c(2,2), mar=c(5, 2, 4, 6) + 0.1)
plot(sum_mindex, main="sum mindex wc1.4 current low resolution")
plot(raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio12.asc"), main="bio12 wc1.2 used in paper niche (high resolution)")
dev.off()#units are more or less similar, altouh the maximum value of sum_mindex is the lower resolution (0.1666667 instead of 0.08333333)

#save work space
save.image("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/rdata/etpt_check_april_2019.RData")
require(raster)