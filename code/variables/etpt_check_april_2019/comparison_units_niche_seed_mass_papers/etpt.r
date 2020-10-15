#code written by rafi for creating evapotranspiration and moisture variable. If you want explanations, ou can see the doce of future moisture index

library(raster)
ta1<-raster("/Users/zimmerma/Data/WClim2/wc2.0_5m_tavg/wc2.0_5m_tavg_01.tif")
sr1<-raster("/Users/zimmerma/Data/WClim2/wc2.0_5m_srad/wc2.0_5m_srad_01.tif")
et1<-ifelse(ta1<=0,0,(0.4 / 30.0 * ((0.0239001 * sr1) + 50.0) * (ta1 / (ta1 + 15.0)) ))
et1<-ifelse(ta1<0,0,(0.4 / 30.0 * ((0.0239001 * sr1) + 50.0) * (ta1 / 10.0) / (ta1 + 15.0)))
et1<-               (0.4 / 30.0 * ((0.0239001 * B32) + 50.0) * A32 / (A32 + 15.0))
et1<-               (0.4 / 30.0 * ((0.0239001 * sr1) + 50.0) * ((ta1/10) / ((ta1/10) + 15.0))  )

et.calc<- function(x,y){
  ifelse(is.na(x),NA,
  ifelse(x>0,(4.0 / 30.0 * ((0.0239001 * y) + 50.0) * ((x) / ((x) + 15.0))),0))
  } 

rst_test <- function(rst_a,rst_B)
{
return(ifelse(rst_a==0 & rst_B==1,2,0))
}





bio.stk<-brick(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19)

bio.pca<-princomp(bio.stk)

b1.spdf<-as(b1, "SpatialPixelsDataFrame")
b1.strt<-spsample(b1.spdf,n=50,zxpr




b1.spdf<-as(small[[1]], "SpatialPixelsDataFrame")


### Rafi

require(raster)
require(foreach)
require(doParallel)

ta1<-raster("/Users/wueest/Dropbox/WClim2/wc2.0_5m_tavg/wc2.0_5m_tavg_01.tif")
sr1<-raster("/Users/wueest/Dropbox/WClim2/wc2.0_5m_srad/wc2.0_5m_srad_01.tif")


(4.0 / 30.0 * ((0.0239001 * y) + 50.0) * ((x) / ((x) + 15.0)))

overfun1 <- function(x, y) {
	(x / (x + 15)) * ((0.0239001 * y) + 50.0) * (4.0 / 30.0)
}

tst1 <- overlay(ta1, sr1, fun = overfun1)

overfun2 <- function(x, y) {
	ifelse(x < 0, 0, (x / (x + 15)) * ((0.0239001 * y) + 50.0) * (4.0 / 30.0))
}

tst2 <- overlay(ta1, sr1, fun = overfun2)

overfun3 <- function(x, y, z) {
	ifelse(z < 0, 0, (x / (x + 15)) * ((0.0239001 * y) + 50.0) * (4.0 / 30.0))
}

tst3 <- overlay(ta1, sr1, ta1, fun = overfun3)
# => that works!

## run across 12 moths
# preps
months <- formatC(1:12, width = 2, flag = '0')
overfun <- Vectorize(function(x, y, z) {
	ifelse(z < 0, 0, (x / (x + 15)) * ((0.0239001 * y) + 50.0) * (4.0 / 30.0))
})
fefun <- function(m){
	ta <- raster(paste('/Users/wueest/Dropbox/WClim2/wc2.0_5m_tavg/wc2.0_5m_tavg_', m, '.tif', sep = ''))
	sr <- raster(paste('/Users/wueest/Dropbox/WClim2/wc2.0_5m_srad/wc2.0_5m_srad_', m, '.tif', sep = ''))
	overlay(ta, sr, ta, fun = overfun)
}
# set up cluster
clust <- makeCluster(12)
registerDoParallel(clust)
# run
etpt <- foreach(i = months, .packages = 'raster', .combine = stack) %dopar% {
	fefun(m = i)
}
names(etpt) <- paste('month', 1:12, sep = '_')
stopCluster(clust)


precfils <- list.files('/Users/wueest/Dropbox/WClim2/wc2.0_5m_prec/', pattern = '.tif$', full.names = TRUE)
prec <- stack(precfils)
days.in.month <- function(month, year = NULL){
	month = as.integer(month)
	if (is.null(year))
		year = as.numeric(format(Sys.Date(), '%Y'))
	dt = as.Date(paste(year, month, '01', sep = '-'))
	dates = seq(dt, by = 'month', length = 2)
	as.numeric(difftime(dates[2], dates[1], units = 'days'))
} # from: http://stackoverflow.com/questions/6243088/find-out-the-number-of-days-of-a-month-in-r
days <- sapply(1:12, days.in.month, year = 2015)
mind <- brick(lapply(1:12, FUN = function(x){
	prec[[x]] - etpt[[x]] / 10 * days[x]
}))
names(mind) <- paste('month', 1:12, sep = '_')

## export
writeRaster(etpt, filename = '/Users/wueest/Dropbox/WClim2/wc2.1_5m_etpt/etpt.grd')
writeRaster(mind, filename = '/Users/wueest/Dropbox/WClim2/wc2.1_5m_mind/mind.grd')

