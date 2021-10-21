#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
    #https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
    #https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



######################################################################
################ PLOTTING FINAL GLOBAL FIGURES #######################
######################################################################

#This script plots the global figures for SDMs.



###################################################
##### DIFFERENCES RESPECT TO PREVIOUS VERSION #####
###################################################

#Respect to version 2:
    #I have changed the path for running the analyses in the laptop of David (msi)

    #I have created the map plots with cairo pdf because the symbol delta was not properly displayed.

    #I HAVE NOT revised the entire script after these changes, only I checked that the global figures between versions (2 and 3) are the same.

    #I had to install raster, for that I had to install several ubuntu libraries
        #sudo apt install libgeos-dev #https://stackoverflow.com/questions/53389181/installing-the-r-package-rgeos-on-linux-geos-config-not-found-or-not-executab
        #sudo apt install libgdal-dev #https://stackoverflow.com/questions/12141422/error-gdal-config-not-found-while-installing-r-dependent-packages-whereas-gdal
        #sudo apt-get install libudunits2-dev

#Respect to version 3
    # I have changed the paths for the new computes
    # I have changed the SD by interquartile range to show dispersion around the global median of range loss and change across Pinus.
    # I have not check anything else of this script between version 3 and 4.


########################
##### BEGIN SCRIPT #####
########################

#set working directory
setwd("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus")

#required packages
require(raster)
require(RColorBrewer)

#environment variable for using it as a background
bio1 = raster("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/datos/finals/bio1.asc")
clay = raster("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/datos/finals/clay.asc")
environment_var = bio1*clay
environment_var[which(getValues(environment_var) >= min(getValues(environment_var),na.rm = TRUE))] <- 0 #set all continent areas as 0

#list species
list_species = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/code/presences/species.txt", sep="\t", header=T)
str(list_species)
summary(list_species)

#extract epithet from species list
epithet_species_list = NULL
for(i in 1:nrow(list_species)){

    #selected species
    selected_species = as.vector(list_species[i,])

    #extract epithet
    epithet_species_list = append(epithet_species_list, strsplit(selected_species, split=" ")[[1]][2])
}
summary(is.na(epithet_species_list)) #all false
#check
if(FALSE){
    require(tidyverse)
    paste("Pinus", epithet_species_list, sep=" ") == str_trim(as.vector(list_species[,1])) #the seperated epithet more Pinus are equal to list_species? We use str_trim from tidyverse to remove end spaces in each element of list_species
}#it is in false because loading tidyverse lead to load several packages that have a function name "extract", and this gives problems with the extract function of raster. If you want to check run these lines manually

#remove tecunumanii, jaliscana y discolor. These species are not used for Niche paper. The two first because we should downloadad gbif data now, so we would mix gbif data form 2016 and 2019. The third one was impossible to differentiate from P. cembroides
epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("tecunumanii", "jaliscana", "discolor"))]
#check
c("tecunumanii", "jaliscana", "discolor") %in% epithet_species_list

##load results from sinthesis_figures_v2 in Rafa-pro
current_suit_stack_sum = raster("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/initial_global_figures/current_suit_stack_sum.asc") 
projected_suit_phylo_stack_sum = raster("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/initial_global_figures/projected_suit_phylo_stack_sum.asc") 
projected_suit_stack_sum = raster("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/initial_global_figures/projected_suit_stack_sum.asc")
sum_distributions_sum = raster("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/initial_global_figures/sum_distributions_sum.asc")
raster_range_calc_stack_sum = raster("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/initial_global_figures/raster_range_calc_stack_sum.asc")

#convert raster_range_calc_stack_sum to zero-one
raster_range_calc_stack_sum[which(getValues(raster_range_calc_stack_sum) > 0)] <- 1
#convert to polygon
polygon_range_calc_stack_sum = rasterToPolygons(raster_range_calc_stack_sum, fun=function(x){x==1}, n=16, dissolve = TRUE) #convert to a polygon

#use that polygon for masking predictions rasters and remove areas outside the buffers calc ranges
current_suit_stack_sum = mask(current_suit_stack_sum, polygon_range_calc_stack_sum)
projected_suit_phylo_stack_sum = mask(projected_suit_phylo_stack_sum, polygon_range_calc_stack_sum)
projected_suit_stack_sum = mask(projected_suit_stack_sum, polygon_range_calc_stack_sum)

#use env_var to remove water bodies inside the calc_range_ buffer. Remember, that for each speices you set as zero all areas of the glob with zero suitaiblity or NA to sum rasters of all species. Previous suitaiblity rasters per species have no suitaiblity in water bodies as they were masked with bio1*clay When we maske with calc_range_, we remove areas otuside the buffer, but water bodeis remain inside.
current_suit_stack_sum = mask(current_suit_stack_sum, environment_var)
projected_suit_phylo_stack_sum = mask(projected_suit_phylo_stack_sum, environment_var)
projected_suit_stack_sum = mask(projected_suit_stack_sum, environment_var)
    #after mask with env_var, sea border are changed. They are more precise. This is becasue env_var has more resolution than the raster from which polygon_range_calc_stack_sum came

#we remove these lines form the script because we don't want to include as zero terrestrial areas where we did not predict pine richness under current conditions. We only want zero in areas where models predict pine richness to be zero. The current pine richness map is a prediction, is not the sum of the actual distribution of pines. We use this to compare future model predictions with model predictions under current conditions (avoding mixing expert maps with models maps)
if(FALSE){
    #remove from environment_var_new, areas with data of current predicted pine richness. This raster will be used to include continental areas with no pines
    environment_var_new = environment_var
    environment_var_new[which(getValues(current_suit_stack_sum) >=0)] <- NA

    #set as zero areas from current_suit_stack_sum that have env data 
    current_suit_stack_sum = current_suit_stack_sum
    current_suit_stack_sum[which(!is.na(getValues(environment_var_new)))] <- 0
}


#calculate difference in suitability with and without phylogenetic correction
differ_suit_phylo = projected_suit_phylo_stack_sum - current_suit_stack_sum
differ_suit_no_phylo = projected_suit_stack_sum - current_suit_stack_sum

#change in future richnnes between applying or not the phylogenetic correction
phylo_difference = projected_suit_phylo_stack_sum - projected_suit_stack_sum

#remove zeros from current_suit_stack_sum 
current_suit_stack_sum_no_zeros = current_suit_stack_sum
current_suit_stack_sum_no_zeros[which(getValues(current_suit_stack_sum_no_zeros) == 0)] <- NA

##################
##### MAPS #######
##################

#set extent of the raster to be plotted
plot_extent = c(-180,180,-10,90) #if you change the extent of the plot, some sea border could change (only difference of 1 cell). For example, when you plot the whole globe, save the plot as pdf and then zoom to the Canary Islands, the shape of the island is not very accurate. If you crop the maps to P. canariensis buffer, then the shape of the islands in the pdf is more similar to the reality.

###GENERAL TIP FOR PLOTTING WITH RASTER
#breaks indicate the number of partitions between colors, whilst "at" indicate the numbers in the legend. The number of colors have to be EQUAL to the number o breaks. Breaks should encompass the RANGE of values of the raster


######################
##### PALLETES #######
######################
require(RColorBrewer)
#color palletes for divergence (gain_loss richness under future) and sequence (current richness).
#We selected from Colorbrewer a a single hue pallete with green and a divergent pallete (blue - red with white in the middle).
mypalette_sequence <-brewer.pal(9,"Greens") 
mypalette_divergence <-brewer.pal(11,"RdBu")
    #Names taken from "http://colorbrewer2.org/#type=sequential&scheme=Greens&n=9"
    #All works for anomalous trychromacy and dychromacy ("http://www.color-blindness.com/coblis-color-blindness-simulator/")
#these palletes are used in colorRampPalette to create a function that can create a great number of colors 
colfunc_sequence <- colorRampPalette(mypalette_sequence)
colfunc_divergence <- colorRampPalette(mypalette_divergence)


###########################################
####### current_future_sui_no_phylo #######
###########################################

##plot
cairo_pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/final_global_figures/current_future_suit_v3.pdf", width=12, height=12) #we have to use cairo_pdf. In Ubuntu 18.04, the delta symbol is not correctly showed.
par(mfcol=c(2,1), mai=c(0,0.4,0,1.1), oma=c(0,0,2,1))

##first pannel
#background
plot(crop(environment_var, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, axes=FALSE, box=FALSE, main="") #higher values in argument start of gray colors lead to brighter gray

#add title of the plot
mtext(text="Predicted pine richness under current conditions", side=3, line=-1, outer=FALSE, cex=2, font = 2)

#sum_distributions
plot(crop(current_suit_stack_sum_no_zeros, plot_extent), add=TRUE, axes=FALSE, box=FALSE, col=colfunc_sequence(281), axis.args=list(at=seq(0,28,4), cex.axis=1.3), legend.shrink=0.8, breaks=seq(0,28,0.1), legend=TRUE, legend.args=list(text='Pine richness', side=4, font=2, line=3.2, cex=1.7)) #rev invert vector of colors. "At" argument in axis.args indicate the numbers in the legend

##second pannel
#background
plot(crop(environment_var, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, axes=FALSE, box=FALSE, main="")

#add title of the plot
mtext(text="Predicted richness differences under future conditions", side=3, line=-1, outer=FALSE, cex=2, font = 2)

#sum suitability across species
plot(crop(differ_suit_phylo, plot_extent), add=TRUE, col=colfunc_divergence(161), breaks=seq(-8,8,0.1), axes=FALSE, box=FALSE, axis.args=list(at=seq(-8,8,2), cex.axis=1.3), legend=TRUE, legend.shrink=0.8, legend.args=list(text=expression(bold(paste(Delta~'pine richness'))), side=4, font=2, line=3.2, cex=1.7)) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"
dev.off()


##plot
cairo_pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/final_global_figures/future_suit_phylo_no_phylo_v3.pdf", width=12, height=12) #we have to use cairo_pdf. In Ubuntu 18.04, the delta symbol is not correctly showed.
par(mfcol=c(2,1), mai=c(0,0.4,0,1.1), oma=c(0,0,2,1))

##first pannel
#background
plot(crop(environment_var, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, axes=FALSE, box=FALSE, main="")

#add title of the plot
mtext(text="Predicted richness differences under future conditions", side=3, line=-1, outer=FALSE, cex=2, font = 2)
mtext(text="without 'phylogenetic correction'", side=3, line=-3.5, outer=FALSE, cex=2, font = 2)


#sum suitability across species
plot(crop(differ_suit_no_phylo, plot_extent), add=TRUE, col=colfunc_divergence(181), breaks=seq(-9,9,0.1), axes=FALSE, box=FALSE, axis.args=list(at=seq(-9,9,2), cex.axis=1.3), legend=TRUE, legend.shrink=0.8, legend.args=list(text=expression(bold(paste(Delta~'pine richness'))), side=4, font=2, line=3.2, cex=1.7)) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"

##second pannel
#background
plot(crop(environment_var, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, axes=FALSE, box=FALSE, main="")

#add title of the plot
mtext(text="Predicted richness differences under future conditions", side=3, line=-1, outer=FALSE, cex=2, font = 2)
mtext(text="with 'phylogenetic correction'", side=3, line=-3.5, outer=FALSE, cex=2, font = 2)

#sum suitability across species
plot(crop(differ_suit_phylo, plot_extent), add=TRUE, col=colfunc_divergence(181), breaks=seq(-9,9,0.1), axes=FALSE, box=FALSE, axis.args=list(at=seq(-9,9,2), cex.axis=1.3), legend=TRUE, legend.shrink=0.8, legend.args=list(text=expression(bold(paste(Delta~'pine richness'))), side=4, font=2, line=3.2, cex=1.7)) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"
dev.off()


###############################
####### phylo__no_phylo #######
###############################

#color pallete like terrain
colfunc_terrain <- colorRampPalette(c("gray93", terrain.colors(24)[20],terrain.colors(24)[12], terrain.colors(24)[1]))

#see palette
barplot(rep(1,24), col = terrain.colors(24))


##plot
cairo_pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/final_global_figures/phylo_no_phylo_difference_v3.pdf", width=20, height=8) #we have to use cairo_pdf. In Ubuntu 18.04, the delta symbol is not correctly showed.
par(mfcol=c(1,1), mai=c(0,1.1,0,1.1), oma=c(0,0,2,1))


##second pannel
#background
plot(crop(environment_var, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, axes=FALSE, box=FALSE, main="")

#set title
mtext(text=expression(bold("Increase of predicted pine richness with the 'Phylogenetic correction'")), side=3, line=-4, outer=FALSE, cex=3, font = 2)

#sum suitability across species
plot(crop(phylo_difference, plot_extent), add=TRUE, axes=FALSE, box=FALSE, axis.args=list(at=seq(0,9,1), cex.axis=2.3), col=colfunc_sequence(91), breaks=seq(0,9,0.1), legend=TRUE, legend.shrink=0.65, legend.args=list(text=expression(bold(paste(Delta~'pine richness'))), side=4, font=2, line=3.7, cex=2.3)) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"

dev.off()



###############################
######### HISTOGRAMS ##########
###############################

#load the table with range loss
suitability_changes = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/initial_global_figures/suitability_changes.csv", sep=",", header=T)
str(suitability_changes)

###stack histogram for range loss
# 1) Define the breaks to use on your Histogram
xrange_range_loss = seq(0,100,10)

# 2) Have your vectors ready
v1_range_loss = suitability_changes$range_loss_no_phylo #no_filo
v2_range_loss = suitability_changes$range_loss_phylo #filo

# 3) subset your vectors to be inside xrange_range_loss
v1_range_loss = subset(v1_range_loss,v1_range_loss<=max(xrange_range_loss) & v1_range_loss>=min(xrange_range_loss)) #subsetea el vector con los valores de suitability loss con y sin corrección filo, cogiendo todos los valores menores iguales y mayores o iguales que el rango que elegido previamente.
v2_range_loss = subset(v2_range_loss,v2_range_loss<=max(xrange_range_loss) & v2_range_loss>=min(xrange_range_loss)) #the xrange_range_loss covers all data included in v1_range_loss and v2_range_loss (this is the idea)

# 4) Now, use hist to compute the counts per interval
hv1_range_loss = hist(v1_range_loss,breaks=xrange_range_loss,plot=F)$counts #hacemos un histograma con los % de perdida de suit sin filo correction, pero no lo ploteamos. De ahí sacamos los counts, o sea, el número de casos de cada % de perdida, es decir, la altura de las barras: Cuantas especies pierden suitability desde -10 (sin incluirlo) hasta 0 (incluyendolo, right close, left open histogram)
length(which(v1_range_loss>0 & v1_range_loss<=10)) == head(hv1_range_loss, 1) #It gives 8, like the first count, which goes from -10 to 0.
length(which(v1_range_loss>90 & v1_range_loss<=100)) == tail(hv1_range_loss, 1) #It gives 46, like the first count, which goes from 90 (not including it) to 100 (including it). 
hv2_range_loss = hist(v2_range_loss,breaks=xrange_range_loss,plot=F)$counts
length(which(v2_range_loss>0 & v2_range_loss<=10)) == head(hv2_range_loss, 1) #It gives 8, like the first count, which goes from -10 to 0.
length(which(v2_range_loss>90 & v2_range_loss<=100)) == tail(hv2_range_loss, 1) #It gives 44, like the first count, which goes from 90 (not including it) to 100 (including it). 

###stack histogram for range change
# 1) Define the breaks to use on your Histogram
xrange_range_change = seq(-100,100,20)

# 2) Have your vectors ready
v1_range_change = suitability_changes$range_change_no_phylo #no_filo
v2_range_change = suitability_changes$range_change_phylo #filo

# 3) subset your vectors to be inside xrange_range_change
v1_range_change = subset(v1_range_change,v1_range_change<=max(xrange_range_change) & v1_range_change>=min(xrange_range_change)) #subsetea el vector con los valores de suitability loss con y sin corrección filo, cogiendo todos los valores menores iguales y mayores o iguales que el rango que elegido previamente.
v2_range_change = subset(v2_range_change,v2_range_change<=max(xrange_range_change) & v2_range_change>=min(xrange_range_change)) #the xrange_range_change covers all data included in v1_range_change and v2_range_change (this is the idea)

# 4) Now, use hist to compute the counts per interval
hv1_range_change = hist(v1_range_change,breaks=xrange_range_change,plot=F)$counts #hacemos un histograma con los % de perdida de suit sin filo correction, pero no lo ploteamos. De ahí sacamos los counts, o sea, el número de casos de cada % de perdida, es decir, la altura de las barras: Cuantas especies pierden suitability desde -10 (sin incluirlo) hasta 0 (incluyendolo, right close, left open histogram)
length(which(v1_range_change>-100 & v1_range_change<=-80)) == head(hv1_range_change, 1) #It gives 8, like the first count, which goes from -10 to 0.
length(which(v1_range_change>80 & v1_range_change<=100)) == tail(hv1_range_change, 1) #It gives 46, like the first count, which goes from 90 (not including it) to 100 (including it). 
hv2_range_change = hist(v2_range_change,breaks=xrange_range_change,plot=F)$counts
length(which(v2_range_change>-100 & v2_range_change<=-80)) == head(hv2_range_change, 1) #It gives 8, like the first count, which goes from -10 to 0.
length(which(v2_range_change>80 & v2_range_change<=100)) == tail(hv2_range_change, 1) #It gives 44, like the first count, which goes from 90 (not including it) to 100 (including it). 

#5a) Generate a Frequency BarPlot with bars parallels
pdf(file="/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/final_global_figures/range_loss_change_histograms.pdf",width=6,height=8)
par(mfcol=c(2,1), mar=c(5.2,5.2,5.2,5.2), mai=c(1,0.75,0.2,0.1), mgp=c(2.5,1.7,0), cex.main=1.7, cex.axis=1.8)

##range loss
#plot bars
xs_range_loss = barplot(c(rbind(hv1_range_loss, hv2_range_loss)), #counts with and without phylo combinated, the first count wihtout correction, then the first count with, the second count without, second count with... It is to say, both vectors are merged alternating indexes
    col=c("gray39", "gray87"), #two different grays
    #names.arg=xrange[-1], #nombres del eje X, son los breaks
    space=c(0, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0), #espacio cero entre las dos barras de una misma categoria, y luego 0.2 entre cateogiras de perdida de idoneidad. 
    las=0, #axis labels paralelas al eje
    xlab="Range loss (%)", ylab="Nº species", cex.lab=1.5, cex.axis=1.5, main=NULL, lwd=1, font.lab=2, ylim=c(0,60), yaxt='n')

#add label breaks of x axis
mtext(seq(0,100,10), side=1, at=c(xs_range_loss[seq(1,20,2)]-0.5,xs_range_loss[length(xs_range_loss)]+0.5), cex=1, line=0.8) #at the beginning of each bar (a little bit move to left), and another label at the end, because a seq from 1 to 22 each 2 skip 22, thus we added that 22 with length(xs)

#add label breaks of y axis
mtext(seq(0,60,10), side=2, at=seq(0,60,10), cex=1, line=0.8)

#add legend
legend(x=11, y=57, legend=c("Without phylo. correction", "With phylo. correction"), fill=c("gray39", "gray87"), cex=1, horiz = FALSE)

#add box
box(lwd=1)


##range change
#plot bars
xs_range_change = barplot(c(rbind(hv1_range_change, hv2_range_change)), #counts with and without phylo combinated, the first count wihtout correction, then the first count with, the second count without, second count with... It is to say, both vectors are merged alternating indexes
    col=c("gray39", "gray87"), #two different grays
    #names.arg=xrange[-1], #nombres del eje X, son los breaks
    space=c(0, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0), #espacio cero entre las dos barras de una misma categoria, y luego 0.2 entre cateogiras de perdida de idoneidad. 
    las=0, #axis labels paralelas al eje
    xlab="Range change (%)", ylab="Nº species", cex.lab=1.5, cex.axis=1.5, main=NULL, lwd=1, font.lab=2, ylim=c(0,70), yaxt='n')

#add label breaks of x axis
mtext(seq(-100,100,20), side=1, at=c(xs_range_change[seq(1,20,2)]-0.5,xs_range_change[length(xs_range_change)]+0.5), cex=1, line=0.8) #at the beginning of each bar (a little bit move to left), and another label at the end, because a seq from 1 to 22 each 2 skip 22, thus we added that 22 with length(xs)

#add label breaks of y axis
mtext(seq(0,70,10), side=2, at=seq(0,70,10), cex=1, line=0.8)

#add legend
#legend(x=9, y=40, legend=c("Without phylo. correction", "With phylo. correction"), fill=c("gray39", "gray87"), cex=1, horiz = FALSE)

#add box
box(lwd=1)
dev.off()


# 5b) Generate a data.frame with counts and labels

#range loss
#for each break (xrange_range_loss) except the last one 
labels_counts_range_loss = NULL
for(i in 1:(length(xrange_range_loss)-1)){

    #create range label
    range_label_range_loss = paste("'", xrange_range_loss[i], "-", xrange_range_loss[i+1], "'", sep="")#we add '' to avoid that excel considers the column names as number and does estrange things

    #save it
    labels_counts_range_loss = append(labels_counts_range_loss, range_label_range_loss)
}
#check that all except the last one are included
length(labels_counts_range_loss) == length(xrange_range_loss)-1

#bind counts without and with phylogenetic correction
counts_text_range_loss = rbind.data.frame(hv1_range_loss, hv2_range_loss)

#set colnames as labels_counts_range_loss
colnames(counts_text_range_loss) <- labels_counts_range_loss

#set row names as no-phylo and phylo
rownames(counts_text_range_loss) <- c("without_phylo", "with_phylo")
str(counts_text_range_loss)
counts_text_range_loss

#save it
write.table(counts_text_range_loss, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/final_global_figures/hist_counts_range_loss.csv", sep=",", row.names = TRUE, col.names = NA) #we set row.names = TRUE and col.names = NA to get as first column of the csv an empty space and then the two row names. The rest of column have as first cell the range label. If we set col.names=TRUE, then the first column with row names has a range label, and the following columns are wrong set. See "https://stackoverflow.com/questions/2478352/write-table-writes-unwanted-leading-empty-column-to-header-when-has-rownames".


#range change
# 5b) Generate a data.frame with counts and labels

#for each break (xrange_range_change) except the last one 
labels_counts_range_change = NULL
for(i in 1:(length(xrange_range_change)-1)){

    #create range label
    range_label_range_change = paste("'", xrange_range_change[i], "-", xrange_range_change[i+1], "'", sep="")#we add '' to avoid that excel considers the column names as number and does estrange things

    #save it
    labels_counts_range_change = append(labels_counts_range_change, range_label_range_change)
}
#check that all except the last one are included
length(labels_counts_range_change) == length(xrange_range_change)-1

#bind counts without and with phylogenetic correction
counts_text_range_change = rbind.data.frame(hv1_range_change, hv2_range_change)

#set colnames as labels_counts_range_change
colnames(counts_text_range_change) <- labels_counts_range_change

#set row names as no-phylo and phylo
rownames(counts_text_range_change) <- c("without_phylo", "with_phylo")
str(counts_text_range_change)
counts_text_range_change

#save it
write.table(counts_text_range_change, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/final_global_figures/hist_counts_range_change.csv", sep=",", row.names = TRUE, col.names = NA) #we set row.names = TRUE and col.names = NA to get as first column of the csv an empty space and then the two row names. The rest of column have as first cell the range label. If we set col.names=TRUE, then the first column with row names has a range label, and the following columns are wrong set. See "https://stackoverflow.com/questions/2478352/write-table-writes-unwanted-leading-empty-column-to-header-when-has-rownames".



### extract percentage for the paper results
#extract the number of species
#list species
list_species = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/code/presences/species.txt", sep="\t", header=T)
str(list_species)
summary(list_species)

#extract epithet from species list
epithet_species_list = NULL
for(i in 1:nrow(list_species)){

    #selected species
    selected_species = as.vector(list_species[i,])

    #extract epithet
    epithet_species_list = append(epithet_species_list, strsplit(selected_species, split=" ")[[1]][2])
}
summary(is.na(epithet_species_list)) #all false
#check
if(FALSE){
    require(tidyverse)
    paste("Pinus", epithet_species_list, sep=" ") == str_trim(as.vector(list_species[,1])) #the seperated epithet more Pinus are equal to list_species? We use str_trim from tidyverse to remove end spaces in each element of list_species
}#it is in false because loading tidyverse lead to load several packages that have a function name "extract", and this gives problems with the extract function of raster. If you want to check run these lines manually

#remove tecunumanii, jaliscana y discolor. These species are not used for Niche paper. The two first because we should downloadad gbif data now, so we would mix gbif data form 2016 and 2019. The third one was impossible to differentiate from P. cembroides
epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("tecunumanii", "jaliscana", "discolor"))]
#check
c("tecunumanii", "jaliscana", "discolor") %in% epithet_species_list

#species number
number_species = length(epithet_species_list)
number_species == 112

##extract some exact percentages
##range loss
#no phylo
#number of species in that situation
high_lost_no_phylo_range_loss = sum(counts_text_range_loss[which(rownames(counts_text_range_loss)=="without_phylo"), which(colnames(counts_text_range_loss) %in% c("'80-90'", "'90-100'"))])
#percentage respect to the total number of species
print((high_lost_no_phylo_range_loss*100)/number_species) #112 species is the 100%, high_lost_no_phylo_range_loss will be X; (100*high_lost_no_phylo_range_loss)/112 will be the percentage

##cases with intermediate (from 11 to 50)
intermediate_lost_no_phylo_range_loss = sum(counts_text_range_loss[which(rownames(counts_text_range_loss)=="without_phylo"), which(colnames(counts_text_range_loss) %in% c("'10-20'", "'20-30'", "'30-40'", "'40-50'"))])
#percentage respect to the total number of species
print((intermediate_lost_no_phylo_range_loss*100)/number_species) #112 species is the 100%, intermediate_lost_no_phylo_range_loss will be X; (100*intermediate_lost_no_phylo_range_loss)/112 will be the percentage

##cases with intermediate more high
intermediate_high_lost_no_phylo_range_loss = sum(counts_text_range_loss[which(rownames(counts_text_range_loss)=="without_phylo"), which(colnames(counts_text_range_loss) %in% c("'10-20'", "'20-30'", "'30-40'", "'40-50'", "'50-60'", "'60-70'", "'70-80'", "'80-90'", "'90-100'"))])
#percentage respect to the total number of species
print((intermediate_high_lost_no_phylo_range_loss*100)/number_species) #112 species is the 100%, intermediate_lost_no_phylo_range_loss will be X; (100*intermediate_lost_no_phylo_range_loss)/112 will be the percentage

##cases with intermediate more high
more20_no_phylo_range_loss = sum(counts_text_range_loss[which(rownames(counts_text_range_loss)=="without_phylo"), which(colnames(counts_text_range_loss) %in% c("'20-30'", "'30-40'", "'40-50'", "'50-60'", "'60-70'", "'70-80'", "'80-90'", "'90-100'"))])
#percentage respect to the total number of species
print((more20_no_phylo_range_loss*100)/number_species) #112 species is the 100%, intermediate_lost_no_phylo_range_loss will be X; (100*intermediate_lost_no_phylo_range_loss)/112 will be the percentage

##cases with low (from 0 to 10)
low_lost_no_phylo_range_loss = sum(counts_text_range_loss[which(rownames(counts_text_range_loss)=="without_phylo"), which(colnames(counts_text_range_loss) %in% c("'0-10'"))])
#percentage respect to the total number of species
print((low_lost_no_phylo_range_loss*100)/number_species) #112 species is the 100%, low_lost_no_phylo_range_loss will be X; (100*low_lost_no_phylo_range_loss)/112 will be the percentage

#with phylo
#number of species in that situation
high_lost_phylo_range_loss = sum(counts_text_range_loss[which(rownames(counts_text_range_loss)=="with_phylo"), which(colnames(counts_text_range_loss) %in% c("'80-90'", "'90-100'"))])
#percentage respect to the total number of species
print((high_lost_phylo_range_loss*100)/number_species) #112 species is the 100%, high_lost_phylo_range_loss will be X; (100*high_lost_phylo_range_loss)/112 will be the percentage

##cases with intermediate (from 11 to 50)
intermediate_lost_phylo_range_loss = sum(counts_text_range_loss[which(rownames(counts_text_range_loss)=="with_phylo"), which(colnames(counts_text_range_loss) %in% c("'10-20'", "'20-30'", "'30-40'", "'40-50'"))])
#percentage respect to the total number of species
print((intermediate_lost_phylo_range_loss*100)/number_species) #112 species is the 100%, intermediate_lost_phylo_range_loss will be X; (100*intermediate_lost_phylo_range_loss)/112 will be the percentage

##cases with intermediate more high
intermediate_high_lost_phylo_range_loss = sum(counts_text_range_loss[which(rownames(counts_text_range_loss)=="with_phylo"), which(colnames(counts_text_range_loss) %in% c("'10-20'", "'20-30'", "'30-40'", "'40-50'", "'50-60'", "'60-70'", "'70-80'", "'80-90'", "'90-100'"))])
#percentage respect to the total number of species
print((intermediate_high_lost_phylo_range_loss*100)/number_species) #112 species is the 100%, intermediate_lost_phylo_range_loss will be X; (100*intermediate_lost_phylo_range_loss)/112 will be the percentage

##cases with intermediate more high
more20_phylo_range_loss = sum(counts_text_range_loss[which(rownames(counts_text_range_loss)=="with_phylo"), which(colnames(counts_text_range_loss) %in% c("'20-30'", "'30-40'", "'40-50'", "'50-60'", "'60-70'", "'70-80'", "'80-90'", "'90-100'"))])
#percentage respect to the total number of species
print((more20_phylo_range_loss*100)/number_species) #112 species is the 100%, intermediate_lost_no_phylo_range_loss will be X; (100*intermediate_lost_no_phylo_range_loss)/112 will be the percentage

##cases with low (from 0 to 10)
low_lost_phylo_range_loss = sum(counts_text_range_loss[which(rownames(counts_text_range_loss)=="with_phylo"), which(colnames(counts_text_range_loss) %in% c("'0-10'"))])
#percentage respect to the total number of species
print((low_lost_phylo_range_loss*100)/number_species) #112 species is the 100%, low_lost_phylo_range_loss will be X; (100*low_lost_phylo_range_loss)/112 will be the percentage


##range_change
##no phylo
#number of species with high loss (less -100 more than -80)
high_lost_no_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="without_phylo"), which(colnames(counts_text_range_change) %in% c("'-100--80'"))])
#percentage respect to the total number of species
print((high_lost_no_phylo_range_change*100)/number_species) #112 species is the 100%, high_lost_no_phylo_range_change will be X; (100*high_lost_no_phylo_range_change)/112 will be the percentage

##cases with intermediate (from -79.9 to -20)
intermediate_lost_no_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="without_phylo"), which(colnames(counts_text_range_change) %in% c("'-100--80'", "'-80--60'", "'-60--40'", "'-40--20'"))])
#percentage respect to the total number of species
print((intermediate_lost_no_phylo_range_change*100)/number_species) #112 species is the 100%, intermediate_lost_no_phylo_range_change will be X; (100*intermediate_lost_no_phylo_range_change)/112 will be the percentage

##cases with low (from -20 to 0)
low_lost_no_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="without_phylo"), which(colnames(counts_text_range_change) %in% c("'-20-0'"))])
#percentage respect to the total number of species
print((low_lost_no_phylo_range_change*100)/number_species) #112 species is the 100%, low_lost_no_phylo_range_change will be X; (100*low_lost_no_phylo_range_change)/112 will be the percentage


#number of species with high gain (to 80.1 to 100)
high_gain_no_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="without_phylo"), which(colnames(counts_text_range_change) %in% c("'80-100'"))])
#percentage respect to the total number of species
print((high_gain_no_phylo_range_change*100)/number_species) #112 species is the 100%, high_gain_no_phylo_range_change will be X; (100*high_gain_no_phylo_range_change)/112 will be the percentage

##cases with intermediate (from 20.1 to 80)
intermediate_gain_no_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="without_phylo"), which(colnames(counts_text_range_change) %in% c("'20-40'", "'40-60'", "'60-80'"))])
#percentage respect to the total number of species
print((intermediate_gain_no_phylo_range_change*100)/number_species) #112 species is the 100%, intermediate_gain_no_phylo_range_change will be X; (100*intermediate_gain_no_phylo_range_change)/112 will be the percentage

##cases with low gain (from 0 to 40)
low_gain_no_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="without_phylo"), which(colnames(counts_text_range_change) %in% c("'0-20'", "'20-40'"))])
#percentage respect to the total number of species
print((low_gain_no_phylo_range_change*100)/number_species) #112 species is the 100%, low_gain_no_phylo_range_change will be X; (100*low_gain_no_phylo_range_change)/112 will be the percentage

##cases with some gain
some_gain_no_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="without_phylo"), which(colnames(counts_text_range_change) %in% c("'0-20'", "'20-40'", "'40-60'", "'60-80'", "'80-100'"))])
#percentage respect to the total number of species
print((some_gain_no_phylo_range_change*100)/number_species) #112 species is the 100%, low_gain_no_phylo_range_change will be X; (100*low_gain_no_phylo_range_change)/112 will be the percentage

##with phylo
#number of species with with high loss (less -100 more than -80)
high_lost_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="with_phylo"), which(colnames(counts_text_range_change) %in% c("'-100--80'"))])
#percentage respect to the total number of species
print((high_lost_phylo_range_change*100)/number_species) #112 species is the 100%, high_lost_phylo_range_change will be X; (100*high_lost_phylo_range_change)/112 will be the percentage

##cases with intermediate (from -79.9 to -20)
intermediate_lost_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="with_phylo"), which(colnames(counts_text_range_change) %in% c("'-80--60'", "'-60--40'", "'-40--20'"))])
#percentage respect to the total number of species
print((intermediate_lost_phylo_range_change*100)/number_species) #112 species is the 100%, intermediate_lost_phylo_range_change will be X; (100*intermediate_lost_phylo_range_change)/112 will be the percentage

##cases with low (from -20 to 0)
low_lost_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="with_phylo"), which(colnames(counts_text_range_change) %in% c("'-20-0'"))])
#percentage respect to the total number of species
print((low_lost_phylo_range_change*100)/number_species) #112 species is the 100%, low_lost_phylo_range_change will be X; (100*low_lost_phylo_range_change)/112 will be the percentage


#number of species with high gain (to 80.1 to 100)
high_gain_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="with_phylo"), which(colnames(counts_text_range_change) %in% c("'80-100'"))])
#percentage respect to the total number of species
print((high_gain_phylo_range_change*100)/number_species) #112 species is the 100%, high_gain_phylo_range_change will be X; (100*high_gain_phylo_range_change)/112 will be the percentage

##cases with intermediate gain (from 20.1 to 80)
intermediate_gain_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="with_phylo"), which(colnames(counts_text_range_change) %in% c("'20-40'", "'40-60'", "'60-80'"))])
#percentage respect to the total number of species
print((intermediate_gain_phylo_range_change*100)/number_species) #112 species is the 100%, intermediate_gain_phylo_range_change will be X; (100*intermediate_gain_phylo_range_change)/112 will be the percentage

##cases with low gain (from 0 to 40)
low_gain_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="with_phylo"), which(colnames(counts_text_range_change) %in% c("'0-20'"))])
#percentage respect to the total number of species
print((low_gain_phylo_range_change*100)/number_species) #112 species is the 100%, low_gain_phylo_range_change will be X; (100*low_gain_phylo_range_change)/112 will be the percentage

##cases with some gain
some_gain_phylo_range_change = sum(counts_text_range_change[which(rownames(counts_text_range_change)=="with_phylo"), which(colnames(counts_text_range_change) %in% c("'0-20'", "'20-40'", "'20-40'", "'40-60'", "'60-80'", "'80-100'"))])
#percentage respect to the total number of species
print((some_gain_phylo_range_change*100)/number_species) #112 species is the 100%, low_gain_no_phylo_range_change will be X; (100*low_gain_no_phylo_range_change)/112 will be the percentage

# 5c) Calculate the differences between phylo - no phylo
differ_percent = data.frame(selected_species=NA, range_loss_no_phylo=NA, range_loss_phylo=NA, differ_range_loss=NA, range_change_no_phylo=NA, range_change_phylo=NA, differ_range_change=NA)
for(i in 1:length(epithet_species_list)){

    #selected species
    selected_species = epithet_species_list[i]

    #extract percentage
    range_loss_no_phylo = suitability_changes[which(suitability_changes$species==selected_species),]$range_loss_no_phylo
    range_loss_phylo = suitability_changes[which(suitability_changes$species==selected_species),]$range_loss_phylo   
    range_change_no_phylo = suitability_changes[which(suitability_changes$species==selected_species),]$range_change_no_phylo
    range_change_phylo = suitability_changes[which(suitability_changes$species==selected_species),]$range_change_phylo  

    #calculate absolute difference
    differ_range_loss = abs(range_loss_phylo-range_loss_no_phylo)
    differ_range_change = abs(range_change_phylo-range_change_no_phylo)

    #save it
    differ_percent = rbind.data.frame(
        differ_percent, 
        cbind.data.frame(selected_species, range_loss_no_phylo, range_loss_phylo, differ_range_loss, range_change_no_phylo, range_change_phylo, differ_range_change))
}
differ_percent = differ_percent[which(!apply(is.na(differ_percent), MARGIN=1, FUN=all)),]#remove the row with all NA taking only rows without all NA. This was made applying all to a data frame with true or false for is.na. Rows with all TRUe would be NA rows
nrow(differ_percent)

#calculate median and the interquartile range
median_values = cbind.data.frame("global median", rbind.data.frame(apply(differ_percent[,-1], 2, median)))
first_quartile_values = cbind.data.frame("global first quartile", rbind.data.frame(apply(differ_percent[,-1], 2, quantile, 0.25)))
third_quartile_values = cbind.data.frame("global third quartile", rbind.data.frame(apply(differ_percent[,-1], 2, quantile, 0.75)))
interquartile_range = cbind.data.frame("global interquartile range", third_quartile_values[-1]-first_quartile_values[-1]) #we have to add [-1] to avoid selecting the first column with the name
#check
interquartile_range[,-1] == third_quartile_values[-1] - first_quartile_values[-1]
#add column names
names(median_values) <- c("selected_species", colnames(differ_percent[,-1]))
names(first_quartile_values) <- c("selected_species", colnames(differ_percent[,-1]))
names(third_quartile_values) <- c("selected_species", colnames(differ_percent[,-1]))
names(interquartile_range) <- c("selected_species", colnames(differ_percent[,-1]))
    #It makes sense to use the mean and the standard deviation as measures of center and spread only for distributions that are reasonably symmetric with a central peak. When outliers are present, the mean and standard deviation are not a good choice (unless you want that these outliers influence the summary statistic). An outlier can decrease/increase the mean so that the mean is too low or too high to be representative of whole sample The outlier can also decrease/increase the standard deviation, which gives the impression of a wide variability in the variable This makes sense because the standard deviation measures the average deviation of the data from the mean. So a point that has a large deviation from the mean will increase the average of the deviations. This is the case for several of the variables presented in this table (see annotated plot line). Therefore, we are going to use the median and the interquartile range.
        #par(mfcol=c(3,2)); for(i in 2:ncol(differ_percent)){plot(differ_percent[,i])}
        #Link very useful and simple: https://courses.lumenlearning.com/wmopen-concepts-statistics/chapter/standard-deviation-4-of-4/
    # I guess it does not make sense to use median +- IQR because the distance between the median and the second quartile could be different from the distance to the third quartile. 
        #For example, consider the following sequence: c(1,2,3,4,5,6,7,8,9,9,9,10). The second quartile would be 3.750, the median would be 6.500, and the third quartile would be 9.000. 9.000-6.500=2.5, while 6.500-3.750=2.75. Therefore, the median is closer to the second quartile.

#bind to the data.frame
differ_percent = rbind.data.frame(median_values, first_quartile_values, third_quartile_values, interquartile_range, differ_percent)
str(differ_percent)

#save it
write.table(differ_percent, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/global_figures/final_global_figures/differ_phylo_inside_v3.csv", col.names = TRUE, row.names = FALSE, sep=",")