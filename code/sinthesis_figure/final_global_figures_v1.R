#set working directory
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus")

#required packages
require(raster)
require(RColorBrewer)

#load bio1 for using it as a background
bio1 = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/bio1.asc")#better bio1 than bioclim variables because some water bodies inside continent have climate data.
albicaulis_distribution = raster("/Volumes/GoogleDrive/My Drive/science/phd/dispersal_heterogeneity/datos/DATOS/MAPS/p_albicaulis_01.img") #load ablicaulis buffer to get resolution of distribution maps
bio1 = resample(bio1, albicaulis_distribution, method="bilinear") #reduce resolution
bio1[which(getValues(bio1) >= min(getValues(bio1),na.rm = TRUE))] <- 0 #set all continent areas as 0

#list species
list_species = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/species.txt", sep="\t", header=T)
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

#drop discolor (problem tazonomy, no diferetiaced from cembriodes)
if("discolor" %in% epithet_species_list){
    epithet_species_list = epithet_species_list[-which(epithet_species_list == "discolor")]
}

#check it
!"discolor" %in% epithet_species_list
length(epithet_species_list) == 112

##load distribution rasters
#number of pines per cell
sum_distributions = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/sum_distributions.asc") 

#number of pines per cell including distribution and migration buffer
sum_buffers_distrib = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/sum_buffers_distrib.asc") 

#number of pines per cell including migration but not distribution buffer
sum_buffers = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/sum_buffers.asc")

##################
##### PLOTS ######
##################

#set extent of the raster to be plotted
plot_extent = c(-180,180,-10,90)

###############################
####### sum suitability #######
###############################
#Sum suitability across pines divided by the maximum suitability per site (nº pines*100)

##load sum suitability raster
sum_suitability = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/suit_change_stack_sum_proport_sum_suitability_phylo_cor_no.asc") #sum of suitability per cell for all pine species and divided by the maximum possible suitability according to the number of species (e.g. if 6 species are present in a cell, 6*100 would be the maximum suitability)

##plot
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/figures_main_text/global_map_sum_suit.pdf", width=6, height=6)
par(mfcol=c(2,1), mai=c(0.5,0.5,0.5,1))

##first pannel
#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="Pine richness") #higher values in argument start of gray colors lead to brighter gray

#sum_distributions
plot(crop(sum_distributions, plot_extent), add=TRUE, col= rev(terrain.colors(41))[-1], legend=TRUE, axis.args=list(at=seq(1,19,3)), legend.shrink=0.75) #rev invert vector of colors. "At" argument in axis.args indicate the numbers in the legend

##second pannel
#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main=expression(bold("Predicted habitat suitability for" ~ bolditalic("Pinus")~ "spp.")))

#sum suitability across species
plot(crop(sum_suitability, plot_extent), add=TRUE, col= brewer.pal(7, "RdYlGn"), legend.shrink=0.75) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"
dev.off() 

#other options of colors
#rev(heat.colors(40))
#colorRampPalette(c("chartreuse3", "yellow", "red"))(100)
#rev(brewer.pal(7, "YlOrRd")

###outside suitability inside larger for talks
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/figures_talks/global_map_sum_suit_alone.pdf", width=12, height=6)
#pannel 1
par(mfcol=c(1,1), oma=c(0,0,1,1))
##second pannel
#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="")

#sum suitability across species
plot(crop(sum_suitability, plot_extent), add=TRUE, col= brewer.pal(7, "RdYlGn"), legend.shrink=0.75) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"

#set title
mtext(text=expression(bold("Predicted global suitability inside pine ranges")), side=3, line=1, outer=FALSE, cex=2, font = 2)

dev.off() 

###richness larger for talks
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/figures_talks/global_richness_alone.pdf", width=12, height=6)
#pannel 1
par(mfcol=c(1,1), oma=c(0,0,1,1))
##second pannel
#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="")

#sum_distributions
plot(crop(sum_distributions, plot_extent), add=TRUE, col= rev(terrain.colors(41))[-1], legend=TRUE, axis.args=list(at=seq(1,19,3)), legend.shrink=0.75) #rev invert vector of colors. "At" argument in axis.args indicate the numbers in the legend

#set title
mtext(text=expression(bold("Pine richness")), side=3, line=1, outer=FALSE, cex=2, font = 2)

dev.off() 


##############################################################
####### sum suitability outside + evol. potential ############
##############################################################
#Sum suitability across pines outside distribution divided by the maximum suitability per site (nº pines*100).

####required raster calculations
## sum_suitability_outside
#sum of suitability per cell OUTSIDE range for all pine species and divided by the maximum possible suitability according to the number of species (e.g. if 6 species are present in a cell, 6*100 would be the maximum suitability)
sum_suitability_outside = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/suit_change_stack_sum_proport_sum_suitability_outside_phylo_cor_no.asc")

#Raster to plot: from sum_suitability_outside create a plot with only 1 on those areas with data to plotting
sum_suitability_outside_plot = raster(extent(sum_suitability_outside), resolution = res(sum_suitability_outside))
sum_suitability_outside_plot[]<-0
sum_suitability_outside_plot[which(!is.na(getValues(sum_suitability_outside)))] <- 1

#Raster for transparency: from sum_suitability_outside creates a plot in which NA values are also 0. This is not a problem, because NA and 0 will be transparent in the plot (0 transparency), perfect. 
sum_suitability_outside_plot_zeros = sum_suitability_outside+(100-max(getValues(sum_suitability_outside), na.rm = TRUE)) #sum the difference to 100 form the maximum value to convert max values into 100. This will decrease transparency in the plots, and it is applied in all the map for equal (and the same for evol potential) so there is no bias, only it is an improvement of the visualization
sum_suitability_outside_plot_zeros[which(is.na(getValues(sum_suitability_outside)))] <- 0 #set NA values as zero
sum_suitability_outside_plot_zeros[which(getValues(sum_suitability_outside)==0)] <- 0 #set as zero values that are zero in sum_suitability_outside. Note that this values were summed the maximum value of the raster, so there were artificially increased. 

##phylo_suit_inout
#sum of phylo suitability per cell (inside and outside of distribution) for all pine species and divided by the maximum possible suitability according to the number of species (e.g. if 6 species are present in a cell, 6 would be the maximum suitability; phylo suit is from 0 to 1)
phylo_suit_inout = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/suit_change_stack_sum_proport_sum_suitability_phylo_inout_phylo_cor_yes.asc")

#Raster to plot: from phylo_suit_inout create a plot with only 1 on those areas with data to plotting
phylo_suit_inout_plot = raster(extent(phylo_suit_inout), resolution = res(phylo_suit_inout))
phylo_suit_inout_plot[]<-0
phylo_suit_inout_plot[which(!is.na(getValues(phylo_suit_inout)))] <- 1

#Raster for transparency: from phylo_suit_inout creates a plot in which NA values are also 0. This is not a problem, because NA and 0 will be transparent in the plot (0 transparency), perfect. 
phylo_suit_inout_plot_zeros = phylo_suit_inout+(100-max(getValues(phylo_suit_inout), na.rm = TRUE)) #we sum to get 100 in areas with the maximum value, then we set as 0 those areas that in the begining were zero (last line)
phylo_suit_inout_plot_zeros[which(is.na(getValues(phylo_suit_inout)))] <- 0 #set NA values as zero
phylo_suit_inout_plot_zeros[which(getValues(phylo_suit_inout)==0)] <- 0 #set as zero values that are zero in phylo_suit. Note that this values were summed the maximum value of the raster, so there were artificially increased.

##combination of suitability future outside and phylo-suit in-out
#combine phylo_suit_inout and sum_suitability_outside
sum_suitability_outside_to_combine = sum_suitability_outside

#set as 0 those areas that are NA in sum_suitability_outside but 0 in phylo_suit_inout. This avoid problems in the sum (NA propagation)
sum_suitability_outside_to_combine[which(is.na(getValues(sum_suitability_outside)) & !is.na(getValues(phylo_suit_inout)))] <- 0 

#sum both rasters
final_uncertainty_raster = sum_suitability_outside_to_combine+phylo_suit_inout

#the raster with more than 100 are set as 100. This is similar to set as suitable those cell suitables accoridng to the SDM and in addition to the phylo range. 
final_uncertainty_raster[which(getValues(final_uncertainty_raster) > 100)] <- 100 #avoid in that way vlaues higher than 100 (1 for alpha argument in plot.raster) 

#from final_uncertainty_raster create a plot with only 1 on those areas with data to plotting
final_uncertainty_raster_plot = raster(extent(final_uncertainty_raster), resolution = res(final_uncertainty_raster))
final_uncertainty_raster_plot[]<-0
final_uncertainty_raster_plot[which(!is.na(getValues(final_uncertainty_raster)))] <- 1

#from final_uncertainty_raster creates a plot in which NA values are also 0. This is not a problem, because NA and 0 will be transparent in the plot (0 transparency), perfect. 
final_uncertainty_raster_plot_zeros = final_uncertainty_raster+(100-max(getValues(final_uncertainty_raster), na.rm = TRUE)) #we sum to get 100 in areas with the maximum value, then we set as 0 those areas that in the begining were zero (last line)
final_uncertainty_raster_plot_zeros[which(is.na(getValues(final_uncertainty_raster)))] <- 0 #set NA values as zero
final_uncertainty_raster_plot_zeros[which(getValues(final_uncertainty_raster)==0)] <- 0 #set as zero values that are zero in phylo_suit. Note that this values were summed the maximum value of the raster, so there were artificially increased.

###outside suitability and evol potential in the same raster 
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/figures_main_text/global_map_migration.pdf", width=12, height=6)
#pannel 1
par(mfcol=c(1,1), oma=c(0,0,1,1))

#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="")

#plot sum_suitability
plot(crop(sum_suitability, plot_extent), add=TRUE, col= (brewer.pal(7, "RdYlGn")), legend.shrink=0.75, legend=TRUE, alpha=1) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"

#plot final_uncertainty_raster (sum of suitaiblity outside and phylo-suit in-out)
plot(crop(final_uncertainty_raster_plot, plot_extent), add=TRUE, col= "#08fc51", alpha=crop(final_uncertainty_raster_plot_zeros/100, plot_extent), legend=FALSE)

#set title
mtext(text=expression(bold("Predicted global suitability for" ~ bolditalic("Pinus") ~ "spp. in 2070")), side=3, line=1, outer=FALSE, cex=2.5, font = 2)
dev.off()

#the same plto with migration and phylo for talks
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/figures_talks/global_map_migration_talk.pdf", width=12, height=6)
#pannel 1
par(mfcol=c(1,1), oma=c(0,0,1,1))

#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="")

#plot sum_suitability
plot(crop(sum_suitability, plot_extent), add=TRUE, col= (brewer.pal(7, "RdYlGn")), legend.shrink=0.75, legend=TRUE, alpha=1) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"

#plot final_uncertainty_raster (sum of suitaiblity outside and phylo-suit in-out)
plot(crop(final_uncertainty_raster_plot, plot_extent), add=TRUE, col= "#08fc51", alpha=crop(final_uncertainty_raster_plot_zeros/100, plot_extent), legend=FALSE)

#set title
mtext(text=expression(bold("Predicted global suitability outside pine ranges + Evol. potential")), side=3, line=1, outer=FALSE, cex=2, font = 2)
dev.off()

###outside suitability and evol potential in the same raster compare to other map with only outside suitability
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/figures_supple/global_map_migration_with_without_phylo.pdf", width=6, height=6)
par(mfcol=c(2,1), mai=c(0.5,0.5,0.5,1))

#panel 1
#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="")

#plot sum_suitability
plot(crop(sum_suitability, plot_extent), add=TRUE, col= (brewer.pal(7, "RdYlGn")), legend.shrink=0.75, legend=TRUE, alpha=1) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"

#sum_suitability_outside_plot
plot(crop(sum_suitability_outside_plot, plot_extent), add=TRUE, col= "#08fc51", alpha=crop(sum_suitability_outside_plot_zeros/100, plot_extent), legend=FALSE)

#pannel 2
#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="")

#plot sum_suitability
plot(crop(sum_suitability, plot_extent), add=TRUE, col= (brewer.pal(7, "RdYlGn")), legend.shrink=0.75, legend=TRUE, alpha=1) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"

#plot final_uncertainty_raster (sum of suitaiblity outside and phylo-suit in-out)
plot(crop(final_uncertainty_raster_plot, plot_extent), add=TRUE, col= "#08fc51", alpha=crop(final_uncertainty_raster_plot_zeros/100, plot_extent), legend=FALSE)
dev.off()


#####both outside suitability and evol separated in different rasters
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/draft_figures/global_map_migration_two_rasters_version.pdf", width=12, height=6)

##pannel 1
par(mfcol=c(1,1), oma=c(0,0,1,1))

#gray background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="")

#sum suitability
plot(crop(sum_suitability, plot_extent), add=TRUE, col= (brewer.pal(7, "RdYlGn")), legend.shrink=0.75, legend=TRUE, alpha=1) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"

#set again par with new=TRUE to avoid the misalignment of overlapped rasters
par(mfcol=c(1,1), oma=c(0,0,1,1), new=TRUE) #overlay raster: "https://stackoverflow.com/questions/24213453/overlay-raster-plot-using-plot-add-t-leads-to-arbitrary-misalignment-of-fin"

#plot sum_suitability_outside
plot(crop(sum_suitability_outside_plot, plot_extent), add=FALSE, col= "#08fc51", alpha=crop(sum_suitability_outside_plot_zeros/100, plot_extent), legend=FALSE)#alpha is the transparency parameter, form 0 to 1, because of this we divided the raster for transparency (0-100) by 100. 

#set again par with new=TRUE to avoid the misalignment of overlapped rasters
par(mfcol=c(1,1), oma=c(0,0,1,1), new=TRUE) #overlay raster: "https://stackoverflow.com/questions/24213453/overlay-raster-plot-using-plot-add-t-leads-to-arbitrary-misalignment-of-fin"

#plot phylo_suit_inout
plot(crop(phylo_suit_inout_plot, plot_extent), add=FALSE, col= "#1609fc", alpha=crop(phylo_suit_inout_plot_zeros/100, plot_extent), legend=FALSE)#alpha is the transparency parameter, form 0 to 1, because of this we divided the raster for transparency (0-100) by 100.

#set title
mtext(text=expression(bold("Predicted global suitability for" ~ bolditalic("Pinus") ~ "spp. in 2070")), side=3, line=1, outer=FALSE, cex=2.5, font = 2)

#set legend
legend("topleft", legend="Future suit outside distribution", fill="#08fc51", cex=0.8)
legend("topright", legend="Evolutionary potential", fill="#1609fc", cex=0.8)
dev.off() #La nueva figura me parece una buena aproximación, pero eso es solo una opinión personal, así que si vosotros preferís un solo raster, lo dejamos así. La leyenda para los colores debería ser algo así como "Lime green and blue denote areas outside extant distribution that can be suitable in the future. Lime green: areas that can be reiched by dispersal; Blue: areas that can be colonized because of the evolutionary potential within pine species" Luego en el mapa solo pondría "Dispersal" y "Adaptation".


#plot with suit otuisde and whiput phylo for talks
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/figures_talks/global_map_migration_no_phylo.pdf", width=12, height=6)

#pannel 1
par(mfcol=c(1,1), oma=c(0,0,1,1))

#gray background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="")

#sum suitability
plot(crop(sum_suitability, plot_extent), add=TRUE, col= (brewer.pal(7, "RdYlGn")), legend.shrink=0.75, legend=TRUE, alpha=1) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"

#plot sum_suitability_outside
plot(crop(sum_suitability_outside_plot, plot_extent), add=TRUE, col= "#08fc51", alpha=crop(sum_suitability_outside_plot_zeros/100, plot_extent), legend=FALSE)#alpha is the transparency parameter, form 0 to 1, because of this we divided the raster for transparency (0-100) by 100. 

#set title
mtext(text=expression(bold("Predicted global suitability outside pine ranges")), side=3, line=1, outer=FALSE, cex=2, font = 2)
dev.off()

###############################
######### Histograms ##########
###############################

#load data for histogram of decreases of suitability
percent_no_phylo = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/percentage_data_phylo_cor_no.csv", sep=",", header=TRUE)
percent_phylo = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/percentage_data_phylo_cor_yes.csv", sep=",", header=TRUE)
percent_no_phylo$species == percent_phylo$species#species order match

###stack histogram 
# 1) Define the breaks to use on your Histogram
xrange = seq(-10,100,10)

# 2) Have your vectors ready
v1 = percent_no_phylo$percent_loss #no_filo
v2 = percent_phylo$percent_loss #filo

# 3) subset your vectors to be inside xrange
v1 = subset(v1,v1<=max(xrange) & v1>=min(xrange)) #subsetea el vector con los valores de suitability loss con y sin corrección filo, cogiendo todos los valores menores iguales y mayores o iguales que el rango que elegido previamente.
v2 = subset(v2,v2<=max(xrange) & v2>=min(xrange)) #the xrange covers all data included in v1 and v2 (this is the idea)

# 4) Now, use hist to compute the counts per interval
hv1 = hist(v1,breaks=xrange,plot=F)$counts #hacemos un histograma con los % de perdida de suit sin filo correction, pero no lo ploteamos. De ahí sacamos los counts, o sea, el número de casos de cada % de perdida, es decir, la altura de las barras: Cuantas especies pierden suitability desde -10 (sin incluirlo) hasta 0 (incluyendolo, right close, left open histogram)
length(which(v1>-10 & v1<=0)) == head(hv1, 1) #It gives 8, like the first count, which goes from -10 to 0.
length(which(v1>90 & v1<=100)) == tail(hv1, 1) #It gives 46, like the first count, which goes from 90 (not including it) to 100 (including it). 
hv2 = hist(v2,breaks=xrange,plot=F)$counts
length(which(v2>-10 & v2<=0)) == head(hv2, 1) #It gives 8, like the first count, which goes from -10 to 0.
length(which(v2>90 & v2<=100)) == tail(hv2, 1) #It gives 44, like the first count, which goes from 90 (not including it) to 100 (including it). 

#5a) Generate a Frequency BarPlot with bars parallels
pdf(file="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/figures_main_text/percent_histogram_final.pdf",width=12,height=10)
par(mar=c(5.2,5.2,5.2,5.2), mai=c(2,1.5,0.8,1), mgp=c(4.5,1.7,0), cex.main=1.7, cex.axis=1.8)

#plot bars
xs = barplot(c(rbind(hv1, hv2)), #counts with and without phylo combinated, the first count wihtout correction, then the first count with, the second count without, second count with... It is to say, both vectors are merged alternating indexes
    col=c("gray39", "gray87"), #two different grays
    #names.arg=xrange[-1], #nombres del eje X, son los breaks
    space=c(0, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0), #espacio cero entre las dos barras de una misma categoria, y luego 0.2 entre cateogiras de perdida de idoneidad. 
    las=0, #axis labels paralelas al eje
    xlab="Suitability loss (%)", ylab="Number of species", cex.lab=3, cex.axis=2.3, main=NULL, lwd=1, font.lab=2, ylim=c(0,50), yaxt='n')

#add label breaks of x axis
mtext(seq(-10,100,10), side=1, at=c(xs[seq(1,22,2)]-0.5,xs[length(xs)]+0.5), cex=2, line=1.5) #at the beginning of each bar (a little bit move to left), and another label at the end, because a seq from 1 to 22 each 2 skip 22, thus we added that 22 with length(xs)

#add label breaks of y axis
mtext(seq(0,50,10), side=2, at=seq(0,50,10), cex=2, line=1.5)

#add legend
legend(x=-0.1, y=48, legend=c("Without phylogenetic correction", "With phylogenetic correction"), fill=c("gray39", "gray87"), cex=2)

#add box
box(lwd=1)
dev.off()

# 5b) Generate a data.frame with counts and labels

#for each break (xrange) except the last one 
labels_counts = NULL
for(i in 1:(length(xrange)-1)){

    #create range label
    range_label = paste("'", xrange[i], "-", xrange[i+1], "'", sep="")#we add '' to avoid that excel considers the column names as number and does estrange things

    #save it
    labels_counts = append(labels_counts, range_label)
}
#check that all except the last one are included
length(labels_counts) == length(xrange)-1

#bind counts without and with phylogenetic correction
counts_text = rbind.data.frame(hv1, hv2)

#set colnames as labels_counts
colnames(counts_text) <- labels_counts

#set row names as no-phylo and phylo
rownames(counts_text) <- c("without_phylo", "with_phylo")
str(counts_text)
counts_text

#save it
write.table(counts_text, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/figures_main_text/hist_counts.csv", sep=",", row.names = TRUE, col.names = NA) #we set row.names = TRUE and col.names = NA to get as first column of the csv an empty space and then the two row names. The rest of column have as first cell the range label. If we set col.names=TRUE, then the first column with row names has a range label, and the following columns are wrong set. See "https://stackoverflow.com/questions/2478352/write-table-writes-unwanted-leading-empty-column-to-header-when-has-rownames".

### extract percentage for the paper results
#extract the number of species
species = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/species.txt", sep="\t", header=TRUE)
nrow(species) == 112
number_species = nrow(species)

##percentage of species that are predicted to lose more than 80% suitable areas inside their range without the phylogenetic correction. We can include 80-90 because this is a right-closed and left open histogram. Thus in 80-90 are not included species that lost 80% of their suitable are, but species that lost exactly  the 90% are included
#number of species in that situation
high_lost_no_phylo = sum(counts_text[which(rownames(counts_text)=="without_phylo"), which(colnames(counts_text) %in% c("'80-90'", "'90-100'"))])
#percentage respect to the total number of species
print((high_lost_no_phylo*100)/number_species) #112 species is the 100%, high_lost_no_phylo will be X; (100*high_lost_no_phylo)/112 will be the percentage


##percentage of species that are predicted to increase suitable areas inside their range without the phylogenetic correction.
#number of species in that situation
increase_suit_no_phylo = sum(counts_text[which(rownames(counts_text)=="without_phylo"), which(colnames(counts_text) %in% c("'-10-0'"))])
#percentage respect to the total number of species
print((increase_suit_no_phylo/ number_species)*100)

##percentage of species that are predicted to lose more than 0% of suitable area until 80% inside their range without the phylogenetic correction. 
rest_cases_no_phylo = sum(counts_text[which(rownames(counts_text)=="without_phylo"), which(!colnames(counts_text) %in% c("'-10-0'", "'80-90'", "'90-100'"))])
#percentage respect to the total number of species
print((rest_cases_no_phylo/ number_species)*100)

##percentage of species that are predicted to lose more than 80% suitable areas inside their range with the phylogenetic correction.
#number of species in that situation
high_lost_yes_phylo = sum(counts_text[which(rownames(counts_text)=="with_phylo"), which(colnames(counts_text) %in% c("'80-90'", "'90-100'"))])
#percentage respect to the total number of species
print((high_lost_yes_phylo/ number_species)*100)

##percentage of species that are predicted to lose more than 0% of suitable area until 80% inside their range with the phylogenetic correction. 
rest_cases_yes_phylo = sum(counts_text[which(rownames(counts_text)=="with_phylo"), which(!colnames(counts_text) %in% c("'-10-0'", "'80-90'", "'90-100'"))])
#percentage respect to the total number of species
print((rest_cases_yes_phylo/ number_species)*100)

# 5c) Calculate the differences between phylo - no phylo
differ_percent = data.frame(selected_species=NA, per_no_phylo=NA, per_phylo=NA, differ=NA)
for(i in 1:length(epithet_species_list)){

    #selected species
    selected_species = epithet_species_list[i]

    #extract percentage
    per_no_phylo = percent_no_phylo[which(percent_no_phylo$species==selected_species),]$percent_loss
    per_phylo = percent_phylo[which(percent_phylo$species==selected_species),]$percent_loss   

    #calculate absolute difference
    differ = abs(per_phylo-per_no_phylo)

    #save it
    differ_percent = rbind.data.frame(
        differ_percent, 
        cbind.data.frame(selected_species, per_no_phylo, per_phylo, differ))

}
differ_percent = differ_percent[which(!apply(is.na(differ_percent), MARGIN=1, FUN=all)),]#remove the row with all NA taking only rows without all NA. This was made applying all to a data frame with true or false for is.na. Rows with all TRUe would be NA rows
nrow(differ_percent)

#reorder differ_percent in basis on the differences
differ_percent = differ_percent[order(differ_percent$differ, decreasing=TRUE),]

#save it
write.table(differ_percent, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/figures_main_text/differ_phylo_inside.csv", col.names = TRUE, row.names = FALSE, sep=",")

# 5d) Generate a Frequency BarPlot that is equivalent to a Stacked histogram
pdf(file="/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/draft_figures/percent_histogram.pdf",width=12,height=12)
par(mar=c(5.2,5.2,5.2,5.2), mai=c(2,1.5,0.8,1), mgp=c(4.5,1.7,0), cex.main=1.7, cex.axis=1.8)

#plot bars
xs = barplot(rbind(hv1,hv2), #juntamos todos los counts como filas
    col=c("gray39", "gray87"), #color de siempre
    #names.arg=xrange[-1], #nombres del eje X, son los breaks
    space=0, #espacio cero entre barras
    las=0, #axis labels paralelas al eje
    xlab="Suitability loss (%)", ylab="Number of species", cex.lab=3, cex.axis=2.3, main=NULL, lwd=1, font.lab=2, ylim=c(0,100), yaxt='n')

#add break labels of x and y axis
mtext(seq(-10,100,10), side=1, at=0:11, cex=2, line=1.5)
mtext(seq(0,100,10), side=2, at=seq(0,100,10), cex=2, line=1.5)

#legend
legend(x=-0.1, y=97, legend=c("Without phylogenetic correction", "With phylogenetic correction"), fill=c("gray39", "gray87"), cex=2)

#add box
box(lwd=1)
dev.off() #more info in http://stackoverflow.com/questions/26612805/r-histogram-with-multiple-populations

#6) Other posibility is a density plot of both data againted.
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/draft_figures/densityplot_histo.pdf",width=8,height=8)

#sum counts without and with phylogenetic correction
Total = hv1 + hv2

#plot
barplot(rbind(hv1/Total,hv2/Total),col=c("gray39", "gray87"),names.arg=xrange[-1],space=0,las=1)

#legend
legend(x=-0.1, y=97, legend=c("Without phylogenetic correction", "With phylogenetic correction"), fill=c("gray39", "gray87"), cex=2)
dev.off()

#########################################
######### From 75 to 25 INside ##########
#########################################

###NO phylo 
#from more than 75 to less than 25% of certainty about suitability without correction
more_75_less_25_no_phylo = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/suit_change_stack_sum_proport_more_75_less_25_phylo_cor_no.asc") 

#from more than 75 to less than 25% of certainty about suitability
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/draft_figures/global_map_loss_suit.pdf", width=6, height=6)
par(mfcol=c(2,1), mai=c(0.5,0.5,0.5,1))

##panel 1 
#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="Pine richness")

#sum distributions
plot(crop(sum_distributions, plot_extent), add=TRUE, col= rev(terrain.colors(41))[-1], legend=TRUE, axis.args=list(at=seq(1,19,3)), legend.shrink=0.75) #rev invert a vector

#panel 2
#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="Percentage of pines with suitability losses")

#more_75_less_25_no_phylo
plot(crop(more_75_less_25_no_phylo, plot_extent), add=TRUE, col= rev(brewer.pal(7, "RdYlGn")), legend.shrink=0.75) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"
dev.off() 


#### YES/NO phylo
#from more than 75 to less than 25% of certainty about suitability 
more_75_less_25_yes_phylo = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/suit_change_stack_sum_proport_more_75_less_25_phylo_cor_yes.asc") 

#from more than 75 to less than 25% of certainty about suitability YES/NO phylo
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/draft_figures/global_map_loss_suit_both_phylo.pdf", width=6, height=6)
par(mfcol=c(2,1), mai=c(0.5,0.5,0.5,1))

##panel 1
#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="Suitability loss")

#more_75_less_25_no_phylo
plot(crop(more_75_less_25_no_phylo, plot_extent), add=TRUE, col= rev(brewer.pal(7, "RdYlGn")), legend.shrink=0.75) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"

##panel 2
#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="")

#more_75_less_25_yes_phylo
plot(crop(more_75_less_25_yes_phylo, plot_extent), add=TRUE, col= rev(brewer.pal(7, "RdYlGn")), legend.shrink=0.75) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"

#draw interest circles
require(plotrix)
draw.circle(-65,50,15,nv=100,border=NULL,col=NA,lty=1,density=NULL,angle=45,lwd=1.5)
draw.circle(-85,35,11,nv=100,border=NULL,col=NA,lty=1,density=NULL,angle=45,lwd=1.5)
draw.circle(12,45,13,nv=100,border=NULL,col=NA,lty=1,density=NULL,angle=45,lwd=1.5)
draw.circle(140,40,11,nv=100,border=NULL,col=NA,lty=1,density=NULL,angle=45,lwd=1.5)
dev.off() 

##########################################################
######### more_75_less_25 + categoric migration ##########
##########################################################

#areas suitables after CC without phylogenetic correction and not including distribution
total_75_migra_no_distrib_phylo = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/suit_change_stack_sum_proport_migration_phylo_cor_no.asc") 
#from that plot create a plot with only 1 on those areas with data to plotting
total_75_migra_no_distrib_phylo_plot = raster(extent(total_75_migra_no_distrib_phylo), resolution = res(total_75_migra_no_distrib_phylo))
total_75_migra_no_distrib_phylo_plot[]<-0
total_75_migra_no_distrib_phylo_plot[which(!is.na(getValues(total_75_migra_no_distrib_phylo)))] <- 1
#from that plot creates a plot in which NA values are also 0. This is not a problem, because NA and 0 will be transparent in the plot, perfect. 
total_75_migra_no_distrib_phylo_zeros = total_75_migra_no_distrib_phylo
total_75_migra_no_distrib_phylo_zeros[which(is.na(getValues(total_75_migra_no_distrib_phylo)))] <- 0

#areas suitables after CC with phylogenetic correction and not including distribution
total_75_migra_yes_distrib_phylo = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/suit_change_stack_sum_proport_migration_phylo_cor_yes.asc") 
#from that plot create a plot with only 1 on those areas with data to plotting
total_75_migra_yes_distrib_phylo_plot = raster(extent(total_75_migra_yes_distrib_phylo), resolution = res(total_75_migra_yes_distrib_phylo))
total_75_migra_yes_distrib_phylo_plot[]<-0
total_75_migra_yes_distrib_phylo_plot[which(!is.na(getValues(total_75_migra_yes_distrib_phylo)))] <- 1
#from that plot creates a plot in which NA values are also 0. This is not a problem, because NA and 0 will be transparent in the plot, perfect. 
total_75_migra_yes_distrib_phylo_zeros = total_75_migra_yes_distrib_phylo
total_75_migra_yes_distrib_phylo_zeros[which(is.na(getValues(total_75_migra_yes_distrib_phylo)))] <- 0

#from more than 75 to less than 25% of certainty about suitability more migration possibility but NO phylo
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/draft_figures/global_map_loss_suit_migration_no_phylo.pdf", width=6, height=6)
par(mfcol=c(2,1), mai=c(0.5,0.5,0.5,1))
##panel 1
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="Suitability losses")
plot(crop(more_75_less_25_no_phylo, plot_extent), add=TRUE, col= rev(brewer.pal(7, "RdYlGn")), legend.shrink=0.75)

#pannel 2
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="Suitability losses + Migration opportunity")
plot(crop(more_75_less_25_no_phylo, plot_extent), add=TRUE, col= rev(brewer.pal(7, "RdYlGn")), legend.shrink=0.75, legend=TRUE) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"
plot(crop(total_75_migra_no_distrib_phylo_plot, plot_extent), add=TRUE, col= "#08fc51", alpha=crop(total_75_migra_no_distrib_phylo_zeros/100, plot_extent), legend=FALSE)
legend("topleft", legend="Migration window", fill="#08fc51", cex=0.8)
dev.off() 

#from more than 75 to less than 25% of certainty about suitability more migration possibility and YES/NO phylo
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/draft_figures/global_map_loss_suit_migration_both_phylo.pdf", width=6, height=6)
par(mfcol=c(2,1), mai=c(0.5,0.6,0.5,1))

##panel 1
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.3), legend=FALSE, main="Suitability loss + Migration window")
#more_75_less_25_no_phylo
plot(crop(more_75_less_25_no_phylo, plot_extent), add=TRUE, col= rev(brewer.pal(7, "RdYlGn")), legend.shrink=0.75, legend=TRUE, alpha=1) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"
#migration categoric without phylo
plot(crop(total_75_migra_no_distrib_phylo_plot, plot_extent), add=TRUE, col= "#08fc51", alpha=crop(total_75_migra_no_distrib_phylo_zeros/100, plot_extent), legend=FALSE)
legend("topleft", legend="Migration window", fill="#08fc51", cex=0.8)

##panel 2
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.3), legend=FALSE, main="")
#more_75_less_25_no_phylo

plot(crop(more_75_less_25_no_phylo, plot_extent), add=TRUE, col= rev(brewer.pal(7, "RdYlGn")), legend.shrink=0.75, legend=TRUE, alpha=1) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"
#plot circles around interesting areas
require(plotrix)
draw.circle(-65,52.5,15,nv=100,border=NULL,col=NA,lty=1,density=NULL,angle=45,lwd=1.5)
draw.circle(-120,50,15,nv=100,border=NULL,col=NA,lty=1,density=NULL,angle=45,lwd=1.5)
draw.circle(10,45,12,nv=100,border=NULL,col=NA,lty=1,density=NULL,angle=45,lwd=1.5)
draw.circle(140,40,11,nv=100,border=NULL,col=NA,lty=1,density=NULL,angle=45,lwd=1.5)
#migration categoric with phylo
plot(crop(total_75_migra_yes_distrib_phylo_plot, plot_extent), add=TRUE, col= "#08fc51", alpha=crop(total_75_migra_yes_distrib_phylo_zeros/100, plot_extent), legend=FALSE)
dev.off() 


#################################################
######### From 75 to 25 INside+outside ##########
#################################################

#areas suitables after CC including distribution areas and without phylogenetic correction 
total_75_migra_no_phylo = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/suit_change_stack_sum_proport_total_more_75_distri_migra_phylo_cor_no.asc") 

#areas suitables after CC including distribution areas and with phylogenetic correction 
total_75_migra_phylo = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/initial_global_figures/suit_change_stack_sum_proport_total_more_75_distri_migra_phylo_cor_yes.asc")

#areas suitables after CC without phylogenetic correction
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/draft_figures/global_map_total_suit_no_phylo.pdf", width=6, height=6)
par(mfcol=c(2,1), mai=c(0.5,0.5,0.5,1))
##panel 1
#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="Pine richness")
#sum migration buffers including distribution
plot(crop(sum_buffers_distrib, plot_extent), add=TRUE, col= rev(terrain.colors(41))[-1], legend=TRUE, axis.args=list(at=seq(1,49,7)), legend.shrink=0.75) #rev invert a vector

##panel 2
#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="Percentage of pines with suitability")
#plot cases of more than 75% suit in the future
plot(crop(total_75_migra_no_phylo, plot_extent), add=TRUE, col= brewer.pal(7, "RdYlGn"), legend.shrink=0.75) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"
dev.off()

#areas suitable after CC with phylogenetic correction
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/global_figures/final_global_figures/draft_figures/global_map_total_suit_yes_phylo.pdf", width=6, height=6)
par(mfcol=c(2,1), mai=c(0.5,0.5,0.5,1))

##panel 1
#background
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="Pine richness")
#sum migration buffers including distribution
plot(crop(sum_buffers_distrib, plot_extent), add=TRUE, col= rev(terrain.colors(41))[-1], legend=TRUE, axis.args=list(at=seq(1,49,7)), legend.shrink=0.75) #rev invert a vector

##panel 2
plot(crop(bio1, plot_extent), col=gray.colors(1, start=0.2), legend=FALSE, main="Percentage of pines with suitability")
#plot cases of more than 75% suit in the future
plot(crop(total_75_migra_phylo, plot_extent), add=TRUE, col= brewer.pal(7, "RdYlGn"), legend.shrink=0.75) #colors of colorbrewer2 can be seen at "http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=3"
dev.off()