#####compare rasters of wordclim (used in niche paper) and chelsa (used in seed mass) to see if the units are similar

####define work directory
setwd("/Volumes/GoogleDrive/My\ Drive/science/phd/nicho_pinus")

#Libraries
require(raster)

#extract raster of bioclimatic variables used under current condicions in pines niche. They are stored in the hard drive "my Book"
stack_current_bioclim = stack(list.files("/Volumes/My Book/niche_evolution/climate/current", pattern="bio", full.names=TRUE))
nlayers(stack_current_bioclim) == 19

#for each raster plot it
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/comparison_units_niche_seed_mass_papers/current_wc_nicho.pdf")
for(i in 1:nlayers(stack_current_bioclim)){

    #select the [i] layer
    selected_raster = stack_current_bioclim[[i]]

    #plot it
    plot(selected_raster, main=names(stack_current_bioclim)[i])
}
dev.off()

#extract raster of bioclimatic variables used under current condicions in pines niche. They are stored in the hard drive "my Book"
stack_future_bioclim_bc26 = stack(list.files("/Volumes/My Book/niche_evolution/climate/future", pattern="bc26", full.names=TRUE))
nlayers(stack_future_bioclim_bc26) == 19

#for each raster plot it
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/comparison_units_niche_seed_mass_papers/future_bc26_wc_nicho.pdf")
for(i in 1:nlayers(stack_future_bioclim_bc26)){

    #select the [i] layer
    selected_raster = stack_future_bioclim_bc26[[i]]

    #plot it
    plot(selected_raster, main=names(stack_future_bioclim_bc26)[i])
}
dev.off()

#compare this rasters with "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/comparison_units_niche_seed_mass_papers/bioclim_chelsa.pdf"

#CONCLUSION: THE SAME UNITS IN ALL CASES.