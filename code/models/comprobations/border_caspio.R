#Code para ver como el problema del borde del caspio no es tal. Comparando el mismo modelo y la misma partición con un raster y otro vemos que la idoneidad es exactmaente igual. Lo único es que hay idoneidad hasta el borde pues en un caso sigue el borde suave y en el tro lo sigue de forma gruesa, pero se refiere a la misma zona y le da la misma idoneidad. 

#require package
require(raster)

#load raster with two types of borders
brutia_no_sierra = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/borde_caspio_comprobations/rasters/rf_brutia_caspio.asc")
brutia_sierra = raster("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/borde_caspio_comprobations/rasters/rf_brutia_caspio_sierra.asc")

#plot
pdf("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/borde_caspio_comprobations/plot/sierra_no_sierra_caspio_brutia.pdf", width=12, height=6)
par(mfcol=c(1,2))
plot(brutia_no_sierra, main="Fine border Caspio")
plot(brutia_sierra, main="Thick border Caspio")
dev.off()
