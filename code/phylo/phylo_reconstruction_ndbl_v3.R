#######code for reconstructing climate niche with the ndbl phylogeny

### differences respect to previous versions

# Respect to version 2
    #I have changed the path and modified the final plot to include the test statistic of the spearman correlation test. I have not checked the script, BUT the result of the plot is the same than the original plot.
    #NOT RUN THE ENTIRE SCRIPT, AVOID PARTS WHEN THERE IS WRITTING OF FILES. These files has to be modified by hand, so we are using those created with the previous version of the script (version 2).
    #I installed phangor (used for other packages) by source (version 2.7.0), becuase the current version (2.8.0) does not accept R.3.6.4, my current version.


#######################################################
########## RECONSTRUCTION OF ANCESTRAL STATE ##########
#######################################################

## load required packages
require(ape)
require(phytools)
require(geiger)
require(diversitree)

#load climate data obtained when the first ancestral reconstruction was performed
climate_medians = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_medians/climate_medians.csv", header=TRUE, sep=",") 
str(climate_medians)

#load sd climate data
climate_sd = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/climate_sd.csv", header=TRUE, sep=",") 
str(climate_sd)
climate_se = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction/climatic_var_ranges/climate_se.csv", header=TRUE, sep=",")
str(climate_se)

####load a prepare tree_fbdl
tree_fbdl<-read.nexus("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/datos/phlyo/phylogeny/FBDl_MCC_commAnc.tre") 
## new species icnldued by bianca in tree_fbdl that we have to drop, and also discolor
species_to_drop = tree_fbdl$tip.label[which(!tree_fbdl$tip.label %in% paste("Pinus_", climate_medians$species, sep=""))]

## prune the tree_fbdl of speceis without seed mass data
tree_fbdl_prunned = drop.tip(tree_fbdl, species_to_drop)

## check
c("Pinus_discolor", "Pinus_jaliscana", "Pinus_tecunumanii") %in% tree_fbdl_prunned$tip.label

####load a prepare tree_ndbl
tree_ndbl<-read.nexus("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/datos/phlyo/phylogeny/NDbl_MCC_commAnc.tre") 
## new species icnldued by bianca in tree_ndbl that we have to drop, and also discolor
species_to_drop = tree_ndbl$tip.label[which(!tree_ndbl$tip.label %in% paste("Pinus_", climate_medians$species, sep=""))]

## prune the tree_ndbl of speceis without seed mass data
tree_ndbl_prunned = drop.tip(tree_ndbl, species_to_drop)

## check
c("Pinus_discolor", "Pinus_jaliscana", "Pinus_tecunumanii") %in% tree_ndbl_prunned$tip.label

####check both tree have the same tips
summary(tree_ndbl_prunned$tip.label == tree_fbdl_prunned$tip.label)

##reorder rows of climate data in basis on tip labels 
climate_medians = climate_medians[match(tree_ndbl_prunned$tip.label, paste("Pinus_", climate_medians$species,sep="")),]
climate_sd = climate_sd[match(tree_ndbl_prunned$tip.label, paste("Pinus_", climate_sd$species,sep="")),]
climate_se = climate_se[match(tree_ndbl_prunned$tip.label, paste("Pinus_", climate_se$species,sep="")),]

## save climatic variables in a vector with species names as names
bio4_vector = climate_medians$median_bio4
bio17_vector = climate_medians$median_bio17

## set names of these variables as species names
names(bio4_vector) <- paste("Pinus_", climate_medians$species, sep="")
names(bio17_vector) <- paste("Pinus_", climate_medians$species, sep="")

##intra variability as SE of all data across distribution
intra_var_bio4 = climate_se$se_bio4
names(intra_var_bio4) <- paste("Pinus_", climate_se$species, sep="")
intra_var_bio17 = climate_se$se_bio17
names(intra_var_bio17) <- paste("Pinus_", climate_se$species, sep="")

##check order
names(bio4_vector) == tree_ndbl_prunned$tip.label
names(bio17_vector) == tree_ndbl_prunned$tip.label
names(intra_var_bio4) == tree_ndbl_prunned$tip.label
names(intra_var_bio17) == tree_ndbl_prunned$tip.label
names(bio4_vector) == tree_fbdl_prunned$tip.label
names(bio17_vector) == tree_fbdl_prunned$tip.label
names(intra_var_bio4) == tree_fbdl_prunned$tip.label
names(intra_var_bio17) == tree_fbdl_prunned$tip.label


###########################
##### BROWNIAN MOTION######
###########################
#ndbl
pgls_recon_bio17_ndbl = ace(x=bio17_vector, phy=tree_ndbl_prunned, type="continuous", method = "GLS", CI = TRUE, model="BM", marginal=FALSE, corStruct = corBrownian(1, phy = tree_ndbl_prunned))
pgls_recon_bio4_ndbl = ace(x=bio4_vector, phy=tree_ndbl_prunned, type="continuous", method = "GLS", CI = TRUE, model="BM", marginal=FALSE, corStruct = corBrownian(1, phy = tree_ndbl_prunned))
#fbdl
pgls_recon_bio17_fbdl = ace(x=bio17_vector, phy=tree_fbdl_prunned, type="continuous", method = "GLS", CI = TRUE, model="BM", marginal=FALSE, corStruct = corBrownian(1, phy = tree_fbdl_prunned))
pgls_recon_bio4_fbdl = ace(x=bio4_vector, phy=tree_fbdl_prunned, type="continuous", method = "GLS", CI = TRUE, model="BM", marginal=FALSE, corStruct = corBrownian(1, phy = tree_fbdl_prunned))


###########################
###### OU #################
###########################

####### OU fbdl in Comare 4.6 #########
#extract species names
species_names = paste("Pinus_", climate_medians$species, sep="")

#write species names
write.table(species_names, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction/compare4.6_data/compare_species_names.txt", sep="\t", col.names = FALSE, row.names =  FALSE, quote=FALSE) #We don't want quotes because these files will be copied to compare. quote: a logical value (‘TRUE’ or ‘FALSE’) or a numeric vector.  If ‘TRUE’, any character or factor columns will be surrounded by double quotes.  If a numeric vector, its elements are taken as the indices of columns to quote.  In both cases, row and column names are quoted if they are written.  If ‘FALSE’, nothing is quoted.

#bind climatic data to SE estimates (0 in our cases for all speices because we have not intraspecies data)
climate_medians$species == climate_se$species
compare_bio4 = paste(climate_medians$median_bio4, "<", climate_se$se_bio4, ">", sep="")
compare_bio17 = paste(climate_medians$median_bio17, "<", climate_se$se_bio17, ">", sep="")

#write variables
write.table(compare_bio4, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction/compare4.6_data/compare_bio4.txt", sep="\t", col.names = FALSE, row.names =  FALSE, quote=FALSE) #YOU HAVE TO DROP " FROM THE FILE
write.table(compare_bio17, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction/compare4.6_data/compare_bio17.txt", sep="\t", col.names = FALSE, row.names =  FALSE, quote=FALSE) #YOU HAVE TO DROP " FROM THE FILE

#save the tree prunned for compare 4.6
write.tree(tree_fbdl_prunned, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction/compare4.6_data/FBDl_MCC_commAnc_prunned.tree") #You have to copy hasta el ";", no etas espacio, sino no funciona


#En "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/software/Comp46bExec" pinchas startForm.class para activar compare. Le das a main window. 

#Copias todos estos archivos creados con al código de justo arriba en Taxon names, Taxon means (rasgos) y Enter Phylogeny. Indicas que son 112 taxa, 1 rasgo y 1 filogenia. Hay que indicar que si queremos usar los SE dentro de especie ó asumir que la variación dentro de especie es desconocida. En este caso incluímos los valores de SE como variabilidad intraespecífica. Seleccionar PGLS-ancestros y exeecute.   

#Correr modelo: 
    #BM: Linear model which is the option by default. Not select specyfing alpha. Nothing else. 
    #OU:Luego Exponential model, specyfing Alpha, pon el valor de alfa de bio4  ó bio17 redondeados a dos decimales (0.076 y 0.036 respectivamente), 100 iteraciones y run (he comprobado el resultado con 1000 interaciones en ambas variablws y sale exactamente lo mismo). Asú se corre un OU en comapre4.6. Esto se ha seguido de Guerrero et al.,... Wiens ., 2013 ("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3710863/")

#Le das a guardar en un archivo llamado "compare_res_bioX_XX.txt". De ese archivo copias la parte de "Trait #1: Ancestral state estimates" y la pegas en un excel, le das a pegar con el importador de datos (se hace en el boton que surge al pegar como el de mantener-quitar formato). Así te separará cada columna. Solo falta añadir a "Adj." el "SE" que queda en la siguiente columna (es SE adjusted) y guardar como .csv.

#load results of OU with SE intraespecífica (alpha = 0.07 for bio4 and 0.03 for bio17)
anc_ou_bio4_fbdl = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction/compare_results/anc_bio4_ou.csv", sep=",", header=TRUE)
str(anc_ou_bio4_fbdl)
anc_ou_bio17_fbdl = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction/compare_results/anc_bio17_ou.csv", sep=",", header=TRUE)
str(anc_ou_bio17_fbdl)

#change names of root by 112+1 (as notation of ape)
#bio4
first_node_fbdl = length(tree_fbdl_prunned$tip.label) + 1
levels(anc_ou_bio4_fbdl$Node) = c(levels(anc_ou_bio4_fbdl$Node), eval(first_node_fbdl))
anc_ou_bio4_fbdl$Node[which(anc_ou_bio4_fbdl$Node == "Root")] <-  length(tree_fbdl_prunned$tip.label) + 1
anc_ou_bio4_fbdl$Node = droplevels(anc_ou_bio4_fbdl$Node)
#bio17
first_node_fbdl = length(tree_fbdl_prunned$tip.label) + 1
levels(anc_ou_bio17_fbdl$Node) = c(levels(anc_ou_bio17_fbdl$Node), eval(first_node_fbdl))
anc_ou_bio17_fbdl$Node[which(anc_ou_bio17_fbdl$Node == "Root")] <-  length(tree_fbdl_prunned$tip.label) + 1
anc_ou_bio17_fbdl$Node = droplevels(anc_ou_bio17_fbdl$Node)

#change the name of Node to nodo1 for mergin with data.frame of node-species numbers
colnames(anc_ou_bio4_fbdl)[which(colnames(anc_ou_bio4_fbdl) == "Node")] <- "nodo1" 
colnames(anc_ou_bio17_fbdl)[which(colnames(anc_ou_bio17_fbdl) == "Node")] <- "nodo1"


###### OU ndbl with geiger
require(geiger)

##bio4
fitContinuous(phy = tree_ndbl_prunned, dat = bio4_vector, SE=0, model="OU", control = list(niter = 1000, CI = 0.95, method = c("subplex","L-BFGS-B")), bounds=list(SE=c(0,0.5)), ncores=2) #alpha = 0.100792 (0.100 for compare)
fitContinuous(phy = tree_ndbl_prunned, dat = bio4_vector, SE=NA, model="OU", control = list(niter = 1000, CI = 0.95, method = c("subplex","L-BFGS-B")), bounds=list(SE=c(0,0.5)), ncores=2) #alpha = 0.100792 (0.100 for compare)
fitContinuous(phy = tree_ndbl_prunned, dat = bio4_vector, SE=intra_var_bio4, model="OU", control = list(niter = 1000, CI = 0.95, method = c("subplex","L-BFGS-B")), bounds=list(SE=c(0,0.5)), ncores=2) #alpha = 0.100795 (0.100 for compare)

##bio17
fitContinuous(phy = tree_ndbl_prunned, dat = bio17_vector, SE=0, model="OU", control = list(niter = 1000, CI = 0.95, method = c("subplex","L-BFGS-B")), bounds=list(SE=c(0,0.5)), ncores=2) #alpha = 0.046440 (0.046 for compare)
fitContinuous(phy = tree_ndbl_prunned, dat = bio17_vector, SE=NA, model="OU", control = list(niter = 1000, CI = 0.95, method = c("subplex","L-BFGS-B")), bounds=list(SE=c(0,0.5)), ncores=2) #alpha = 0.046439 (0.046 for compare)
fitContinuous(phy = tree_ndbl_prunned, dat = bio17_vector, SE=intra_var_bio4, model="OU", control = list(niter = 1000, CI = 0.95, method = c("subplex","L-BFGS-B")), bounds=list(SE=c(0,0.5)), ncores=2) #alpha = 0.046410 (0.046 for compare)

####### OU nbdl in Comare 4.6 #########
#extract species names
species_names = paste("Pinus_", climate_medians$species, sep="")

#write species names
write.table(species_names, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction_ndbl/compare4.6_data/compare_species_names.txt", sep="\t", col.names = FALSE, row.names =  FALSE, quote=FALSE) #We don't want quotes because these files will be copied to compare. quote: a logical value (‘TRUE’ or ‘FALSE’) or a numeric vector.  If ‘TRUE’, any character or factor columns will be surrounded by double quotes.  If a numeric vector, its elements are taken as the indices of columns to quote.  In both cases, row and column names are quoted if they are written.  If ‘FALSE’, nothing is quoted.

#bind climatic data to SE estimates (0 in our cases for all speices because we have not intraspecies data)
climate_medians$species == climate_se$species
compare_bio4 = paste(climate_medians$median_bio4, "<", climate_se$se_bio4, ">", sep="")
compare_bio17 = paste(climate_medians$median_bio17, "<", climate_se$se_bio17, ">", sep="")

#write variables
write.table(compare_bio4, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction_ndbl/compare4.6_data/compare_bio4.txt", sep="\t", col.names = FALSE, row.names =  FALSE, quote=FALSE) #YOU HAVE TO DROP " FROM THE FILE
write.table(compare_bio17, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction_ndbl/compare4.6_data/compare_bio17.txt", sep="\t", col.names = FALSE, row.names =  FALSE, quote=FALSE) #YOU HAVE TO DROP " FROM THE FILE

#save the tree prunned for compare 4.6
write.tree(tree_ndbl_prunned, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction_ndbl/compare4.6_data/NDbl_MCC_commAnc_prunned.tree") #You have to copy hasta el ";", no etas espacio, sino no funciona


#En "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/software/Comp46bExec" pinchas startForm.class para activar compare. Le das a main window. 

#Copias todos estos archivos creados con al código de justo arriba en Taxon names (nombres), Taxon means (rasgos) y Enter Phylogeny. Indicas que son 112 taxa, 1 rasgo y 1 filogenia. Hay que indicar que si queremos usar los SE dentro de especie ó asumir que la variación dentro de especie es desconocida. En este caso incluímos los valores de SE como variabilidad intraespecífica. Seleccionar PGLS-ancestros y exeecute.   

#Correr modelo: 
    #BM: Linear model which is the option by default. Not select specyfing alpha. Nothing else. 
    #OU:Luego Exponential model, specyfing Alpha, pon el valor de alfa de bio4  ó bio17 redondeados a dos decimales (0.08 y 0.03 respectivamente), 100 iteraciones y run (he comprobado el resultado con 1000 interaciones en ambas variablws y sale exactamente lo mismo). Asú se corre un OU en comapre4.6. Esto se ha seguido de Guerrero et al.,... Wiens ., 2013 ("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3710863/")

#Le das a guardar en un archivo llamado "compare_res_bioX_XX.txt". De ese archivo copias la parte de "Trait #1: Ancestral state estimates" y la pegas en un excel, le das a pegar con el importador de datos (se hace en el boton que surge al pegar como el de mantener-quitar formato). Así te separará cada columna. Solo falta añadir a "Adj." el "SE" que queda en la siguiente columna (es SE adjusted) y guardar como .csv.

#load results of OU with SE intraespecífica (alpha = 0.1 for bio4 and 0.046 for bio17)
anc_ou_bio4_ndbl = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction_ndbl/compare_results/anc_bio4_ou_ndbl.csv", sep=",", header=TRUE)
str(anc_ou_bio4_ndbl)
anc_ou_bio17_ndbl = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction_ndbl/compare_results/anc_bio17_ou_ndbl.csv", sep=",", header=TRUE)
str(anc_ou_bio17_ndbl)


#change names of root by 112+1 (ass notation of ape)
#bio4
first_node_ndbl = length(tree_ndbl_prunned$tip.label) + 1
levels(anc_ou_bio4_ndbl$Node) = c(levels(anc_ou_bio4_ndbl$Node), eval(first_node_ndbl))
anc_ou_bio4_ndbl$Node[which(anc_ou_bio4_ndbl$Node == "Root")] <-  length(tree_ndbl_prunned$tip.label) + 1
anc_ou_bio4_ndbl$Node = droplevels(anc_ou_bio4_ndbl$Node)
#bio17
first_node_ndbl = length(tree_ndbl_prunned$tip.label) + 1
levels(anc_ou_bio17_ndbl$Node) = c(levels(anc_ou_bio17_ndbl$Node), eval(first_node_ndbl))
anc_ou_bio17_ndbl$Node[which(anc_ou_bio17_ndbl$Node == "Root")] <-  length(tree_ndbl_prunned$tip.label) + 1
anc_ou_bio17_ndbl$Node = droplevels(anc_ou_bio17_ndbl$Node)

#change the name of Node to nodo1 for mergin with data.frame of node-species numbers
colnames(anc_ou_bio4_ndbl)[which(colnames(anc_ou_bio4_ndbl) == "Node")] <- "nodo1" 
colnames(anc_ou_bio17_ndbl)[which(colnames(anc_ou_bio17_ndbl) == "Node")] <- "nodo1"


###plot ancestral states of bio17 and bio4 under both phylgenies
#correlation tests ou
test_bio17_ou = cor.test(anc_ou_bio17_ndbl$State, anc_ou_bio17_fbdl$State, method="spearman")
test_bio4_ou = cor.test(anc_ou_bio4_ndbl$State, anc_ou_bio4_fbdl$State, method="spearman")
#correlation test bm
test_bio17_bm = cor.test(pgls_recon_bio17_ndbl$ace, pgls_recon_bio17_fbdl$ace, method="spearman")
test_bio4_bm = cor.test(pgls_recon_bio4_ndbl$ace, pgls_recon_bio4_fbdl$ace, method="spearman")

#figure
pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/phd/nicho_pinus/results/phylo_reconstruction_ndbl/cors_bio4_bio17_ndb_fbdl_v2.pdf")
par(mfrow=c(2,2))

#bio4
plot(pgls_recon_bio4_ndbl$ace~pgls_recon_bio4_fbdl$ace, xlab="BIO4 anc across FBDl phylogeny", ylab="BIO4 anc across NDbl phylogeny")
estimate_bio4_bm = bquote(italic(rho) == .(format(round(test_bio4_bm$estimate,2))))
text(x=9500, y=3600, labels = estimate_bio4_bm, cex=1)
if(test_bio4_bm$p.value < 2.2e-16){
    p_bio4_bm = bquote(italic(p.value) < .(format(2.2e-16)))
}else{
    p_bio4_bm = bquote(italic(p.value) == .(format(round(test_bio4_bm$p.value,4))))
}
text(x=9500, y=2800, labels = p_bio4_bm, cex=1)
s_bio4_bm = bquote(S == .(format(round(test_bio4_bm$statistic,2))))
text(x=9500, y=2000, labels = s_bio4_bm, cex=1)

#title two first pannels
title(main=paste("Ancestral state under Brownian motion", sep=""), outer=TRUE, cex.main=2, font.main= 2, line=-2.5)

#bio17
plot(pgls_recon_bio17_ndbl$ace~pgls_recon_bio17_fbdl$ace, xlab="BIO17 anc across FBDl phylogeny", ylab="BIO17 anc across NDbl phylogeny")
estimate_bio17_bm = bquote(italic(rho) == .(format(round(test_bio17_bm$estimate,2))))
text(x=-250, y=-600, labels = estimate_bio17_bm, cex=1)
if(test_bio17_bm$p.value < 2.2e-16){
    p_bio17_bm = bquote(italic(p.value) < .(format(2.2e-16)))
}else{
    p_bio17_bm = bquote(italic(p.value) == .(format(round(test_bio17_bm$p.value,4))))
}
text(x=-250, y=-650, labels = p_bio17_bm, cex=1)
s_bio17_bm = bquote(S == .(format(round(test_bio17_bm$statistic,2))))
text(x=-250, y=-700, labels = s_bio17_bm, cex=1)

#bio4
plot(anc_ou_bio4_ndbl$State~anc_ou_bio4_fbdl$State, xlab="BIO4 anc across FBDl phylogeny", ylab="BIO4 anc across NDbl phylogeny")
estimate_bio4_ou = bquote(italic(rho) == .(format(round(test_bio4_ou$estimate,2))))
text(x=7000, y=5770, labels = estimate_bio4_ou, cex=1)
if(test_bio4_ou$p.value < 2.2e-16){
    p_bio4_ou = bquote(italic(p.value) < .(format(2.2e-16)))
}else{
    p_bio4_ou = bquote(italic(p.value) == .(format(round(test_bio4_ou$p.value,4))))
}
text(x=7000, y=5500, labels = p_bio4_ou, cex=1)
s_bio4_ou = bquote(S == .(format(round(test_bio4_ou$statistic,2))))
text(x=7000, y=5280, labels = s_bio4_ou, cex=1)

#bio17
plot(anc_ou_bio17_ndbl$State~anc_ou_bio17_fbdl$State, xlab="BIO17 anc across FBDl phylogeny", ylab="BIO17 anc across NDbl phylogeny")
estimate_bio17_ou = bquote(italic(rho) == .(format(round(test_bio17_ou$estimate,2))))
text(x=-400, y=-600, labels = estimate_bio17_ou, cex=1)
if(test_bio17_ou$p.value < 2.2e-16){
    p_bio17_ou = bquote(italic(p.value) < .(format(2.2e-16)))
}else{
    p_bio17_ou = bquote(italic(p.value) == .(format(round(test_bio17_ou$p.value,4))))
}
text(x=-400, y=-635, labels = p_bio17_ou, cex=1)
s_bio17_ou = bquote(S == .(format(round(test_bio17_ou$statistic,2))))
text(x=-400, y=-670, labels = s_bio17_ou, cex=1)

#title two second pannels
title(main=paste("Ancestral state under Ornstein-Uhlenbeck", sep=""), outer=TRUE, cex.main=2, font.main= 2, line = -23.5)
dev.off()