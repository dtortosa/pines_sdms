## load workspace
load("/contMap_error.RData")

## load required packages
require(phytools)

## climate data loaded
str(climate_medians)

## tree loaded
summary(tree)

## species without climatic data in the tree
species_to_drop = tree$tip.label[which(!tree$tip.label %in% paste("Pinus_", climate_medians$species, sep=""))]

## prune the tree of speceis without seed mass data
tree_prunned = drop.tip(tree, species_to_drop)

## save bio17 data in a vector with species names as names
bio17_vector = climate_medians$median_bio17

## set names of these variables as species names
names(bio17_vector) <- paste("Pinus_", climate_medians$species, sep="")

## indicate the same order than tip labels
bio17_vector = bio17_vector[tree_prunned$tip.label]

## create contMap object
obj_bm_bio17 = contMap(tree_prunned, bio17_vector, plot=TRUE, lwd = 2) 

## plot error bars with errorbar.contMap
errorbar.contMap(obj=obj_bm_bio17, scale.by.ci=TRUE) #Error obtained: "Error in obj$cols[ii:jj] : only 0's may be mixed with negative subscripts"

## modify function to add error bars of ancestral state (errorbar.contMap), because of a problem in the calculation of ii and jj. There are problems when the lower limit of CI95 and the lower limit of the range of current values of the trait are both negatives. The function calculates the difference, but the resulting number is negative, so I have added abs()
#Additional modification to add ancestral states previously reconstructed. In "ancestral.states" include an object with the state as ace and the confidence intervals as CI95. Typical ace object. 
errorbar_contMap_modified = function (obj, user=FALSE, anc.states=NULL, ...){
    if (hasArg(x)){
        x <- list(...)$x
    } else{
        x <- setNames(sapply(1:Ntip(obj$tree), function(x, obj) {
        ii <- which(obj$tree$edge[, 2] == x)
        ss <- names(obj$tree$maps[[ii]][length(obj$tree$maps[[ii]])])
        obj$lims[1] + as.numeric(ss)/(length(obj$cols) - 1) * 
            diff(obj$lims)
        }, obj = obj), obj$tree$tip.label)
    }
    if (hasArg(scale.by.ci)) {
        scale.by.ci <- list(...)$scale.by.ci
    } else {
        scale.by.ci <- TRUE
    }
    if (hasArg(lwd)){
        lwd <- list(...)$lwd
    } else {
        lwd <- 14
    }
    tree <- obj$tree
    if(user==FALSE){ #MODIFIED LINE
        aa <- fastAnc(tree, x, CI = TRUE)
    } else {
        aa = anc.states
    }
    xlim <- range(aa$CI95)
    if (xlim[2] > obj$lims[2] || xlim[1] < obj$lims[1]) {
        cat(paste("  -----\n  The range of the contMap object, presently (", 
            round(obj$lims[1], 4), ",", round(obj$lims[2], 4), 
            "), should be equal to\n  or greater than the range of the CIs on ancestral states: (", 
            round(xlim[1], 4), ",", round(xlim[2], 4), ").\n", 
            sep = ""))
        cat(paste("  To ensure that your error bars are correctly plotted, please recompute your\n", 
            "  contMap object and increase lims.\n  -----\n", 
            sep = ""))
    }
    d <- diff(obj$lims)
    if (scale.by.ci) {
        v <- aa$CI95[, 2] - aa$CI95[, 1]
        v <- v/max(v)
    } else {
        v <- rep(0.5, tree$Nnode)
    }    
    n <- length(obj$cols) - 1
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    h <- max(nodeHeights(tree))
    for (i in 1:tree$Nnode) {
        ii <- round((abs(aa$CI95[i, 1] - obj$lims[1]))/d * n) #MODIFIED LINE
        jj <- round((abs(aa$CI95[i, 2] - obj$lims[1]))/d * (n + 1)) #MODIFIED LINE
        cols <- obj$cols[ii:jj]
        add.color.bar(leg = 0.1 * h * v[i], cols = cols, prompt = FALSE, 
            x = lastPP$xx[i + Ntip(tree)] - 0.05 * h * v[i], 
            y = lastPP$yy[i + Ntip(tree)], title = "", subtitle = "", 
            lims = NULL, lwd = lwd)
    }
}
errorbar_contMap_modified(obj=obj_bm_bio17, scale.by.ci=TRUE)