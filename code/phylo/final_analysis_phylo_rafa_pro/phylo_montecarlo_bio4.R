#load climate data
climate_medians = read.table("/Users/dsalazar/nicho_pinus/data/climate_medians/climate_medians.csv", header=TRUE, sep=",") 
str(climate_medians)

## cargamos el arbol
require(ape)
pine_tree<-read.nexus("/Users/dsalazar/nicho_pinus/data/phylogeny/FBDl_MCC_commAnc.tre") 

## new species icnldued by bianca in tree that we have to drop
species_to_drop = pine_tree$tip.label[which(!pine_tree$tip.label %in% paste("Pinus_", climate_medians$species, sep=""))]

## prune the tree of speceis without seed mass data
tree_prunned = drop.tip(pine_tree, species_to_drop)

## save climatic variables in a vector with species names as names
bio4_vector = climate_medians$median_bio4
bio17_vector = climate_medians$median_bio17

## set names of these variables as species names
names(bio4_vector) <- paste("Pinus_", climate_medians$species, sep="")
names(bio17_vector) <- paste("Pinus_", climate_medians$species, sep="")

## indicate the same order than tip labels
bio4_vector = bio4_vector[tree_prunned$tip.label]
bio17_vector = bio17_vector[tree_prunned$tip.label]

#check order names
tree_prunned$tip.label == names(bio4_vector)
tree_prunned$tip.label == names(bio17_vector)

################################################
############  PHYLOGENETIC MONTECARLO ##########
################################################

#Package for phylogenetic montecarlo
require(pmc) #Con este paquete vamos a correr un "Phylogenetic Monte Carlo" para la selección entre dos modelos: ModeloA más simple y modeloB más complejo. Primero, los parámetros de ambos modelos son estimados usando los datos originales. Entonces, se simulan n datasets (1000 en nuestro caso como hacen en el paper del paquete) siguiendo la evolución que dicta cada modelo con los parámetros que hemos obtenido previamente. Así obtrendemos un rasgo que evoluciona en nuestra filogenia bajo un modelo BM con el mismo sigma que el de WP, y otro rasgo que evoluciona en nuestra filogenia con el mismo signa, theta y alpha que el de WP. Con cada rasgo se reestiman los paramatroes de los DOS modelos (BM y OU), es decir, el rasgo que evoluciona segun BM se usa para ajustar un BM y un OU, mientras que el rasgo que evoluciona bajo OU se usa para ajustar otro BM y otro OU. Se hace un likelihood ratio test entre ambos modelos para cada rasgo, es decir obtendríamos dos valores de LRT en cada simulación: El valor null, ó hipotesis nula, que sería la diferencia de likelihood entre modelo BM y OU asjutados con un rasgo simulado bajo BM; y el valor test ó nuestro test de interés, que sería la diferencia de likelihood entre modelo BM y OU ambos ajustados con un rasgo que sigue evolución OU. Cuanto mayor sea LRT, más apoyo para el segundo modelo, más complejo, OU en nuestro caso. Esperaríamos, que si el rasgo evouciona por OU y no hay un bias del árbol a favor del modelo más complejo, el LRT del rasgo que evoluciona bajo BM será más bajo que el del rasgo que evoluciona bajo OU, es decir, el modelo OU no es mejor para el caso del rasgo que evoluciona bajo BM. Al final tendremos una distribución de LRT bajo BM y OU que podremos comparar. Ojo al detalle que este procesi implica 4 ajusted por maxima verosimiltud (maximum likelihood), mientras que para un AIC solo usas dos. La figura 2 de "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS" explica muy bien esto. 

############################
#### Geiger - OU.1 BIO4 ####
############################
#OU con un optimo, en este caso no tenemos estado depedencia
require(geiger)

### extract WP data with species as row.names
dat = climate_medians[,which(colnames(climate_medians) %in% c("species", "median_bio4"))]
row.names(dat) <- paste("Pinus_", dat$species, sep="")
dat[,which(colnames(dat) == "species")] <- NULL
str(dat)

### bind tree and data into tmp
tmp = treedata(tree_prunned, dat)

### extract phylogeny
phy = tmp$phy

### extract WP data with species as row.names
datos = tmp$data

### run the phylo montecarlo
simulations_ou.1_bio4 = pmc(phy, datos, "BM", "OU", nboot=2000, mc.cores=4)

### save simulations
save(simulations_ou.1_bio4, file="//Users/dsalazar/nicho_pinus/results/rdata/bio_4_geiger_BM_OU.1_nboot_2000.rda")

### load it
#load("/Users/dsalazar/nicho_pinus/results/rdata/bio_4_geiger_BM_OU.1_nboot_2000.rda")


### plot likelihood ratio test between BM (modelA) and OU.1 (modelB) models fitted with data simulated under BM and OU.1 respectively
pdf("/Users/dsalazar/nicho_pinus/results/pmc/figures/bio_4_geiger_BM_OU_1_nboot_2000.pdf")
par(mfrow=c(2,2))

## alpha values for data simulated under OU.1
#comparison BB (OU.1 model fitted with OU.1-simulated data) and the value of alpha (selection strength)
lr_ou = simulations_ou.1_bio4$test

#calculate density distribution
density_lr_OU = density(lr_ou)

#plot
plot(density_lr_OU, xlim=c(-0.5,75), ylim=c(0,0.9), main="Likelihood ratio test")

#add color to the full area under the curve
x1 <- min(which(density_lr_OU$x >= 0))  
x2 <- max(which(density_lr_OU$x <  75))
with(density_lr_OU, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("black",alpha.f=0.65)))

#plot 95CI
#quantiles_ou = quantile(lr_ou, prob=c(0.025, 0.975))
#abline(v=c(quantiles_ou[1], quantiles_ou[2]), col="black", lty=2) #No los añado porque podría haber problemas para calcularlos así y no encuentro la función de pmc para ello: "Given the noisy nature of parameters estimated from phylogenies, we recommend that confidence interval should routinely be reported, and to facilitate this, have implemented this as pmc::confidenceIntervals.pow. Confidence intervals could also be estimated from the curvature of the likelihood surface, but these can be unreliable and problematic to compute." From "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS". 


## alpha values for data simulated under BM
#comparison AB (OU.1 model fitted with BM-simulated data) and the value of alpha (selection strength)
lr_bm = simulations_ou.1_bio4$null

#calculate density distribution
density_lr_BM = density(lr_bm)

#add density plot to the previous plot
lines(density_lr_BM)

#add color to the full area under the curve
x1 <- min(which(density_lr_BM$x >= -0.8))  
x2 <- max(which(density_lr_BM$x <  13))
with(density_lr_BM, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("gray",alpha.f=0.65)))

#plot 95CI
#quantiles_BM = quantile(lr_bm, prob=c(0.025, 0.975))
#abline(v=c(quantiles_BM[1], quantiles_BM[2]), col="black", lty=2)

#add legen
#legend(53, 0.85, legend=c("BM", "OU.1"), fill=c("gray", adjustcolor("black",alpha.f=0.65)))

#add 95% LR test
abline(v=c(simulations_ou.1_bio4$lr), col="black", lty=5, lwd=2) #Este valor tiene que estar por encima del percentil 95 del LRT del modeloA (null) y por debajo del percentil 95 del LRT del modeloB (test). Mira la sección Methods-Model selection en "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS".
abline(v=c(quantile(lr_ou, probs=0.95)), col="black", lty=1)
abline(v=c(quantile(lr_bm, probs=0.95)), col="black", lty=1) #La proporción de valores simulados mayore que LRT de los datos reales nos da una especie de P.value para el test de selección de modelos. La probabilidad de que la diferencia de likelihood entre modelos (LRT) observada se de bajo el modelo cero. Si cogemos el valor de LRT que es mayor que el 95% de los LRT bajo modeloA, si el LRT observado es mayor que ese threshold, entonces es poco probable que el valor obtenido se haya dado bajo el modeloA (más simple). También podemos calcular el poder del test, la probabilidad de acertar rechazando el modeloA porque los datos vienen del B. Para eso tenemos que usar la distribución de LRTs bajo modelo 1. Si como antes cogemos el valor de LRT mayor que el 95% de LRT simulados bajo modeloB, la cantidad de distribución que queda a la izquierda de ese threshold se aproxima a la probabiidad de rechaza el mdoeloA cuando los datos son producidos por el modeloB. 


### plot alpha between data simulated under BM (modelA) and OU.1 (modelB)

## alpha values for data simulated under OU.1
#comparison BB (OU.1 model fitted with OU.1-simulated data) and the value of alpha (selection strength)
alpha_ou = simulations_ou.1_bio4$par_dists$value[which(simulations_ou.1_bio4$par_dists$comparison == "BB" & simulations_ou.1_bio4$par_dists$parameter == "alpha")] 

#calculate density distribution
density_OU = density(alpha_ou)

#plot
plot(density_OU, xlim=c(0,0.31), ylim=c(0,120), , main="Alpha")

#add color to the full area under the curve
x1 <- min(which(density_OU$x >= 0))  
x2 <- max(which(density_OU$x <  0.2))
with(density_OU, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("black",alpha.f=0.65)))

#plot 95CI
#quantiles_ou = quantile(alpha_ou, prob=c(0.025, 0.975))
#abline(v=c(quantiles_ou[1], quantiles_ou[2]), col="black", lty=2) #No los añado porque podría haber problemas para calcularlos así y no encuentro la función de pmc para ello: "Given the noisy nature of parameters estimated from phylogenies, we recommend that confidence interval should routinely be reported, and to facilitate this, have implemented this as pmc::confidenceIntervals.pow. Confidence intervals could also be estimated from the curvature of the likelihood surface, but these can be unreliable and problematic to compute." From "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS". 

## alpha values for data simulated under BM
#comparison AB (OU.1 model fitted with BM-simulated data) and the value of alpha (selection strength)
alpha_BM = simulations_ou.1_bio4$par_dists$value[which(simulations_ou.1_bio4$par_dists$comparison == "AB" & simulations_ou.1_bio4$par_dists$parameter == "alpha")]

#calculate density distribution
density_BM = density(alpha_BM)

#add density plot to the previous plot
lines(density_BM)

#add color to the full area under the curve
x1 <- min(which(density_BM$x >= -0.2))  
x2 <- max(which(density_BM$x <  0.2))
with(density_BM, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=adjustcolor("gray",alpha.f=0.65)))

#plot 95CI
#quantiles_BM = quantile(alpha_bm, prob=c(0.025, 0.975))
#abline(v=c(quantiles_BM[1], quantiles_BM[2]), col="black", lty=2)

#add legen
legend(0.205, 117, legend=c("BM", "OU.1"), fill=c("gray", adjustcolor("black",alpha.f=0.65)))

#add 95% LR test
abline(v=c(simulations_ou.1_bio4$B$opt$alpha), col="black", lty=5, lwd=2) #Este valor tiene que estar por encima del percentil 95 del LRT del modeloA (null) y por debajo del percentil 95 del LRT del modeloA (test). Mira la sección Methods-Model selection en "IS YOUR PHYLOGENY INFORMATIVE? MEASURING THE POWER OF COMPARATIVE METHODS".
#abline(v=c(quantile(alpha_ou, probs=0.95)), col="black", lty=1)
#abline(v=c(quantile(alpha_BM, probs=0.95)), col="black", lty=1) #Quito estas lineas porque no se si tiene sentido aplicar el approach de poner el 95% de cada parametros en cada modelo, el paper de coop usan este approach con el likelihood ratio test nada más. 
dev.off()