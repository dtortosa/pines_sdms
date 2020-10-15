#SEWAL. Code for modelling and project into the future. It is prepared for run in sewall.  

###################################
#ESTABLECE EL DIRECTORIO DE TRABAJO
###################################
#DIRECTORIO DE TRABAJO
setwd("/Users/dsalazar/nicho_pinus/")


#function for binarize the current projections  
#I have separated the binarization from the fitting of modelling because I have problems in gam when I load ecospat package, which is used in caculates of threshold.      
binarize_curent_projections = function(species){

    #required libraries
    #library(raster)
    #library(dismo)
    #library(ecospat) #for threshold
    #library(gtools) #for mixedsort function, which order elements of a list by numbers. Used in pre-projections

    #begin the species
    print(paste("begin",species))

    #load models 
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/models", paste(species, "glm_model.rda", sep="_"), sep="/"))
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/models", paste(species, "gam_model.rda", sep="_"), sep="/"))
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/models", paste(species, "rf_model.rda", sep="_"), sep="/"))  

    #load evaluations
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", paste(species, "glm_evaluation.rda", sep="_"), sep="/"))
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", paste(species, "gam_evaluation.rda", sep="_"), sep="/"))
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", paste(species, "rf_evaluation.rda", sep="_"), sep="/")) 

    #load evaluation data
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", paste(species, "glm_evaluation_data.rda", sep="_"), sep="/"))
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", paste(species, "gam_evaluation_data.rda", sep="_"), sep="/"))
    load(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", paste(species, "rf_evaluation_data.rda", sep="_"), sep="/")) 

    #load continuous predictions of glm 
    glm_predict = stack(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/continuous_predictions/continuous_predictions_glm", paste(species, "tif", sep="."), sep="_"))

    #load continuous predictions of gam     
    gam_predict = stack(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/continuous_predictions/continuous_predictions_gam", paste(species, "tif", sep="."), sep="_"))

    #load binary predictions of Random forest (already finished)
    rf_predict = stack(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_predictions/binary_predictions_rf", paste(species, "tif", sep="."), sep="_"))

    #lists for saving thresholds
    glm_threshold = list()
    gam_threshold = list()

    #lists for saving binary projections
    glm_prediction_bin = list()
    gam_prediction_bin = list()
    rf_prediction_bin = list()

    for (k in 1:12){

        ###################################
        ############Threshold #############
        ################################### 

        #calculate the best threshold according to kappa
        glm_kappa = ecospat.max.kappa(glm_evaluation_predict[[k]]$predict_glm, glm_evaluation_predict[[k]]$presence)[[2]]
        gam_kappa = ecospat.max.kappa(gam_evaluation_predict[[k]]$predict_gam, gam_evaluation_predict[[k]]$presence)[[2]] #we don't have to binarize the predictions of random forest model because is a clasificaction tree, the predictions are binarized. 

        #calculate the best threshold according to TSS
        glm_tss = ecospat.max.tss(glm_evaluation_predict[[k]]$predict_glm, glm_evaluation_predict[[k]]$presence)[[2]]
        gam_tss = ecospat.max.tss(gam_evaluation_predict[[k]]$predict_gam, gam_evaluation_predict[[k]]$presence)[[2]]

        #save all 
        glm_optim_threshold = as.data.frame(cbind(rbind(as.numeric(glm_kappa[,2][1]),as.numeric(glm_kappa[,2][2])), rbind(as.numeric(glm_tss[,2][1]),as.numeric(glm_tss[,2][2]))))
        names(glm_optim_threshold) = c("kappa", "TSS")
        row.names(glm_optim_threshold) = c("maximum value", "threshold")

        gam_optim_threshold = as.data.frame(cbind(rbind(as.numeric(gam_kappa[,2][1]),as.numeric(gam_kappa[,2][2])), rbind(as.numeric(gam_tss[,2][1]),as.numeric(gam_tss[,2][2]))))
        names(gam_optim_threshold) = c("kappa", "TSS")
        row.names(gam_optim_threshold) = c("maximum value", "threshold")

        glm_threshold[[k]] = glm_optim_threshold
        gam_threshold[[k]] = gam_optim_threshold

        ###################################
        ######Binarize predictions#########
        ###################################

        #convert each projection in binary using the threshold calculated with fit_eval_models (TSS)
        #glm
        glm_prediction_bin[[k]] = glm_predict[[k]] #copy the raster 
        glm_prediction_bin[[k]][glm_prediction_bin[[k]]>=(glm_threshold[[k]][2,2]),] <- 1 #give 1 to the pixels with a predicted value higher or equal than the threshold
        glm_prediction_bin[[k]][glm_prediction_bin[[k]]<(glm_threshold[[k]][2,2]),] <- 0 #give 0 to the pixels with a predicted value lower than the threshold

        #gam
        gam_prediction_bin[[k]] = gam_predict[[k]] #copy the raster 
        gam_prediction_bin[[k]][gam_prediction_bin[[k]]>=(gam_threshold[[k]][2,2]),] <- 1 #give 1 to the pixels with a predicted value higher or equal than the threshold
        gam_prediction_bin[[k]][gam_prediction_bin[[k]]<(gam_threshold[[k]][2,2]),] <- 0 #give 0 to the pixels with a predicted value lower than the threshold

        #rf
        rf_prediction_bin[[k]] = rf_predict[[k]] #we don't have to binarize the predictions of this model because is a clasificaction tree, the predictions are not continuous probabilities, on the contrary, they are 0 or 1 (binary) 
    }

    #stack all binary predictions 
    binary_predictions = stack(stack(glm_prediction_bin), stack(gam_prediction_bin), stack(rf_prediction_bin))    

    ##create the ensamble
    ensamble_predictions_bin = calc(binary_predictions, function(x) (sum(x)*100)/nlayers(binary_predictions)) #calculate the percentage of models for which a pixel is suitable

    #save threshold
    save(glm_threshold, file=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/threshold", paste(species, "glm_threshold.rda", sep="_"), sep="/"))
    save(gam_threshold, file=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/threshold", paste(species, "gam_threshold.rda", sep="_"), sep="/"))

    #save the binary predictions of glm and gam (rf make them before)
    writeRaster(stack(glm_prediction_bin), filename=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_predictions/binary_predictions_glm", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
    writeRaster(stack(gam_prediction_bin), filename=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_predictions/binary_predictions_gam", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)

    #save ensamble
    writeRaster(ensamble_predictions_bin, filename=paste("/Users/dsalazar/nicho_pinus/results/final_analyses/ensamble_predictions_bin", paste("ensamble_predictions_bin", paste(species, "tif", sep="."), sep="_"), sep="/"), overwrite=TRUE)

    #########################################
    ######ZIP AND DELETE EVALUATIONS ########
    ######################################### 
    ####FOR RELEASE SPACE IN SEAWELL
    #list continuous projections created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files
        evaluations = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", pattern=glob2rx(paste("*", species, "*", "rda", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        evaluations = evaluations[!grepl("pseudostrobus", evaluations)]
    } else {
        evaluations = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations", pattern=glob2rx(paste("*", species, "*", "rda", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE) 
    }

    #zip all continuous projections for lack of space in disk
    zip(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/evaluations/evaluations", paste(species, "zip", sep="."), sep="_"), evaluations, flags="-j") #j indicate that you don't want all the directory structure

    #delete continuous projections
    file.remove(evaluations) 

    ########################################################
    ######ZIP AND DELETE BIN AND CONTINUOUS PREDICTIONS ####
    ######################################################## 
    #list continuous predictions created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files    
        continuous_predictions = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/continuous_predictions", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        continuous_predictions = continuous_predictions[!grepl("pseudostrobus", continuous_predictions)]
    } else {
        continuous_predictions = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/continuous_predictions", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
    } 

    #zip all continuous predictions for lack of space in disk
    zip(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/continuous_predictions/continuous_predictions", paste(species, "zip", sep="."), sep="_"), continuous_predictions, flags="-j") #j indicate that you don't want all the directory structure

    #delete continuous predictions
    file.remove(continuous_predictions)

    #list binary predictions created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files    
        binary_predictions = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_predictions", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        binary_predictions = binary_predictions[!grepl("pseudostrobus", binary_predictions)]        
    } else {
        binary_predictions = list.files("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_predictions", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)        
    }

    #zip all binary predictions for lack of space in disk
    zip(paste("/Users/dsalazar/nicho_pinus/results/final_analyses/binary_predictions/binary_predictions", paste(species, "zip", sep="."), sep="_"), binary_predictions, flags="-j") #j indicate that you don't want all the directory structure

    #delete binary predictions
    file.remove(binary_predictions)

    #name of species
    print(paste(species, "ended")) 
}

########Paralelize the process######
require(foreach)
require(doParallel) #for parallel


#load species names
list_species = read.table("/Users/dsalazar/nicho_pinus/data/list_species.txt", sep="\t", header=TRUE)

#extract epithet from species list
epithet_species_list = NULL
for(i in 1:nrow(list_species)){

    #selected species
    selected_species = as.vector(list_species[i,])

    #extract epithet
    epithet_species_list = append(epithet_species_list, strsplit(selected_species, split=" ")[[1]][2])
}
summary(is.na(epithet_species_list)) #all false

# set up cluster
clust <- makeCluster(6) 
registerDoParallel(clust)


###########################
#########BINARIZE##########
###########################
#binarize current predictions for non-problematic species
foreach(i = epithet_species_list, .packages=c("raster", "dismo", "ecospat", "gtools")) %dopar% { 
    binarize_curent_projections(species = i)
} 


#stop the cluster 
stopCluster(clust)
