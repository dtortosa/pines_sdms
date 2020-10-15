###Code for compare the performance of models with and without soil variables using metrics of evaluation. 

###definimos el directorio de trabajo
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus")

#Librerias
require(raster) #for work with rasters
require(dismo) #for reading evaluations

#load list of species
list_species = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/presences/species.txt", sep="\t", header=T)

#extract epithet from species list
epithet_species_list = NULL
for(i in 1:nrow(list_species)){

    #selected species
    selected_species = as.vector(list_species[i,])

    #extract epithet
    epithet_species_list = append(epithet_species_list, strsplit(selected_species, split=" ")[[1]][2])
}
summary(is.na(epithet_species_list)) #all false

list_results = list()
#loop for model comparison
for(i in 1:length(epithet_species_list)){ #for each species

    #select the [i] species
    selected_species = epithet_species_list[i]

    ### load evaluation with soil
    ## unzip
    path_unzipped_eval_with_soil = unzip(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_final/evaluations/evaluations_", selected_species, ".zip", sep=""), list=FALSE, exdir="/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_final/evaluations/")
    path_unzipped_thres_with_soil = unzip(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_final/thresholds/threshold_", selected_species, ".zip", sep=""), list=FALSE, exdir="/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_final/thresholds/")


    ## load all files
    for(k in 1:length(path_unzipped_eval_with_soil)){ #for each file of evaluations.zip
        load(path_unzipped_eval_with_soil[k], verbose=TRUE) #load the [k] file
    }
    for(k in 1:length(path_unzipped_thres_with_soil)){ #for each file of evaluations.zip
        load(path_unzipped_thres_with_soil[k], verbose=TRUE) #load the [k] file
    }

    ## change names
    gam_evaluation_predict_with_soil <- gam_evaluation_predict
    gam_evaluation_with_soil <- gam_evaluation
    glm_evaluation_predict_with_soil <- glm_evaluation_predict
    glm_evaluation_with_soil <- glm_evaluation
    rf_evaluation_predict_with_soil <- rf_evaluation_predict
    rf_evaluation_with_soil <- rf_evaluation  
    gam_threshold_with_soil <- gam_threshold
    glm_threshold_with_soil <- glm_threshold


    ### load evaluation without soil
    ## unzip
    path_unzipped_eval_without_soil = unzip(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_without_soil/evaluations_without_soil/evaluations_", selected_species, ".zip", sep=""), list=FALSE, exdir="/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_without_soil/evaluations_without_soil/")
    path_unzipped_thres_without_soil = unzip(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_without_soil/thresholds_without_soil/threshold_", selected_species, ".zip", sep=""), list=FALSE, exdir="/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_without_soil/thresholds_without_soil/")


    ## load all files
    for(k in 1:length(path_unzipped_eval_without_soil)){ #for each file of evaluations.zip
        load(path_unzipped_eval_without_soil[k], verbose=TRUE) #load the [k] file
    }
    for(k in 1:length(path_unzipped_thres_without_soil)){ #for each file of evaluations.zip
        load(path_unzipped_thres_without_soil[k], verbose=TRUE) #load the [k] file
    }

    ## change names
    gam_evaluation_predict_without_soil <- gam_evaluation_predict
    gam_evaluation_without_soil <- gam_evaluation
    glm_evaluation_predict_without_soil <- glm_evaluation_predict
    glm_evaluation_without_soil <- glm_evaluation
    rf_evaluation_predict_without_soil <- rf_evaluation_predict
    rf_evaluation_without_soil <- rf_evaluation  
    gam_threshold_without_soil <- gam_threshold
    glm_threshold_without_soil <- glm_threshold


    ### calculate median values 
    ## median values of statistic across 12 data partition. We use sapply, which applies a function across a list and generates a data.frame (lapply generates another list). For evaluations we have a list of slots (extracted with @), whilst for thresholds we have data.frames. En each case indicate the type of object and the statistic to extract. In the case of threshold we obtain the value of maximum TSS and kappa but also the threshold itself, because of this we select the first row of the sapply output (TSS an kappa values, not the threshold). https://stackoverflow.com/questions/30131542/lapply-with-function

    #without soil
    median_gam_auc_without_soil = median(sapply(gam_evaluation_without_soil, 'slot', "auc"))
    median_glm_auc_without_soil = median(sapply(glm_evaluation_without_soil, 'slot', "auc"))
    median_rf_auc_without_soil = median(sapply(rf_evaluation_without_soil, 'slot', "auc"))
    median_gam_tss_without_soil = median(sapply(gam_threshold_without_soil, '$.data.frame', "TSS")[1,])
    median_gam_kappa_without_soil = median(sapply(gam_threshold_without_soil, '$.data.frame', "kappa")[1,])
    median_glm_tss_without_soil = median(sapply(glm_threshold_without_soil, '$.data.frame', "TSS")[1,])
    median_glm_kappa_without_soil = median(sapply(glm_threshold_without_soil, '$.data.frame', "kappa")[1,])

    #with soil
    median_gam_auc_with_soil = median(sapply(gam_evaluation_with_soil, 'slot', "auc"))
    median_glm_auc_with_soil = median(sapply(glm_evaluation_with_soil, 'slot', "auc"))
    median_rf_auc_with_soil = median(sapply(rf_evaluation_with_soil, 'slot', "auc"))
    median_gam_tss_with_soil = median(sapply(gam_threshold_with_soil, '$.data.frame', "TSS")[1,])
    median_gam_kappa_with_soil = median(sapply(gam_threshold_with_soil, '$.data.frame', "kappa")[1,])
    median_glm_tss_with_soil = median(sapply(glm_threshold_with_soil, '$.data.frame', "TSS")[1,])
    median_glm_kappa_with_soil = median(sapply(glm_threshold_with_soil, '$.data.frame', "kappa")[1,])

    #bind both results and also calculate the difference between the,
    results = cbind.data.frame(
        rbind.data.frame(
            median_gam_auc_without_soil,
            median_glm_auc_without_soil,
            median_rf_auc_without_soil,
            median_gam_tss_without_soil,
            median_gam_kappa_without_soil,
            median_glm_tss_without_soil,
            median_glm_kappa_without_soil),   
        rbind.data.frame(
            median_gam_auc_with_soil,
            median_glm_auc_with_soil,
            median_rf_auc_with_soil,
            median_gam_tss_with_soil,
            median_gam_kappa_with_soil,
            median_glm_tss_with_soil,
            median_glm_kappa_with_soil), 
        rbind.data.frame(
            median_gam_auc_with_soil - median_gam_auc_without_soil,
            median_glm_auc_with_soil - median_glm_auc_without_soil,
            median_rf_auc_with_soil - median_rf_auc_without_soil,
            median_gam_tss_with_soil - median_gam_tss_without_soil,
            median_gam_kappa_with_soil - median_gam_kappa_without_soil,
            median_glm_tss_with_soil - median_glm_tss_without_soil,
            median_glm_kappa_with_soil - median_glm_kappa_without_soil))

    #set col and row names
    colnames(results) <- c("without_soil", "with_soil", "with_soil-without_soil")
    rownames(results) <- c("gam_auc", "glm_auc", "rf_auc", "gam_tss", "gam_kappa", "glm_tss", "glm_kappa")
       
    #remove the file extracted   
    file.remove(path_unzipped_eval_without_soil, path_unzipped_thres_without_soil)
    file.remove(path_unzipped_eval_with_soil, path_unzipped_thres_with_soil)

    #save results in the list    
    list_results[[i]] <- results 

    #set the name of [i] species as the name of this element of the list
    names(list_results)[i] <- selected_species
}
length(list_results) 
summary(names(list_results) == epithet_species_list)


#extract the median value of the differences between soil and no soil for all statistics using sapply. WE extract from a data.frame, the data included in the column with_soil, and select the row corresponding to gam_auc, glm_auc....
final_results = cbind.data.frame(
    row.names(list_results[[1]]),
    rbind.data.frame(
        median(sapply(list_results, "$.data.frame", "with_soil-without_soil")[which(row.names(list_results[[1]])=="gam_auc"),]),
        median(sapply(list_results, "$.data.frame", "with_soil-without_soil")[which(row.names(list_results[[1]])=="glm_auc"),]),
        median(sapply(list_results, "$.data.frame", "with_soil-without_soil")[which(row.names(list_results[[1]])=="rf_auc"),]),
        median(sapply(list_results, "$.data.frame", "with_soil-without_soil")[which(row.names(list_results[[1]])=="gam_tss"),]),
        median(sapply(list_results, "$.data.frame", "with_soil-without_soil")[which(row.names(list_results[[1]])=="gam_kappa"),]),
        median(sapply(list_results, "$.data.frame", "with_soil-without_soil")[which(row.names(list_results[[1]])=="glm_tss"),]),
        median(sapply(list_results, "$.data.frame", "with_soil-without_soil")[which(row.names(list_results[[1]])=="glm_kappa"),]))) 

colnames(final_results) <- c("parameters", "median_difference_with-without_soil")
final_results

#save final results
write.table(final_results, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/soil_comprobations/model_evaluation/eval_comparison.csv", sep=",", row.names = FALSE)