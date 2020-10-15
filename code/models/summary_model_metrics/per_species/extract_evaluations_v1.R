#set working directory
setwd("/Users/dsalazar/nicho_pinus/results/final_analyses/")

#required packages
require(dismo)
require(raster)
require(randomForest)

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

#empty data frame to save results 
medians_evaluations = as.data.frame(matrix(NA, nrow=length(epithet_species_list), ncol=17))
colnames(medians_evaluations) <- c("species", "glm_kappa_median", "glm_kappa_sd", "glm_tss_median", "glm_tss_sd", "glm_auc_median", "glm_auc_sd", "gam_kappa_median", "gam_kappa_sd", "gam_tss_median", "gam_tss_sd", "gam_auc_median", "gam_auc_sd", "rf_oob_median", "rf_oob_sd", "rf_auc_median", "rf_auc_sd")
str(medians_evaluations)

#loop for obtaining medians of evaluations across species
for(i in 1:length(epithet_species_list)){

    #select the [i] species
    selected_species = epithet_species_list[i]
    
    #save the name of the [i] species
    medians_evaluations[i, "species"] <- selected_species

    #unzip data of threshold
    unzip(paste("threshold/threshold_", selected_species, ".zip", sep=""), exdir="threshold")

    #unzip data of AUC
    unzip(paste("evaluations/evaluations_", selected_species, ".zip", sep=""), exdir="evaluations")

    #unzip data of modeling (only for RF, calculate OOB)
    unzip(paste("models/models_", selected_species, ".zip", sep=""), exdir="models")

    ####################################################
    ###### ALL evaluation metrics for GLM and GAM ######
    ####################################################
    models_to_extract = c("glm", "gam")
    for(m in 1:length(models_to_extract)){

        #selected the [m] model
        selected_model = models_to_extract[m]

        #### TSS and Kappa ###
        #load threshold data for the [m] model
        load(paste("threshold/", selected_species, "_", selected_model, "_threshold.rda", sep=""))

        #put into an object kappa and tss data
        kappa_tss_data = eval(parse(text=paste(selected_model, "_threshold", sep="")))

        #extract values of kappa and TSS for all data partitions
        tss_values = sapply(kappa_tss_data, '$.data.frame', "TSS")[1,]       
        kappa_values = sapply(kappa_tss_data, '$.data.frame', "kappa")[1,]       

        #save the median and sd of kappa  across partitions of the data
        medians_evaluations[i, paste(selected_model, "_kappa_median", sep="")] <- median(kappa_values)
        medians_evaluations[i, paste(selected_model, "_kappa_sd", sep="")] <- sd(kappa_values)
            
        #save the median and sd of TSS   across partitions of the data              
        medians_evaluations[i, paste(selected_model, "_tss_median", sep="")] <- median(tss_values)
        medians_evaluations[i, paste(selected_model, "_tss_sd", sep="")] <- sd(tss_values)

        #### AUC ####
        #load AUC data for the [m] model 
        load(paste("evaluations/", selected_species, "_", selected_model, "_evaluation.rda", sep=""))

        #put into an object AUC data of the [m] model
        evaluation_model = eval(parse(text=paste(selected_model, "_evaluation", sep="")))

        #extract AUC data for each data partition
        auc_values = sapply(evaluation_model, "slot", "auc")

        #save the median and sd of AUC for [m] model across partitions of the data
        medians_evaluations[i, paste(selected_model, "_auc_median", sep="")] <- median(auc_values)
        medians_evaluations[i, paste(selected_model, "_auc_sd", sep="")] <- sd(auc_values)     
    }

    ###########################################
    ###### ALL evaluation metrics for RF ######
    ###########################################
    
    ### AUC ###
    #load the data of AUC for RF
    load(paste("evaluations/", selected_species, "_rf_evaluation.rda", sep=""))

    #save into an object AUC data of RF
    evaluation_rf = eval(parse(text=paste("rf_evaluation", sep="")))

    #extract values of AUC across data partitions
    auc_values = sapply(evaluation_rf, "slot", "auc")

    #save median and sd of AUC across partitions of the data
    medians_evaluations[i, paste("rf_auc_median", sep="")] <- median(auc_values)
    medians_evaluations[i, paste("rf_auc_sd", sep="")] <- sd(auc_values)   

    #### OOB ####
    #load the RF model and the training data used for fitting it
    load(paste("models/", selected_species, "_rf_model.rda", sep=""))    
    load(paste("models/", selected_species, "_training_data.rda", sep=""))  

    #calculate global OOB for all the trees fitted in each data partition. We want OOB value for each partition. 
    oob_values = NULL
    #for each partition of the data
    for(j in 1:length(rf_resample)){

        #select the [j] partition of the data
        selected_partition = rf_resample[[j]]

        #select the confusion matrix but without the class.error column
        confusion_matrix = selected_partition$confusion[,-3]

        #from the confusion matrix select those numbers in the upper or the lower triangle of the matrix, which are the fails (real 0 indicated as 1 or real 1 indicated as 0).
        fails = confusion_matrix[which(upper.tri(confusion_matrix) | lower.tri(confusion_matrix))]

        #calculate global OOB for the model as the total number of fails during predictions and the total number of data for the [j] training partition. 
        oob_values = append(oob_values, (sum(fails) / nrow(training_data[[j]]$training_data))*100) #more info about calculations in "http://www.blopig.com/blog/2017/04/a-very-basic-introduction-to-random-forests-using-r/"
    }

    #save the median and sd of OOB values across partitions
    medians_evaluations[i, paste("rf_oob_median", sep="")] <- median(oob_values)
    medians_evaluations[i, paste("rf_oob_sd", sep="")] <- sd(oob_values)   

    #remove unzipped files
    file.remove(list.files("threshold", pattern=".rda", full.names = TRUE))
    file.remove(list.files("evaluations", pattern=".rda", full.names = TRUE))
    file.remove(list.files("models", pattern=".rda", full.names = TRUE))
}

#save the results
write.table(medians_evaluations, "medians_evaluations/medians_evaluations.csv", sep=",", col.names=TRUE, row.names=FALSE)