#changes respect previous versions
    #we use median and sd in previous versions but this is not correct. You have to use median and interquartile difference or meand and sd. 

#set working directory
setwd("/Users/dsalazar/nicho_pinus/results/final_analyses/")

#required packages
require(dismo)
require(raster)

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
#no NAs
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


#empty data frame to save results 
means_evaluations = as.data.frame(matrix(NA, nrow=length(epithet_species_list), ncol=17))
colnames(means_evaluations) <- c("species", "glm_kappa_mean", "glm_kappa_sd", "glm_tss_mean", "glm_tss_sd", "glm_auc_mean", "glm_auc_sd", "gam_kappa_mean", "gam_kappa_sd", "gam_tss_mean", "gam_tss_sd", "gam_auc_mean", "gam_auc_sd", "rf_kappa_mean", "rf_kappa_sd", "rf_tss_mean", "rf_tss_sd")
str(means_evaluations)

#loop for obtaining means of evaluations across species
for(i in 1:length(epithet_species_list)){

    #select the [i] species
    selected_species = epithet_species_list[i]
    
    #print species name
    print(selected_species)

    #save the name of the [i] species
    means_evaluations[i, "species"] <- selected_species

    #unzip data of threshold
    unzip(paste("threshold/threshold_", selected_species, ".zip", sep=""), exdir="threshold")

    #unzip data of AUC
    unzip(paste("evaluations/evaluations_", selected_species, ".zip", sep=""), exdir="evaluations")

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

        #save the mean and sd of kappa  across partitions of the data
        means_evaluations[i, paste(selected_model, "_kappa_mean", sep="")] <- mean(kappa_values)
        means_evaluations[i, paste(selected_model, "_kappa_sd", sep="")] <- sd(kappa_values)
            
        #save the mean and sd of TSS   across partitions of the data              
        means_evaluations[i, paste(selected_model, "_tss_mean", sep="")] <- mean(tss_values)
        means_evaluations[i, paste(selected_model, "_tss_sd", sep="")] <- sd(tss_values)

        #### AUC ####
        #load AUC data for the [m] model 
        load(paste("evaluations/", selected_species, "_", selected_model, "_evaluation.rda", sep=""))

        #put into an object AUC data of the [m] model
        evaluation_model = eval(parse(text=paste(selected_model, "_evaluation", sep="")))

        #extract AUC data for each data partition
        auc_values = sapply(evaluation_model, "slot", "auc")

        #save the mean and sd of AUC for [m] model across partitions of the data
        means_evaluations[i, paste(selected_model, "_auc_mean", sep="")] <- mean(auc_values)
        means_evaluations[i, paste(selected_model, "_auc_sd", sep="")] <- sd(auc_values)     
    }

    ###########################################
    ###### ALL evaluation metrics for RF ######
    ###########################################
    
    #load the evaluation data for RF
    load(paste("evaluations/", selected_species, "_rf_evaluation.rda", sep=""))

    #save into an object AUC data of RF
    evaluation_rf = eval(parse(text=paste("rf_evaluation", sep="")))

    ### AUC ###
    #AUC is calculated using different thresholds. For each threshold (a value separating presences from absences), you calculate the True positive rate and the False positive rate. As the threshold increase, more true positives you have, but also more false positive. AUC provides an aggregate measure of performance across all possible classification thresholds. One way of interpreting AUC is as the probability that the model ranks a random positive example more highly than a random negative example. Is an integrative measurement ACROSS THRESHOLDS, therefore, if you have very few, it is not a good idea. For example, en random forest, we have a threshold in -0.0001 and another in 1.0001. In the first case, points with a probability higher than -0.0001 will be considered as presences (i.e. all), whilst in the second case only those points with a probability higher than 1 will be considered presences (i.e. no one because we have 0 and 1 values). The rest of thresholds (0.999 and 1) will considered as presences those points with 1 and as absences those points with 0. This is basically the separation in 0 and 1 of a binary random forest.
    evaluation_rf[[1]]@t
    #IN CONCLUSION, you cannot calculate AUC from a classification tree. For further details see "https://developers.google.com/machine-learning/crash-course/classification/roc-and-auc" and "https://stats.stackexchange.com/questions/372236/what-is-the-formula-to-calculate-the-area-under-the-roc-curve-from-a-contingency".

    #### TSS and Kappa####
    #We calculate the TSS metric for Random forest. This metric (true skill statistic) can be used for calculating the optimal threshold when your are binarizing continuos probability of presence. Therefore you calculate the optimal value of probability from which you establish a cell as suitable, i.e. the value that best separate absences from presences. However, this value can be also calculated from a confusion matrix obtained from binary probabilities. Therefore, we can calcualte TSS from evaluation data using random forest!!  

    #TSS = sensitivity + specificity – 1. Sensitivity = nº true presences / (nº true presences + nº false absences); What proportion of total presences are well established, because of this the denominator is the number of true presences more the number of false absences, i.e. presences bad defined as absences. specificity = nº true absences / (nº true absences + nº false presences); What proportion of total absences are well established, because of this the denominator is the number of true absences more the number of false presences, i.e. absences bad defined as presences. See "Assessing the accuracy of species distribution models: prevalence, kappa and the true skill statistic (TSS)" for details about TSS and kapppa calculations

    #The same goes for Kappa, this metric is used to calculate the best threshold, but you only need a confusion matrix to calculate it, so we can obtained from binary probabilities. The formula is very long, you can see it at "Assessing the accuracy of species distribution models: prevalence, kappa and the true skill statistic (TSS)"

    #for each partition of the data we extact confusion matrix stored in evaluate objects of dismo
    tss_values = NULL
    kappa_values = NULL  
    for(j in 1:length(evaluation_rf)){

        #You have to bear in mind that the evaluate function of dismo takes your presences and absences for evaluation, along with your model and the evaluate, but consider that your probabilities are continuous. Because of this, you have 4-5 confusion matrix, each one calculated with a different threshold, a different probability form which you consider suitable. In the case of our RF, which is a clasification tree, probabilities are 0 and 1, because of this, the first confusion matrix has zero false negative and zero true negatives, this matrix is calculated with a threshold = -0.0001, therefore all points are considered as presences, because all have at least a probabiity of zero. The result is that you select well all the true presences, but of course the rest of points are false presences (absences bad defined), as you don't have any absence. The same goes for the last correlation matrix but in the opposite sense. The threshold used for that confusion maitrx was 1.0001, therefore, all point are considered as absences (no point has a probaility of presence higher than 1). In that way you get all the true absences, but the rest of points as false absences (presences bad defined). Intermediate thesholds, i.e. any threshold that differentiates points with 0 and 1 of probability will have the higher rate of success, because our prediciton maps of RF as binary!! Because of this, we will select the first of the threshold that have the maximum true positive rate + true negative rate. All the confuion matrix with the highest TNR+TPR are ok. 

        #sum TPR (true positive rate) and TNR (true negative rate) for each threshold evaluation_rf[[j]]@TPR and evaluation_rf[[j]]@TNR are vectors with the values of TPR and TNR for all dfferent confusion matrix calculated with the different thresholds
        tpr_tnr_different_thresholds = evaluation_rf[[j]]@TPR + evaluation_rf[[j]]@TNR
            #TRP is calculated as true positives / (true positives + false negatives). Therefore is equal to sensitivity.
            #TNR is calculated as true negatives / (true negatives + false positives). Therefore is equal to specificity 

        #which are those with the highest sum? select the first one
        selected_row = which(tpr_tnr_different_thresholds == max(tpr_tnr_different_thresholds))[1]

        #from the confusion matrix selected (selected_row; confusin matrix are presented in rows, each in a row), select the different elements: true positive (true presence), false positive (false presence), false negative (false absence) and true negative (true presence)
        tp = as.vector(evaluation_rf[[j]]@confusion[selected_row, "tp"])
        fp = as.vector(evaluation_rf[[j]]@confusion[selected_row, "fp"])
        fn = as.vector(evaluation_rf[[j]]@confusion[selected_row, "fn"])
        tn = as.vector(evaluation_rf[[j]]@confusion[selected_row, "tn"])

        #calculate the total number of presences and absences (this will be used for kappa calculation)
        n = tp + fp + fn + tn

        #calculate sensitivity
        sensitivity = tp / (tp + fn)#true presences / the sum of true presences more false negtive, i.e. presences bad defined as absences

        #compare with the TPR from the evaluate object (TPR is equal to sensitivity)
        print(paste("check sensitivity:", round(sensitivity - evaluation_rf[[j]]@TPR[selected_row],4)))

        #calculate specificity
        specificity = tn / (tn + fp)#true absences / the sum of true absences more false presences, i.e. absences bad defined as presences
        
        #compare with the TNR from the evaluate object (TNR is equal to specificity)
        print(paste("check specificity:", round(specificity - evaluation_rf[[j]]@TNR[selected_row],4)))

        #calculate TSS
        TSS = sensitivity + specificity - 1

        #other form more direct to calculate tss
        TSS2 = ( (tp*tn) - (fp*fn) ) / ( (tp+fn) * (fp+tn))

        #check both are equal
        check_tss = round(TSS - TSS2, 4)
        print(paste("check TSS:", check_tss))

        #calculate kappa
        kappa_1 = ( ( (tp + tn) / n ) - ( ( (tp+fp) * (tp+fn) + (fn+tn) * (tn+fp) ) / (n^2) ) ) / (1 - ( ( (tp+fp) * (tp+fn) + (fn+tn) * (tn+fp) ) / (n^2) ) )

        #check that the kappa calculate by me is the same than the kappa obtained with the "evaluate" function
        check_kappa = round(kappa_1 - evaluation_rf[[j]]@kappa[selected_row],4)
        print(paste("check kappa:", check_kappa))#we use the kappa value obtained for one of the best thresholds, which is indicated with selected_row (see above)

        #save it in case all is ok
        tss_values = append(tss_values, ifelse(check_tss==0, TSS, NA))
        kappa_values = append(kappa_values, ifelse(check_kappa==0, kappa_1, NA))
    }

    #save the mean and sd of TSS and kappa values across partitions
    means_evaluations[i, paste("rf_tss_mean", sep="")] <- mean(tss_values)
    means_evaluations[i, paste("rf_tss_sd", sep="")] <- sd(tss_values)      
    means_evaluations[i, paste("rf_kappa_mean", sep="")] <- mean(kappa_values)
    means_evaluations[i, paste("rf_kappa_sd", sep="")] <- sd(kappa_values) 

    #remove unzipped files
    file.remove(list.files("threshold", pattern=".rda", full.names = TRUE))
    file.remove(list.files("evaluations", pattern=".rda", full.names = TRUE))
}

#save the results
write.table(means_evaluations, "medians_evaluations/medians_evaluations_v4.csv", sep=",", col.names=TRUE, row.names=FALSE)