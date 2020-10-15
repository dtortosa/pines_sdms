## Code for checking that weight of points is correct. Problem of high and low precision occurrences in random forest and problem of pseudoabsences in GLM and GAM. See above for further details. 

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

#####################################################
####### PROBLEM HIGH VS. LOW PRECISION IN RF ########
#####################################################
#For random forest, the weighted regression has to be made by means of an stratified sampling. We split the presence/pseudoabsence data into three strata: A = High precision occurrences, B = Low precision occurrences and C = pseudoabsences. Then we established the number of points that will be selected from each strata in each tree fitted. We needed that the number of rows selected for strata=A to be double than the number of rows selected with strata = B, because in GLM and GAM the weight of high precision occurrences is the double (1 vs. 0.5). Moreover, rows of A plus rows of B must be equal to the number of rows of C, this would equivalent to sum weight of all occurrences and divided by the number of PAs, which is much higher than the number of occurrences. In that way, the probability of sampling a given PAs in the fit is lower than the probability to select a occurrence.In the case of species with only low precision occurrences, we have to select the same number of low precision occurrences than PAs for sharing the weight of all occurrences in all PAs, giving a lower probability for each PA. 

#Rafi proposed to me an equation system to solve this problem. He considered that A+B+C must be equal than the total number of rows:  
        #Therefore: A+B+C=Nrow, and A=2*B, and C=A+B. 
        #Therefore: 2B + B + 2B + B = Nrow
        #Therefore: 6B = Nrow
        #Therefore: B=Nrow/6, and we round to the floor. 

#The first equation is incompatible with the other two. First, if Nrow = A+B+C, then we are selecting all points and there is no sampling. Second, the number of PAs is much higher compare to occurrences (C very high), thus whether A+B+C=Nrow, Nrow would very higher. Therefore, Nrow/6 would much higher than the actual number of low precision occurrences that we have. The same occur with A, if B (Nrow/6) is very high, also A (2*(Nrow/6)) would be higher than the actual number of high precisio points. In the analyses I fixed this problem indicating if the number of ocurrences (low or high) is lower than needed (Nrow/6 and 2*(Nrow/6)), we select all occurrences we have, then the sum of A+B = C. In this way, the weight of occurrences is shared across the PAs (A+B=C), but the problem is that A always is lower than B, because we select all point with high and low precision, being the latter more common in all species, since they come from Critfield maps.

#the following code calculate the ratio between the total number of high and low precision points, which were used as the sample size of stratum A and B in randpom forest. This is problematic only for species with high precision points and low precision points. If there is only one type of occurrences, all point of this type will be selected and the same number of PAs will be selected, thus A/B=C. This is exactly that GLM and GAM do, if for example we have only low precision points, all occurrences will have 0.5 of weight (weight = 1), we sum these weights and divided by the number of PAs. The weight of PAs is established in basis on low precision points only, it does not matter the number given to low precision points, the key is the relative difference of weights between the data. I have checked this, adding 0.5 to the weight of all low precision points and the sum and divided by the number of PAs and the results are the same. 

#empty data.frame to save results 
results_high_vs_low_occurrences = data.frame(species=character(), nrow_divided_2=numeric(), n_low_precision_points=character(), high_low_precision=character(), n_high_precision_points=character(), high_low_ratio=character())

#for each species
for(i in 1:length(epithet_species_list)){

    #select [i] species
    selected_species = epithet_species_list[i]

    #read data used in modelling for [i] species 
    data = read.table(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ocurrences/", selected_species, "_complete.presences.csv", sep=""), sep=",", header=TRUE)
    
    #calculate the corresponding value of A or B when only there is one type as calculate the wrong code for modelling
    nrow_divided_2 = nrow(data)/2

    #check if there data of low and high precision (the problematic cases, see above)
    high_low_precision = 1 %in% data$precision_weight & 0.5 %in% data$precision_weight

    #number of low precision points
    n_low_precision_points = length(which(data$precision_weight == 0.5))

    #if there are both types of occurrences
    enough_high_precision_points = NULL 
    if(high_low_precision == TRUE){

        #calculate the number of high precision occrurences
        n_high_precision_points = length(which(data$precision_weight == 1))
        
        #calculate the ratio between them
        high_low_ratio = n_high_precision_points/n_low_precision_points

    } else { #if there are not high and low precision occurrences set as N
        n_high_precision_points = NA
        high_low_ratio = NA
    }    

    #save results in the data.frame
    results_high_vs_low_occurrences = rbind.data.frame(results_high_vs_low_occurrences, cbind.data.frame(selected_species, nrow_divided_2, n_low_precision_points, high_low_precision, n_high_precision_points, high_low_ratio))
}

#select those species with high and low precision points
results_high_vs_low_occurrences[which(results_high_vs_low_occurrences$high_low_precision==TRUE),] #canariensis, halepensis, mugo, nigra, pinaster, pinea, sylvestris. In any case there is exactly the double of high precision occurrences than low precision, in all cases there are more low precision occurrences, except in canariensis and pinaster, in which there is 3.85 and 1.01 higher high precision occurrences. This is not correct, we need exactly the double of high precision occurrences than low precision. 

#SOLUTION: The solution will recalculate the sample size of each strata in basis on the number of points of each stratum. If GLM consider high precision points as 1 and low precision points as 0.5, means that the probability of take into account a low precision point is the half than a high precision point. Therefore, if we calculate the number of B points in basis on A we are not taking into account the number of B points and hence the probability of take one of them. Example: There are 50 high precision points, 500 with low precision and 1000 PAs. We need to take all high precision, thus A = 50, so the probability of take into account a high precision in a tree is 100. If we calculated B as A/2, 25 low precision points should be taken. 25 es the half of 50, but this entails a 0.5 of probability for a low precision point being taken? NO, the probability for a given low precision point to be taken is 25 divided the total of low precision points, i.e. 25/500 = 0.05, ten times lower than 0.5. This is wrong, we have to select half of low precision points in each tree, 500/2 = 250. Now, the probability of a low precision point to be taken in 250/500=0.5. Therefore, for each tree the probability to select a high and low precision point is 1 and 0.5 respectively. 

#For each of these 7 species we have to check that the final number of occurrences is enough for the number of varaibles used. The number of variables was established using only occurrences (see "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/variable_selection_inside_clusters/data_number_variables.csv"). In addition, we have to bear in mind the number of high precision points, if we have only 1 high precision point is worthy to divided the low precision by 2? Actually we don't have data of high precision, we should take low precision as the maximum probability an use all points. We should establish a threshold for that. 

#Following the same logic, C = A + B. Therefore, C = 50+250 = 300 points. The probability of a PAs to be taken is 300 divided the total of PAs (1000), i.e. 300/1000=0.3. In GLM the result should be the same: weight of high precision points is 1*50 and weight of the low precision points is 0.5*500, they sum 300, which divided by 1000 is equal to 0.3. It does not matter that the low precision points are labeled as 0.5 for weight, the key is the relative difference to the other points. If there is no points with 1, 0.5 will be the maximum probability and all points will be considered. Indeed, changing weight of low precision points to 1 in species without high precision points does not change anything, the model is the same. 

#This is a problem for the 7 species with high and low precision points. The rest of species with only one type of occurrences (actually the rest has only low precisiion points) are ok. Because in these species never Nrow/2 is lower than the number of low precision points, so all the low precision points are considered to be included in each tree (maximum probability), and the same number of PAs was selected, giving a much lower probability for a PA to be taken. Remember that in the wrong code we made: 
                #Therefore: B+C=Nrow, and C=B. 
                #Therefore: B + B = Nrow
                #Therefore: 2B = Nrow
                #Therefore: B=Nrow/2, and we round to the floor. 

#species for which the number of low precision points is higher than total_points/2
length(which(results_high_vs_low_occurrences$n_low_precision_points > results_high_vs_low_occurrences$nrow_divided_2)) #Any species, therefore in all cases the total number of low precision points were used for each tree, and thus there is no problem with previous analysis of species with only high precision points. Used the same number of occurrences that the used to calculate the number of variables per species. 

#write
write.table(results_high_vs_low_occurrences, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/error_comprobations/weight_data_comprobations/results_high_vs_low_occurrences.txt", sep="\t", col.names=TRUE, row.names=FALSE)

#####################################################
############## PROBLEM WEIGHT PAs ###################
#####################################################
#The weight of each PA was calculated as the sum of weights of all occurrences (low and high precision) divided by the number of PAs. The problem is that the number of PAs used was the initial obtained form ecospat. From these PAs we performed a stratified sampling with more effort on larger strata (combination of environmental variables). Therefero, the final number of PAs was lower than the used for weighting PA, thus we are giving lower weight for PAs that corresponding.

#The following code calculate the weight of each PA according to the number of PAs used in the analyses. In addition is compared to the weight calculated previously. In addition, the weight for PAs finally used in a training data set is calculated to compare with the weight for PAs correctly calculate but using all the data, which was used to discuss with Rafa.  

#empty data.frame
results_weight_PA = data.frame(species=character(), sum_weight_occurrences=character(), n_PA=character(), correct_weight_PA=character(), correct_weight_PA_training=character(), incorrect_weight_PA=character(), correct_incorrect_ratio=character())

#loop
set.seed(56756)
for(i in 1:length(epithet_species_list)){ #for each species

    #select the [i] species
    selected_species = epithet_species_list[i]

    #load data of the [i] species
    data = read.table(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ocurrences/", selected_species, "_complete.presences.csv", sep=""), sep=",", header=TRUE)
    #sum weight of all occruences
    sum_weight_occurrences = sum(data$precision_weight[which(data$precision_weight %in% c(1,0.5))])

    #number of PAs in the final data
    n_PA = nrow(data[which(!data$precision_weight %in% c(1,0.5)),])

    #calculate the correct weight 
    correct_weight_PA = sum_weight_occurrences/n_PA

    #calculate the correct weight of PAs for a subset of training
    presences = data[data$presence==1,] #subset the presences of the final data
    pseudo_absences = data[data$presence==0,] #subset the psuedo-absences
    index_evaluating_presences = sample(1:nrow(presences), round((nrow(presences)*30)/100)) #create index for selecting randmoly the 30% of presences for evaluation
    index_evaluating_absences = sample(1:nrow(pseudo_absences), round((nrow(pseudo_absences)*30)/100)) #create index for selecting randmoly the 70% of presences for training
    train_presen = presences[-index_evaluating_presences, ]
    train_pseudo_absences = pseudo_absences[-index_evaluating_absences,] #the same but wit the rest of ocurrences for obtaining the training data set. 
    training = rbind(train_presen, train_pseudo_absences) #bind presences and PAs for training

    #calculate the correct_weight_PA_training
    correct_weight_PA_training = sum(training[which(training$presence==1),]$precision_weight)/nrow(training[which(training$presence==0),])

    #incorrect weight of PAs
    incorrect_weight_PA = unique(data$precision_weight[which(!data$precision_weight %in% c(1,0.5))])

    #ratio between correct and incorrect weights for PAs
    correct_incorrect_ratio = correct_weight_PA/incorrect_weight_PA

    #save
    results_weight_PA = rbind.data.frame(results_weight_PA, cbind.data.frame(selected_species, sum_weight_occurrences, n_PA, correct_weight_PA, correct_weight_PA_training, incorrect_weight_PA, correct_incorrect_ratio))
}

#max, min, median
max(results_weight_PA$correct_incorrect_ratio)
min(results_weight_PA$correct_incorrect_ratio)
median(results_weight_PA$correct_incorrect_ratio)

#species with ratio higher than 3
results_weight_PA$selected_species[which(results_weight_PA$correct_incorrect_ratio > 3)]

#correlation between weight calculated for the full dataset and the weight calculate for a training subset
par(mfcol=c(2,1))
plot(results_weight_PA$correct_weight_PA, results_weight_PA$correct_weight_PA_training)
plot(density(results_weight_PA$correct_weight_PA-results_weight_PA$correct_weight_PA_training)) #almost the same, so the results discussed with rafa are valid. IMPORTANT: The exact value of the PA weight can change because the random split of the dataset produces different number of high and low precision points between each run, and hence different sum of weight that will be dividided by the total number of PAs. 

#write
write.table(results_weight_PA, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/error_comprobations/weight_data_comprobations/results_weight_PA.txt", sep="\t", col.names=TRUE, row.names=FALSE)

########## compare ensambles of current suitaiblity with correct and incorrect wieght for PAs ############

#load data about problem with weight of PAs
results_weight_PA = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/error_comprobations/weight_data_comprobations/results_weight_PA.txt", sep="\t", header=TRUE)

#species for which correct weight is 3 times or more higher than the incorrect weight
weight_problem_species = results_weight_PA$selected_species[which(results_weight_PA$correct_incorrect_ratio >= 3)]
weight_problem_species = as.vector(weight_problem_species)


if(FALSE){
#comparison for all species
for(i in 1:length(weight_problem_species)){ #for each problematic species
    
    #select the [i] species
    selected_species = weight_problem_species[i]

    #load the ensamble of current suitability with correct PA weight
    correct_weight = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/error_comprobations/error_comprobations_current/ensamble_predictions_bin_correct_weight/ensamble_predictions_bin_", selected_species, ".tif", sep=""))

    #load the ensamble of current suitability with incorrect PA weight
    incorrect_weight = raster(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/error_comprobations/error_comprobations_current/ensamble_predictions_bin_incorrect_weight/ensamble_predictions_bin_", selected_species, ".tif", sep=""))

    #plot both of them
    pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/error_comprobations/weigth_PAs/comparison/comparison_PA_weight_", selected_species, ".pdf", sep=""), width=12, height = 6)
    par(mfcol=c(1,2))
    par(oma=c(0,0,2.7,2))
    plot(correct_weight, main="Correct PA weight")
    mtext(text=bquote(italic('Pinus') ~italic(.(selected_species))), side=3, line=3, adj=1.65, outer=FALSE, cex=1.4, font = 2)
    plot(incorrect_weight, , main="Incorrect PA weight")
    dev.off()
}

#bind all pdf generated into one single pdf file
system("rm /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/error_comprobations/error_comprobations_current/full_comparison/comparison_PA_weight_full.pdf;

    cd /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/error_comprobations/error_comprobations_current/comparison; 

    pdftk *.pdf cat output /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/error_comprobations/error_comprobations_current/full_comparison/comparison_PA_weight_full.pdf") #first delete the prevous full file in full directory; next change the drectory to the that with separate plots; finally run pdftk to bind all pdf in that directory and save the result in full directory
}

#none of these 32 species have high and low precision points
results_weight_PA$selected_species[which(results_weight_PA$correct_incorrect_ratio > 3)] %in% results_high_vs_low_occurrences[which(results_high_vs_low_occurrences$high_low_precision==TRUE),]$selected_species #therefore differences are only related to the weight of the PAs, and the random partition of the data (in a lesser extent).

### From these 32 species the most problematic species according to visual analsys
species_with_highest_differences = c("bhutanica", "bungeana", "culminicola", "dalatensis", "densata", "fenzeliana", "fragilisima", "kesiya", "krempfii", "latteri", "massoniana", "merkusii", "morrisonicola", "muricata", "parviflora",  "patula", "pringlei", "squamata", "tabuliformis", "taiwanensis", "yunnanensis")

#median of correct_incorrect ratio for all species
median(results_weight_PA$correct_incorrect_ratio)

#median of correct_incorrect ratio for species with that ratio higher than 3
median(results_weight_PA[which(results_weight_PA$correct_incorrect_ratio > 3),]$correct_incorrect_ratio)

#median of correct_incorrect ratio for species_with_highest_differences
median(results_weight_PA[which(results_weight_PA$selected_species %in% species_with_highest_differences),]$correct_incorrect_ratio)


########## compare ensambles of future suitaiblity with correct and incorrect wieght for PAs ############

#load data about problem with weight of PAs
results_weight_PA = read.table("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/error_comprobations/weight_data_comprobations/results_weight_PA.txt", sep="\t", header=TRUE)

#species for which correct weight is 3 times or more higher than the incorrect weight
weight_problem_species = results_weight_PA$selected_species[which(results_weight_PA$correct_incorrect_ratio >= 3)]
weight_problem_species = as.vector(weight_problem_species)


if(FALSE){
#comparison for all species
for(i in 1:length(weight_problem_species)){ #for each problematic species
    
    #select the [i] species
    selected_species = weight_problem_species[i]

    #load the ensamble of future suitability with correct PA weight
    correct_weight = raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_final/ensamble_projections_bin/ensamble_projections_bin_", selected_species, ".tif", sep=""))

    #load the ensamble of future suitability with incorrect PA weight
    incorrect_weight = raster(paste("/Users/diegosalazar/phd_big_documents/pines_niche/ensambles_ancient/ensamble_projections_bin_ancient/ensamble_projections_bin_", selected_species, ".tif", sep=""))

    #plot both of them
    pdf(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/error_comprobations/error_comprobations_future/comparison/comparison_PA_weight_", selected_species, ".pdf", sep=""), width=12, height = 6)
    par(mfcol=c(1,2))
    par(oma=c(0,0,2.7,2))
    plot(correct_weight, main="Correct PA weight")
    mtext(text=bquote(italic('Pinus') ~italic(.(selected_species))), side=3, line=3, adj=1.65, outer=FALSE, cex=1.4, font = 2)
    plot(incorrect_weight, , main="Incorrect PA weight")
    dev.off()
}

#bind all pdf generated into one single pdf file
system("rm /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/error_comprobations/error_comprobations_future/full_comparison/comparison_PA_weight_full.pdf;

    cd /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/error_comprobations/error_comprobations_future/comparison; 

    pdftk *.pdf cat output /Volumes/GoogleDrive/My\\ Drive/science/phd/nicho_pinus/results/error_comprobations/error_comprobations_future/full_comparison/comparison_PA_weight_full.pdf") #first delete the prevous full file in full directory; next change the drectory to the that with separate plots; finally run pdftk to bind all pdf in that directory and save the result in full directory
}