#code for comparison the two criterias of variable selection: The difference between these two criterias is the number of variables selectec, in the new criteria we have a fixed number of variables (5) and a fixed ratio soil/climate (2/3). In the ancient, can be selected more than 5 variables, and the ratio soil/climate  is variable (bias potential response to climate of species). The selection with both types of criteria has been made in variable_selection_inside_clusters.R

###definimos el directorio de trabajo
setwd("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus")

#Librerias
require(raster)
require(foreach) #for repeat a process several times 
require(doParallel) #for parallel

###list ocurrences
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

##Load variables
list_variables = list.files("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals", full.names=TRUE, pattern=".asc")
variables = stack(list_variables)

##Load the group variable and the ranks 
group_variables = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/species_clustering/_tables/complete_2_10_g_5.csv", header=TRUE)  
length(unique(group_variables$groups)) #we have 5 groups


####################################################
########### Ancient criteria #########################
####################################################

####Create data frames with presence and environmental variables for each species
calculate_d2 = function(species){

    #requireed functions
    Dsquared = function(model, adjust = TRUE) {
      # version 1.1 (13 Aug 2013)
      # calculates the explained deviance of a GLM
      # model: a model object of class "glm"
      # adjust: logical, whether or not to use the adjusted deviance taking into acount the nr of observations and parameters (Weisberg 1980; Guisan & Zimmermann 2000)
      d2 <- (model$null.deviance - model$deviance) / model$null.deviance
      if (adjust) {
        n <- length(model$fitted.values)
        p <- length(model$coefficients)
        d2 <- 1 - ((n - 1) / (n - p)) * (1 - d2)
      }
      return(d2)
    } #D2. If you have problem with this you can use Dsquared function of modEvA (http://modeva.r-forge.r-project.org/). 

    #select the group of the corresponding species
    variables_cluster = group_variables[which(group_variables$species == species),]$groups 

    #select the number cluster of this species 
    rasters_list = list.files("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals", pattern=".asc", full.names=TRUE) #list the corresponding group of variables
    variables_stack = stack(rasters_list) #stack them 

    #load the number of ocurrences for each species 
    number_ocurrences = read.csv("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/ocurrences/ocurrences_per_species.csv", header=TRUE)

    #extract the number ocurrences of the corresponding species
    n_ocurrences = number_ocurrences[number_ocurrences$species == species,]$number_ocurrences
   
    #load variables of low number ocurrences species according to ancient criteria
    load("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/final_variables_low_number_ocurrence_species_ancient.rda")
    
    #load variables of low number ocurrences species according to new criteria
    load("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/datos/finals/final_variables_low_number_ocurrence_species.rda")

    #load selected variables according to the ancient criteria (more than 5 variables)
    load("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/variables/final_variables_MAL/list_selected_variables_ancient_criteria.rda") 
    ultimate_variables_ancient = ultimate_variables

    #load selected variables according to the new criteria (5 variables: 3 climatic and 2 of soil)
    load("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/final_variables/list_selected_variables.rda") #this change the ultimate variable object
    ultimate_variables_new = ultimate_variables


    #if else for determining if the species has low or higher number of ocurrences
    if (species %in% names(final_variables_low_number_ocurrence_species_new)){ #if the species is included in the list of low number of ocurrences 

        selected_variables_new = final_variables_low_number_ocurrence_species_new[[species]] #the variables will be selected from the list of low number

    } else { #if not

        selected_variables_new = ultimate_variables_new[[variables_cluster]] #the variables will be selected from the list of high number of ocurrences

    }

    #the same for variabls selected with ancient criteria 
    if (species %in% names(final_variables_low_number_ocurrence_species_ancient)){

        selected_variables_ancient = final_variables_low_number_ocurrence_species_ancient[[species]]
    } else { 

        selected_variables_ancient = ultimate_variables_ancient[[variables_cluster]] #select the variables of this cluster

    }

    #select variables from the stack
    variables_stack_ancient = variables_stack[[selected_variables_ancient]] 
    variables_stack_new = variables_stack[[selected_variables_new]] 

    ##Extract values of variables 
    presences = read.csv(paste("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/results/pseudo_absences", paste(species, "complete.presences.csv", sep="_"), sep="/")) #read the data final data with presences and pseudoabsences
    variables_ancient = extract(variables_stack_ancient, presences[, c("longitude", "latitude")]) #extract the value of the variable in these points
    variables_new = extract(variables_stack_new, presences[, c("longitude", "latitude")]) #extract the value of the variable in these points
    data_ancient = cbind(presences, variables_ancient) #bind the presence data and the variable data in one data frame 
    data_new = cbind(presences, variables_new) #bind the presence data and the variable data in one data frame 

    #create the formula 
    formula.regresion_poly_ancient = as.formula(paste("presence ~ poly(", paste(names(data_ancient[,-c(1:4)]), collapse=", 2) + poly("), ", 2)", collapse=""))
    formula.regresion_poly_new = as.formula(paste("presence ~ poly(", paste(names(data_new[,-c(1:4)]), collapse=", 2) + poly("), ", 2)", collapse=""))

    #run the glm
    glm_pseudo_absence_ancient<-glm(formula.regresion_poly_ancient, family=binomial(link=logit), weights=precision_weight, data=data_ancient)     
    glm_pseudo_absence_new<-glm(formula.regresion_poly_new, family=binomial(link=logit), weights=precision_weight, data=data_new) 
    
    #calcualte D2
    d2_ancient = Dsquared(glm_pseudo_absence_ancient)
    d2_new = Dsquared(glm_pseudo_absence_new)

    #bind both values
    d2 = cbind(d2_ancient, d2_new)

    #save them
    d2_data_frame = data.frame()
    d2_data_frame = rbind(d2_data_frame, d2)
    d2_data_frame = cbind(species, d2_data_frame)
    d2_data_frame
}

#create a vector with species of cluster number 4 
species = as.vector(group_variables[group_variables$groups==4,]$species) #we include as.vector with the purpose of avoid that species object is created as a column of group_varialbes (factor), instead of as a vector. 

# set up cluster
clust <- makeCluster(1) 
registerDoParallel(clust)

# run for all species
d2.final = foreach(species = species, .packages="raster", .combine="rbind") %dopar% { 
    calculate_d2(species = species)
} 

#stop the cluster 
stopCluster(clust)

#check if all is correct
str(d2.final)

#write the rank
write.csv(d2.final, "/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/variable_selection_inside_clusters/comparison_selection_criterias/comparison_d2_two_types_criteria.csv")

##see all data
d2.final

#see cases in which d2_ancient is higher than 2_new
d2.final[d2.final$d2_new<d2.final$d2_ancient,] #in general theres is few differences

#see cases in which d2_ancient is higher than 2_new
d2.final[d2.final$d2_new>d2.final$d2_ancient,] #there is a little bit increase of the d2 en some cases. The three extrange cases are due to problems in models, which don't converge. 

#plot
png("/Volumes/GoogleDrive/My Drive/science/phd/nicho_pinus/code/variables/variable_selection/variable_selection_inside_clusters/comparison_selection_criterias/comparison_d2_two_types_criteria.png")
plot(d2.final$d2_ancient~d2.final$d2_new)
dev.off() #In general the results are very similar, and the cases with more differences are cases with mode d2 for new criteria variables, because the models with more variables don't converge. The new criteria of variable selection works perfect. 

