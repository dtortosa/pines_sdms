###################################
#ESTABLECE EL DIRECTORIO DE TRABAJO
###################################
#DIRECTORIO DE TRABAJO
setwd("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus")

####################################################
########### Loop preparing data #########################
####################################################
##RUN only ONE time. Already runned. 

#function for preparing data
library(raster)
library(dismo)
prepar_data = function(species) {

    #load the group of species according to cluster
    group_species = read.csv("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/code/variables_presences/species_clustering/_tables/complete_2_10_g_5.csv", header=TRUE)  

    #select the group of the corresponding species
    variables_cluster = group_species[which(group_species$species == species),]$groups 

    #select the number cluster of this species 
    rasters_list = list.files("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/finals", pattern=".asc", full.names=TRUE) #list the corresponding group of variables
    variables_stack = stack(rasters_list) #stack them 

    #load names of selected variables
    load("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/final_variables/list_selected_variables.rda") 

    #take selected variables from the list of variables
    selected_variables = ultimate_variables[[variables_cluster]]

    #calculate the number of ocurrences for each species 
    number_ocurrences = read.csv("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/ocurrences/ocurrences_per_species.csv", header=TRUE)
    n_ocurrence = number_ocurrences[number_ocurrences$species==species,]$number_ocurrences

    if (n_ocurrence<length(selected_variables)*10){ #if there is not 10 ocurrences for each variable 
        
        #load the names of selected variables for los number ocurrences species
        load("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/finals/final_variables_low_number_ocurrence_species.rda")

        #select the selecte variables for low number ocurrences species
        variables_stack = variables_stack[[final_variables_low_number_ocurrence_species_new[[species]]]]    
    
    } else {

        variables_stack = variables_stack[[selected_variables]]
    
    }     
 
    ##Extract values of variables 
    presences = read.csv(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/pseudo_absences", paste(species, "complete.presences.csv", sep="_"), sep="/")) #read the data final data with presences and pseudoabsences
    variables = extract(variables_stack, presences[, c("longitude", "latitude")])

    #extract the value of the variable in these points
    data = cbind(presences, variables) #bind the presence data and the variable data in one data frame 

    #change the names of variables if the number of selected variables is 1 (it is to say, only 5 columns in data)
    if (ncol(data)<6){
        colnames(data)[5] = names(variables_stack)
    }

    #write the final.data file
    write.csv(data, paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/models", paste(species, "csv", sep="."), sep="/"), row.names=FALSE)    
}

####################################################
########### Loop of fitting ########################
####################################################

#In each step we will create a dataset for training and evaluation and make both process, then we will repeat with other partition of training and evaluation. 
#We have to create a loop with for each species and before create empty lists (glm_species...) for include all the models for each species. In these lists we will sabe XX_resample lists


#function for fitting models 
library(raster)
library(dismo)
library(gam) #OTRA LIBRERIA PARA GAM
library(randomForest) #RANDOM FOREST
fit_eval_models = function(species){ #for the corresponding species 
    
    #begin the species
    print(paste("begin",species))

    #load cluster number for each species
    group_species = read.csv("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/code/variables_presences/species_clustering/_tables/complete_2_10_g_5.csv", header=TRUE)[,c("species", "groups")] 

    #select the group of the corresponding species
    variables_cluster = group_species[which(group_species$species == species),]$groups 

    #Load data with presences and values of environmental variables
    data = read.csv(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/models", paste(species, "csv", sep="."), sep="/"), header=TRUE)

    #select the number cluster of this species 
    rasters_list = list.files("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/finals", pattern=".asc", full.names=TRUE) #list the corresponding group of variables
    variables_stack = stack(rasters_list) #stack them 

    #calculate the number of ocurrences for each species 
    number_ocurrences = read.csv("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/ocurrences/ocurrences_per_species.csv", header=TRUE)
    n_ocurrence = number_ocurrences[number_ocurrences$species==species,]$number_ocurrences

    #load the final variables selected for each cluster 
    load("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/final_variables/list_selected_variables.rda") #the object is ultimate_variables

    #select the selected variable for the corresponding cluster
    selected_variables = ultimate_variables[[variables_cluster]]
    
    #if else for deal with species with low and high number of ocurrences
    if (n_ocurrence<length(selected_variables)*10){ #if there is not 10 ocurrences for each variable of the cluster 
        
        #load the number of ocurrences
        load("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/finals/final_variables_low_number_ocurrence_species.rda") #the object is final_variables_low_number_ocurrence_species.rda

        #select the variables for this species from the list of reduced number in basis on raking
        reduced_variables = final_variables_low_number_ocurrence_species_new[[species]] #select element "species" of the list

        #select these variable from the stack of current variables 
        variables_stack = variables_stack[[reduced_variables]]    
    
    } else { #if not, and thus there are a lot of points, we can use all the variables for the corresponding cluster

        variables_stack = variables_stack[[selected_variables]]
    
    }         
 
    #load raster of PA buffer
    raster_PA_buffer = raster(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/pseudo_absences", paste(species, "PA_buffer.asc", sep="_"), sep="/"))

    #crop variables using distribution buffer
    variables_stack = crop(variables_stack, raster_PA_buffer)

    #lists for saving models 
    glm_resample = list() 
    gam_resample = list() 
    rf_resample = list()

    #lists for saving predictions 
    glm_predict = list() 
    gam_predict = list() 
    rf_predict = list()

    #lists for saving predicted values on evaluation points 
    glm_evaluation_predict = list() 
    gam_evaluation_predict = list() 
    rf_evaluation_predict = list()

    #lists for saving evaluation of models 
    glm_evaluation = list()
    gam_evaluation = list()
    rf_evaluation = list()

    #For 12 times 
    for (k in 1:12){

        ######################################################
        ######separate data for training and evaluation#######
        ######################################################
        presences = data[data$presence==1,] #subset the presences of the final data
        pseudo_absences = data[data$presence==0,] #subset the psuedo-absences
        index_evaluating_presences = sample(1:nrow(presences), round((nrow(presences)*30)/100)) #create index for selecting randmoly the 30% of presences for evaluation
        index_evaluating_absences = sample(1:nrow(pseudo_absences), round((nrow(pseudo_absences)*30)/100)) #create index for selecting randmoly the 70% of presences for training
        evaluate_presen = presences[index_evaluating_presences, ] #subset presences for evaluation
        evaluate_pseudo_absences = pseudo_absences[index_evaluating_absences,] #subset PAs for evaluation
        train_presen = presences[-index_evaluating_presences, ]
        train_pseudo_absences = pseudo_absences[-index_evaluating_absences,] #the same but wit the rest of ocurrences for obtaining the training data set. 
        training = rbind(train_presen, train_pseudo_absences) #bind presences and PAs for training
        evaluation = rbind(evaluate_presen, evaluate_pseudo_absences) #bind presences and PAs for evaluation

        ###################################
        ############Fit models ############
        ###################################

        ################################
        #GENERALIZED LINEAR MODELS (GLM)
        ################################
        formula.regresion.poly = as.formula(paste("presence ~ poly(", paste(names(training)[-c(1:4)], collapse=", 2) + poly("), ", 2)", collapse=""))
        glm_pseudo_absence<-glm(formula.regresion.poly, family=binomial(link=logit), weights=precision_weight, data=training) #usa presencia como variable respuesta y como predictores solo dos variable. La familia e sbinomial porque nuestros datos de presencia tienen 1 y 0, y la curva que voy a usar es logística (polínomo de primer grado). Scamos los datos de la tabla presencia.pseudoausencia.entrenamiento, que tiene todas las presencias menos las de evaluacion, y todas las pseudoausencias. "quasibinomial" family bacause of the presences are weighted by precision_weight, and thus there is numbers with decimals, not only 0 and 1. We will use polynomial the level 2 with the purpose that the response curve can differ from the lineality, but no more than 2, because we want a simple model. 

        #make stepwise
        glm_pseudo_absence = step(glm_pseudo_absence, direction="both", trace=0)

        #####################################
        #MODELOS ADITIVOS GENERALIZADOS (GAM)
        #####################################
        formula.gam<-as.formula(paste("presence ~ s(", paste(names(training)[-c(1:4)], collapse=",4) + s("), ",4)", collapse="")) #cada varaible está envuelta en una funcion que suaviza (smooth), gam traaja sobre valore suavizados de las variables. Por eso es s(). 
        if (is.element("package:mgcv", search())) detach("package:mgcv", force=TRUE) #make sure the mgcv package is not loaded to avoid conflicts between packages.  
        gam_pseudo_absence<-gam(formula.gam, family=binomial(link=logit), data=training, weights=precision_weight)

        #make stepwise
        #genereate a scope for step.gam
        #Given a data.frame as an argument, generate a scope list for use in step.gam, each element of which gives the candidates for that term.
        gam_scope = gam.scope(frame=training[,-c(1:3)], response=1, smoother="s", arg="df=4") #frame is the data: response variable and explicativa variables; response indicate the column of the response; smoother is the type of smooth that we use; arg: additional argument, like for example the freedom degree of smooth. 
        gam_pseudo_absence = step.gam(gam_pseudo_absence, scope=gam_scope, direction="both", trace=0, data=training) #we include the list of terms in scope. 


        ##############
        #RANDOM FOREST
        ##############

        #create a strata variable with diferent value for each precision_weight
        if(nrow(training[training$precision_weight==1,])>0){ #if there are points with high precision
            if(nrow(training[training$precision_weight==0.5,])>0){ #if there is low precision ocurrences

                training$strata = factor(NA, levels=c("A", "B", "C")) #create an empty factor 
                training[training$precision_weight==1,]$strata <- "A" #if precision_weight==1 strata is equal to "A"
                training[training$precision_weight==0.5,]$strata <- "B" #if precision_weight==0.5 strata is equal to "B"
                training[!(training$precision_weight==1 | training$precision_weight==0.5),]$strata <- "C" #if precision_weight is different of 1 and 0.5, strata is equal to "C"
            } else { #if there is not low precision ocurrences
                training$strata = factor(NA, levels=c("A", "C")) #create an empty factor 
                training[training$precision_weight==1,]$strata <- "A" #if precision_weight==0.5 strata is equal to "B"
                training[!(training$precision_weight==1 | training$precision_weight==0.5),]$strata <- "C" #if precision_weight is different of 1 and 0.5, strata is equal to "C"
            }
        } else { #if not 
            training$strata = factor(NA, levels=c("B", "C")) #create an empty factor 
            training[training$precision_weight==0.5,]$strata <- "B" #if precision_weight==0.5 strata is equal to "B"
            training[!(training$precision_weight==1 | training$precision_weight==0.5),]$strata <- "C" #if precision_weight is different of 1 and 0.5, strata is equal to "C"
        }

        #comprobation
        #table(training[training$strata=="A",]$precision_weight)
        #table(training[training$strata=="B",]$precision_weight)
        #table(training[training$strata=="C",]$precision_weight)

        #save the new factor as a vector 
        strata_vector = training$strata 

        #delete the factor from the data.frame         
        training = training[,-which(names(training)=="strata")] 

        #create the N sample size
        #this vector have the same length than unique(strata_vector), and indicate the number of rows that will be taken for each strata (category=high precision, low precision and PAs).
        #we need that the number of rows selected with strata=A to be higher than the number of rows selected with strata = B. Moreover, rows of A plus rows of B must be equal to the number of rows of C. In addition, A+B+C must be equal than the total number of rows. 
        #Therefore: A+B+C=Nrow, and A=2*B, and C=A+B. 
        #Therefore: 2B + B + 2B + B = Nrow
        #Therefore: 6B = Nrow
        #Therefore: B=Nrow/6, and we round to the floor. 
 
        if (nrow(training[training$precision_weight==1,])>0){ #if there are points with high precision (A higher than 0):
            if(nrow(training[training$precision_weight==0.5,])>0){ #if there is low precision ocurrences
                b = ifelse(floor(nrow(training)/6) > table(strata_vector)["B"], floor(table(strata_vector)["B"]), floor((nrow(training))/6)) #we divided the number of rows by 6 to calculate B, but we apply a ifelse. If Nrow/6 surpasses the total number of points in the stata B, the number will be this limits. If not the number will be Nrow/6
                a = ifelse(2 * b > table(strata_vector)["A"], floor(table(strata_vector)["A"]), 2 * b) #A=2*B, and we apply the ifelse, like in the last case, but using the number of points in the strata A
                c = a + b #calcualte C

                sampsize_vector = c("A"=as.vector(a), "B"=as.vector(b), "C"=as.vector(c)) #bind all numbers in a unique vector of length equal to 3.  
            } else { #if there is not low precision ocurrences
                    #Therefore: A+C=Nrow, and C=A. 
                    #Therefore: A + A = Nrow
                    #Therefore: 2A = Nrow
                    #Therefore: A=Nrow/2, and we round to the floor. 
                a = ifelse(floor(nrow(training)/2) > table(strata_vector)["A"], floor(table(strata_vector)["A"]), floor((nrow(training))/2)) #we divided the number of rows by 6 to calculate B, but we apply a ifelse. If Nrow/6 surpasses the total number of points in the stata B, the number will be this limits. If not the number will be Nrow/6
                c = a #calcualte C

                sampsize_vector = c("A"=as.vector(a), "C"=as.vector(c)) #bind all numbers in a unique vector of length equal to 2 (there is no ocurrences with high precision, A=0)  

            }
        } else { #if not, it is to say, there is not high precision ocurrences (A=0)
                #Therefore: B+C=Nrow, and C=B. 
                #Therefore: B + B = Nrow
                #Therefore: 2B = Nrow
                #Therefore: B=Nrow/2, and we round to the floor. 
            b = ifelse(floor(nrow(training)/2) > table(strata_vector)["B"], floor(table(strata_vector)["B"]), floor((nrow(training))/2)) #we divided the number of rows by 6 to calculate B, but we apply a ifelse. If Nrow/6 surpasses the total number of points in the stata B, the number will be this limits. If not the number will be Nrow/6
            c = b #calcualte C

            sampsize_vector = c("B"=as.vector(b), "C"=as.vector(c)) #bind all numbers in a unique vector of length equal to 2 (there is no ocurrences with high precision, A=0)  
        }

        #write regresion formula
        formula.regresion<-as.formula(paste("as.factor(presence) ~ ", paste(names(training)[-c(1:4)], collapse="+"), collapse=""))

        #run the model
        rf_pseudo_absence<-randomForest(formula.regresion, data=training, importance=TRUE, ntree=500, strata=strata_vector, sampsize=sampsize_vector) #It is a classification tree, we select this becasue our variable is a factor of two levels. See "http://www.simafore.com/blog/bid/62482/2-main-differences-between-classification-and-regression-trees" for further information. 

        ###################################
        ############Save models ###########
        ###################################
        glm_resample[[k]] = glm_pseudo_absence
        gam_resample[[k]] = gam_pseudo_absence
        rf_resample[[k]] = rf_pseudo_absence

        ###################################
        ############Predictions############
        ################################### 

        glm_predict[[k]] = predict(variables_stack, glm_resample[[k]], type="response")
        gam_predict[[k]] = predict(variables_stack, gam_resample[[k]], type="response")
        rf_predict[[k]] = predict(variables_stack, rf_resample[[k]], type="response") #the prediction of rf is binary because we are using classification forest. therefore, we don't have to binarize the predictions of this model. 


        ###################################
        ############Evaluation ############
        ################################### 

        #extract presence probability on evaluation points
        predict_glm = extract(glm_predict[[k]], evaluation[,c("longitude", "latitude")])
        predict_gam = extract(gam_predict[[k]], evaluation[,c("longitude", "latitude")])
        predict_rf = extract(rf_predict[[k]], evaluation[,c("longitude", "latitude")])


        #bind predictions and cordinates of evaluation points
        glm_evaluation_predict[[k]] = cbind(evaluation[,c("longitude", "latitude", "presence")], predict_glm)
        gam_evaluation_predict[[k]] = cbind(evaluation[,c("longitude", "latitude", "presence")], predict_gam)
        rf_evaluation_predict[[k]] = cbind(evaluation[,c("longitude", "latitude", "presence")], predict_rf)     

        #make the evaluation
        glm_evaluation[[k]] = evaluate(p=glm_evaluation_predict[[k]][glm_evaluation_predict[[k]]$presence==1, "predict_glm"], a=glm_evaluation_predict[[k]][glm_evaluation_predict[[k]]$presence==0, "predict_glm"])
        gam_evaluation[[k]] = evaluate(p=gam_evaluation_predict[[k]][gam_evaluation_predict[[k]]$presence==1, "predict_gam"], a=gam_evaluation_predict[[k]][gam_evaluation_predict[[k]]$presence==0, "predict_gam"])
        rf_evaluation[[k]] = evaluate(p=rf_evaluation_predict[[k]][rf_evaluation_predict[[k]]$presence==1, "predict_rf"], a=rf_evaluation_predict[[k]][rf_evaluation_predict[[k]]$presence==0, "predict_rf"])
    }

    #stack all continuous predictions 
    continuous_predictions_glm= stack(glm_predict)
    continuous_predictions_gam= stack(gam_predict)

    #strack predictions of random forest (binary)
    binary_predictions_rf = stack(rf_predict)   

    #save models 
    save(glm_resample, file=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/models", paste(species, "glm_model.rda", sep="_"), sep="/"))
    save(gam_resample, file=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/models", paste(species, "gam_model.rda", sep="_"), sep="/"))
    save(rf_resample, file=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/models", paste(species, "rf_model.rda", sep="_"), sep="/"))  

    #save the continuous predictions
    writeRaster(continuous_predictions_glm, filename=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/predictions/continuous_predictions_glm", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
    writeRaster(continuous_predictions_gam, filename=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/predictions/continuous_predictions_gam", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)

    #save the predictions of random forest (binary)
    writeRaster(binary_predictions_rf, filename=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/predictions/binary_predictions_rf", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
    file.remove(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/predictions/binary_predictions_rf", paste(species, "tif.aux.xml", sep="."), sep="_"))

    #save data for evaluation
    save(glm_evaluation_predict, file=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/evaluations", paste(species, "glm_evaluation_data.rda", sep="_"), sep="/"))
    save(gam_evaluation_predict, file=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/evaluations", paste(species, "gam_evaluation_data.rda", sep="_"), sep="/"))
    save(rf_evaluation_predict, file=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/evaluations", paste(species, "rf_evaluation_data.rda", sep="_"), sep="/")) 

    #save evaluations
    save(glm_evaluation, file=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/evaluations", paste(species, "glm_evaluation.rda", sep="_"), sep="/"))
    save(gam_evaluation, file=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/evaluations", paste(species, "gam_evaluation.rda", sep="_"), sep="/"))
    save(rf_evaluation, file=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/evaluations", paste(species, "rf_evaluation.rda", sep="_"), sep="/"))

    #end the species
    print(paste("end",species))  
}

#function for binarize the current projections  
#I have separated the binarization from the fitting of modelling because I have problems in gam when I load ecospat package, which is used in caculates of threshold.      
library(raster)
library(dismo)
library(ecospat) #for threshold
binarize_curent_projections = function(species){

    #begin the species
    print(paste("begin",species))

    #load models 
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/models", paste(species, "glm_model.rda", sep="_"), sep="/"))
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/models", paste(species, "gam_model.rda", sep="_"), sep="/"))
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/models", paste(species, "rf_model.rda", sep="_"), sep="/"))  

    #load evaluations
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/evaluations", paste(species, "glm_evaluation.rda", sep="_"), sep="/"))
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/evaluations", paste(species, "gam_evaluation.rda", sep="_"), sep="/"))
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/evaluations", paste(species, "rf_evaluation.rda", sep="_"), sep="/")) 

    #load evaluation data
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/evaluations", paste(species, "glm_evaluation_data.rda", sep="_"), sep="/"))
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/evaluations", paste(species, "gam_evaluation_data.rda", sep="_"), sep="/"))
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/evaluations", paste(species, "rf_evaluation_data.rda", sep="_"), sep="/")) 

    #load continuous predictions of glm and gam
    glm_predict = stack(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/predictions/continuous_predictions_glm", paste(species, "tif", sep="."), sep="_"))
    gam_predict = stack(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/predictions/continuous_predictions_gam", paste(species, "tif", sep="."), sep="_"))

    #load binary predictions of Random forest (already finished)
    rf_predict = stack(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/predictions/binary_predictions_rf", paste(species, "tif", sep="."), sep="_"))

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

        #convert each projection in binary using the threshold calculated with fit_eval_models
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
    save(glm_threshold, file=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/threshold", paste(species, "glm_threshold.rda", sep="_"), sep="/"))
    save(gam_threshold, file=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/threshold", paste(species, "gam_threshold.rda", sep="_"), sep="/"))

    #save the binary predictions
    writeRaster(binary_predictions, filename=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/predictions/binary_predictions", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)

    #save ensamble
    writeRaster(ensamble_predictions_bin, filename=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/predictions", paste("ensamble_predictions_bin", paste(species, "tif", sep="."), sep="_"), sep="/"), overwrite=TRUE)

    #end the species
    print(paste("end",species))  
}


#function for make future projections
library(raster)
library(dismo)
library(randomForest) #for predict rf models
library(gam) #for predict gam model
projections = function(species){

    #begin the species
    print(paste("begin",species))

    #load cluster number for each species
    group_species = read.csv("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/code/variables_presences/species_clustering/_tables/complete_2_10_g_5.csv", header=TRUE)[,c("species", "groups")] 

    #select the group of the corresponding species
    variables_cluster = group_species[which(group_species$species == species),]$groups 

    #Load data with presences and values of environmental variables
    data = read.csv(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/models", paste(species, "csv", sep="."), sep="/"), header=TRUE)

    #select the number cluster of this species 
    rasters_list = list.files("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/finals", pattern=".asc", full.names=TRUE) #list the corresponding group of variables
    variables_stack = stack(rasters_list) #stack them 

    #calculate the number of ocurrences for each species 
    number_ocurrences = read.csv("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/ocurrences/ocurrences_per_species.csv", header=TRUE)
    n_ocurrence = number_ocurrences[number_ocurrences$species==species,]$number_ocurrences

    #load the final variables selected for each cluster 
    load("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/final_variables/list_selected_variables.rda") #the object is ultimate_variables

    #select the selected variable for the corresponding cluster
    selected_variables = ultimate_variables[[variables_cluster]]
    
    #if else for deal with species with low and high number of ocurrences
    if (n_ocurrence<length(selected_variables)*10){ #if there is not 10 ocurrences for each variable of the cluster 
        
        #load the number of ocurrences
        load("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/finals/final_variables_low_number_ocurrence_species.rda") #the object is final_variables_low_number_ocurrence_species.rda

        #select the variables for this species from the list of reduced number in basis on raking
        reduced_variables = final_variables_low_number_ocurrence_species_new[[species]] #select element "species" of the list

        #select these variable from the stack of current variables 
        variables_stack = variables_stack[[reduced_variables]]    
    
    } else { #if not, and thus there are a lot of points, we can use all the variables for the corresponding cluster

        variables_stack = variables_stack[[selected_variables]]
    
    }         
 
    #extract names of the current variables for select soil variable in the loop
    names_variables_stack = names(variables_stack)

    #load raster of PA buffer
    raster_PA_buffer = raster(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/pseudo_absences", paste(species, "PA_buffer.asc", sep="_"), sep="/"))

    #crop variables using distribution buffer
    variables_stack = crop(variables_stack, raster_PA_buffer)

    #create a vector with all bio variables
    bio = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

    #select from selected variables only bio variables, which will be replaced by future variables
    future_selected = selected_variables[selected_variables %in% bio]
    
    #list of continuos projections and binary projections for all scenarios and climatic models
    climatic_scenarios = c("bc26", "bc45", "bc60", "bc85", "cc26", "cc45", "cc60", "cc85", "gs26", "gs45", "gs60", "gs85", "he26", "he45", "he60", "he85", "ip26", "ip45", "ip60", "ip85", "mr26", "mr45", "mr60", "mr85", "mg26", "mg45", "mg60", "mg85")
    
    #load models 
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/models", paste(species, "glm_model.rda", sep="_"), sep="/"))
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/models", paste(species, "gam_model.rda", sep="_"), sep="/"))
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/models", paste(species, "rf_model.rda", sep="_"), sep="/"))

    #load thresholds
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/threshold", paste(species, "glm_threshold.rda", sep="_"), sep="/"))
    load(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/threshold", paste(species, "gam_threshold.rda", sep="_"), sep="/"))

    #loop for run 12 glms, 12 gams and 12 rf for each climatic_model*IPCC_scenario
    for (i in climatic_scenarios){

        #load future climatic variables 
        list_future_climate = list.files("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/datos/future_climate/bioclim_moisture", pattern=i, full.names=TRUE) #list all raster of the [i] modelo*scenario 

        #stack all of them 
        future_climate = stack(list_future_climate)

        #crop the future variable with the PA buffer (only predict in the interesting area)
        future_climate = crop(future_climate, raster_PA_buffer)

        #resample the future climatic variables 
        future_climate = resample(future_climate, variables_stack[[1]], method="bilinear")
        
        #stack future climatic variables and the current soil variables
        final_future_climate = stack(
            future_climate[[paste(future_selected, i, sep="_")]], #select the selected variables using i (model*scenario) and "bioXX"
            variables_stack[[names_variables_stack[names_variables_stack %in% c("ph", "cec", "carbon", "depth", "sand", "silt", "clay")]]]) #select from the stack the only the SOIL selected variables, becuase of this we use a vector with all soil variables names
           
        #loop for change the names of variables and match with variable names in the models   
        require(stringr) #require for nchar (calculate the number of characters) and str_split_fixed

        final_names = NULL 
        for (k in names(final_future_climate)){ #for each name of the variables in the final stack

            name_splitted = str_split_fixed(k, "_", 2)
           
            if (nchar(name_splitted[2])>0){ #if the number has more than 0 characters
               
                final_names =  append(final_names, name_splitted[1]) #select the first part of the split (bio name)

            } else { #if not and thus is a soil variable 

                 final_names = append(final_names, k) #save exactly the same name 
            }
        }

        #change the names of the variables in the stack 
        names(final_future_climate) = final_names

        #list of continuos projections for each model 
        glm_projections = list()
        gam_projections = list()
        rf_projections = list()

        #list of binary projections for each model 
        glm_projection_bin = list()
        gam_projection_bin = list()
        rf_projection_bin = list()

        for (k in 1:12){

            glm_projections[[k]] = predict(final_future_climate, glm_resample[[k]], type="response")
            gam_projections[[k]] = predict(final_future_climate, gam_resample[[k]], type="response")
            rf_projections[[k]] = predict(final_future_climate, rf_resample[[k]], type="response")

            #convert each projection in binary using the threshold calculated with fit_eval_models. 
            #glm
            glm_projection_bin[[k]] = glm_projections[[k]] #copy the raster 
            glm_projection_bin[[k]][glm_projection_bin[[k]]>=(glm_threshold[[k]][2,2]),] <- 1 #give 1 to the pixels with a predicted value higher or equal than the threshold
            glm_projection_bin[[k]][glm_projection_bin[[k]]<(glm_threshold[[k]][2,2]),] <- 0 #give 0 to the pixels with a predicted value lower than the threshold

            #gam
            gam_projection_bin[[k]] = gam_projections[[k]] #copy the raster 
            gam_projection_bin[[k]][gam_projection_bin[[k]]>=(gam_threshold[[k]][2,2]),] <- 1 #give 1 to the pixels with a predicted value higher or equal than the threshold
            gam_projection_bin[[k]][gam_projection_bin[[k]]<(gam_threshold[[k]][2,2]),] <- 0 #give 0 to the pixels with a predicted value lower than the threshold

            #rf
            rf_projection_bin[[k]] = rf_projections[[k]] #It does not change  because it is already a binary projection 
        }

        #bind glm and gam projections of this scenario in a stack
        scenarios_projections_only_glm_gam = stack(stack(glm_projections), stack(gam_projections))

        #stack all binary projections
        scenarios_projections_bin_only_glm_gam  = stack(stack(glm_projection_bin), stack(gam_projection_bin))
        scenarios_projections_bin_only_rf = stack(rf_projection_bin)

        #save continous projections of glm and gam
        writeRaster(scenarios_projections_only_glm_gam, filename=paste(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/continous_projections", paste("continuous", paste("projection_only_glm_gam", i, sep="_"), sep="_"), sep="/"), paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)

        #save binary projections of gam and glm
        writeRaster(scenarios_projections_bin_only_glm_gam, filename=paste(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/binary_projections", paste("binary", paste(paste("projection", "only_glm_gam", sep="_"), i, sep="_"), sep="_"), sep="/"), paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)

        #save binary projectiosn of rf
        writeRaster(scenarios_projections_bin_only_rf, filename=paste(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/binary_projections", paste("binary", paste(paste("projection", "only_rf", sep="_"), i, sep="_"), sep="_"), sep="/"), paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
        file.remove(paste(paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/binary_projections", paste("binary", paste(paste("projection", "only_rf", sep="_"), i, sep="_"), sep="_"), sep="/"), paste(species, "tif.aux.xml", sep="."), sep="_")) #delte the .xlm file produced because in the next step I load all the data with pinaster name and .tf, in this situation theses .xlm files would give problems. 
    }

    ##reate the ensamble
    #load all projections
    list_projections = list.files("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/binary_projections", pattern=paste(species, "tif", sep="."), full.names=TRUE)
    all_projections = stack(list_projections)

    #calculate the percentage of models for which a pixel is suitable
    ensamble_projections_bin = calc(all_projections, function(x) (sum(x)*100)/nlayers(all_projections))

    #save ensamble
    writeRaster(ensamble_projections_bin, filename=paste("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/results/modelos/ensamble_projections_bin", paste("ensamble_projections_bin", paste(species, "tif", sep="."), sep="_"), sep="/"), overwrite=TRUE)

    #end the species
    print(paste("end",species))  
}

########Paralelize the process######
require(foreach) #for repeat a process several times 
require(doParallel) #for parallel

#create a vector with two species to test the functions
list_species = read.table("/Users/diegosalazar/Google Drive/academico/tesis/nicho_pinus/code/search/species.txt", sep="", header=T)
species = list_species$specific_epithet

# set up cluster
clust <- makeCluster(2) 
registerDoParallel(clust)

##run each process separately

#preparing data: ONLY ONE TIME 
foreach(species = species, .packages=c("raster", "dismo")) %dopar% { 
    prepar_data(species = species)
} 

#fitting and evaluation
foreach(species = species, .packages=c("raster", "dismo", "gam", "randomForest")) %dopar% { 
    fit_eval_models(species = species)
} 

#binarize current predictions
foreach(species = species, .packages=c("raster", "dismo", "ecospat")) %dopar% { 
    binarize_curent_projections(species = species)
} 

#project to the future
foreach(species = species, .packages=c("raster", "dismo", "randomForest", "gam")) %dopar% { 
    projections(species = species)
} 


#stop the cluster 
stopCluster(clust)
