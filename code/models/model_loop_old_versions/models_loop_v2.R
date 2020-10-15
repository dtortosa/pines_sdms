#SEWAL. Code for modelling and project into the future. It is prepared for run in sewall.  

###################################
#ESTABLECE EL DIRECTORIO DE TRABAJO
###################################
#DIRECTORIO DE TRABAJO
setwd("/home/dsalazar/")


####################################################
########### Loop preparing data #########################
####################################################
##RUN only ONE time. Already runned. 

#function for preparing data (extract data of variables in ocurrences points)
prepar_data = function(species) {

    #required libraries
    #library(raster)
    #library(dismo)

    #load the group of species according to cluster
    group_species = read.csv("/home/dsalazar/data/complete_2_10_g_5.csv", header=TRUE)  
    #select the group of the corresponding species
    variables_cluster = group_species[which(group_species$species == species),]$groups 

    #select the number cluster of this species 
    rasters_list = list.files("/home/dsalazar/climate/current", pattern=".asc", full.names=TRUE) #list the corresponding group of variables
    rasters_names = list.files("/home/dsalazar/climate/current", pattern=".asc", full.names=FALSE) #list names of raster
    names_variables = NULL #loop for separate raster names from extension ".asc"
    for (i in rasters_names){
        names_variables = append(names_variables, strsplit(i, split=".asc")[[1]])
    }
    variables_stack = stack(rasters_list) #stack them 
    names(variables_stack) = names_variables #give the names to layers of the stack

    #load names of selected variables  
    load("/home/dsalazar/data/list_selected_variables.rda") 

    #take selected variables from the list of variables
    selected_variables = ultimate_variables[[variables_cluster]]

    #calculate the number of ocurrences for each species 
    number_ocurrences = read.csv("/home/dsalazar/data/ocurrences_per_species.csv", header=TRUE)
    n_ocurrence = number_ocurrences[number_ocurrences$species==species,]$number_ocurrences

    #load names of selected variables 
    if (n_ocurrence<length(selected_variables)*10){ #if there is not 10 ocurrences for each variable 
        
        #load the names of selected variables of low number ocurrences species
        load("/home/dsalazar/data/final_variables_low_number_ocurrence_species.rda")

        #select the selecte variables for low number ocurrences species
        variables_stack = variables_stack[[final_variables_low_number_ocurrence_species_new[[species]]]]    
    
    } else {

        variables_stack = variables_stack[[selected_variables]]
    
    }     
 
    ##Extract values of variables 
    presences = read.csv(paste("/home/dsalazar/data/ocurrences", paste(species, "complete.presences.csv", sep="_"), sep="/")) #read the data final data with presences and pseudoabsences
    variables = extract(variables_stack, presences[, c("longitude", "latitude")]) #extract the value of the variable in these points
    data = cbind(presences, variables) #bind the presence data and the variable data in one data frame 

    #change the names of variables if the number of selected variables is 1 (it is to say, only 5 columns in data)
    if (ncol(data)<6){
        colnames(data)[5] = names(variables_stack)
    }

    #write the final.data file
    write.csv(data, paste("/home/dsalazar/data/data_prepare_modelling", paste(species, "csv", sep="."), sep="/"), row.names=FALSE) 

    #name of species
    print(paste(species, "ended"))   
}


####################################################
########### Loop of fitting ########################
####################################################

#In each step we will create a dataset for training and evaluation and make both process, then we will repeat with other partition of training and evaluation. 
#We have to create a loop with for each species and before create empty lists (glm_species...) for include all the models for each species. In these lists we will sabe XX_resample lists


#function for fitting models of current habitat suitability
fit_eval_models = function(species){ #for the corresponding species 
    
    #required libraries
    #library(raster)
    #library(dismo)
    #library(gam) #OTRA LIBRERIA PARA GAM
    #library(randomForest) #RANDOM FOREST

    #begin the species
    print(paste("begin",species))

    #load cluster number for each species
    group_species = read.csv("/home/dsalazar/data/complete_2_10_g_5.csv", header=TRUE)[,c("species", "groups")] 

    #select the group of the corresponding species
    variables_cluster = group_species[which(group_species$species == species),]$groups 

    #Load data with presences and values of environmental variables
    data = read.csv(paste("/home/dsalazar/data/data_prepare_modelling", paste(species, "csv", sep="."), sep="/"), header=TRUE)

    #select the number cluster of this species 
    rasters_list = list.files("/home/dsalazar/climate/current", pattern=".asc", full.names=TRUE) #list the corresponding group of variables
    rasters_names = list.files("/home/dsalazar/climate/current", pattern=".asc", full.names=FALSE) #list names of raster
    names_variables = NULL #loop for separate raster names from extension ".asc"
    for (i in rasters_names){
        names_variables = append(names_variables, strsplit(i, split=".asc")[[1]])
    }
    variables_stack = stack(rasters_list) #stack them 
    names(variables_stack) = names_variables #give the names to layers of the stack

    #load names of selected variables  
    load("/home/dsalazar/data/list_selected_variables.rda") 

    #take selected variables from the list of variables
    selected_variables = ultimate_variables[[variables_cluster]]

    #calculate the number of ocurrences for each species 
    number_ocurrences = read.csv("/home/dsalazar/data/ocurrences_per_species.csv", header=TRUE)
    n_ocurrence = number_ocurrences[number_ocurrences$species==species,]$number_ocurrences

    if (n_ocurrence<length(selected_variables)*10){ #if there is not 10 ocurrences for each variable 
        
        #load the names of selected variables for los number ocurrences species
        load("/home/dsalazar/data/final_variables_low_number_ocurrence_species.rda")

        #select the selecte variables for low number ocurrences species
        variables_stack = variables_stack[[final_variables_low_number_ocurrence_species_new[[species]]]]    
    
    } else {

        variables_stack = variables_stack[[selected_variables]]
    
    }       
 
    #load raster of PA buffer
    raster_PA_buffer = raster(paste("/home/dsalazar/data/pa_buffers", paste(species, "PA_buffer.asc", sep="_"), sep="/"))

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

    #list for save training data
    training_data = list()

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
        glm_pseudo_absence = step(glm_pseudo_absence, direction="both", trace=0) #direction: the mode the mode of stepwise search, can be one of ‘"both"’, ‘"backward"’, or ‘"forward"’, with a default of ‘"both"’. 
            #trace: if positive, information is printed during the running of ‘step’. Larger values may give more detailed information. 

        #####################################
        #MODELOS ADITIVOS GENERALIZADOS (GAM)
        #####################################
        formula.gam<-as.formula(paste("presence ~ s(", paste(names(training)[-c(1:4)], collapse=",4) + s("), ",4)", collapse="")) #cada varaible está envuelta en una funcion que suaviza (smooth), gam traaja sobre valore suavizados de las variables. Por eso es s(). 
        if (is.element("package:mgcv", search())) detach("package:mgcv", force=TRUE) #make sure the mgcv package is not loaded to avoid conflicts with gam package.   
        gam_pseudo_absence<-gam(formula.gam, family=binomial(link=logit), data=training, weights=precision_weight)

        #make stepwise
        #genereate a scope for step.gam
        #Given a data.frame as an argument, generate a scope list for use in step.gam, each element of which gives the candidates for that term.
        gam_scope = gam.scope(frame=training[,-c(1:3)], response=1, smoother="s", arg="df=4") #frame is the data: response variable and explicativa variables; response indicate the column of the response; smoother is the type of smooth that we use; arg: additional argument, like for example the freedom degree of smooth. 
        gam_pseudo_absence = step.gam(gam_pseudo_absence, scope=gam_scope, direction="both", trace=0, data=training) #we include the list of terms in scope. 


        ##############
        #RANDOM FOREST
        ##############

        ###WEIGHT the random forest analysis###
        ##Strata variable
        #create a strata variable with diferent value for each precision_weight. This variables let us increase the probability of sample of data in relation to precision weight. Points with higher precision weight were more proably sampled than those with lower precision weight. We will give a letter to each point according to the precision weight (A = high precision presences; B = low precision presences; C = pseudoabsences).  
        if(nrow(training[training$precision_weight==1,])>0){ #if there are points with high precision (we have to consider the posibility that some species don't have high precision points)
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
        } else { #if there is not high precision points 
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

        ##N sample size
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
        ############Save training data ####
        ###################################
        training_data[[k]]=training

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

    #save the training data
    save(training_data, file=paste("/home/dsalazar/modelos/models", paste(species, "training_data.rda", sep="_"), sep="/"))

    #save models 
    save(glm_resample, file=paste("/home/dsalazar/modelos/models", paste(species, "glm_model.rda", sep="_"), sep="/"))
    save(gam_resample, file=paste("/home/dsalazar/modelos/models", paste(species, "gam_model.rda", sep="_"), sep="/"))
    save(rf_resample, file=paste("/home/dsalazar/modelos/models", paste(species, "rf_model.rda", sep="_"), sep="/"))  

    #save the continuous predictions
    writeRaster(continuous_predictions_glm, filename=paste("/home/dsalazar/modelos/continuous_predictions/continuous_predictions_glm", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
    writeRaster(continuous_predictions_gam, filename=paste("/home/dsalazar/modelos/continuous_predictions/continuous_predictions_gam", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)

    #save the predictions of random forest (binary)
    writeRaster(binary_predictions_rf, filename=paste("/home/dsalazar/modelos/binary_predictions/binary_predictions_rf", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)

    #save data for evaluation
    save(glm_evaluation_predict, file=paste("/home/dsalazar/modelos/evaluations", paste(species, "glm_evaluation_data.rda", sep="_"), sep="/"))
    save(gam_evaluation_predict, file=paste("/home/dsalazar/modelos/evaluations", paste(species, "gam_evaluation_data.rda", sep="_"), sep="/"))
    save(rf_evaluation_predict, file=paste("/home/dsalazar/modelos/evaluations", paste(species, "rf_evaluation_data.rda", sep="_"), sep="/")) 

    #save evaluations
    save(glm_evaluation, file=paste("/home/dsalazar/modelos/evaluations", paste(species, "glm_evaluation.rda", sep="_"), sep="/"))
    save(gam_evaluation, file=paste("/home/dsalazar/modelos/evaluations", paste(species, "gam_evaluation.rda", sep="_"), sep="/"))
    save(rf_evaluation, file=paste("/home/dsalazar/modelos/evaluations", paste(species, "rf_evaluation.rda", sep="_"), sep="/"))

    #name of species
    print(paste(species, "ended"))  
}

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
    load(paste("/home/dsalazar/modelos/models", paste(species, "glm_model.rda", sep="_"), sep="/"))
    load(paste("/home/dsalazar/modelos/models", paste(species, "gam_model.rda", sep="_"), sep="/"))
    load(paste("/home/dsalazar/modelos/models", paste(species, "rf_model.rda", sep="_"), sep="/"))  

    #load evaluations
    load(paste("/home/dsalazar/modelos/evaluations", paste(species, "glm_evaluation.rda", sep="_"), sep="/"))
    load(paste("/home/dsalazar/modelos/evaluations", paste(species, "gam_evaluation.rda", sep="_"), sep="/"))
    load(paste("/home/dsalazar/modelos/evaluations", paste(species, "rf_evaluation.rda", sep="_"), sep="/")) 

    #load evaluation data
    load(paste("/home/dsalazar/modelos/evaluations", paste(species, "glm_evaluation_data.rda", sep="_"), sep="/"))
    load(paste("/home/dsalazar/modelos/evaluations", paste(species, "gam_evaluation_data.rda", sep="_"), sep="/"))
    load(paste("/home/dsalazar/modelos/evaluations", paste(species, "rf_evaluation_data.rda", sep="_"), sep="/")) 

    #load continuous predictions of glm 
    glm_predict = stack(paste("/home/dsalazar/modelos/continuous_predictions/continuous_predictions_glm", paste(species, "tif", sep="."), sep="_"))

    #load continuous predictions of gam     
    gam_predict = stack(paste("/home/dsalazar/modelos/continuous_predictions/continuous_predictions_gam", paste(species, "tif", sep="."), sep="_"))

    #load binary predictions of Random forest (already finished)
    rf_predict = stack(paste("/home/dsalazar/modelos/binary_predictions/binary_predictions_rf", paste(species, "tif", sep="."), sep="_"))

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
    save(glm_threshold, file=paste("/home/dsalazar/modelos/threshold", paste(species, "glm_threshold.rda", sep="_"), sep="/"))
    save(gam_threshold, file=paste("/home/dsalazar/modelos/threshold", paste(species, "gam_threshold.rda", sep="_"), sep="/"))

    #save the binary predictions of glm and gam (rf make them before)
    writeRaster(stack(glm_prediction_bin), filename=paste("/home/dsalazar/modelos/binary_predictions/binary_predictions_glm", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
    writeRaster(stack(gam_prediction_bin), filename=paste("/home/dsalazar/modelos/binary_predictions/binary_predictions_gam", paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)

    #save ensamble
    writeRaster(ensamble_predictions_bin, filename=paste("/home/dsalazar/modelos/ensamble_predictions_bin", paste("ensamble_predictions_bin", paste(species, "tif", sep="."), sep="_"), sep="/"), overwrite=TRUE)


    #########################################
    ######ZIP AND DELETE EVALUATIONS ########
    ######################################### 
    ####FOR RELEASE SPACE IN SEAWELL
    #list continuous projections created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files
        evaluations = list.files("/home/dsalazar/modelos/evaluations", pattern=glob2rx(paste("*", species, "*", "rda", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        evaluations = evaluations[!grepl("pseudostrobus", evaluations)]
    } else {
        evaluations = list.files("/home/dsalazar/modelos/evaluations", pattern=glob2rx(paste("*", species, "*", "rda", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE) 
    }

    #zip all continuous projections for lack of space in disk
    zip(paste("/home/dsalazar/modelos/evaluations/evaluations", paste(species, "zip", sep="."), sep="_"), evaluations, flags="-j") #j indicate that you don't want all the directory structure

    #delete continuous projections
    file.remove(evaluations)  

    ########################################################
    ######ZIP AND DELETE BIN AND CONTINUOUS PREDICTIONS ####
    ######################################################## 
    #list continuous predictions created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files    
        continuous_predictions = list.files("/home/dsalazar/modelos/continuous_predictions", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        continuous_predictions = continuous_predictions[!grepl("pseudostrobus", continuous_predictions)]
    } else {
        continuous_predictions = list.files("/home/dsalazar/modelos/continuous_predictions", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
    } 

    #zip all continuous predictions for lack of space in disk
    zip(paste("/home/dsalazar/modelos/continuous_predictions/continuous_predictions", paste(species, "zip", sep="."), sep="_"), continuous_predictions, flags="-j") #j indicate that you don't want all the directory structure

    #delete continuous predictions
    file.remove(continuous_predictions)

    #list binary predictions created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files    
        binary_predictions = list.files("/home/dsalazar/modelos/binary_predictions", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        binary_predictions = binary_predictions[!grepl("pseudostrobus", binary_predictions)]        
    } else {
        binary_predictions = list.files("/home/dsalazar/modelos/binary_predictions", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)        
    }
    #zip all binary predictions for lack of space in disk
    zip(paste("/home/dsalazar/modelos/binary_predictions/binary_predictions", paste(species, "zip", sep="."), sep="_"), binary_predictions, flags="-j") #j indicate that you don't want all the directory structure

    #delete binary predictions
    file.remove(binary_predictions)
 
    #name of species
    print(paste(species, "ended")) 
}


#function for make future projections
projections = function(species){

    #required libraries
    #library(raster)
    #library(dismo)
    #library(randomForest) #for predict rf models
    #library(gam) #for predict gam model
    #library(rgdal) #for problem in ensamble predictions

    #begin the species
    print(paste("begin",species))

    #load cluster number for each species
    group_species = read.csv("/home/dsalazar/data/complete_2_10_g_5.csv", header=TRUE)[,c("species", "groups")] 

    #select the group of the corresponding species
    variables_cluster = group_species[which(group_species$species == species),]$groups 

    #Load data with presences and values of environmental variables
    data = read.csv(paste("/home/dsalazar/data/data_prepare_modelling", paste(species, "csv", sep="."), sep="/"), header=TRUE)

    #select the number cluster of this species 
    rasters_list = list.files("/home/dsalazar/climate/current", pattern=".asc", full.names=TRUE) #list the corresponding group of variables
    rasters_names = list.files("/home/dsalazar/climate/current", pattern=".asc", full.names=FALSE) #list names of raster
    names_variables = NULL #loop for separate raster names from extension ".asc"
    for (i in rasters_names){
        names_variables = append(names_variables, strsplit(i, split=".asc")[[1]])
    }
    variables_stack = stack(rasters_list) #stack them 
    names(variables_stack) = names_variables #give the names to layers of the stack

    #load names of selected variables  
    load("/home/dsalazar/data/list_selected_variables.rda") 

    #take selected variables from the list of variables
    selected_variables = ultimate_variables[[variables_cluster]]

    #calculate the number of ocurrences for each species 
    number_ocurrences = read.csv("/home/dsalazar/data/ocurrences_per_species.csv", header=TRUE)
    n_ocurrence = number_ocurrences[number_ocurrences$species==species,]$number_ocurrences

    if (n_ocurrence<length(selected_variables)*10){ #if there is not 10 ocurrences for each variable 
        
        #load the names of selected variables for los number ocurrences species
        load("/home/dsalazar/data/final_variables_low_number_ocurrence_species.rda")

        #select the selecte variables for low number ocurrences species
        variables_stack = variables_stack[[final_variables_low_number_ocurrence_species_new[[species]]]]    
    
    } else {

        variables_stack = variables_stack[[selected_variables]]
    
    }       
 
    #extract names of the current variables for select soil variable in the loop
    names_variables_stack = names(variables_stack)

    #load raster of PA buffer
    raster_PA_buffer = raster(paste("/home/dsalazar/data/pa_buffers", paste(species, "PA_buffer.asc", sep="_"), sep="/"))

    #crop variables using distribution buffer
    variables_stack = crop(variables_stack, raster_PA_buffer)

    #create a vector with all bio variables
    bio = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

    #select from selected variables only bio variables, which will be replaced by future variables
    if (n_ocurrence<length(selected_variables)*10){
        future_selected = final_variables_low_number_ocurrence_species_new[[species]][final_variables_low_number_ocurrence_species_new[[species]] %in% bio]        
    } else {
        future_selected = selected_variables[selected_variables %in% bio]
    }
    #list of continuos projections and binary projections for all scenarios and climatic models
    climatic_scenarios = c("bc26", "bc45", "bc60", "bc85", "cc26", "cc45", "cc60", "cc85", "gs26", "gs45", "gs60", "gs85", "he26", "he45", "he60", "he85", "ip26", "ip45", "ip60", "ip85", "mg26", "mg45", "mg60", "mg85", "mr26", "mr45", "mr60", "mr85")
    
    #load models 
    load(paste("/home/dsalazar/modelos/models", paste(species, "glm_model.rda", sep="_"), sep="/"))
    load(paste("//home/dsalazar/modelos/models", paste(species, "gam_model.rda", sep="_"), sep="/"))
    load(paste("/home/dsalazar/modelos/models", paste(species, "rf_model.rda", sep="_"), sep="/"))

    #load thresholds
    load(paste("/home/dsalazar/modelos/threshold", paste(species, "glm_threshold.rda", sep="_"), sep="/"))
    load(paste("/home/dsalazar/modelos/threshold", paste(species, "gam_threshold.rda", sep="_"), sep="/"))

    #loop for run 12 glms, 12 gams and 12 rf for each climatic_model*IPCC_scenario
    for (i in climatic_scenarios){

        #load future climatic variables 
        list_future_climate = list.files("/home/dsalazar/climate/future", pattern=i, full.names=TRUE) #list all raster of the [i] modelo*scenario 
        list_names_climatic_variables = list.files("/home/dsalazar/climate/future", pattern=i) #list names 
        names_climatic_variables = NULL #drop .asc extension
        for (k in list_names_climatic_variables){
            names_climatic_variables = append(names_climatic_variables, strsplit(k, split=".asc")[[1]])
        }

        #stack all of them 
        future_climate = stack(list_future_climate)
        names(future_climate) = names_climatic_variables #give names to the layers

        #crop the future variable with the PA buffer (only predict in the interesting area)
        future_climate = crop(future_climate, raster_PA_buffer)

        #resample the future climatic variables 
        future_climate = resample(future_climate, variables_stack[[1]], method="bilinear")
        
        #stack future climatic variables and the current soil variables
        if(nlayers(variables_stack)>1){
            final_future_climate = stack(
                future_climate[[paste(future_selected, i, sep="_")]], #select the selected variables using i (model*scenario) and "bioXX"
                variables_stack[[names_variables_stack[names_variables_stack %in% c("ph", "cec", "carbon", "depth", "sand", "silt", "clay")]]]) #select from the stack the only the SOIL selected variables, becuase of this we use a vector with all soil variables names
        } else {
            final_future_climate = stack(future_climate[[paste(future_selected, i, sep="_")]])
        }   

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
        writeRaster(stack(glm_projections), filename=paste(paste("/home/dsalazar/modelos/continuous_projections", paste("continuous", paste("projection_glm", i, sep="_"), sep="_"), sep="/"), paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
        writeRaster(stack(gam_projections), filename=paste(paste("/home/dsalazar/modelos/continuous_projections", paste("continuous", paste("projection_gam", i, sep="_"), sep="_"), sep="/"), paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)        

        #save binary projections of gam and glm
        writeRaster(stack(glm_projection_bin), filename=paste(paste("/home/dsalazar/modelos/binary_projections", paste("binary", paste(paste("projection", "glm", sep="_"), i, sep="_"), sep="_"), sep="/"), paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
        writeRaster(stack(gam_projection_bin), filename=paste(paste("/home/dsalazar/modelos/binary_projections", paste("binary", paste(paste("projection", "gam", sep="_"), i, sep="_"), sep="_"), sep="/"), paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)

        #save binary projectiosn of rf
        writeRaster(stack(rf_projection_bin), filename=paste(paste("/home/dsalazar/modelos/binary_projections", paste("binary", paste(paste("projection", "rf", sep="_"), i, sep="_"), sep="_"), sep="/"), paste(species, "tif", sep="."), sep="_"), options="INTERLEAVE=BAND", overwrite=TRUE)
    }

    ##reate the ensamble
    #load all projections
    list_projections = list.files("/home/dsalazar/modelos/binary_projections", pattern=glob2rx(paste("*", species, "*", "tif", sep="")), full.names=TRUE)
    all_projections = stack(list_projections)

    #calculate the percentage of models for which a pixel is suitable
    ensamble_projections_bin = calc(all_projections, function(x) (sum(x)*100)/nlayers(all_projections))

    #save ensamble
    writeRaster(ensamble_projections_bin, filename=paste("/home/dsalazar/modelos/ensamble_projections_bin", paste("ensamble_projections_bin", paste(species, "tif", sep="."), sep="_"), sep="/"), overwrite=TRUE)

    #########################################
    ######ZIP AND DELETE MODELS #############
    ######################################### 
    #list continuous projections created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files    
        models = list.files("/home/dsalazar/modelos/models", pattern=glob2rx(paste("*", species, "*", "rda", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        models = models[!grepl("pseudostrobus", models)]
    } else {
        models = list.files("/home/dsalazar/modelos/models", pattern=glob2rx(paste("*", species, "*", "rda", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)        
    } 

    #zip all continuous projections for lack of space in disk
    zip(paste("/home/dsalazar/modelos/models/models", paste(species, "zip", sep="."), sep="_"), models, flags="-j") #j indicate that you don't want all the directory structure

    #delete continuous projections
    file.remove(models) #binary will be used in the ensamble, because of this we don't delete them now    

    #########################################
    ######ZIP AND DELETE THRESHOLDS #########
    ######################################### 
    #list continuous projections created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files        
        threshold = list.files("/home/dsalazar/modelos/threshold", pattern=glob2rx(paste("*", species, "*", "rda", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        threshold = threshold[!grepl("pseudostrobus", threshold)]
    } else {
        threshold = list.files("/home/dsalazar/modelos/threshold", pattern=glob2rx(paste("*", species, "*", "rda", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
    }

    #zip all continuous projections for lack of space in disk
    zip(paste("/home/dsalazar/modelos/threshold/threshold", paste(species, "zip", sep="."), sep="_"), threshold, flags="-j") #j indicate that you don't want all the directory structure

    #delete continuous projections
    file.remove(threshold)  

    ########################################################
    ######ZIP AND DELETE BIN AND CONTINUOUS PROJECTIONS ####
    ######################################################## 
    #list continuous projections created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files       
        continuous_projections = list.files("/home/dsalazar/modelos/continuous_projections", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        continuous_projections = continuous_projections[!grepl("pseudostrobus", continuous_projections)]
    } else {
        continuous_projections = list.files("/home/dsalazar/modelos/continuous_projections", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)        
    }
    #zip all continuous projections for lack of space in disk
    zip(paste("/home/dsalazar/modelos/continuous_projections/continuous_projections", paste(species, "zip", sep="."), sep="_"), continuous_projections, flags="-j") #j indicate that you don't want all the directory structure

    #delete continuous projections
    file.remove(continuous_projections) #binary will be used in the ensamble, because of this we don't delete them now    

    #list binary projections created 
    if (species=="strobus"){ #conditional for avoiding problems with strobus and pseudostrobus listing files           
        binary_projections = list.files("/home/dsalazar/modelos/binary_projections", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)
        binary_projections = binary_projections[!grepl("pseudostrobus", binary_projections)]
    } else {
        binary_projections = list.files("/home/dsalazar/modelos/binary_projections", pattern=glob2rx(paste("*", species, "*", "tif", sep=""), trim.head=TRUE, trim.tail=FALSE), full.names=TRUE)        
    }
    #zip all binary projections for lack of space in disk
    zip(paste("/home/dsalazar/modelos/binary_projections/binary_projections", paste(species, "zip", sep="."), sep="_"), binary_projections, flags="-j") #j indicate that you don't want all the directory structure

    #delete binary projections
    file.remove(binary_projections)

    #name of species
    print(paste(species, "ended"))  
}

########Paralelize the process######
require(foreach)
require(doParallel) #for parallel

#create a vector with species names
list_species = read.table("/home/dsalazar/data/species.txt", sep="", header=T)
species = as.vector(list_species$specific_epithet)

#load problematic species
list_problematic_species = read.csv("/home/dsalazar/data/problematic_species.csv", header=TRUE)
problematic_species = as.vector(list_problematic_species$specific_epithet)

#select non_problematic species from the pool of species
non_problematic_species = species[!species %in% problematic_species]

# set up cluster
clust <- makeCluster(6) 
registerDoParallel(clust)

#####################
###PREPARING DATA####
#####################
#preparing data for non-problematic species
foreach(i =non_problematic_species, .packages=c("raster", "dismo")) %dopar% { 
    prepar_data(species = i)
} 

#preparing data for problematic species
foreach(i = problematic_species, .packages=c("raster", "dismo")) %dopar% { 
    prepar_data(species = i)
} 

###########################
###FITTING & EVALUATION####
###########################
#fitting and evaluation for non-problematic species
foreach(i = non_problematic_species, .packages=c("raster", "dismo", "gam", "randomForest")) %dopar% { 
    fit_eval_models(species = i)
} 

#fitting and evaluation for problematic species
foreach(i = problematic_species, .packages=c("raster", "dismo", "gam", "randomForest")) %dopar% { 
    fit_eval_models(species = i)
} 

###########################
#########BINARIZE##########
###########################
#binarize current predictions for non-problematic species
foreach(i = non_problematic_species, .packages=c("raster", "dismo", "ecospat", "gtools")) %dopar% { 
    binarize_curent_projections(species = i)
} 

#binarize current predictions for problematic species
foreach(i = problematic_species, .packages=c("raster", "dismo", "ecospat", "gtools")) %dopar% { 
    binarize_curent_projections(species = i)
} 

###########################
#########PROJECTIONS#######
###########################
#project to the future for non_problematic species
foreach(i = non_problematic_species, .packages=c("raster", "dismo", "randomForest", "gam", "rgdal")) %dopar% { 
    projections(species = i)
} 


#project to the future for problematic species
foreach(i = problematic_species, .packages=c("raster", "dismo", "randomForest", "gam", "rgdal")) %dopar% { 
    projections(species = i)
} 

 
#stop the cluster 
stopCluster(clust)
