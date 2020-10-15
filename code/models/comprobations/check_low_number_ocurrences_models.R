#SEWALL. code for check the models, evaluations and threshold of species with a low number of ocurrences in sewall. We extract files from .zips wihout unzip. 

#############
###library###
#############
require(randomForest)
require(gam) #we load parckages of models used for avoding problems during checking. Without RF packages, the model is not recognaised, and is showed a interminable list of numbers. 
require(dismo) #for seeing evaluations

#select the species
species = "sylvestris"

##############
### models ###
##############

#check the files included in the zip of the models
unzip(paste("/home/dsalazar/modelos/models/models_", species, ".zip", sep=""), files = NULL, list = TRUE, overwrite = TRUE, junkpaths = FALSE, exdir = ".", unzip = "internal", setTimes = FALSE)

#extract the rda files of models and training data
#when you use this code several times (open 16 connections), R fails, because yo can't opne more than 16. We run closeAllConnections() for avoiding this
closeAllConnections()
load(unz(paste("/home/dsalazar/modelos/models/models_", species, ".zip", sep=""), paste(species, "_glm_model.rda", sep="")))
load(unz(paste("/home/dsalazar/modelos/models/models_", species, ".zip", sep=""), paste(species, "_gam_model.rda", sep="")))
load(unz(paste("/home/dsalazar/modelos/models/models_", species, ".zip", sep=""), paste(species, "_rf_model.rda", sep="")))
load(unz(paste("/home/dsalazar/modelos/models/models_", species, ".zip", sep=""), paste(species, "_training_data.rda", sep="")))

#open the models
eval(parse(text="glm_resample")) #a way to call a variable with its name as string (http://stackoverflow.com/questions/9057006/getting-strings-recognized-as-variable-names-in-r?lq=1)
eval(parse(text="gam_resample"))
eval(parse(text="rf_resample")) #Info about OBB estimate of error rate: Random forests technique involves sampling of the input data with replacement (bootstrap sampling). In this sampling, about one thrird of the data is not used for training and can be used to testing.These are called the out of bag samples. Error estimated on these out of bag samples is the out of bag error. See (https://www.quora.com/What-is-the-out-of-bag-error-in-Random-Forests and http://stats.stackexchange.com/questions/109232/low-explained-variance-in-random-forest-r-randomforest/144375#144375). Random forest also select ranomd variables, thus each arbol is different in relation to the data and the variables used. When the process is finished, it is calculated for each point the probability of presence according to all the trees, and we take the mode. In this way we can avoid the overeffect of certain points, because the existente one single point can change the better line that separates presences of absences (see video 10 of 2015 SDM course of blas, min 13.)
    #The confusion matrix indicate the cases of error. The columns indicate the referencia and the rows the predictions. This, 0/0 and 1/1 are presences and pseudoabsences well predicted, whilst 0/1 is a presence predicted as absence, and 1/0 is a absence precited as presence. The total of points included in the matrix is equal to the rows of the training data (see http://stats.stackexchange.com/questions/90919/confusion-matrix-calculation-in-random-forest-classifier-in-r)

#open the training data
str(eval(parse(text="training_data")))


##################
### evaluation ###
##################

#check the files included in the zip of the evaluations
unzip(paste("/home/dsalazar/modelos/evaluations/evaluations_", species, ".zip", sep=""), files = NULL, list = TRUE, overwrite = TRUE, junkpaths = FALSE, exdir = ".", unzip = "internal", setTimes = FALSE)

#extract the rda files of the evaluations
load(unz(paste("/home/dsalazar/modelos/evaluations/evaluations_", species, ".zip", sep=""), paste(species, "_glm_evaluation.rda", sep="")))
load(unz(paste("/home/dsalazar/modelos/evaluations/evaluations_", species, ".zip", sep=""), paste(species, "_gam_evaluation.rda", sep="")))
load(unz(paste("/home/dsalazar/modelos/evaluations/evaluations_", species, ".zip", sep=""), paste(species, "_rf_evaluation.rda", sep="")))

#open the evaluations
eval(parse(text="glm_evaluation"))
eval(parse(text="gam_evaluation"))
eval(parse(text="rf_evaluation"))


##################
### thresholds ###
##################

#check the files included in the zip of the threshold
unzip(paste("/home/dsalazar/modelos/threshold/threshold_", species, ".zip", sep=""), files = NULL, list = TRUE, overwrite = TRUE, junkpaths = FALSE, exdir = ".", unzip = "internal", setTimes = FALSE)

#check the files included in the zip of the models
load(unz(paste("/home/dsalazar/modelos/threshold/threshold_", species, ".zip", sep=""), paste(species, "_glm_threshold.rda", sep="")))
load(unz(paste("/home/dsalazar/modelos/threshold/threshold_", species, ".zip", sep=""), paste(species, "_gam_threshold.rda", sep="")))

#open threshold
eval(parse(text="glm_threshold"))
eval(parse(text="gam_threshold"))

####Results
#I have checked the fitting of models of some species randomly selected and the models of all species with lower than 50 ocurrences (16 species), which are species restricted to small areas. In general, models of species with enough ocurrences (ie more than 50) have a good fitting: Maximum values of kappa higher than 0.7 (even above 0.9 in some cases; glm-gam) and OOB estimates of error rate below 3% (random forest). In contrast, models fitted with lower ocurrences exhibit worse perfomance, not in all species, but in some of them: maximum kappa below 0.6 and OOB above 5%. These species are being problematic, and we will have to discuss if it is a good idea use them. 




