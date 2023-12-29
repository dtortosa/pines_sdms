#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
    #https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
    #https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



#########################################################
####### loop that cleans ocurrences and save them for modelling ####### #########################################################



###################################################
##### DIFFERENCES RESPECT TO PREVIOUS VERSION #####
###################################################

#Respect to version 1:



########################
##### BEGIN SCRIPT #####
########################

##occurrences resampling

#species="radiata"
#exsitu="yes"
#moisture="yes"
#phylo="yes"
master_pipeline=function(species, exsitu, moisture, phylo){

    system(paste("mkdir -p ./results/pipeline_2024/", species, sep=""))

    system(paste("Rscript ./code/pipeline_2024/01_occurrences_v1.R ", species, " ", exsitu, sep=""))

}


master_pipeline(species="halepensis", exsitu="yes", moisture="yes", phylo="yes")
