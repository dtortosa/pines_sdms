#!/usr/bin/env Rscript

#we use the shebang "#!/usr/bin/env Rscript" to facilitate portability between machines. "env" offers more portability because it searches for the first occurrence of the R interpreter in the environment where we are running. It should thus look for R version installed in the container
    #https://www.r-bloggers.com/2019/11/r-scripts-as-command-line-tools/
    #https://www.baeldung.com/linux/bash-shebang-lines

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
    #https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
    #https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



##################################
####### CREATE SLURM FILES ####### 
##################################

#create slurm and bash scripts



###################################################
##### DIFFERENCES RESPECT TO PREVIOUS VERSION #####
###################################################

#Respect to version 1:



########################
##### BEGIN SCRIPT #####
########################

#set the seed for reproducibility
set.seed(56756)

#create some folders
system("mkdir -p ./code/phylo/recipes/01_global_test_phylo_ubuntu_20_04_v1_slurm_files")

#load species names
list_species = read.table("code/presences/species.txt", sep="\t", header=TRUE)

#extract epithet from species list
epithet_species_list = NULL
for(i in 1:nrow(list_species)){

    #selected species
    selected_species = as.vector(list_species[i,])

    #extract epithet
    epithet_species_list = append(epithet_species_list, strsplit(selected_species, split=" ")[[1]][2])
}
#check there is no NA
summary(!is.na(epithet_species_list))
#check
if(FALSE){
    require(tidyverse)
    paste("Pinus", epithet_species_list, sep=" ") == str_trim(as.vector(list_species[,1])) #the seperated epithet more Pinus are equal to list_species? We use str_trim from tidyverse to remove end spaces in each element of list_species
}#it is in false because loading tidyverse lead to load several packages that have a function name "extract", and this gives problems with the extract function of raster. If you want to check run these lines manually

#remove tecunumanii, jaliscana y discolor. These species are not used for Niche paper. The two first because we should downloadad gbif data now, so we would mix gbif data form 2016 and 2019. The third one was impossible to differentiate from P. cembroides
epithet_species_list = epithet_species_list[which(!epithet_species_list %in% c("tecunumanii", "jaliscana", "discolor"))]
#check these species are not present
!c("tecunumanii", "jaliscana", "discolor") %in% epithet_species_list

#extract occurrence data
#we are using data from Perret et al 2018
    #the dataset of Perret is carefully curated, obtaining data from herbarium specimens and accounts in the literature. They only consider as naturalized those records with clear evidence of self-sustainability and geographical coordinates. This solves a very important limitation compared to just use our GBIF data, as we do not include occurrences like from gardens... This is also a great independent validation, as the data does not come from GBIF, it is completely independent data. Different regions and different data sources.
    #Naturalized distributions show that climatic disequilibrium is structured by niche size in pines (Pinus L.)

#extract and load the naturalized occurrences from zip file
if(!file.exists("./results/global_test_phylo_current/exsitu_occurrences/pinus_occurrences_fordryad_exoticonly.csv")){
    system(paste(" \\
        unzip \\
            -o \\
                ./datos/phlyo/method_validation/doi_10_5061_dryad_1hr1n52__v20181213.zip \\
                pinus_occurrences_fordryad_exoticonly.csv \\
            -d ./results/global_test_phylo_current/exsitu_occurrences/", sep="")) 
}
naturalized_occurrences=read.csv("./results/global_test_phylo_current/exsitu_occurrences/pinus_occurrences_fordryad_exoticonly.csv", header=TRUE)
unique_species=unique(naturalized_occurrences$species)
summary(naturalized_occurrences)
#check
if(nrow(naturalized_occurrences)!=597){
    stop("ERROR! FALSE! WE DO NOT HAVE ALL NATURALIZED OCCURRENCES WE SHOULD HAVE ACCORDING TO THE PERRET'S MANUSCRIPT")
}
if(sum(unique_species %in% epithet_species_list)!= length(unique_species)){
    stop("ERROR! FALSE! WE DO NOT HAVE THE SAME SPECIES NAMES IN PERRET DATA")
}




##################################
###### CHECK OUTPUTS SPECIES #####
##################################
#define function
#species=unique_species[1]
check_outputs_species=function(species){

    #output file path
    path_output_file=paste("./results/global_test_phylo_current/species_output_files/", species, ".txt", sep="")


    ##check warnings about suitability bin sample size
    #check whether we have the warning about suitability bin size
    count_non_phylo_error_bin = system(paste(" \\
        grep \\
            --count \\
            'ERROR! FALSE! WE HAVE AT LEAST 1 NON-PHYLO SUITABILITY BIN WITH LESS THAN 60 DATA.POINTS FOR SPECIES ", species, "' \\
            ", path_output_file, sep=""), intern=TRUE)
        #--count: print only a count of selected lines

    #if the check appears 1 or several times
    #we can have several times the same error because this is checked for each partition
    if(count_non_phylo_error_bin[1]>0){

        #get the line number where eval_phylo starts
        first_line_eval_phylo=as.numeric(system(paste(" \\
            grep \\
                --line-number \\
                'STARTING predict_eval_phylo FOR  ", species, "' \\
                ", path_output_file, "| \\
            cut \\
                --delimiter : \\
                --fields 1", sep=""), intern=TRUE))
            #get line number with grep and, from the output, get the first element before :, which is the line number
            #https://stackoverflow.com/a/21412398

        #now get the next line
        second_line_eval_phylo=system(paste("\\
            awk \\
                '{if(NR==", first_line_eval_phylo+1, "){print $0}}' \\
                ", path_output_file, sep=""), intern=TRUE)

        #the next line should be the ENDING because no analyses should be done when the ERROR previously indicated appears, if not, we have a problem
        if(!grepl(paste("ENDING predict_eval_phylo FOR  ", species, sep=""), second_line_eval_phylo, fixed=TRUE)){
            stop(paste("ERROR! FALSE! PHYLO EVAL SHOULD NOT HAVE BEEN RUN FOR ", species, sep=""))
        }
    }


    #count cases with error about suitability bin sample size
    count_phylo_error_bin=system(paste("\\
        grep \\
            --count \\
            'ERROR! FALSE! WE HAVE AT LEAST 1 PHYLO SUITABILITY BIN WITH LESS THAN 60 DATA.POINTS FOR PHYLO MODEL ' \\
        ", path_output_file, sep=""), intern=TRUE)

    #if appears more than 0 times
    #it can appear several times as this is checked by partition and phylo model
    if(count_phylo_error_bin[1]>0){
        stop(paste("ERROR! FALSE! WE HAVE A PROBLEM WITH PHYLO SUITABILITY BIN SAMPLE SIZE FOR ONE OF THE PHYLO MODELS FOR SPECIES ", species, "CHECK THIS IS NOT AFFECTING THE PHYLO MODEL WE ARE USING THE PAPER", sep=""))
    }


    ##check correlation between proportion and non-proportion phylo APPROACHES
    #look for a line whit the warning
    warning_proporition_cor=system(paste("\\
        grep \\
            'WARNING! PROPORTION AND NON-PROPORTION ARE NOT SIMILAR AS EXPECTED. CORRELATION IS' \\
            ", path_output_file, sep=""), intern=TRUE)

    #if we indeed have the line
    if(!is.na(warning_proporition_cor[1]) & warning_proporition_cor[1]<0.9){

        #stop and check the correlation
        stop(paste("CHECK THE CORRELATION BETWEEN PROPORTION AND NON-PROPORTION PHYLO APPROACHES FOR ", species, ". IT IS TOO LOW: ", strsplit(warning_proporition_cor, split=" ")[[1]][13], sep=""))

    }

    #print the finish
    return("FINISH")
}

#run on one species
#check_outputs_species("banksiana")

#run across species
output_results=sapply(unique_species, check_outputs_species)

#check we have run all species
if(sum(output_results=="FINISH"))!=length(unique_species)){
    stop("ERROR! FALSE! WE HAVE A PROBLEM CHECKING THE OUTPUT FILE OF EACH SPECIES")
}




##################################
###### CHECK OUTPUTS BATCHES #####
##################################

##get all batches names
#list batches
list_paths_batches=list.files(path="./code/phylo/recipes/01_global_test_phylo_ubuntu_20_04_v1_slurm_files/", pattern="slurm_global_test_batch_*", full.names=TRUE)
    #path: a character vector of full path names
    #pattern: an optional regular expression.
    #full.names: a logical value.  If 'TRUE', the directory path is prepended to the file names to give a relative file path.

#get the names
#x=list_paths_batches[1]
batch_names_check_output=sapply(list_paths_batches, {function(x) paste("batch_", strsplit(x, split="\\.|_")[[1]][15], sep="")})
    #for each path, split using "." (scape as "\\.") and "_". Then get the 15th element which is the bath number. Paste that number with "batch_" to have the batch number
if(FALSE %in% (names(batch_names_check_output)==list_paths_batches)){
    stop("ERROR! FALSE! PROBLEM CHECKING OUTPUT BATCHES")
}
    #the names of the vector should be the paths
names(batch_names_check_output)=NULL
    #remove the names


##por aquii

##define function
#batch=batch_names_check_output[1]
check_outputs_species=function(batch){

    path_batch_output=paste("./code/phylo/global_test_phylo_current_v1_", batch, ".Rout", sep="")

    n_workers=as.numeric(system(paste(" \\
        grep \\
            'starting worker' \\
            --only-matching \\
            ", path_batch_output, " |
        awk \\
            'END{print NR}'", sep=""), intern=TRUE))
        #we are going to use this numbre to calculate the accepted number of "false" we can have in the output file, as we have 1 false for each new worker started plus 1 when the general session is opened 

    n_error_false=as.numeric(system(paste(" \\
        grep \\
            'error|false' \\
            --extended-regexp \\
            --ignore-case \\
            --only-matching \\
            ", path_batch_output, " | \\
        awk 'END{print NR}'" ,sep=""), intern=TRUE))
        #we can have falses when running each independent worker, so we can only use "error" as signal of problems


    if(n_error_false!=(n_workers+1)){
        stop(paste("ERROR! FALSE! WE HAVE AN ERROR IN THE OUTPUT OF ", batch, sep=""))
    }



    check_finish=as.numeric(system(paste(" \\
        grep \\
            '## FINISH ##' \\
            --count \\
            ", path_batch_output, sep=""), intern=TRUE))

    if(check_finish!=1){
        stop(paste("ERROR! FALSE! THERE IS NO FINISH IN THE SCRIPT OF ", batch, sep=""))
    }


}







#########################################################
###### PROCESS FILES WITH THE NUMBER OF DATA.POINTS #####
#########################################################

##create complete table
#list the tables of each batch
list_paths_n_table=list.files(path="./results/global_test_phylo_current/n_points_before_resampling/", pattern="n_points_before_resampling_batch_*", full.names=TRUE)
    #path: a character vector of full path names
    #pattern: an optional regular expression.
    #full.names: a logical value.  If 'TRUE', the directory path is prepended to the file names to give a relative file path.

#get the name 
batch_names=sapply(strsplit(list_paths_n_table, split="_|\\."), {function(x) paste(x[12], "_", x[13], sep="")})
    #for each path
        #split the path using "_" and "."
            #we have to use "\\" to use "."
        #for the resulting split, get the 13th and 14th element and past them together, so you can have the batch name

list_tables=lapply(list_paths_n_table, read.table, sep="\t", header=TRUE)

names(list_tables)=batch_names

require(dplyr)


n_points=bind_rows(list_tables, .id = "batch_name")


n_points$percent_points_inside=(n_points$points_inside/n_points$total_points)*100

n_points

if(length(which(n_points$percent_points_inside>10))>6){
    stop("ERROR! FALSE! WE HAVE MORE THAN 5 SPECIES WITH MORE THAN 10% OF THE NATURALIZED OCCURRENCES INSIDE THE PA BUFFER")
}


#compare median boyce vs number of occurrences
    #you should get the median boyce across partitions just wihtout phylo and check correlation with number of occurrences

##SEND AGAIN NECESARRY FILES



##halepensis has 97 otuptus of modEvA::Boyce, it should be 12 partitions * 3 algorithms * 2 dup vs non-dup
    #CHECK THIS