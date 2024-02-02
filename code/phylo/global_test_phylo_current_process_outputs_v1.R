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
                ./datos/phylo/method_validation/doi_10_5061_dryad_1hr1n52__v20181213.zip \\
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
            --extended-regexp \\
            'WARNING! PROPORTION AND NON-PROPORTION ARE NOT SIMILAR AS EXPECTED. CORRELATION IS .*' \\
            ", path_output_file, sep=""), intern=TRUE)

    #get the correlation value
    cor_value=strsplit(warning_proporition_cor, split=" ")[[1]][12]

    #stop only if we get a correlation value below 0.9
    if(cor_value!="NA" & cor_value<0.9){

        #stop and check the correlation
        stop(paste("CHECK THE CORRELATION BETWEEN PROPORTION AND NON-PROPORTION PHYLO APPROACHES FOR ", species, ". IT IS TOO LOW: ", , sep=""))
    }

        ###por aquiii
            #problem here, what if other species besides elliotti has NA? we are not counting them...

    ##check we have run the boyce function the expected number of times
    #count number of times we have output for this function
    n_boyce=as.numeric(system(paste(" \\
        grep \\
            'Data include .* presences and .* absences\\.' \\
            --extended-regexp \\
            --count \\
            ", path_output_file, sep=""), intern=TRUE))
        #use regular expression to look for lines having Data include... and any character for the number of presences and absences except a new line character (.*). The string ends with an actual dot, so we have to scape "."
            #https://stackoverflow.com/a/2912904

    #check the number of outputs is the correct
    if(n_boyce!=((((3*2)+2)*12)+(3*12*8))){
        stop("ERROR! FALSE! WE DO NOT HAVE THE EXPECTED NUMBER OF OUTPUTS FROM BOYCE FUNCTION")
    }
        #in non-phylo: ((3*2)+2)*12
            #we run boyce for 3 algorithms twice in each case, i.e., removing or not duplicates. We also run two times the internal function used by modeva::boyce to calculate the input presence/absence file
            #this is repeated across the 12 partitions
        #in phylo: 3*12*8
            #we run boyce for 3 algorithms across 12 partitions and 8 phylo-approaches


    ##print the finish
    return("FINISH")
}

#run on one species
#check_outputs_species("banksiana")

#run across species
output_results=sapply(unique_species, check_outputs_species)

#check we have run all species
if(sum(output_results=="FINISH")!=length(unique_species)){
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
#the names of the vector should be the paths
if(FALSE %in% (names(batch_names_check_output)==list_paths_batches)){
    stop("ERROR! FALSE! PROBLEM CHECKING OUTPUT BATCHES")
}
#remove the names
names(batch_names_check_output)=NULL


##define function
#batch=batch_names_check_output[1]
check_outputs_batch=function(batch){

    #select the output file of the selected batch
    path_batch_output=paste("./code/phylo/global_test_phylo_current_v1_", batch, ".Rout", sep="")

    #calculate the number of workers indicated in the output file
    n_workers=as.numeric(system(paste(" \\
        grep \\
            'starting worker' \\
            --only-matching \\
            ", path_batch_output, " |
        awk \\
            'END{print NR}'", sep=""), intern=TRUE))
        #--only-matching: show only nonempty parts of lines that match
        #count the number of times 'starting worker' appears
        #we are going to use this number to calculate the accepted number of "false" we can have in the output file, as we have 1 false for each new worker started plus 1 when the general session is opened 

    #count number of error OR falses
    n_error_false=as.numeric(system(paste(" \\
        grep \\
            'error|false' \\
            --extended-regexp \\
            --ignore-case \\
            --only-matching \\
            ", path_batch_output, " | \\
        awk 'END{print NR}'" ,sep=""), intern=TRUE))
        #look for "error" OR "false" using regular expression and ignoring case (e.g., error or ERROR is considered). Then show only nonempty parts of lines that match (--only-matching) and count the number of these cases.

    #the number of error|false should be equal to the number of workers we have plus 1, i.e., we should have a FALSE after starting each R session, the general session and the parallel workers
    if(n_error_false!=(n_workers+1)){
        stop(paste("ERROR! FALSE! WE HAVE MORE THAN EXPECTED ERROR|FALSE IN THE OUTPUT OF ", batch, sep=""))
    }

    #FINISH should appear 1 time, indicating the end of the script
    check_finish=as.numeric(system(paste(" \\
        grep \\
            '## FINISH ##' \\
            --count \\
            ", path_batch_output, sep=""), intern=TRUE))
        #--count: print only a count of selected lines
    if(check_finish!=1){
        stop(paste("ERROR! FALSE! THERE IS NO FINISH IN THE SCRIPT OF ", batch, sep=""))
    }

    #end this
    return("FINISH")
}

#run on one batch
#check_outputs_batch("batch_3")

#run across species
output_batch_results=sapply(batch_names_check_output, check_outputs_batch)

#check we have run all species
if(sum(output_batch_results=="FINISH")!=length(batch_names_check_output)){
    stop("ERROR! FALSE! WE HAVE A PROBLEM CHECKING THE OUTPUT FILE OF EACH BATCH")
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

#get the name of each table, i.e., the batch number
#x=strsplit(list_paths_n_table[1], split="_|\\.")[[1]]
batch_names=sapply(strsplit(list_paths_n_table, split="_|\\."), {function(x) paste(x[12], "_", x[13], sep="")})
    #split each path using "_" and "."
            #we have to use "\\" to use "."
    #for each element in the resulting list, i.e., a vector with the split elements
        #get the 13th and 14th element and past them together, so you can have the batch name

#get a list with all the tables
list_tables=lapply(list_paths_n_table, read.table, sep="\t", header=TRUE)
if(length(list_tables)!=length(list_paths_n_table)){
    stop("ERROR! FALSE! WE HAVE A PROBLEM CREATING THE TABLE WITH NUMBER OF POINTS")
}

#add the names
names(list_tables)=batch_names

#bind all tables together as they have the same columns
require(dplyr)
n_points=bind_rows(list_tables, .id = "batch_name")
    #Bind any number of data frames by row, making a longer result. This is similar to 'do.call(rbind, dfs)', but the output will contain all columns that appear in any of the inputs

#calculate the percentage of points inside
n_points$percent_points_inside=(n_points$points_inside/n_points$total_points)*100

#stop if we have many species with a percentage above 10%
if(length(which(n_points$percent_points_inside>10))>6){
    stop("ERROR! FALSE! WE HAVE MORE THAN 6 SPECIES WITH MORE THAN 10% OF THE NATURALIZED OCCURRENCES INSIDE THE PA BUFFER")
}


##extract boyce non-phylo per species
#species=unique(naturalized_occurrences$species)[1]
get_boyce=function(species){

    #check whether file exists
    boyce_exists=file.exists(paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/", species, "_boyce_table.tsv.gz", sep=""))
    
    #if boyce does exists
    if(boyce_exists){

        #read the boyce table
        boyce_table=read.table(
            paste("./results/global_test_phylo_current/predict_eval_no_phylo/", species, "/boyce_index/", species, "_boyce_table.tsv.gz", sep=""),
            sep="\t",
            header=TRUE)

        #open empty DF to save results
        boyce_rows=data.frame(species=NA, algorithm=NA, percentile_2.5=NA, percentile_50=NA, percentile_97.5=NA)

        #for each algorithm
        #algorithm="glm"
        for(algorithm in c("glm", "gam", "rf")){

            #extract the percentiles of the selected model and then transponse and convert to DF
            boyce_row=data.frame(t(boyce_table[which(boyce_table$partition %in% c("percentile_2.5", "percentile_50", "percentile_97.5")), which(colnames(boyce_table)==paste(algorithm, "_eval", sep=""))]))

            #set the column names using the percentile names
            colnames(boyce_row)=boyce_table[which(boyce_table$partition %in% c("percentile_2.5", "percentile_50", "percentile_97.5")), which(colnames(boyce_table)=="partition")]
                #we get the names of the percentile from the partition column to follow the same order we have in boyce_row

            #bind the percentile to the algorithm and species names
            boyce_row=cbind.data.frame(species=species, algorithm=algorithm, boyce_row)

            #bind to the previous data.frame
            boyce_rows=rbind.data.frame(boyce_rows, boyce_row) 
        }

        #remove the first row with all NANs
        boyce_rows=boyce_rows[which(apply(is.na(boyce_rows), 1, sum)!=ncol(boyce_rows)),]
            #select only those rows for which the sum of is.na() is NOT equal to the total number of columns, i.e., rows without all NA

        #update row names
        rownames(boyce_rows)=1:nrow(boyce_rows)

        #return the table
        return(boyce_rows)
    }
}

#apply the function to one species
#get_boyce("halepensis")

#apply the function across all species
list_results_boyce=lapply(unique(naturalized_occurrences$species), get_boyce)
    #get a list with the cleaned boyce table per species

#check we have the correct rows and columns in each table
#x=list_results_boyce[[1]]
check_rows_boyce=FALSE %in% (sapply(list_results_boyce, {function(x) if(!is.null(rownames(x))){identical(rownames(x), c("1","2","3"))}}))
check_cols_boyce=FALSE %in% (sapply(list_results_boyce, {function(x) if(!is.null(colnames(x))){identical(colnames(x), c("species", "algorithm", "percentile_2.5", "percentile_50", "percentile_97.5"))}}))
if(check_rows_boyce | check_cols_boyce){
    stop("ERROR! FALSE! WE HAVE A PROBLEM WITH THE EXTRACTION OF BOYCE")
}

#bind all tables as the have the same columns
results_boyce=bind_rows(list_results_boyce, .id=NULL)


##plot boyce against the number of occurrences after resampling
#algorithms
algorithms=c("glm", "gam", "rf")

#open empty plot
pdf(paste("./results/global_test_phylo_current/n_points_before_resampling/boyce_vs_n_occurrences.pdf", sep=""), width=12, height=8)
plot(1, type="n", xlab="", ylab="", xlim=c(0, 100), ylim=c(-1, 1))
    #the x axis is the number of points, we do not have species with more than 70-80 occurrences
    #the y axis is the boyce index, which can go from -1 to 1

#for each algorithm
#algorithm_index=2
for(algorithm_index in 1:length(algorithms)){

    #select the algorithm
    selected_algorithm=algorithms[algorithm_index]

    #select rows of the selected algorithm
    results_boyce_subset=results_boyce[which(results_boyce$algorithm==selected_algorithm), ]
    
    #merge with the number of occurrences
    results_boyce_subset=merge(results_boyce_subset, n_points, by="species")

    #create a sequence from 1 to the number of species for which we have data
    species_shape=seq(1,length(results_boyce_subset$species),1)
        #this will be used to get a different shape for each specie

    #plot points
    points(x=results_boyce_subset$points_outside_after_resampling, y=results_boyce_subset$percentile_50, col=algorithm_index, pch=species_shape)
        #occurrences after resampling against the median boyce
        #color 1,2,3... based on the index of the algorithm, i.e., 1 for all points of glm, 2 for all points of gam, and 3 for RF
        #the shape will be different across species, as we are plotting here all species for a given model, we need a vector with all the shapes for these species

    #add 95CI as error bars
    arrows(
        x0=results_boyce_subset$points_outside_after_resampling, 
        y0=results_boyce_subset$percentile_2.5, 
        x1=results_boyce_subset$points_outside_after_resampling, 
        y1=results_boyce_subset$percentile_97.5, 
        code=3, 
        angle=90, 
        length=0.02,
        lwd=0.4,
        col=algorithm_index)
        #the bar is in the same X position, i.e., the number of occurrences
        #the bar moves across the Y axis, from the percentile 2.5 to percentile 97.5
        #select the type of arrow and the dimensions
        #add the color of the algorithm
}

#add the legend
legend(x="topright", legend=algorithms, fill=1:length(algorithms))
    #using the name of the models and the same color for each algorithm as used in the points, i.e., from 1 to the total number of algorithms.
dev.off()
    #CHECK THE PLOT:
        #the number of occurrences is biasing our results? If so, we should see a positive correlation between boyce and the number of occurrences. If we see that boyce decreases with the number of occurrences, then it is possible that we do not have ability to evaluate species with a low number of occurrences?




###########################################
###### PLOT BOYCE PHYLO VS. NON-PHYLO #####
###########################################

##define function to extract the phylo-boyce
#species=unique(naturalized_occurrences$species)[1]
get_boyce_phylo=function(species){

    #check whether file exists
    boyce_exists=file.exists(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/", species, "_boyce_table_non_phylo_vs_phylo.tsv.gz", sep=""))
    
    #if boyce does exists
    if(boyce_exists){

        #read the boyce table
        boyce_table=read.table(
            paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/", species, "_boyce_table_non_phylo_vs_phylo.tsv.gz", sep=""),
            sep="\t",
            header=TRUE)

        #open empty DF to save results
        boyce_rows=data.frame(species=NA, algorithm=NA, no_phylo_percentile_2.5=NA, no_phylo_percentile_50=NA, no_phylo_percentile_97.5=NA, phylo_percentile_2.5=NA, phylo_percentile_50=NA, phylo_percentile_97.5=NA)

        #for each algorithm
        #algorithm="glm"
        for(algorithm in c("glm", "gam", "rf")){

            ##extract the percentiles of the selected model without phylo
            #get the row
            boyce_row_no_phylo=boyce_table[
                which(boyce_table$model==paste(algorithm, "_eval", sep="")), 
                which(colnames(boyce_table) %in% c("model", "percentile_2.5", "percentile_50", "percentile_97.5"))]

            #check we have selected the correct model
            if(boyce_row_no_phylo$model!=paste(algorithm, "_eval", sep="")){
                stop(paste("ERROR! FALSE! PROBLEM EXTRACTING BOYCE PHYLO AND NON-PHYLO FOR ", species, sep=""))
            } else { #if so, remove the model column
                boyce_row_no_phylo$model=NULL
            }

            #add no-phylo to the column names
            colnames(boyce_row_no_phylo)=paste("no_phylo_", colnames(boyce_row_no_phylo), sep="")


            ##extract the percentiles of the selected model with phylo
            #get the row
            boyce_row_phylo=boyce_table[
                which(boyce_table$model==paste("phylo_rasters_subset_inter_", algorithm, "_boyce", sep="")), 
                which(colnames(boyce_table) %in% c("model", "percentile_2.5", "percentile_50", "percentile_97.5"))]

            #check we have selected the correct model
            if(boyce_row_phylo$model!=paste("phylo_rasters_subset_inter_", algorithm, "_boyce", sep="")){
                stop(paste("ERROR! FALSE! PROBLEM EXTRACTING BOYCE PHYLO AND NON-PHYLO FOR ", species, sep=""))
            } else { #if so, remove the model column
                boyce_row_phylo$model=NULL
            }

            #add phylo to the column names
            colnames(boyce_row_phylo)=paste("phylo_", colnames(boyce_row_phylo), sep="")

            #bind the phylo and no-phylo percentile to the algorithm and species names
            boyce_row=cbind.data.frame(species=species, algorithm=algorithm, boyce_row_no_phylo, boyce_row_phylo)

            #bind to the previous data.frame
            boyce_rows=rbind.data.frame(boyce_rows, boyce_row) 
        }

        #remove the first row with all NANs
        boyce_rows=boyce_rows[which(apply(is.na(boyce_rows), 1, sum)!=ncol(boyce_rows)),]
            #select only those rows for which the sum of is.na() is NOT equal to the total number of columns, i.e., rows without all NA

        #update row names
        rownames(boyce_rows)=1:nrow(boyce_rows)

        #return the table
        return(boyce_rows)
    }
}

#apply the function to one species
#get_boyce_phylo("banksiana")

#apply the function across all species
list_results_boyce_phylo=lapply(unique(naturalized_occurrences$species), get_boyce_phylo)
    #get a list with the cleaned boyce table per species

#check we have the correct rows and columns in each table
#x=list_results_boyce_phylo[[1]]
check_rows_boyce_phylo=FALSE %in% (sapply(list_results_boyce_phylo, {function(x) if(!is.null(rownames(x))){identical(rownames(x), c("1","2","3"))}}))
check_cols_boyce_phylo=FALSE %in% (sapply(list_results_boyce_phylo, {function(x) if(!is.null(colnames(x))){identical(colnames(x), c("species", "algorithm", "no_phylo_percentile_2.5", "no_phylo_percentile_50", "no_phylo_percentile_97.5", "phylo_percentile_2.5", "phylo_percentile_50", "phylo_percentile_97.5"))}}))
if(check_rows_boyce_phylo | check_cols_boyce_phylo){
    stop("ERROR! FALSE! WE HAVE A PROBLEM WITH THE EXTRACTION OF PHYLO AND NON-PHYLO BOYCE")
}

#bind all tables as the have the same columns
results_boyce_phylo=bind_rows(list_results_boyce_phylo, .id=NULL)

#check we have the same non-phylo results than in the previous boyce table
#get the percentile columns from both tables
subset_boyce_no_phylo_1=results_boyce[,which(colnames(results_boyce) %in% c("percentile_2.5", "percentile_50", "percentile_97.5"))]
subset_boyce_no_phylo_2=results_boyce_phylo[,which(grepl("no_phylo", colnames(results_boyce_phylo), fixed=TRUE))]
#compare each column between the previous boyce table and the non-phylo data from the new table
#i=1
for(i in 1:ncol(subset_boyce_no_phylo_1)){

    #if the difference of any algorithm*species combination between both tables is higher than <1e-15, we have a problem
    if(FALSE %in% (subset_boyce_no_phylo_1[,i]-subset_boyce_no_phylo_2[,i])<1e-15){
        stop("ERROR! FALSE! PROBLEM EXTRACTING BOYCE RESULTS")
    }
}



#now plot phylo vs non-phylo





##SEND AGAIN NECESARRY FILES

