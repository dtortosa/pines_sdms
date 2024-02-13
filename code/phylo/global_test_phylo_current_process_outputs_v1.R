#!/usr/bin/env Rscript

#we use the shebang "#!/usr/bin/env Rscript" to facilitate portability between machines. "env" offers more portability because it searches for the first occurrence of the R interpreter in the environment where we are running. It should thus look for R version installed in the container
    #https://www.r-bloggers.com/2019/11/r-scripts-as-command-line-tools/
    #https://www.baeldung.com/linux/bash-shebang-lines

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
    #https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
    #https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file

#Note about terminal
    #you cannot run more than 4096 bytes of characters from the script to the terminal. It is not a problem of sublime, but a limit imposed in R. This is hard-coded and I should have to modify the C-code an re-compile to solve this.
        #https://stackoverflow.com/questions/13216480/paste-character-limit




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

    #check whether the output file exists
    output_file_exists=file.exists(path_output_file)

    #if the file exits do stuff
    if(output_file_exists){

        ##get lines of start and end eval phylo
        #get the line number where eval_phylo starts and then the next one, so we can know whether phylo_eval has been run or not
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

            #the next line after STARTING PHYLO should be the ENDING because no analyses should be done when the ERROR previously indicated appears, if not, we have a problem
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

        #get the correlation value for each case
        #x=warning_proporition_cor[[1]]
        cor_values=sapply(warning_proporition_cor, {function(x) strsplit(x, split=" ")[[1]][13]})
        names(cor_values)=NULL

        #if we have at least one value of correlation obtained from the warning
        if(length(cor_values)!=0){

            #low cor if at least 1 is NA
            if("NA" %in% cor_values){
                cor_flag="low_cor"
            } else {

                #if not NA, then calculate the min cor
                cor_value=min(as.numeric(cor_values))

                #if the min cor is lower than 0.9
                if(cor_value<0.9){
                    cor_flag="low_cor"     
                } else {
                    cor_flag="high_cor"
                }
            }
        } else { #we do not have warning

            #if the next line after STARTING phylo is NOT END
            if(!grepl("ENDING predict_eval_phylo FOR", second_line_eval_phylo, fixed=TRUE)){

                #then the phylo analyses have been run and we still do not have warning, so correlation is high
                cor_flag="high_cor"
            } else {

                #if not, then the phylo analyses have not been run
                cor_flag="cor_not_done"
            }
        }


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

        #calculate the expected number of boyce outputs
        expected_boyce=((((3*2)+6)*12)+(3*12*8))
            #in non-phylo: ((3*2)+2)*12
                #we run boyce for 3 algorithms twice in each case, i.e., removing or not duplicates. We also run 6 times the internal function used by modeva::boyce to calculate the input presence/absence file
                #this is repeated across the 12 partitions
            #in phylo: 3*12*8
                #we run boyce for 3 algorithms across 12 partitions and 8 phylo-approaches

        #get the line number where eval_phylo starts
        first_line_no_eval_phylo=as.numeric(system(paste(" \\
            grep \\
                --line-number \\
                'STARTING predict_eval_no_phylo FOR  ", species, "' \\
                ", path_output_file, "| \\
            cut \\
                --delimiter : \\
                --fields 1", sep=""), intern=TRUE))
            #get line number with grep and, from the output, get the first element before :, which is the line number
            #https://stackoverflow.com/a/21412398

        #now get the fourth next line
        four_lines_after_no_eval_phylo=system(paste("\\
            awk \\
                '{if(NR==", first_line_no_eval_phylo+4, "){print $0}}' \\
                ", path_output_file, sep=""), intern=TRUE)

        #if the fourth line after non phylo eval is the END
        if(grepl(paste("ENDING  ", species, sep=""), four_lines_after_no_eval_phylo, fixed=TRUE)){

            #means the species has not run and we should have zero boyce outputs
            check_boyce=n_boyce==0
        } else {

            #if not, the specie has run and we should have the expected number of boyce outputs
            check_boyce=n_boyce==expected_boyce
        }

        #check the number of outputs is the correct
        if(!check_boyce){
            stop(paste("ERROR! FALSE! WE DO NOT HAVE THE EXPECTED NUMBER OF OUTPUTS FROM BOYCE FUNCTION ", species, sep=""))
        }


        ##return the cor flag and FINISH
        return(list(cor_flag, "FINISH"))
    }
}

#run on one species
#check_outputs_species("banksiana")

#run across species
output_results=sapply(unique_species, check_outputs_species)

#check we have run all species and the number of low_correlation between proportion and non-proportion is low
#x=output_results[,1]
    #each species is a column
output_results_count_finish=sum(apply(output_results, 2, {function(x) if(!is.null(x[2]) & x[2]=="FINISH") 1 else 0}))
output_results_count_low_cor_cases=sum(apply(output_results, 2, {function(x) if(!is.null(x[1]) & x[1]=="low_cor") 1 else 0}))
output_results_count_low_high_cor_cases=sum(apply(output_results, 2, {function(x) if(!is.null(x[1]) & x[1]%in%c("low_cor", "high_cor")) 1 else 0}))
    #one liner ifelse
        #https://stackoverflow.com/questions/15586566/if-statement-in-r-can-only-have-one-line
if(output_results_count_finish!=length(unique_species)){
    stop("ERROR! FALSE! WE HAVE A PROBLEM CHECKING THE OUTPUT FILE OF EACH SPECIES")
}

#count the number of cases with low correlation between the main phylogenetic approaches
percent_low_cor=(output_results_count_low_cor_cases/output_results_count_low_high_cor_cases)*100
if(percent_low_cor>50){
    stop("ERROR! FALSE! WE HAVE MORE THAN 50% OF SPECIES WITH A LOW CORRELATION BETWEEN PHYLO PROPORTION AND NON-PROPORTION")
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

#stop if we have many species with a percentage above 30%
if(length(which(n_points$percent_points_inside>30))>6){
    stop("ERROR! FALSE! WE HAVE MORE THAN 6 SPECIES WITH MORE THAN 30% OF THE NATURALIZED OCCURRENCES INSIDE THE PA BUFFER")
}

#save table
write.table(n_points, "./results/global_test_phylo_current/n_points_before_resampling/n_points_total.tsv", sep="\t", col.names=TRUE, row.names=FALSE)


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
plot(1, type="n", xlab="Number of occurrences after resampling", ylab="Boyce index", xaxt="n", yaxt="n", 
    xlim=c(0, max(n_points$points_outside_after_resampling)+10), 
    ylim=c(-1, 1))
axis(1, at=seq(0, max(n_points$points_outside_after_resampling)+10, 5))
axis(2, at=seq(-1, 1, 0.2))
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
    species_shape=seq(1,length(unique(results_boyce_subset$species)),1)
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
legend(x="bottomright", legend=unique(results_boyce$species), pch=seq(1,length(unique(results_boyce$species)),1))
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
                which(boyce_table$model==paste("phylo_rasters_proportion_subset_inter_", algorithm, "_boyce", sep="")), 
                which(colnames(boyce_table) %in% c("model", "percentile_2.5", "percentile_50", "percentile_97.5"))]

            #check we have selected the correct model
            if(boyce_row_phylo$model!=paste("phylo_rasters_proportion_subset_inter_", algorithm, "_boyce", sep="")){
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


##plot boyce against the number of occurrences after resampling
#algorithms
algorithms=c("glm", "gam", "rf")

#open empty plot
pdf(paste("./results/global_test_phylo_current/boyce_non_phylo_vs_phylo.pdf", sep=""), width=8, height=8)
    #we need the same height and width in order to properly compare the confidence intervals of X and Y. If the shape of the plot is not cuadrangular, we cannot visualize the difference of the CI between the two variable plots.
plot(1, type="n", xlab="Boyce index non-phylo", ylab="Boyce index phylo", xaxt="n", yaxt="n", 
    xlim=c(-1, 1), 
    ylim=c(-1, 1))
axis(1, at=seq(-1, 1, 0.2))
axis(2, at=seq(-1, 1, 0.2))
    #the x and y axis are the boyce index, which can go from -1 to 1

#for each algorithm
#algorithm_index=1
for(algorithm_index in 1:length(algorithms)){

    #select the algorithm
    selected_algorithm=algorithms[algorithm_index]

    #select rows of the selected algorithm
    results_boyce_phylo_subset=results_boyce_phylo[which(results_boyce_phylo$algorithm==selected_algorithm), ]

    #create a sequence from 1 to the number of species for which we have data
    species_shape=seq(1,length(results_boyce_phylo_subset$species),1)
        #this will be used to get a different shape for each specie

    #plot points
    points(x=results_boyce_phylo_subset$no_phylo_percentile_50, y=results_boyce_phylo_subset$phylo_percentile_50, col=algorithm_index, pch=species_shape)
        #occurrences after resampling against the median boyce
        #color 1,2,3... based on the index of the algorithm, i.e., 1 for all points of glm, 2 for all points of gam, and 3 for RF
        #the shape will be different across species, as we are plotting here all species for a given model, we need a vector with all the shapes for these species

    #add 95CI as error bars
    #Y error bars, i.e, 95CI of the phylo boyce
    arrows(
        x0=results_boyce_phylo_subset$no_phylo_percentile_50, 
        y0=results_boyce_phylo_subset$phylo_percentile_2.5, 
        x1=results_boyce_phylo_subset$no_phylo_percentile_50, 
        y1=results_boyce_phylo_subset$phylo_percentile_97.5, 
        code=3, 
        angle=90, 
        length=0.02,
        lwd=0.4,
        col=algorithm_index)
        #the bar is in the same X position, i.e., the non-phylo boyce index
        #the bar moves across the Y axis, i.e., the phylo boyce index, from the percentile 2.5 to percentile 97.5
        #select the type of arrow and the dimensions
        #add the color of the algorithm
    #X error bars, i.e, 95CI of the non-phylo boyce
    arrows(
        x0=results_boyce_phylo_subset$no_phylo_percentile_2.5, 
        y0=results_boyce_phylo_subset$phylo_percentile_50, 
        x1=results_boyce_phylo_subset$no_phylo_percentile_97.5, 
        y1=results_boyce_phylo_subset$phylo_percentile_50, 
        code=3, 
        angle=90, 
        length=0.02,
        lwd=0.4,
        col=algorithm_index)
        #the bar is in the same Y position, i.e., the phylo boyce index
        #the bar moves across the X axis, i.e., the non-phylo boyce index, from the percentile 2.5 to percentile 97.5
}

#add the legend
legend(x="topleft", legend=algorithms, fill=1:length(algorithms))
    #using the name of the models and the same color for each algorithm as used in the points, i.e., from 1 to the total number of algorithms.
legend(x="bottomright", legend=unique(results_boyce_phylo$species), pch=seq(1,length(unique(results_boyce_phylo$species)),1))

#add diagonal line
abline(coef = c(0,1))
dev.off()
    #CHECK THE PLOT:
        #if they are the same, the points will follow the diagonal. Also the 95CI bars will be the same between phylo and non-phylo
        #if phylo (Y) is lower than non-phylo (X), points will be below diagonal
        #if phylo (Y) is higher than non-phylo (X), points will be above diagonal




###################################
###### TABLE PHYLO DIFF PAPER #####
###################################

##define function to extract the phylo-boyce
#species=unique(naturalized_occurrences$species)[1]
get_boyce_phylo_diff=function(species){

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
        boyce_columns=data.frame(species=species)

        #for each algorithm
        #algorithm="glm"
        for(algorithm in c("glm", "gam", "rf")){

            ##extract the percentiles of the selected model
            #get the row
            boyce_column_phylo_diff=boyce_table[
                which(boyce_table$model==paste("phylo_rasters_proportion_subset_inter_", algorithm, "_diff", sep="")), 
                which(colnames(boyce_table) %in% c("model", "percentile_2.5", "percentile_50", "percentile_97.5"))]

            #check we have selected the correct model
            if(boyce_column_phylo_diff$model!=paste("phylo_rasters_proportion_subset_inter_", algorithm, "_diff", sep="")){
                stop(paste("ERROR! FALSE! PROBLEM EXTRACTING BOYCE PHYLO AND NON-PHYLO FOR ", species, sep=""))
            } else { #if so, remove the model column
                boyce_column_phylo_diff$model=NULL
            }

            #add phylo to the column names
            colnames(boyce_column_phylo_diff)=paste(algorithm, "_phylo_diff_", colnames(boyce_column_phylo_diff), sep="")

            #bind the phylo and no-phylo percentile to the algorithm and species names
            boyce_column_phylo_diff_species=cbind.data.frame(species=species, boyce_column_phylo_diff)

            #bind to the previous data.frame
            boyce_columns=merge(boyce_columns, boyce_column_phylo_diff_species) 
        }

        #check we have the correct column names
        #get all possible combinations of percentiles and algortihms as a DF
        check_columns_raw=expand.grid(c("phylo_diff_percentile_2.5", "phylo_diff_percentile_50", "phylo_diff_percentile_97.5"), c("glm", "gam", "rf"))
            #the first vector is used to sort the DF, so we need to use percentile first, to have 2.5,5,97.5...
        #paste algorithm and percentiles        
        check_columns=paste(check_columns_raw$Var2, check_columns_raw$Var1, sep="_")
        #add species
        check_columns=c("species", check_columns)
        #check identical
        if(!identical(colnames(boyce_columns), check_columns)){
            stop(paste("ERROR! FALSE! WE HAVE A PROBLEM WITH THE EXTRACTION OF BOYCE PHYLO DIFF FOR SPECIES ", species, sep=""))
        }

        #return the table
        return(boyce_columns)
    }
}

#apply the function to one species
#get_boyce_phylo_diff("halepensis")

#apply the function across all species
list_results_boyce_phylo_diff=lapply(unique(naturalized_occurrences$species), get_boyce_phylo_diff)
    #get a list with the cleaned boyce table per species

#check we have the correct rows and columns in each table
#x=list_results_boyce_phylo_diff[[1]]
#x=list_results_boyce_phylo_diff[[17]]
check_rows_boyce_phylo_diff=FALSE %in% (sapply(list_results_boyce_phylo_diff, {function(x) if(!is.null(x)){nrow(x)==1}}))
check_cols_boyce_phylo_diff=FALSE %in% (sapply(list_results_boyce_phylo_diff, {function(x) if(!is.null(colnames(x))){identical(colnames(x), c("species", "glm_phylo_diff_percentile_2.5", "glm_phylo_diff_percentile_50", "glm_phylo_diff_percentile_97.5", "gam_phylo_diff_percentile_2.5", "gam_phylo_diff_percentile_50", "gam_phylo_diff_percentile_97.5", "rf_phylo_diff_percentile_2.5", "rf_phylo_diff_percentile_50", "rf_phylo_diff_percentile_97.5"))}}))
if(check_rows_boyce_phylo_diff | check_cols_boyce_phylo_diff){
    stop("ERROR! FALSE! WE HAVE A PROBLEM WITH THE EXTRACTION OF PHYLO AND NON-PHYLO BOYCE")
}

#bind all tables as the have the same columns
results_boyce_phylo_diff=bind_rows(list_results_boyce_phylo_diff, .id=NULL)

#add columns with the median(95CI) for each algorithm
for(algorithm in c("glm", "gam", "rf")){

    #paste the rounded median and 95CI of the selected algorithm
    new_column=paste(
        signif(results_boyce_phylo_diff[,
            which(colnames(results_boyce_phylo_diff)==paste(algorithm, "_phylo_diff_percentile_50", sep=""))], 3), "(", 
        signif(results_boyce_phylo_diff[,
            which(colnames(results_boyce_phylo_diff)==paste(algorithm, "_phylo_diff_percentile_2.5", sep=""))], 3), ",", 
        signif(results_boyce_phylo_diff[,
            which(colnames(results_boyce_phylo_diff)==paste(algorithm, "_phylo_diff_percentile_97.5", sep=""))], 3), ")", sep="")
        #we use signif to maintain the scientific notation while reducing the decimals
            #https://stackoverflow.com/a/45106054

    #add the new column to the DF
    results_boyce_phylo_diff[,paste(algorithm, "_phylo_diff_median_95CI", sep="")]=new_column
}

#check the creation of the new variable went well
#algorithm="glm"
for(algorithm in c("glm", "gam", "rf")){

    #extract the column indexes of each percentile of the selected algorithm
    per_2.5_col_index=which(colnames(results_boyce_phylo_diff)==paste(algorithm, "_phylo_diff_percentile_2.5",sep=""))
    median_col_index=which(colnames(results_boyce_phylo_diff)==paste(algorithm, "_phylo_diff_percentile_50",sep=""))
    per_97.5_col_index=which(colnames(results_boyce_phylo_diff)==paste(algorithm, "_phylo_diff_percentile_97.5",sep=""))

    #for each row of the DF apply function
    #x=results_boyce_phylo_diff[1,]
    phylo_diff_median_approach_2=apply(results_boyce_phylo_diff, 1, {function(x) paste(signif(as.numeric(x[median_col_index]),3), "(", signif(as.numeric(x[per_2.5_col_index]),3), ",", signif(as.numeric(x[per_97.5_col_index]),3), ")", sep="")})
        #get the percentiles, round maintaining the scientific notation and paste them
        #we do this for all rows

    #check whether the newly created vector is exactly the same than the column with median - 95CI
    if(!identical(phylo_diff_median_approach_2, results_boyce_phylo_diff[,paste(algorithm, "_phylo_diff_median_95CI", sep="")])){
        stop(paste("ERROR! FALSE! PROBLEM CREATING THE TABLE WITH MEDIAN PHYLO FOR ALGORITHM ", algorithm, sep=""))
    }
}

#save the table
write.table(results_boyce_phylo_diff, paste("./results/global_test_phylo_current/results_boyce_non_phylo_vs_phylo_diff.tsv", sep=""), sep="\t", row.names=FALSE)
    #alternative presentation
        #barplot with phylo-diff, 3 bars per species, being positive or negative
        #we will see the great negative impact in strobus.

    #CHECK THE TABLE
        #if the 95CI overlaps with zero, it means that the difference between non-phylo and phylo is positive for some data partitions and negative for others. In other words, we cannot say there are differences in the boyce index after the application of the phylogenetic correction.
        #we need to see 95CIs consistently below 0 for us to say that phylo has a higher boyce index than non-phylo (we do non-phylo vs phylo)
        #silvestris and patula (but patula overlaps with zero)

    #CHECK THE PLOTS ABOUT BOYCE NO PHYLO
        #We need to check whether the P/E ratio vs suitability plots show strange things like in radiata.
        #If one of the partitions is like partition 1 GAM of radiata, then the other ones should have lower boyce indexes reducing the median, if not, take a look in detail.
            #if all partitions show consistently the same situation, and we have a median Boyce very high with just two bins, take a look.
        #also check that the median boyce index correlates well with the ensemble showing predicted suitability outside PA buffer and the naturalized occurrences. Higher boyce should correlated with more matches.
        #check differences between removing or not duplicated presences.




############################
###### GLMM PHYLO DIFF #####
############################

##define function to create table for glmm
#species=unique(naturalized_occurrences$species)[1]
get_boyce_phylo_diff_glmm=function(species){

    #check whether file exists
    boyce_exists=file.exists(paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/", species, "_boyce_table_non_phylo_vs_phylo.tsv.gz", sep=""))
    
    #if boyce does exists
    if(boyce_exists){

        #read the boyce table
        boyce_table=read.table(
            paste("./results/global_test_phylo_current/predict_eval_phylo/", species, "/boyce_index/", species, "_boyce_table_non_phylo_vs_phylo.tsv.gz", sep=""),
            sep="\t",
            header=TRUE)

        #remove percentile columns
        boyce_table=boyce_table[,which(!colnames(boyce_table)%in%c("percentile_2.5", "percentile_50", "percentile_97.5"))]

        #open empty DF to save results
        boyce_glmm_table=data.frame(species=NA, partition=NA, model=NA, non_phylo=NA, phylo=NA, phylo_diff=NA)

        #for each algorithm
        #algorithm="glm"
        for(algorithm in c("glm", "gam", "rf")){

            #across partitions, extract the value of boyce for the selected algorithm without phylo, with phylo and the difference
            non_phylo_algo=boyce_table[
                which(boyce_table$model==paste(algorithm, "_eval", sep="")),
                which(colnames(boyce_table)!="model")]
            phylo_algo=boyce_table[
                which(boyce_table$model==paste("phylo_rasters_proportion_subset_inter_", algorithm, "_boyce", sep="")),
                which(colnames(boyce_table)!="model")]
            phylo_diff_algo=boyce_table[
                which(boyce_table$model==paste("phylo_rasters_proportion_subset_inter_", algorithm, "_diff", sep="")),]
            
            #check we have selected the correct model
            if(phylo_diff_algo$model!=paste("phylo_rasters_proportion_subset_inter_", algorithm, "_diff", sep="")){
                stop(paste("ERROR! FALSE! PROBLEM EXTRACTING BOYCE PHYLO AND NON-PHYLO FOR ", species, sep=""))
            } else { #if so, remove the model column
                phylo_diff_algo$model=NULL
            }

            #check the difference was correctly calculated
            if(TRUE %in% ((non_phylo_algo-phylo_algo)-phylo_diff_algo>1e-15)){
                stop(paste("ERROR! FALSE! WE HAVE A PROBLEM WITH THE CALCULATION OF GLMM PHYLO DIFF TABLE FOR SPECIES ", species, sep=""))
            }

            #convert the phylo diff row to a column
            phylo_diff_column=t(phylo_diff_algo)
            colnames(phylo_diff_column)="phylo_diff"
            rownames(phylo_diff_column)=1:length(phylo_diff_column)

            #do the same with phylo and non-phylo
            non_phylo_column=t(non_phylo_algo)
            phylo_column=t(phylo_algo)
            colnames(non_phylo_column)="non_phylo"
            colnames(phylo_column)="phylo"
            rownames(non_phylo_column)=1:length(non_phylo_column)
            rownames(phylo_column)=1:length(phylo_column)


            #create a vector with the name of the partitions and species
            #x=colnames(phylo_diff_algo)[1]
            partition_names=sapply(colnames(phylo_diff_algo), {function(x) paste(species, "_", strsplit(x,split="X")[[1]][2], sep="")})
            if(FALSE %in% (names(partition_names)==colnames(phylo_diff_algo))){
                stop(paste("ERROR! FALSE! WE HAVE A PROBLEM CALCUALTING THE GLMM TABLE FOR SPECIES ", species, sep=""))
            } else {
                names(partition_names)=NULL
            }

            #combine everything in a DF
            boyce_glmm_table_algo=data.frame(species=species, partition=partition_names, model=algorithm, non_phylo=non_phylo_column, phylo=phylo_column, phylo_diff=phylo_diff_column)

            #bind to the empty data.frame
            boyce_glmm_table=rbind.data.frame(boyce_glmm_table, boyce_glmm_table_algo)
        }

        #remove first row with NAs
        boyce_glmm_table=boyce_glmm_table[which(apply(is.na(boyce_glmm_table),1,sum)!=ncol(boyce_glmm_table)),]

        #update row names
        rownames(boyce_glmm_table)=1:nrow(boyce_glmm_table)

        #make several checks
        check_1=nrow(boyce_glmm_table)!=12*3
            #12 partitions times 3 algorithms
        check_2=FALSE %in% (colnames(boyce_glmm_table)==c("species", "partition", "model", "non_phylo", "phylo", "phylo_diff"))
            #correct column names
        check_3=!identical(boyce_glmm_table$species, rep(species,12*3))
            #species column should be the species name repeated 12*3 times
        check_4=!identical(boyce_glmm_table$partition, rep(paste(species, "_", 1:12, sep=""),3))
            #partition column should be the species name pasted with 1:12 repeated 3 times
        check_5=!identical(boyce_glmm_table$model, c(rep("glm", 12), rep("gam", 12), rep("rf", 12)))
            #model should be the three algorithms repeated 12 times each

        #if any check is TRUE, stop
        if(check_1 | check_2 | check_3 | check_4 | check_5){
            stop(paste("ERROR! FALSE! WE HAVE A PROBLEM WHEN CREATING THE PHYLO DIFF GLMM TABLE FOR SPECIES ", species, sep=""))
        }

        #return the table
        return(boyce_glmm_table)
    }
}

#apply the function to one species
#get_boyce_phylo_diff_glmm("halepensis")

#apply the function across all species
list_results_boyce_phylo_glmm=lapply(unique(naturalized_occurrences$species), get_boyce_phylo_diff_glmm)
    #get a list with the cleaned boyce table per species

#check we have the correct rows and columns in each table
#x=list_results_boyce_phylo_glmm[[1]]
#x=list_results_boyce_phylo_glmm[[17]]
check_rows_boyce_phylo_diff_glmm=FALSE %in% (sapply(list_results_boyce_phylo_glmm, {function(x) if(!is.null(x)){nrow(x)==12*3}}))
    #each table should have 36 rows: 12 partitions times 3 algorithms 
check_cols_boyce_phylo_diff_glmm=FALSE %in% (sapply(list_results_boyce_phylo_glmm, {function(x) if(!is.null(colnames(x))){identical(colnames(x), c("species", "partition", "model", "non_phylo", "phylo", "phylo_diff"))}}))
if(check_rows_boyce_phylo_diff_glmm | check_cols_boyce_phylo_diff_glmm){
    stop("ERROR! FALSE! WE HAVE A PROBLEM WITH THE EXTRACTION OF PHYLO AND NON-PHYLO BOYCE")
}

#bind all tables as the have the same columns
results_boyce_phylo_diff_glmm=bind_rows(list_results_boyce_phylo_glmm, .id=NULL)

#save the table
write.table(results_boyce_phylo_diff_glmm, paste("./results/global_test_phylo_current/results_boyce_phylo_diff_glmm.tsv", sep=""), sep="\t", row.names=FALSE)


##check for significant differences
#info for GLMM that we finally do not used
    #wilcoxon vs GLMM
        #Note that we are just interested in knowing whether the Boyce index is different between applying or not the phylogenetic correction. Using a GLMM would show that the impact differs between species for it is more difficult to show that the phylogenetic correction causes a change over zero for a given model an species, because you have to set a species/model as reference. Putting sylvestris and strobus as references show that they are indeed different from the rest, but again, we obtain many p-values and it is more difficult to present.
        #The wilcoxon test let us to just show if the median Boyce index of with phylo and non-phylo across the 12 partitions of a species/model combination is indeed different.
        #We do not need to add the partition or models as factors as we are calculating the p-value within the 12 partitions of each species/model. We also say to the test that boyce with and without phylo from partition 1 are paired, i.e., they are not independent. See below.
    #Response: 
        #the difference in boyce index between phylo and non-phylo.
        #slope different from zero for factors means that the correction has an impact
    #Main effects
        #species
        #algorithm
        #invaded region?
            #maybe separate NA, EU and AUS? maybe phylo works different depending on the breath of climatic conditions of the continent?
            #like in "Climatic Niche Shifts Are Rare Among Terrestrial Plant Invaders"
            #but many species have occurrences an South Africa, South America and AUS, you would have to do a level factor with these regions so not sure this is meaningful...
    #interactions
        #species*algorithms
            #the key here is the interaction, phylo diff tends to be different depending on the species and the algorithm? in other words, the impact of the phylogenetic correction is dependent on the species and the algorithm used?
            #maybe an overall effect is not visible, but we can see impact in specific cases.
    #random factor: partition of the species
        #partition_1_halepensis, partition_2_halepensis!!!
        #the first partition of halepensis is not the same than the first partition of sylvestris!!!
        #we need to control for this factor but we are not interested to know its effect. Maybe a partition is more difficult to predict than other. We need to let know the model that the boyce index from glm, gam and RF for halepensis_1 are all obtained from the same data!! and they are different from halepensis_2 even though they belong to the same species.
            #this can have an effect on thet variance of the response variable
            #but we do not want to lose degrees of freedom because of this, because we have maaany levels, halepensis_1, halepensis_2.... radiata_1..... better a random factor.
        #we cannot test the interaction with other factors because, again, partition 1 in halepensis is not present in sylvestris, is not fully crossed
        #In cayuelas course, he used a random block with 10 levels and 3 observations in each one, so we should be ok using partition as a random factor. We have three models for each partition_species, i.e., 3 observations in each partition. We also many levels (12 per each species), so we can use as a random factor.
            #we can use partition as a random factor testing the impact on the average value of the response, but we cannot test how influence the impact of a specific factor like model, because we only have 1 boyce value per model and partition, see below!
            #/home/dftortosa/diego_docs/science/formacion/cursos/courses_during_phd/investigacion/R/Modelos mixtos/Modelos mixtos - Luis Cayuelas/5-Modelos lineales mixtos en R.pdf
        #possible structures
            #1|partition:
                #depending on the partition, the boyce index tends to be higher or lower respect to the average (intercep)
            #species|partition
                #all species are not included in each partition, so we cannot test this
            #algorithm|partition
                #we cannot test this because, yes, we have 3 observations per partition and all models are inside each partition, but we only have 1 observation per model and partition!, not 3!
                #indeed, we get the following error testing this
                    #Error: number of observations (=72) <= number of random effects (=72) for term (model | partition); the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
                #but all algorithms are, i.e., we have boyce index for GLM, GAM and RF in each specific partition, so we can test whether the impact of the algorithm depends on the partition

#define a function to apply the wilcoxon test
require(exactRankTests)
    #for wilcox.exact()
#species_algorithm="radiata_glm"
wilcoxon_signed_calc=function(species_algorithm){

    #extract the species and the algorithm
    species=strsplit(species_algorithm, split="_")[[1]][1]
    algorithm=strsplit(species_algorithm, split="_")[[1]][2]

    #select the rows of the selected species and model
    subset_boyce=results_boyce_phylo_diff_glmm[which(results_boyce_phylo_diff_glmm$species==species & results_boyce_phylo_diff_glmm$model==algorithm),]
    
    #check
    if(FALSE %in% (unique(subset_boyce$species)==species) | FALSE %in% (unique(subset_boyce$model)==algorithm)){
        stop(paste("ERROR! FALSE! WE HAVE A PROBLEM WITH THE SUBSET OF DATA FOR WILCOXON FOR ", species))
    }
    if(FALSE %in% (((subset_boyce$non_phylo-subset_boyce$phylo)-subset_boyce$phylo_diff)<1e-15)){
        stop(paste("ERROR! FALSE! WE HAVE A PROBLEM WITH THE CALCULATION OF PHYLO DIFF FOR WILCOXON FOR ", species))
    }

    #calculate the median values of boyce
    median_non_phylo=median(subset_boyce$non_phylo)
    median_phylo=median(subset_boyce$phylo)
    median_phylo_diff=median(subset_boyce$phylo_diff)

    #calculate the test only if boyce phylo and non-phylo are different
    if(!identical(subset_boyce$phylo, subset_boyce$non_phylo)){

        #run the wilcoxon test obtaining exact p-values always, even if ties
        wilcoxon_signed_raw=wilcox.exact(
                x=subset_boyce$phylo, 
                y=subset_boyce$non_phylo, 
                alternative="two.sided",
                paired=TRUE)
            #we are using a R function (wilcox.exact) that computes the paired rank wilcoxon test with exact p-values even in the case of ties, i.e., if two pairs of observations have the difference between phylo and non-phylo so it is no so simple to rank them. For example, partition 1 and partition 2 have both a difference of 0.1 between phylo and non-phylo.
                #https://stats.stackexchange.com/q/597782
            #x and y are the values to be compared
            #alternative="two.sided": test whether x is greater or less than y
            #paired=TRUE: paired test
            #If both 'x' and 'y' are given and 'paired' is 'TRUE', a Wilcoxon signed rank test is performed. The null hypothesis is that the distribution of 'x - y' is symmetric about 'mu'. "mu" is 0 by default, thus we are testing whether the differences are symmetric around zero.
                #The Wilcoxon Signed Rank Test is the non-parametric version of the paired t-test. It is used to test whether or not there is a significant difference between two population medians.
                #We compare pairs of values that are not independent. For example, weight of a patient before and after an intervention. In our case, the boyce index of a model in a partition before and after and applying the phylogenetic correction.
                    #https://www.statology.org/wilcoxon-signed-rank-test/
                    #https://www.statology.org/wilcoxon-signed-rank-test-r/
                #Wilcoxon Signed Rank Test vs T test
                    #Hypothesis: Student’s t-test is a test comparing means, while Wilcoxon’s tests the ordering of the data. For example, if you are analyzing data with many outliers such as individual wealth (where few billionaires can greatly influence the result), Wilcoxon’s test may be more appropriate.
                    #Interpretation: Although confidence intervals can also be computed for Wilcoxon’s test, it may seem more natural to argue about the confidence interval of the mean in the t-test than the pseudomedian for Wilcoxon’s test.
                    #Fulfillment of assumptions: The assumptions of Student’s t-test may not be met for small sample sizes. In this case, it is often safer to select a non-parametric test. However, if the assumptions of the t-test are met, it has greater statistical power than Wilcoxon’s test.
                    #https://www.datascienceblog.net/post/statistical_test/signed_wilcox_rank_test/
                #Given we have only 12 data points for each test, I will go for the wilcoxon test. T-test can deal with violations of normality but having enough sample size, which is not our case. Also, we should check test by test across the 42 tests whether the differences in boyce index are or not indeed normal distributed. Wilcoxon is the best option for our case, as we avoid problems with normality.
                #The hypotheses tested with the wilcoxon signed test are the following
                    #H0: the median difference between the two groups is zero.
                    #H1: The median difference is positive or negative, i.e., the phylogenetic correction makes the boyce index to go higher or lower.

        #get the statistics
        wilcoxon_signed_statistic=wilcoxon_signed_raw$statistic
        names(wilcoxon_signed_statistic)=NULL
            #check this link if you want to know how to calculate the statistic by hand, it is pretty straightforward.
            #It is just calculate the difference between each pair in absolute value. Then rank cases from rank 1 to higher rank (i.e., higher absolute difference). Then separate cases where the difference is negative or positive. Then sum negative and positive ranks, separately. The smaller value of the two sums (in absolute value) is the statistic.
            #For example, if for all partitions phylo-non_phylo is negative, all ranks are negative, and the sum of positive ranks is zero. Then zero is the value of the statistic. 
            #This would be a extreme case of clear differences between groups, hence as the statistic is smaller, more differences between the groups, because this means that differences tend to have the same sign, meaning that it is always the same group the one having a greater value, indicating significant differences.
                #https://www.statology.org/wilcoxon-signed-rank-test/
        wilcoxon_signed_p=wilcoxon_signed_raw$p.value 

        #run also the original wilcoxon function and save the warning in case it exists
        #define a function to get the warnings and the output in a one single object
        #input is code, a function with its arguments
        mytryCatch <- function(expr){

            #try to run the expression
            tryCatch(
                expr=expr,
                warning=function(w){
                    warning_results=list()
                    warning_results$warning=w
                    warning_results$value=expr
                    return(warning_results)
                }
            )
                #we run function within try catch
                #provide a function for warnings, that return the warning and the result of running the expression
        }
            #https://www.statology.org/r-trycatch/
        #run original wilcoxon function
        wilcoxon_signed_raw1=mytryCatch(
            wilcox.test(
                x=subset_boyce$phylo, 
                y=subset_boyce$non_phylo, 
                alternative="two.sided",
                paired=TRUE
            )
        )

        #save the p_value and statistics, along with the variable indicating whether we got or not the warning about the ties
        #if the warning is NULL we get FALSE, if it is not NULL but it does not include the message of "ties", also FALSE. We also get TRUE if there is a warning and it is related to the ties.
        check_ties=ifelse(is.null(wilcoxon_signed_raw1$warning), FALSE, grepl("cannot compute exact p-value with zeroes", wilcoxon_signed_raw1$warning, fixed=TRUE))
        #if TRUE, i.e., we have a warning about the ties
        if(check_ties){

            #we have the warning about ties, so the original wilcoxon function does not calculate exact p-values, it uses instead a normal approximation. Thus, we can only compare the statistic
            if(wilcoxon_signed_statistic!=wilcoxon_signed_raw1$value$statistic){
                stop(paste("ERROR! FALSE! WE HAVE A PROBLEM WITH THE CALCULATION OF THE WILCOXON TEST FOR ", species, sep=""))
            } 
        } else { #we do not have warning about ties, so the p-values should be also equal, we can compare them
            if(
                wilcoxon_signed_statistic!=wilcoxon_signed_raw1$statistic |
                round(wilcoxon_signed_p,7)!=round(wilcoxon_signed_raw1$p.value,7)
            ){
                stop(paste("ERROR! FALSE! WE HAVE A PROBLEM WITH THE CALCULATION OF THE WILCOXON TEST FOR ", species, sep=""))
            }  
        }
            #ties means that when you calculate phylo boyce - non_phylo boyce for each partition, and you try to make a rank of partitions, two or more partitions have exactly the same value of phylo_diff, thus they cannot be sorted by ranking. This makes the function unable to calculate a exact p-value and it has to use a normal approximation.
            #I am not sure if this is a problem when not having normal data, but just in case, I am going to be careful to consider as significant results with this warning, so I have to know it.
            #https://stackoverflow.com/questions/60112609/warnings-in-wilcoxon-rank-sum-test-in-r
    } else {
        wilcoxon_signed_statistic=NaN
        wilcoxon_signed_p=NaN
    }

    return(cbind.data.frame(species, algorithm, median_non_phylo, median_phylo, median_phylo_diff, wilcoxon_signed_statistic, wilcoxon_signed_p))
}

#run the function in one species
#wilcoxon_signed_calc(species="banksiana_glm")

#get a vector with all species/models combinations
species_algorithms=NULL
for(species in unique(results_boyce_phylo_diff_glmm$species)) for(algorithm in unique(results_boyce_phylo_diff_glmm$model)) species_algorithms=append(species_algorithms, paste(species, "_", algorithm, sep=""))
    #for each species and for each model, paste both
if(length(species_algorithms)!=length(unique(results_boyce_phylo_diff_glmm$species))*length(unique(results_boyce_phylo_diff_glmm$model))){
    stop("ERROR! FALSE! WE HAVE A PROBLEM ")
}

#apply the function across all species and model combinations
list_wilcoxon_signed=lapply(species_algorithms, wilcoxon_signed_calc)

#bind all rows into a DF
require(dplyr)
wilcoxon_signed_results=bind_rows(list_wilcoxon_signed, .id=NULL)

#use species/model as row names for easier visualization in the terminal
rownames(wilcoxon_signed_results)=paste(wilcoxon_signed_results$species, "_", wilcoxon_signed_results$algorithm, sep="")

#see species with p<0.05
significant_nominal_p_value=wilcoxon_signed_results[which(wilcoxon_signed_results$wilcoxon_signed_p<0.05), ]
print("significant cases according to nominal p-value")
print(significant_nominal_p_value)

#calculate the FDR across all p-values
wilcoxon_signed_results$wilcoxon_signed_fdr = p.adjust(wilcoxon_signed_results$wilcoxon_signed_p, method="BH")
    #correction for BH, it is valid use this correction when markerse are correlated directly, if two markers are correlated, a higher significance in one, will entail higher significance in the other. In the case the correlation was negative, we should used other correction, I think remember that Benjamini & Yekutieli (2001), but check. For further details see: "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_ucp_cv/p_value_correction/pvalue_correction.pdf"

#see species with FDR<0.05
significant_fdr=wilcoxon_signed_results[which(wilcoxon_signed_results$wilcoxon_signed_fdr<0.1),]
print("significant cases according to FDR<0.1")
print(significant_fdr)

#multiple test correction
    #we have many non-independent tests
        #all Boyce values calculated for sylvestris have in common the species, i.e., the same partitions are used for glm, gam and rf.
        #all Boyce values obtained with glm have in common that have been calculated with the same algorithm. Indeed, algorithm as a factor had a significant effect when running an anova on GLM of phylo-diff (you can easily run it if required to justify the lack of multiple-test correction during revision).
        #Therefore, we have correlation between tests of the same species and tests of the same model, decreasing a lot the number of independent test. If we assume that values of the same species are correlated and values of the same model are correlated, the 42 independent tests divided by the 14 species and the 3 models gives 1!
            #you have 3 values for sylvestris across models, all calculated with the same data!
            #you have 14 values for GLM across species, all calculated with the same model!
    #Because of this, we are not doing so many tests. We should not use Bonferroni and, indeed, we should not need to apply multiple test correction at all, the number of independent tests is really low.

#Results:
    #IMPORTANT: we are discussing results using the percentage respect to the non-phylo value, but this could be misleading because 0.95 vs 0.98 is going to be a lower difference in percentage than 0.25 vs 0.28.
        #maybe we can calculate a percentage with respect to the range of Boyce index, i.e., from -1 to 1, making a length of 2. 
        #an increase of 0.5 would be 0.5/2 would be a 25% increase.
    #Table XXX show that for 13, at least one model has a Boyce index above 0, i.e., performing better than the random expectation. 11 of these species show a high Boyce index (above 0.5), indicating a good ability to predict naturalized occurrences on unseen data. Note that, from the original 25 species of Perret's dataset, we could only perform the complete evaluation of the models for 14 species due to a reduced number of presences after data cleaning. Therefore, our models had a better performance predicting new occurrences than the random expectation for majority of the species analyzed (13 out 14 species). This suggests that, overall, these models could be useful to predict suitability of pines species. Note, however, that we did not perform an exhaustive independent evaluation of the whole Pinus genus as we were bounded by the availability of naturalized occurrences, limiting our ability to determine to what exent our models work well for the whole genus. Neverthelss, we found a good performance for most of the species studied, so it is plausible that our models could also work well for the rest of the species of the genus.
    #As shown in table XXX, from all species analyzed, the phylogenetic has a significant impact on the Boyce index for 4 species (using nominal p-value; 3 after multiple test correction). For P. patula the phylogenetic correction had a slight impact. GLM shows in increase in the Boyce index from 0.427 to 0.437 (2.34% increase; median values across partitions hereafter), while GAM shows an increase from 0.205 to 0.208 (1.46% increase). For P. sylvestris, GAM shows a slight decrease in the Boyce index from 0.914 to 0.906 (0.87% decrease), but still maintaining a very high performance. In this species, RF shows an increase of peformance from 0.525 to 0.541 (3% increase). Finally, for P. strobus, GLM shows an increase from -0.446 to -0.169 (38% increase), while GAM shows an increase from -0.242 to -0.127 (52% increase). Note, however, that despite the great increase in performance due to the phylogenetic correction (Boyce index goes from -1 to 1), the models still perform worse than the random expectation.
    #It is important to note that our ability to evaluate the performance of the phylogenetic correction depends on the availability of naturalized occurrences. We need independent data from outside of the current ranges in order to evaluate the ability of the correction to include parts of the climatic fundamental niche that are not present in current ranges. Therefore, this evaluation approach is limited by the amount of truly naturalized presences of pines. As previously explained, we could perform the full independent evaluation for only 14 pine species (n ocurrences>=5), with only 11 of them having more than 10 occurrences. 
    #Despite this limitation, we found several instances of significant improvements in performence when the phylogenetic correction was applied ranging from 1.46% to 52% increases in the Boyce index. The only case for significant and negative impact of the phylogenetic correction was a model already showing a very higher performance, i.e., Boyce index>0.9, and after the application of the correction the performance still remained above 0.9. Therefore, it does not to make worse the good-performing models. This is congruent with the conservative approach we follow to apply the phylogenetic correction. We limite its application to areas showing uncertainty across SDMs predictions in order to avoid a great influence of phylogenetic incertainty on robust models. Note that this could also explain the limited impact of the phylogenetic correction on the models. We explored more liberal applications of the phylogenetic on the 2 species showing the greatest impact, namely P. sylvestris and P. strobus. In P. sylvestris, the application of the correaction across all areas (independently of the uncertainty) improves the performance of GLM from 0.741 to 0.801 (%%% PVALUE), of GAM from 0.915 to 0.920 (XXX) and finally, the performance of RF increased from 0.525 to 0.547 (XXX). These are greater improvements compared to the more conervative phylogenetic approach used in the manuscript. The improvements with the liberal pylogenetic approach are even more marked in P. strobus. These species showed great improvements of performance whith the conservative phylogenetic correction, but being still below the random expectation. The more liberal approach makes positive the Boyce index for the three models, i.e., all perform above the random expectation, increasing the performance more than 100% (from -0.447 to 0.183, form -0.243 to 0.460 and from -0.047 to 0.324 for GLM, GAM and RF, respectively). Therefore, a more liberal approach was able to improve relatively good-performing models (P. sylvestris) along with models performing much worse than the random expectation (P. strobus). We did not further develop these analyses given the risk of overfitting. The selection of a different phylogenetic approach based on the performance observes in Perret's dataset would lead to the need of an additional independent dataset where test the best performing phylogenetic approach. In other words, if we tune the phylogenetic approach using the test set, this set is no longer unseen or independent as it has been used for fine tunning the models. Neverthelss, these results shows the potential of the phylogenetic correction to improving model performance and supports the interest in further explore it in the future using more liberal approximations.

#save table
write.table(wilcoxon_signed_results, "./results/global_test_phylo_current/wilcoxon_test_phylo.tsv", sep="\t", col.names=TRUE, row.names=FALSE)




##################
##### FINISH #####
##################

#finish the script
print("## FINISH ##")

#singularity exec 01_global_test_phylo_ubuntu_20_04_v1.sif ./code/phylo/global_test_phylo_current_process_outputs_v1.R 2>&1 ./code/phylo/global_test_phylo_current_process_outputs_v1.Rout




######################
##### NEXT STEPS #####
######################
