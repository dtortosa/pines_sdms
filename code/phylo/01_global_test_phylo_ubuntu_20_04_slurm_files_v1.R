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




###########################
##### PREPARE BATCHES #####
###########################

##define the interval for each batch
intervals_batches=c(seq(1,length(unique_species),1), length(unique_species))
    #sequence from 1 to the total number of species, but the end is not included
    #so we have to add the last species (total number of species) as the end


##run loop across batches
#empty list
batches=list()

#loop across the intervals
#i=1
for(i in 1:length(intervals_batches)){

    #do stuff except for the last interval
    if(i!=length(intervals_batches)){

        #the left end is the current interval extreme
        left_end=intervals_batches[i]

        #the right end is the next one
        right_end=intervals_batches[i+1]

        #if this is not the penultimate interval
        if(i!=(length(intervals_batches)-1)){

            #remove 1 to the right end, so we avoid including it in the next interval
            right_end=right_end-1
                #in the penultimate interval there will be NO more steps because the last interval is not run, so we can take the last species
        }

        #select the species belonging to the selected interval
        batches[[i]]=unique_species[left_end:right_end]
    }
}

#add the names of each batch
names(batches)=paste("batch_", 1:length(batches), sep="")


##check no overlap between batches
#iterate across batches
#i=1
for(i in 1:length(batches)){

    #select the [i] batch
    selected_batch=batches[[i]]
    
    #if the selected batch is not the last one
    if(i!=length(batches)){

        #get the next batches together
        next_batches=unlist(batches[(i+1):length(batches)], use.names=FALSE, recursive=FALSE)
            #unlist losing names and avoiding doing the unlisting recursively as the elements of the list are not nested lists

        #if species of the current batch are included in the next ones
        if(TRUE %in% (selected_batch %in% next_batches)){
            stop("ERROR! FALSE! WE HAVE OVERLAP BETWEEN THE BATCHES")
        }
    }
}


##check we have all species in the batches
flatten_batches=unlist(batches, use.names=FALSE, recursive=FALSE)
    #flatten the list of batches removing names
if(!identical(flatten_batches, unique_species)){
    stop("ERROR! FALSE! WE HAVE LOST SPECIES WHEN CREATING THE BATCHES")
}




#################################
##### PREPARE SLURM SCRIPTS #####
#################################

#batch_index=1
for(batch_index in 1:length(batches)){

    #select the [i] batch
    selected_batch=batches[batch_index]
        #we use "[]" instead of "[[]]" to retain the name of the batch

    #open a connection to the new file
    fileConn<-file(paste("./code/phylo/recipes/01_global_test_phylo_ubuntu_20_04_v1_slurm_files/slurm_global_test_", names(selected_batch), ".slurm", sep=""))
        #https://stackoverflow.com/a/2470277
    
    #create a string with the script
    slurm_script=paste(
        "#!/bin/bash/
        # --------------------------------------------------------------
        ### PART 1: Requests resources to run your job.
        # --------------------------------------------------------------
        ###info about slurm commands: https://slurm.schedmd.com/pdfs/summary.pdf
        ### Optional. Set the job name
        #SBATCH --job-name=global_test_phylo_", names(selected_batch), "
        ### Optional. Set the output filename.
        ### SLURM reads %x as the job name and %j as the job ID
        #SBATCH --output=%x-%j.out
        #SBATCH --error=%x-%j.err
        ### Optional. Request email when job begins and ends
        ### SBATCH --mail-type=ALL
        ### Optional. Specify email address to use for notification
        ### SBATCH --mail-user=dftortosa@gmail.com
        ### REQUIRED. Set the partition for your job.
        #SBATCH --partition=albaicin
        ### REQUIRED. Set the number of cores and nodes that will be used for this job. It seems that each node has 56 cores, so you have to adjust accordingly
        #SBATCH --nodes=1
        #SBATCH --ntasks=", length(selected_batch[[1]]), "
        #SBATCH --ntasks-per-node=", length(selected_batch[[1]]), "
        ### OPTIONAL. You can set the number of cores per task. This can be useful for scripts in which inside the parallelized process there are other subprocesses. For example, you run 22 independent processes for each chromosome and then you also parallelize the calculations inside each chromosome (https://stackoverflow.com/a/51141287)
        #SBATCH --cpus-per-task=1
        ### REQUIRED. Set the memory required for this job. I will set 40GC per each of the 100 cores=4000GB; https://public.confluence.arizona.edu/display/UAHPC/Allocation+and+Limits
        ###SBATCH --mem=400gb
        ### REQUIRED. Set the gb per core. YOU HAVE TO SELECT --mem or --mem-per-cpu but NOT BOTH. If you get a .core file, this usually means that the program fails because it asked for too much memory, so it creates a record of the working memory at the time that can be used for debugging. MPI jobs will usually create a core file for each task. You should increase memory limits (https://researchcomputing.princeton.edu/support/knowledge-base/memory).
        #SBATCH --mem-per-cpu=11gb
        ### set the constraint for high memory nodes in case you use a lot of memory per node. Normal nodes have a 512Gb limit.
        ###SBATCH --constraint=hi_mem
        ### REQUIRED. Specify the time required for this job, hhh:mm:ss
        #SBATCH --time=05:00:00

         
        # --------------------------------------------------------------
        ### PART 2: Executes bash commands to run your job
        # --------------------------------------------------------------
        module load singularity
        ### change to your script’s directory
        cd /home/UGR002/dsalazar/phd/nicho_pinus
            #/home is the stable directory, while scratch is where results of analyses can be stored temporary, as stuff gets removed after 20 days
        ### Run your work
        singularity exec 01_global_test_phylo_ubuntu_20_04_v1.sif ./code/phylo/global_test_phylo_current_v1.R --species='", paste(selected_batch[[1]], collapse=","), "' --batch='", names(selected_batch), "' > ./code/phylo/global_test_phylo_current_v1_", names(selected_batch), ".Rout 2>&1", sep="")
    
    #remove the tabs at the beginning of each line
    slurm_script=gsub(pattern="        ", replacement="", slurm_script, fixed=TRUE)
        #The two '*sub' functions differ only in that 'sub' replaces only the first occurrence of a 'pattern' whereas 'gsub' replaces all occurrences.
            #pattern: character string containing a regular expression (or  character string for 'fixed = TRUE')
            #replacement: a replacement for matched pattern in 'sub' and 'gsub'
            #If 'TRUE', 'pattern' is a string to be matched as is.
            #https://stackoverflow.com/a/11936385
    
    #add the slurm script to the file
    writeLines(text=slurm_script, con=fileConn, sep="\n")
        #text: A character vector
        #con: A connection object or a character string.
        #sep: character string.  A string to be written to the connection after each line of text.

    #close the file
    close(fileConn)   
}




###############################
##### PREPARE BASH SCRIPT #####
###############################

##set path to the script
bash_path="./code/phylo/recipes/01_global_test_phylo_ubuntu_20_04_v1_slurm_files/00_master_bash.sh"

##remove the bash files if it previously exists
if(file.exists(bash_path)){
    system(paste("rm ", bash_path, sep=""))  
}


##open a connection
fileConn<-file(description=bash_path, open="a")
    #description is a path to the file to be opened
    #open: A description of how to open the connection
        #we use append to add new lines
    #https://stackoverflow.com/a/75459017


##add lines to the bash script
#open bash script
writeLines(text="#!/bin/bash", con=fileConn, sep="\n")

#give rights to run the R script
writeLines(text="chmod +x ../../global_test_phylo_current_v1.R", con=fileConn, sep="\n")

#send slurm jobs
#batch_name="batch_1"
for(batch_name in names(batches)){
    writeLines(text=paste("sbatch slurm_global_test_", batch_name, ".slurm; ", sep=""), con=fileConn, sep="\n")   
}

#count the number of jobs on the queue
writeLines(text="\nn_jobs=$(squeue -u dsalazar | awk '{if(NR!=1){count++}}END{print count}')", con=fileConn, sep="\n")
writeLines(text="echo 'WE HAVE ' $n_jobs ' jobs'", con=fileConn, sep="\n")
writeLines(text="#to stop all jobs #squeue -u dsalazar | awk '{if(NR!=1){print $1}}' | xargs scancel", con=fileConn, sep="\n")
writeLines(text="#chmod +x 00_master_bash.sh; ./00_master_bash.sh > 00_master_bash.out 2>&1", con=fileConn, sep="\n")    
    #text: A character vector
    #con: A connection object or a character string.
    #sep: character string.  A string to be written to the connection after each line of text.

#close the file
close(fileConn)  




###########################
# check we have all files #
###########################

#we should have 1 file per batch (3) plus the master bash script
n_files=as.numeric(system("\\
    cd ./code/phylo/recipes/01_global_test_phylo_ubuntu_20_04_v1_slurm_files; \\
    n_files=$(ls | awk 'END{print NR}'); \\
    echo $n_files", intern=TRUE))
if(n_files!=length(batches)+1){
    stop("ERROR! FALSE! WE HAVE NOT OBTAINED JOB FILES FOR ALL BATCHES")
}

#check we have all batches in the master bash script
sbatch_appearance = system(paste("\
    grep \\
        'sbatch' \\
        --only-matching \\
        ", bash_path, sep=""), intern=TRUE)
    #--only-matching: show only nonempty parts of lines that match
        #we get only "sbatch" as many times as it appears even if it is repeated in the same row
        #https://stackoverflow.com/a/3249761
if(length(sbatch_appearance)!=length(batches)){
    stop("ERROR! FALSE! WE HAVE NOT OBTAINED JOB FILES FOR ALL BATCHES")
}




##########
# FINISH #
##########
print("FINISH")
#chmod +x ./code/phylo/01_global_test_phylo_ubuntu_20_04_slurm_files_v1.R; ./code/phylo/01_global_test_phylo_ubuntu_20_04_slurm_files_v1.R > ./code/phylo/01_global_test_phylo_ubuntu_20_04_slurm_files_v1.Rout 2>&1
