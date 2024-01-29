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

#create slurm files



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
intervals_batches=c(seq(1,length(unique_species),10), length(unique_species))
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


#i=1
for(i in 1:length(batches)){
    selected_batch=batches[[i]]
    if(i!=length(batches)){
        if(TRUE %in% (batches[[i]] %in% unlist(batches[(i+1):length(batches)]))){
            stop("ERROR! FALSE! WE HAVE OVERLAP BETWEEN THE BATCHES")
        }
    }
}



if(!identical(unlist(batches), unique_species)){
    stop("ERROR! FALSE! WE HAVE LOST SPECIES WHEN CREATING THE BATCHES")
}



#check no overlap between batches


#batch_index=1
for(batch_index in 1:length(batches)){
    batch=batches[batch_index]
    fileConn<-file(paste("./code/phylo/recipes/01_global_test_phylo_ubuntu_20_04_v1_slurm_files/slurm_global_test_", names(batch), ".slurm", sep=""))
        #https://stackoverflow.com/a/2470277
    slurm_script=paste0(
        "#!/bin/bash/
        # --------------------------------------------------------------
        ### PART 1: Requests resources to run your job.
        # --------------------------------------------------------------
        ###info about slurm commands: https://slurm.schedmd.com/pdfs/summary.pdf
        ### Optional. Set the job name
        #SBATCH --job-name=global_test_phylo_", names(batch), "
        ### Optional. Set the output filename.
        ### SLURM reads %x as the job name and %j as the job ID
        #SBATCH --output=%x-%j.out
        #SBATCH --error=%x-%j.err
        ### Optional. Request email when job begins and ends
        ### SBATCH --mail-type=ALL
        ### Optional. Specify email address to use for notification
        ### SBATCH --mail-user=dftortosa@email.arizona.edu
        ### REQUIRED. Set the partition for your job.
        #SBATCH --partition=albaicin
        ### REQUIRED. Set the number of cores and nodes that will be used for this job. It seems that each node has 56 cores, so you have to adjust accordingly
        #SBATCH --nodes=1
        #SBATCH --ntasks=", length(batch[[1]]), "
        #SBATCH --ntasks-per-node=", length(batch[[1]]), "
        ### OPTIONAL. You can set the number of cores per task. This can be useful for scripts in which inside the parallelized process there are other subprocesses. For example, you run 22 independent processes for each chromosome and then you also parallelize the calculations inside each chromosome (https://stackoverflow.com/a/51141287)
        #SBATCH --cpus-per-task=1
        ### REQUIRED. Set the memory required for this job. I will set 40GC per each of the 100 cores=4000GB; https://public.confluence.arizona.edu/display/UAHPC/Allocation+and+Limits
        ###SBATCH --mem=400gb
        ### REQUIRED. Set the gb per core. YOU HAVE TO SELECT --mem or --mem-per-cpu but NOT BOTH. If you get a .core file, this usually means that the program fails because it asked for too much memory, so it creates a record of the working memory at the time that can be used for debugging. MPI jobs will usually create a core file for each task. You should increase memory limits (https://researchcomputing.princeton.edu/support/knowledge-base/memory).
        #SBATCH --mem-per-cpu=12gb
        ### set the constraint for high memory nodes in case you use a lot of memory per node. Normal nodes have a 512Gb limit.
        ###SBATCH --constraint=hi_mem
        ### REQUIRED. Specify the time required for this job, hhh:mm:ss
        #SBATCH --time=24:00:00

         
        # --------------------------------------------------------------
        ### PART 2: Executes bash commands to run your job
        # --------------------------------------------------------------
        module load singularity
        ### change to your script’s directory
        cd /home/UGR002/dsalazar/phd/nicho_pinus
            #/home is the stable directory, while scratch is where results of analyses can be stored temporary, as stuff gets removed after 20 days
        ### Run your work
        singularity exec 01_global_test_phylo_ubuntu_20_04_v1.sif ./scripts/phylo/global_test_phylo_current_v1.R --species='", paste(batch[[1]], collapse=","), "' --batch='", names(batch), "' > ./scripts/phylo/global_test_phylo_current_v1_", names(batch), ".Rout 2>&1", sep="")
    
    #remove the tabs at the beginning
    slurm_script=gsub("        ", "", slurm_script, fixed=TRUE)
        #https://stackoverflow.com/a/11936385
    writeLines(slurm_script, fileConn)


    close(fileConn)   
}
