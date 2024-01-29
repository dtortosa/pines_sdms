#!/bin/bash/
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
###info about slurm commands: https://slurm.schedmd.com/pdfs/summary.pdf
### Optional. Set the job name
#SBATCH --job-name=global_test_phylo_batch_1
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
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=10
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
singularity exec 01_global_test_phylo_ubuntu_20_04_v1.sif ./scripts/phylo/global_test_phylo_current_v1.R --species='banksiana,canariensis,caribaea,clausa,contorta,elliottii,halepensis,jeffreyi,mugo,muricata' --batch='batch_1' > ./scripts/phylo/global_test_phylo_current_v1_batch_1.Rout 2>&1
