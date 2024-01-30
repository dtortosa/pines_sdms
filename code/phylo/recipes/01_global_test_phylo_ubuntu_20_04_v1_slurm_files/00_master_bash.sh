#!/bin/bash
chmod +x ../global_test_phylo_current_v1.R
qsub slurm_global_test_batch_1.slurm; 
qsub slurm_global_test_batch_2.slurm; 
qsub slurm_global_test_batch_3.slurm; 

n_jobs=$(squeue -u dsalazar | awk '{if(NR!=1){count++}}END{print count}')
echo 'WE HAVE ' $n_jobs ' jobs'
#to stop all jobs #squeue -u dsalazar | awk '{if(NR!=1){print $1}}' | xargs qdel
#chmod +x 00_master_bash.sh; ./00_master_bash.sh > 00_master_bash.out 2>&1
