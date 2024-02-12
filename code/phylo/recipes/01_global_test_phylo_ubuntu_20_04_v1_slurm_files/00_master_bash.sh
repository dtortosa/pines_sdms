#!/bin/bash
chmod +x ../../global_test_phylo_current_v1.R
sbatch slurm_global_test_batch_1.slurm; 
sbatch slurm_global_test_batch_2.slurm; 
sbatch slurm_global_test_batch_3.slurm; 
sbatch slurm_global_test_batch_4.slurm; 
sbatch slurm_global_test_batch_5.slurm; 
sbatch slurm_global_test_batch_6.slurm; 
sbatch slurm_global_test_batch_7.slurm; 
sbatch slurm_global_test_batch_8.slurm; 
sbatch slurm_global_test_batch_9.slurm; 
sbatch slurm_global_test_batch_10.slurm; 
sbatch slurm_global_test_batch_11.slurm; 
sbatch slurm_global_test_batch_12.slurm; 
sbatch slurm_global_test_batch_13.slurm; 
sbatch slurm_global_test_batch_14.slurm; 
sbatch slurm_global_test_batch_15.slurm; 
sbatch slurm_global_test_batch_16.slurm; 
sbatch slurm_global_test_batch_17.slurm; 
sbatch slurm_global_test_batch_18.slurm; 
sbatch slurm_global_test_batch_19.slurm; 
sbatch slurm_global_test_batch_20.slurm; 
sbatch slurm_global_test_batch_21.slurm; 
sbatch slurm_global_test_batch_22.slurm; 
sbatch slurm_global_test_batch_23.slurm; 
sbatch slurm_global_test_batch_24.slurm; 
sbatch slurm_global_test_batch_25.slurm; 

n_jobs=$(squeue -u dsalazar | awk '{if(NR!=1 && $4 ~/^global_test_phylo/){count++}}END{print count}')
echo 'WE HAVE ' $n_jobs ' jobs'
#to stop all jobs #squeue -u dsalazar | awk '{if(NR!=1){print $1}}' | xargs scancel
#chmod +x 00_master_bash.sh; ./00_master_bash.sh > 00_master_bash.out 2>&1
