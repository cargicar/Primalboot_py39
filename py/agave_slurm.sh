#!/bin/bash                                                                              

#SBATCH -o RF_playing.%u.out # file to save job's STDOUT (%j = JobId)                    
#SBATCH -e RF_playing.%u.err # file to save job's STDERR (%j = JobId)                    
#SBATCH -p serial # Whatever_the_name_of_your_SLURM_partition_is                         
#SBATCH -c 10 # cores                                                                    
#SBATCH -t 2-12 # 2*24+12 hours wall time                                                
#SBATCH --export=NONE   # Purge the job-submitting shell environment                    \
                                                                                         

module purge
module load conda/rolling
source activate primal39

python run.py
