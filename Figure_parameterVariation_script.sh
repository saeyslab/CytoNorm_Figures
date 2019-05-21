#!/bin/bash
#$ -l h_vmem=10G
#$ -t 1-25
#$ -N CytofNorm
module unload gcc
export LD_LIBRARY_PATH="/software/shared/apps/x86_64/gcc/4.8.0/lib/:/software/shared/apps/x86_64/gcc/4.8.0/lib64:$LD_LIBRARY_PATH"
module load R/x86_64/3.5.1
cd /group/irc/personal/sofievg/Stanford
Rscript Figure_parameterVariation.R $SGE_TASK_ID