#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J cclasso0_1_PGE_NOPGE_sp_MAN
### -- ask for number of cores (default: 1) -- 
#BSUB -n 1
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 100GB of memory per core/slot -- 
#BSUB -R "rusage[mem=100GB]"
### -- specify that we want the job to get killed if it exceeds 100 GB per core/slot -- 
#BSUB -M 100GB
### -- set walltime limit: hh:mm -- 
#BSUB -W 120:00
### Select CPU model
#BSUB -R "select[avx512]" 
### -- set the email address -- 
#BSUB -u asaru@dtu.dk
### -- send notification at start -- 
#BSUB -B 
### -- send notification at completion -- 
#BSUB -N 
### -- Specify the output and error file. %J is the job-id -- 
### -- -o and -e mean append, -oo and -eo mean overwrite -- 
#BSUB -o /work3/asaru/Airbiome_project/Bash_scripts/PGE_NOPGE_sp/log/Output_%J.out 
#BSUB -e /work3/asaru/Airbiome_project/Bash_scripts/PGE_NOPGE_sp/log/Output_%J.err 

### Load R module
module load R/4.3.3-mkl2024

### Set up user library path for R packages
export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.3

# Ensure the directory exists
mkdir -p $R_LIBS_USER

### setting the working directory
PROJECT_DIR=/work3/asaru/Airbiome_project

### Run script
Rscript ${PROJECT_DIR}/Build_microb_net.R --threshold 0.1 --network cclasso --biome PGE_NOPGE_sp