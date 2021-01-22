#!/bin/bash
#SBATCH --job-name=1_HAEC_0_JOB
#SBATCH --ntasks=1                              # Number of PROCESSES
#SBATCH --cpus-per-task=1                       # Number of PROCESSES
#SBATCH --mem-per-cpu=35000                      # Memory specified for each core used ($
#SBATCH -t 6-06:00:00                           # Runtime in D-HH:MM:SS
#SBATCH --share
#SBATCH --partition=long                      # express(2h), short(12h), medium(2d2h), 
#
#SBATCH --mail-user=pepinme@uab.edu
#SBATCH --mail-type=ALL                         # BEGIN, END, ERROR, ALL
#
#SBATCH --error=/data/scratch/pepinme/Napoli/Logs/RRBS_bwameth/Napoli_bwameth.err.txt        # File to which STDERR will be written
#SBATCH --output=/data/scratch/pepinme/Napoli/Logs/RRBS_bwameth/Napoli_bwameth.out.txt	# File to which STDOUT will be written
#
# Mimimum memory required per allocated  CPU  in  MegaBytes.
#SBATCH --mem-per-cpu=35000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=pepinme@uab.edu             #### Modify to your email address.
##
###########################
## Insert your code here ##
###########################

source activate /home/pepinme/.conda/envs/RRBS_env

bwameth.py index /data/scratch/pepinme/Napoli/Input/Genome/GRCh38.p12.genome.fa

###########################
## End of your code here ##
###########################

srun hostname
srun sleep 30                                   # If running many short jobs, look into$

srun /data/scratch/pepinme/Napoli/Scripts/RRBS_bwa/Index.sh
