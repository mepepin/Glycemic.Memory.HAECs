#Script written by Mark E. Pepin, PhD | MD-PhD Trainee
# University of Alabama at Birmingham
# April 2019

#Create relative path to folders
CUR_DIR=`pwd`
SCRIPT_DIR="$(dirname $CUR_DIR)"
BASE_DIR="$(dirname $SCRIPT_DIR)"
JOB_DIR=$BASE_DIR/Jobs
LOG_DIR=$BASE_DIR/Logs
RESULTS_DIR=$BASE_DIR/Results
INPUT_DIR=$BASE_DIR/Input
GENOME_DIR=$INPUT_DIR/Genome
mkdir -p $RESULTS_DIR/
##Assemble the toolbox
#Define the study parameters
STUDY_NAME="Napoli"
VAR_LIST="$INPUT_DIR/fastq_names.txt"

#######################################################################################
#### Any Additional Variables you would like to declare (that are common to all jobs)

mkdir -p $JOB_DIR/RRBS_bwa $LOG_DIR/RRBS_bwa $RESULTS_DIR/RRBS_bwa
mkdir -p $RESULTS_DIR/RRBS_bwa
mkdir -p $RESULTS_DIR/RRBS_bwa/methylation_extractor
mkdir -p $RESULTS_DIR/RRBS_bwa/deduplicate_bismark
mkdir -p $INPUT_DIR/fastq_trimmed

#######################################################################################

# iterate and do work
while read VAR; do
        STAR="$INPUT_DIR/fastq/$VAR"
        jobfile=${JOB_DIR}/RRBS_bwa/${STUDY_NAME}_${VAR}.sh
        logfile=${LOG_DIR}/RRBS_bwa/${STUDY_NAME}_${VAR}.log.txt
        DATE=$(date)
        cat > $jobfile <<EOF
#!/bin/bash
# auto-generated job file
# generated from $PWD/$0
# on ${DATE}
#SBATCH --job-name=${VAR}_JOB
#SBATCH --ntasks=1                              # Number of PROCESSES
#SBATCH --cpus-per-task=8                       # Number of PROCESSES
#SBATCH --mem-per-cpu=16000                     # Memory specified for each core used ($
#SBATCH -t 2-02:00:00                           # Runtime in D-HH:MM:SS
#SBATCH --partition=medium                      # express(2h), short(12h), medium(2d2h), 
#
#SBATCH --mail-user=${USER}@uab.edu
#SBATCH --mail-type=ALL                         # BEGIN, END, ERROR, ALL
#
#SBATCH --error=${LOG_DIR}/RRBS_bwa/${STUDY_NAME}_${VAR}.err.txt        # File to which STDERR will be written
#SBATCH --output=${LOG_DIR}/RRBS_bwa/${STUDY_NAME}_${VAR}.out.txt	# File to which STDOUT will be written
#
# Mimimum memory required per allocated  CPU  in  MegaBytes.
#SBATCH --mem-per-cpu=16000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=${USER}@uab.edu             #### Modify to your email address.
##
###########################
## Insert your code here ##
###########################

source activate /home/pepinme/.conda/envs/RRBS_env

#Trim off the adapter sequences (adapter auto-detected)
#trim_galore \
#-o $INPUT_DIR/fastq_trimmed/ \
#--paired --rrbs --non_directional --length 20 --fastqc \
#$INPUT_DIR/fastq/${VAR}_R1_001.fastq.gz $INPUT_DIR/fastq/${VAR}_R2_001.fastq.gz

##
#bwameth.py index $GENOME_DIR/GRCh38.p12.genome.fa
bwameth.py --threads 8 \
--reference $GENOME_DIR/GRCh38.p12.genome.fa \
$INPUT_DIR/fastq_trimmed/${VAR}_R1_001_val_1.fq.gz $INPUT_DIR/fastq_trimmed/${VAR}_R2_001_val_2.fq.gz \
> $RESULTS_DIR/RRBS_bwa/${VAR}.sam

#Convert .sam to .bam
samtools view -S -b $RESULTS_DIR/RRBS_bwa/${VAR}.sam > $RESULTS_DIR/RRBS_bwa/${VAR}.bam
#Sort using samtools
samtools sort $RESULTS_DIR/RRBS_bwa/${VAR}.bam -o $RESULTS_DIR/RRBS_bwa/${VAR}.sorted.bam
#create an index
samtools index $RESULTS_DIR/RRBS_bwa/${VAR}.sorted.bam
#MethylDackel
MethylDackel extract $GENOME_DIR/GRCh38.p12.genome.fa $RESULTS_DIR/RRBS_bwa/${VAR}.sorted.bam -o $RESULTS_DIR/RRBS_bwa/${VAR}.counted --methylKit

###########################
## End of your code here ##
###########################

srun sleep 30

EOF
   	chmod 770 $jobfile
        echo "sbatch --mail-user=pepinme@uab.edu --job-name=${STUDY_NAME}_${VAR}_JOB $jobfile"
        sbatch --mail-user=pepinme@uab.edu --job-name=${STUDY_NAME}_${VAR}_JOB $jobfile
done < "$VAR_LIST"
