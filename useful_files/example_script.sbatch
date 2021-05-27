#!/bin/bash
#SBATCH --job-name=wf                           # Job name
#SBATCH --mail-type=ALL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tajo5912@colorado.edu                # Where to send mail
#SBATCH --nodes=1                                        # Number of cores job will run on
#SBATCH --ntasks=1                                       # Number of CPU (processers, tasks)
#SBATCH --time=00:10:00                                  # Time limit hrs:min:sec
#SBATCH --partition short                                # Job queue
#SBATCH --mem=2gb                                        # Memory limit
#SBATCH --output=/scratch/Users/tajo5912/rbg/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/tajo5912/rbg/eofiles/%x_%j.err
                           
module purge
module load python/3.6.3
module load bedtools/2.25.0
module load samtools/1.8
module load meme/5.0.3

BIN=/Users/tajo5912/rbg/rbg

GENOME=/scratch/Shares/dowell/genomes/hg38/hg38.fa
MOTIF=/scratch/Users/tajo5912/rbg/test_files/TP53.meme
NAME=pelican

##run the base_content
mkdir -p /Users/tajo5912/rbg/${NAME}

python3 ${BIN} \
-g ${GENOME} \
-a /scratch/Users/tajo5912/rbg/test_files/levandow_nutlin_sig25_3442peaks.bed \
-o /Users/tajo5912/rbg/${NAME} \
-s ${NAME} \
-c 1 \
-m ${MOTIF} \
-n 3442 \
-w 1500 \
-i 6 \
-b /scratch/Users/tajo5912/rbg/test_files/enhancer_background_for_fimo \
-t 0.0001 \
-e \
-v
