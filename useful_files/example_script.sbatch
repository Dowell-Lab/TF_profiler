#!/bin/bash
#SBATCH --job-name=wf                           # Job name
#SBATCH --mail-type=ALL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tajo5912@colorado.edu                # Where to send mail
#SBATCH --nodes=1                                        # Number of cores job will run on
#SBATCH --ntasks=12                                       # Number of CPU (processers, tasks)
#SBATCH --time=23:55:00                                  # Time limit hrs:min:sec
#SBATCH --partition short                                # Job queue
#SBATCH --mem=120gb                                        # Memory limit
#SBATCH --output=/scratch/Users/tajo5912/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/tajo5912/eofiles/%x_%j.err
                           
module purge
module load python/3.6.3
module load bedtools/2.25.0
module load samtools/1.8
module load meme/5.0.3

BIN=/Users/tajo5912/rbg/rbg

GENOME=/scratch/Shares/dowell/genomes/hg38/hg38.fa
BACKGROUND=/scratch/Users/tajo5912/processed_published/scripts/enhancer_background_flat
MOTIFS=/scratch/Users/tajo5912/processed_published/scripts/best_curated_Human.meme

CELL=$1
TX=$2

mkdir -p /scratch/Users/tajo5912/processed_published/cell_type_specific_nutlin/${CELL}${TX}

python3 ${BIN} \
-g ${GENOME} \
-a /scratch/Users/tajo5912/processed_published/cell_type_specific_nutlin/mumerge/${CELL}${TX}_enhancers.bed \
-o /scratch/Users/tajo5912/processed_published/cell_type_specific_nutlin/${CELL}${TX} \
-s ${CELL}${TX} \
-e 75438 \
-w 1500 \
-v \
-c 12 \
-d True \
-p /scratch/Shares/dowell/tajo/motif_scan_flat/motifs \
-m ${MOTIFS} \
-b ${BACKGROUND}
