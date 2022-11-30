#!/bin/bash
#SBATCH --job-name=rbg                          # Job name
#SBATCH --mail-type=ALL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tajo5912@colorado.edu                # Where to send mail
#SBATCH --nodes=1                                        # Number of cores job will run on
#SBATCH --ntasks=12                                       # Number of CPU (processers, tasks)
#SBATCH --time=06:00:00                                  # Time limit hrs:min:sec
#SBATCH --partition short                                # Job queue
#SBATCH --mem=150gb                                        # Memory limit
#SBATCH --output=/scratch/Users/tajo5912/bapplication_v2/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/tajo5912/bapplication_v2/eofiles/%x_%j.err
                           
module purge
module load python/3.6.3
module load bedtools/2.25.0
module load samtools/1.8
module load meme/5.0.3

BIN=/Users/tajo5912/rbg/rbg

GENOME=/scratch/Shares/dowell/genomes/hg38/hg38.fa
BACKGROUND=/scratch/Users/tajo5912/dbapplication/scripts/enhancer_background_flat

###HOCO
MOTIFS=/scratch/Users/tajo5912/dbapplication/scripts/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
PRESCAN_MOTIFS=/scratch/Shares/dowell/tajo/hocomoco_flat/motifs
SIM_DISTANCES=/scratch/Users/tajo5912/bapplication_v2/hoco_background_distances

###BEST CURATED
#MOTIFS=/scratch/Users/tajo5912/dbapplication/scripts/best_curated_Human.meme
#PRESCAN_MOTIFS=/scratch/Shares/dowell/tajo/best_curated_flat/motifs
#SIM_DISTANCES=/scratch/Users/tajo5912/bapplication_v2/basal_background_distances

NAME=$1

mkdir -p /scratch/Users/tajo5912/bapplication_v2/rbg_hoco/${NAME}

python3 ${BIN} \
-g ${GENOME} \
-a /scratch/Users/tajo5912/bapplication_v2/mumerge_outputs/${NAME}.bed \
-o /scratch/Users/tajo5912/bapplication_v2/rbg_hoco/${NAME} \
-s ${NAME} \
-v \
-c 12 \
-e 742 \
-l False \
-d False \
-x False \
-q ${PRESCAN_MOTIFS} \
-k ${SIM_DISTANCES} \
-m ${MOTIFS} \
-b ${BACKGROUND}
