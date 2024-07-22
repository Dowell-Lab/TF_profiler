#!/bin/bash
#SBATCH --job-name=rbg                          # Job name
#SBATCH --mail-type=ALL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremail@colorado.edu                # Where to send mail
#SBATCH --nodes=1                                        # Number of cores job will run on
#SBATCH --ntasks=12                                       # Number of CPU (processers, tasks)
#SBATCH --time=06:00:00                                  # Time limit hrs:min:sec
#SBATCH --partition short                                # Job queue
#SBATCH --mem=150gb                                        # Memory limit
#SBATCH --output=/scratch/Users/yourUser/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/yourUser/eofiles/%x_%j.err
                           
module purge
module load python/3.6.3
module load bedtools/2.25.0
module load samtools/1.8
module load meme/5.0.3

BIN=/Users/yourUser/TF_Profiler

GENOME=/scratch/Shares/dowell/genomes/hg38/hg38.fa
BACKGROUND=/Users/yourUser/TF_Profiler/assets/enhancer_background_flat

###HOCO
MOTIFS=/Users/yourUser/TF_Profiler/assets/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
PRESCAN_MOTIFS=/scratch/Shares/dowell/tajo/hocomoco_flat/motifs
SIM_DISTANCES=/scratch/Shares/dowell/tajo/hoco_background_distances

### If you have access to fiji using the prescan motifs and pre-calculated simulated distances will be way faster and greatly reduce resouce costs. To use these you can use the -q ${PRESCAN_MOTIFS} and -k ${SIM_DISTANCES} commands. Otherwise you will need to scan for experimental motif hits (-x True) and generate simulated sequences and scan them for motif hits (-d True).

NAME=HCT116  #most output files will use the name as the root. It is wise to name this something meaningful, such as the cell type used.

mkdir -p /scratch/Users/yourUser/TF_Profiler/${NAME}

python3 ${BIN} \
-g ${GENOME} \
-a /scratch/Users/yourUser/mumerge_outputs/${NAME}_MUMERGE.bed \
-o /scratch/Users/yourUser/TF_Profiler/${NAME} \
-s ${NAME} \
-m ${MOTIFS} \
-b ${BACKGROUND} \
-v \
-c 12 \
-e 742 \
-l False \
-d True \
-x True \
#-q ${PRESCAN_MOTIFS} \
#-k ${SIM_DISTANCES} \
