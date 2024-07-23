#!/bin/bash
#SBATCH --job-name=test                          # Job name
#SBATCH --mail-type=ALL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tajo5912@colorado.edu                # Where to send mail
#SBATCH --nodes=1                                        # Number of cores job will run on
#SBATCH --ntasks=1                                       # Number of CPU (processers, tasks)
#SBATCH --time=00:20:00                                  # Time limit hrs:min:sec
#SBATCH --partition short                                # Job queue
#SBATCH --mem=15gb                                        # Memory limit
#SBATCH --output=/scratch/Users/tajo5912/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/tajo5912/eofiles/%x_%j.err
                           
module purge
module load python/3.6.3
module load bedtools/2.25.0
module load samtools/1.8
module load meme/5.0.3


BIN=/Users/tajo5912/rbg/gitTEST/TF_profiler
GENOME=/scratch/Shares/dowell/genomes/hg38/hg38.fa

OUTDIR=/scratch/Users/tajo5912/test

mkdir -p ${OUTDIR}

source /Users/tajo5912/rbg/gitTEST/venv_tfprofiler/bin/activate

python ${BIN} \
-g ${GENOME} \
-a ${BIN}/test/MCF7_test_MUMERGE.bed \
-m ${BIN}/test/p53_estrogen_IFN_mini.meme \
-b ${BIN}/assets/enhancer_background_flat \
-o ${OUTDIR} \
-s test_sim \
-v \
-c 1 \
-e 118 \
-l False \
-d True \
-x False \
-n 20 \
-i 1
