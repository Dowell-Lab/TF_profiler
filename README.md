# TF Profiler #
## TF Profiler Overview ##
TF Profiler is a stand alone program that predicts the basally active TFs in a given cellular context. This is not intended to be a differential analysis (please see [TFEA](https://github.com/Dowell-Lab/TFEA)), but rather an assessment of TFs active in the control/baseline conditions. This is achieved by comparing motif co-occurences with bidirectional transcription in experimental control conditions (observed- WT/Vehicle) to those generated from a set of simulated sequences that reflect position specific nucleotide biases (expectation).




This program is associated the publication [A transcription factor (TF) inference method that broadly measures TF activity and identifies mechanistically distinct TF networks](https://www.biorxiv.org/content/10.1101/2024.03.15.585303v1).

Citation:
```
Jones, Taylor, et al. "A transcription factor (TF) inference method that broadly measures TF activity and identifies mechanistically distinct TF networks." bioRxiv (2024): 2024-03.
```

## Installation and Requirements ##
Please see the requirements.txt file. This program was written in python 3.6.3.
It is highly recommended to run this program in a virtual environment.

The python requirements can be installed following these steps:
1. Copy the repository.
```
git clone https://github.com/Dowell-Lab/TF_profiler.git
```
2. Verify that you are running python 3.6.3.
If on the FIJI compute cluster you can load in python 3.6.3, otherwise make sure this is the python version you are using.
```
module load python/3.6.3
which python3
```

3. Create and activate your virtual environment.
```
python3 -m venv venv_TF_Profiler
source venv_TF_Profiler/bin/activate
which python
```
Now your active python should be associated with the virtual environment that you just created.

4. Install package requirements using the requirements.txt file.
```
pip install -r <path/to/TF_Profiler>/requirements.txt
```
Verify that each requirement was successfully downloaded. See requirements here:
```
cycler==0.11.0
joblib==1.1.1
kiwisolver==1.3.1
matplotlib==3.3.4
numpy==1.19.5
pandas==1.1.5
Pillow==8.4.0
pyparsing==3.0.7
python-dateutil==2.8.2
pytz==2022.6
scikit-learn==0.24.2
scipy==1.5.4
six==1.16.0
threadpoolctl==3.1.0
```

Additional dependencies:
If you are on the FIJI compute cluster you may load the following modules:
```
module load python/3.6.3
module load bedtools/2.25.0
module load samtools/1.8
module load meme/5.0.3 
```

Otherwise you can install these programs here:

[bedtools download](https://bedtools.readthedocs.io/en/latest/content/installation.html)

[samtools download](https://www.htslib.org/download/)

[meme suite download](https://meme-suite.org/meme/doc/install.html?man_type=web)

## Testing TF Profiler ##



## Running TF Profiler ##
### Required input files ###
There are three main input files required to run TF profiler. These are:
1. Genome fasta file (.fa)  -We used hg38 availabe [here](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) in our study.
2. Motif file (.meme) -We used the [HOCOMOCO](https://hocomoco11.autosome.org/) core set in our study. This meme file can be found in the assets folder.

And most importantly-

3. Bidirectional annotation file (.bed)
This bidirectional annotation file is derived from the experimental data. For this reason we highly recommend using this [Bidrectional-Flow pipeline](https://github.com/Dowell-Lab/Bidirectional-Flow) to run [Tfit](https://github.com/Dowell-Lab/Tfit). The precise annotation of the center of the bidirectional regions is essential for the success of the MD-score metric (used in TF Profiler).

Suggested flags for the [Bidrectional-Flow pipeline](https://github.com/Dowell-Lab/Bidirectional-Flow):
```
--tfit_prelim \
--prelim_process \
--tfit_split_model \
--savebidirs
```

Along with the requirements for the genome fasta (-g/--genome), the motif file (-m/--motifs) and the bidirecitonal annotation file (-a/--annotation), TF Profiler also requires a path to output the results (-o/--outdir) and a root name for all output files (-s/--sample).

### Assets ###
#For MD-score pre-gen and split (-sr flag?)


This folder contains MD-scores for all TFs in HOCOMOCO v11.
Each file contains the MD-scores based on the relative percentage of promoters from 10-50% for every 2% steps. The concentration of promoters impacts the simulated MD-scores (ie more promoters, higher MD-score for GC rich motifs, lower for AT rich motifs). For this reason simulations within 2% of the promoter content of the experimental/observed data can be used as an appropriate null hypothesis.
The files have names such as: promoter0.5_seed118_md_score_one_hit_per_region_quality; this is a 0.5 (or 50%) promoter containing experiment, the seed used to generate the MD-scores was 118 (same as the publication). A seed is needed as 1 million promoter and 1 million enhancer sequences were generated. These were subset down to 1 million total with relative proportions of enhancers and promoters. The sequences within these sets that were selected by seed 118 for reproducibility.
Finally the score was calculated by only keeping ONE motif instance per region. The region used for scoring was the highest QUALITY region as defined by the fimo output.

### Example run ###

### Additional Run Notes ###
For lower memory usage and speed use mononucleotide simulaiton
Default simulation is n=nregions provided (not a million)


Quickest run time and lowest mem usage will be to use prescan motifs and precalc distances
-q HOCOMOCOv11_full_HUMAN_genomic_motif_hits
-k 


## Help Message ##

```
[-h] [-v [True/False]] -a annotation.bed -o /full/path/to/output -s
         name_of_sample [-c int] [-r [True/False]] [-l [True/False]]
         [-d [True/False]] -g genome.fa [-e int] [-w int] [-n int] [-i int]
         [-q /full/path/to/pre-scanned/motifs] [-x [True/False]]
         [-k /full/path/to/pre-scanned/motifs] -m motif_database.meme
         [-t float] [-b background.csv] [-u [True/False]] [-p float]
         [-y [significance/all/<tf_name>/none]]

This work flow functions to assess active transcription factors in a ATAC-seq or PRO-seq data set. First, it calculates base content per position over a 3000bp window (given by a bed file). Subsequently, it generates sequences using a first or second order MM to obtain a simulated background sequence. Next, it can run FIMO to generate motif calls within the generated sequence (and the original sequence provided) or take a set of pre-scanned sequences. Next, it generates an MD score file for both the original sequence and the simulated sequence. Lastly, it outputs a text file with TF activation information along with various plots displaying this information.

optional arguments:
  -h, --help            show this help message and exit

General Arguments:
  -v [True/False], --verbose [True/False]
                        will generate verbose output file
  -a annotation.bed, --annotation annotation.bed
                        input bed file of bidirectionals or ATAC peaks, ends with .bed or .sorted.bed
  -o /full/path/to/output, --outdir /full/path/to/output
                        directory for output
  -s name_of_sample, --sample name_of_sample
                        name of the sample to be run (str)
  -c int, --cpus int    number of CPUs for multiprocessing. Default=1
  -r [True/False], --continue_run [True/False]
                        if a run was incomplete, add this flag to the original script to pick up where it left off.

Sequence Generation Arguments:
  -l [True/False], --mononucleotide_generation [True/False]
                        If False mononucleotide simulated sequences will not be generated. By default these sequences will also be scanned in FIMO unless the simulated_pre_scan flag is used pointing to a prescanned directory (Default: False)
  -d [True/False], --dinucleotide_generation [True/False]
                        If False dinucleotide simulated sequences will not be generated. By default these sequences will also be scanned in FIMO unless the simulated_pre_scan flag is used pointing to a prescanned directory (Default: True)
  -g genome.fa, --genome genome.fa
                        reference genome in fasta format. Genome index must be in the same directory (genome.fa.fai).
  -e int, --seed int    seed for initializing the random number generator. To set a specific seed -e int. The default is "None" so the seed is based on the clock.
  -w int, --window int  window to extract sequences. Default=1500
  -n int, --sequence_num int
                        number of simulated sequences to be generated. Default: set equal to number of peaks called in the annotation file provided
  -i int, --chrom_num int
                        number of simulated chromosomes. Default: approximately ~11.5 million bases per chromosome (10 percent of the average size of a human chromosome).

FIMO Arguments:
  -q /full/path/to/pre-scanned/motifs, --pre_scan /full/path/to/pre-scanned/motifs
                        directory containing pre-scanned motif hits in bed format over the whole genome obtained using the same fimo parameters (background and threshold) as the simulated dataset. If path is set experimental genome scan will be skipped and pre-scanned motif hits will be used instead.
  -x [True/False], --experimental_fimo [True/False]
                        will run fimo over only the annotated regions from the experimental dataset. True will increase run time. Recommended if you are only looking at a single dataset. If False, provide destination of the pre-scanned genome using the "--pre_scan" flag. Default: False.
  -k /full/path/to/pre-scanned/motifs, --simulated_pre_scan /full/path/to/pre-scanned/motifs
                        directory containing pre-scanned motif hits in bed format over a simulated dataset obtained using the same fimo parameters (background and threshold) as the experimental dataset. If path is set simulated scan will be skipped and pre-scanned motif hits will be used instead. NOTE: Make sure there are sufficient simulated sequences to match the total number of experimental sequences.
  -m motif_database.meme, --motifs motif_database.meme
                        meme file for TFs
  -t float, --threshold_fimo float
                        threshold for motifs called by fimo. Default: 1e-5
  -b background.csv, --background_file background.csv
                        background base composition of a given genome. This flag is HIGHLY recommended. See background options in useful files.

Scoring and Statistics Arguments:
  -u [True/False], --traditional_md [True/False]
                        will calculate md score using traditional method, number of hits within the small window (10 percent of the large window) divided by total hits in the large windiow. Default: True.
  -p float, --pval_cutoff float
                        Cutoff for defining enrichment/depetion. Default: 0.05
  -y [significance/all/<tf_name>/none], --plot_barcode [significance/all/<tf_name>/none]
                        This argument runs the barcode plotting module. There are 4 options for this flag. 1, significance or True plots barcodes only for significant TFs as defined by the p-value cutoff. 2 or all plots all barcodes from the MD-score file. 3,only functions if you provide <tf_name> and plots only the barcode of the TF specified. Note: if the string is ambigous then it will plot ALL TFs that fit the <tf_name> parameters within the md-score file. If no TF matches the string provided then nothing will be plotted. 4, none or False results in no barcodes plotted. Default: significance
```



## Example Output ##
Describe directory structure
Describe key outputs

## Contact Information ##

[Taylor Jones](tcaron360@gmail.com)
