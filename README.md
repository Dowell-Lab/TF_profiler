# TF Profiler Overview #


Difference from paper- does not internally split enhancers and promoters (adding now)

# Installation and Requirements #
module load python/3.6.3
module load bedtools/2.25.0 (https://bedtools.readthedocs.io/en/latest/content/installation.html)
module load samtools/1.8 (https://www.htslib.org/download/)
module load meme/5.0.3 (https://meme-suite.org/meme/doc/install.html?man_type=web)

# Testing TF Profiler #



# Running TF Profiler #
## Required input files ##
Genome fasta file (.fa)
Motif file (.meme)
Bidirectional annotation file (.bed) -- describe how to do this (https://github.com/Dowell-Lab/Bidirectional-Flow built, tfit run https://github.com/Dowell-Lab/Tfit) --tfit_prelim \
--prelim_process \
--tfit_split_model \
--savebidirs \

## Required flags ##

## Example run ##

## Additional Run Notes ##
For lower memory usage and speed use mononucleotide simulaiton
Default simulation is n=nregions provided (not a million)


Quickest run time and lowest mem usage will be to use prescan motifs and precalc distances
-q HOCOMOCOv11_full_HUMAN_genomic_motif_hits
-k 


# Help Message #

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



# Example Output #
Describe directory structure
Describe key outputs

# Contact Information #

[Taylor Jones](tcaron360@gmail.com)
