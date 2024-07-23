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
This bidirectional annotation file is derived from the experimental data. For this reason we highly recommend using this [Bidirectional-Flow pipeline](https://github.com/Dowell-Lab/Bidirectional-Flow) to run [Tfit](https://github.com/Dowell-Lab/Tfit). The precise annotation of the center of the bidirectional regions is essential for the success of the MD-score metric (used in TF Profiler).

Suggested flags for the [Bidirectional-Flow pipeline](https://github.com/Dowell-Lab/Bidirectional-Flow):
```
--tfit_prelim \
--prelim_process \
--tfit_split_model \
--savebidirs
```
For biological replicates we recommend using [muMerge](https://github.com/Dowell-Lab/mumerge) to generate a master annotation file before running through TF Profiler.


Along with the requirements for the genome fasta (-g/--genome), the motif file (-m/--motifs) and the bidirecitonal annotation file (-a/--annotation), TF Profiler also requires a path to output the results (-o/--outdir) and a root name for all output files (-s/--sample).

### Example run ###
A run command may look something like this:
```
python /path/to/TF_Profiler \
-g hg38.fa \
-a bidirectional_annotation.bed \
-o /path/to/output \
-s meaningful_name \
-v \
-c 12 \
-m HOCOMOCOv11_core_HUMAN_mono_meme_format.meme \
-b enhancer_background_flat \
-d True \
-x True
```
See run_rbg_hoco_example.sbatch in assets for addtional information regarding the sbatch set-up.

There is a --continue_run option if a submitted job fails.

Motif scanning is time consuming, and distance calculations are memory intensive. There are multiple flags to minimize recsource requirements:

1. For sequence simulation, it is recommended to use position specific dinucleotide frequencies (-d True). This indicates that the **next** base simulated is informed by the **previous** base in the sequence. Biological sequences (and TF motifs) tend to have higher order patterns rather than being completely independent of the surrounding bases (eg CpG islands). However, this simulation is more memory intensive than a mononucleotide position specific sequence (ie ignorant of the surrounding bases). Thus, the -l flag was created as a simplistic version of the model to reduce resource costs. To use this version of the expectation model set -d False and -l True.
   
2. The default number of simulated sequences is set to the number of regions in the bidirectional annotation file provided. For example, if you provided a bed file of 22,654 regions, exactly 22,654 regions will be simulated based on the position specific dinucleotide frequencies of that annotation file. The -n flag enables the user to vary the number of simulated sequences, with more sequences providing a more stable expectation value. As -n increases however, so do the resource costs. In the publication 1M sequences were simulated, this is resource intensive and likely unnecessary for most applications. Two additional related flags are -s (seed), where you can set the simulation seed for consistency and -i (chromosome number) where you can vary the number of simulated chromosomes fed into the fimo scan. The -i flag defaults to approximately ~11.5M bases per simulated chromosome. When the number of bases per simulated chromosome gets too small (or too large) fimo scan greatly slows.
   
3. The -q flag allows you to use motif hits that are pre-scanned genome wide. This flag enables the time consuming fimo motif scan to be performed **once** and can subsequently be re-applied to multiple datasets. To use this flag provide a complete path to a directory containing bed files that annotate motif hits genome wide. An example of how to generate these files can be found [here](https://github.com/Dowell-Lab/motif_scanning_and_distances). The bed format is partially derived from the fimo output text file, thus is slightly irregular. The columns follow this pattern: ['chr','start', 'stop', 'score', 'strand', 'motif_id'], where the motif_id must be **unique** to every region within the bed file. The motif pre-scans must match the TF PSSM names within the MEME file. Note: On FIJI these are found here: /scratch/Shares/dowell/tajo/hoco_flat/motifs, or can be found on [Zenodo](https://zenodo.org/records/12797230) in HOCOMOCOv11_full_HUMAN_genomic_motif_hits.zip.
   
4. The -k flag is similar to the -q flag but for simulated data sets. This flag enables the user to provide pre-calculated distance calculations for the expectation data. This allows both the time consuming motif scanning step and the memory intensive distance calculations to be performed **once** and reapplied to mutliple datasets. To use this flag provide a complete path to a directory containing distance.txt files from simulated data. A way to generate this data is run TF_Profiler once, and on subsequent runs use the -k flag to the /output/distances/simulated directory (see outputs explained below). In this case, be sure to check that the base composition bias is similar across the runs to ensure an appropriate expectation model. The columns of the distance table follow this pattern: ['region_id','motif_id','distance','distance_rank','quality_rank']. Note: On FIJI these are found here: /scratch/Shares/dowell/tajo/hoco_background_distances, or can be found on [Zenodo](https://zenodo.org/records/12797230) in distance_tables_simulated.zip. This version is simulated from separate promoter and enhancer populations, and combined (as done in the publication).

### Assets ###
The assets folder contains useful files for running TF Profiler.

In this folder various HOCOMOCO motif files (meme) can be found. Motif files are a requirement for TF Profiler. The files here are HOCOMOCO v11 in human (HOCOMOCOv11_full_HUMAN_mono_meme_format.meme) and in mouse (HOCOMOCOv11_full_MOUSE_mono_meme_format.meme), as well as core in human (HOCOMOCOv11_core_HUMAN_mono_meme_format.meme). The core file contains the highest confidence identified human motifs and was used for this publication. Each motif file has an associated *motif_gc_percentages.txt file. This file contains the calculate GC% from the related motif file.

There are two "transcription start site" (TSS) files, one for human (hg38_refseq_merge_1000bp_TSSs.bed) and one for mouse (mm10_refseq_merge_1000bp_TSSs.bed). In brief, these files are based on the RefSeq annotations where the annotated start site is windowed 300bp upstream and 700bp downstream in a strand specific manner. This gives a 1000bp window in which the start of an annotated transcript (NM/NR identifiers by RefSeq) is called a transcript associated "TSS". These files are used to designate a promoter (within 1000bp of a transcript associated TSS) or an enhancer (outside of 1000bp of a transcript associated TSS).

It's important to note that while the mouse data was not presented in the publication, TF Profiler has been run successfully on PRO/GRO-seq from mouse cells.

enhancer_background_flat is a particularly useful file for designating the expected background for the fimo scan (part of meme suite). In this file the background is "flat" ie equal proportions of A/T/C/G. This file is used in the fimo_scan command. It also is a template for altering the flat background expectation if desired.

md_scores.zip contains simulated MD-scores generated from 1M sequences. These MD-scores were generated using conditional probabilities to reflect the dinucleotide preferences within biological sequences. The probabilities for both enhancers and promoters can be found in the probabilities folder in assets and are depicted in [Supplemental Figure 4](https://www.biorxiv.org/content/10.1101/2024.03.15.585303v1), shown here:
![condProbs](https://github.com/user-attachments/assets/bf489eb8-f1b2-4892-99df-0876b36035cb)
These probabilities were derived from hg38_promoters.bed and hg38_enhancers.bed which can be found on [Zenodo](https://zenodo.org/records/12797230).

Each file in md_scores.zip contains the MD-scores calculated for all TFs in HOCOMOCO v11 (n=732). The difference across the files is the relative amount of promoters used within the sequence simulation. The concentration of promoters impacts the simulated MD-scores (ie more promoters, higher MD-score for GC rich motifs, lower for AT rich motifs). For this reason simulations within 2% of the experimentally observed promoter content can be used as an appropriate null hypothesis. The files have names such as: promoter0.22_seed118_md_score_one_hit_per_region_quality; this is a 0.22 (or 22%) promoter containing experiment, the seed used to generate the MD-scores was 118 (same as the publication). A seed is needed as 1 million promoter and 1 million enhancer sequences were generated. These were subset down to 1 million total with set proportions of enhancers and promoters. The sequences within these sets that were selected by seed 118 for reproducibility. Finally, the score was calculated by only keeping one motif instance per region. The region used for scoring was the nearest to the **center of a given bidirectional region**.

Lastly, run_rbg_hoco_example.sbatch is an example sbatch script for submission of FIJI.

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
  -k /full/path/to/pre-calc/distances, --simulated_pre_scan /full/path/to/pre-calc/distances
                        directory containing pre-calculated distances from a simulated dataset obtained using the same fimo parameters (background and threshold) as the experimental dataset. If path is set simulated scan will be skipped and pre-scanned motif hits will be used instead. NOTE: Make sure there are sufficient simulated sequences.
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

## Output Description ##
You will specify the output directory- from there multple subdirectories will be generated. "Sample" will be substituted for whatever rootname is specified by the required -s flag.

In the output you will find:

**1. annotations**
   - sample_experimental_prescan_windowed.bed - if using the -k flag. Windows and adds identifier to bidirectional annotation file.
   - sample_experimental_centered.bed and sample_experimental_window.bed - if using -x flag. Centers/windows and adds identifier to annotation bidirectional annotation file.
   - sample_promoters.bed - the bidirectional annotations that fall within promoter regions.
   - sample_simulated_centered.bed and sample_simulated_window.bed - if using the -d flag. Creates annotation files for simulated data after it has been formated onto chromosomes for the fimo scan.

**2. distances**
   - Contains two folders, one for each simulated and experimental data- within these folders there are distance tables. These are tables that contain all motif hits within +/-window (1500bp) of the center of all annotated bidirections.
   - Simulated will be missing if using the -k flag
Example distance output:

| region_id      | motif_id     | distance | distance_rank | quality_rank |
| -------------- | ------------ | -------- | ------------- | ------------ |
| chr1;region_16 | chr1;motif_1 | 56       | 1             | 1            |
| chr1;region_20 | chr1;motif_4 | 1303     | 3             | 1            | 
| chr1;region_20 | chr1;motif_5 | -742     | 1             | 3            |
| chr1;region_20 | chr1;motif_6 | -751     | 2             | 2            |

   - The region_id matches the ids in the provided annotation file.
   - The motif_id matches the 6th column of the motif bedfiles for whichever TF you are looking at from fimo_scan.
   - The distance is the distance from the center of the annotated region (ie the bidir) to the center of the motif for whatever TF you're looking at. All motif distances within +/-1500bp of the center of the bidirectional were calculated. A negative distance indicates upstream of the center of the bidirectional, whereas a positive distance indicates downstream of the bidirectional. Typically, we consider a motif hit within +/-150bp of the center of the bidirectional as "active."
   - Distance rank is for the case where there are 2+ motif hits within ONE bidirectional. The example here is chr1;region_20... For the distance rank 1 means closest to the center.
   - Quality rank is for the case where there are 2+ motif hits within ONE bidirectional. The example here is chr1;region_20... For the quality rank 1 means the highest confidence fimo call.

**3. generated_sequences**
   - sample_conditional_probabilites_givenX.tsv - the conditional probablities calculated from the underlying sequences of the provided bidirectional annotation file. One for each A/T/C/G.
   - sanple_position1_dinucleotide_probabilities.tsv - probabilities for the 16 possible dinucleotide combinations in positions 1 and 2 to set the seed for the remaining base generation.
   - sample_mononucleotide_probabilites.tsv - the position specific mononucleotide probabilities calculated from the underlying sequences of the provided bidirectional annotation file. See sample_single_position_BaseDistribution.png and sample_single_position_SmoothedBaseDistribution.png in plots for the plots related to this data. There should be a sharp uptick in GC composition around position 0.
   - dinucleotide_mononucleotide_probabilites.tsv - the position specific mononucleotide probablitilies back calculated from the simulated sequences (should be roughly equal to sample_mononucleotide_probabilites.tsv). See dinucleotide_single_position_BaseDistribution.png and dinucleotide_single_position_SmoothedBaseDistribution.png in plots for the plots related to this data.
   - sample_simulated.chrom.sizes - chromosome size file for simulated sequences.
   - sample_simulated.fa - fasta file of simulated sequences formated onto chromosomes for motif scanning.
   - sample_simulated.fa.fai - index for simulated fasta file.

**4. motifs**
   - Contains two folders, one for each simulated and experimental data- within these folders there are motif bed files. These are annotations of every motif hit within the windowed input regions.
   - This directory will be missing if using the -k and -q flags are used.
Example motif bed file:

| chr      | start     | stop | score | strand | motif_id |
| -------------- | ------------ | -------- | ------------- | ------------ |------------ |
| chr1 | 22299 | 22318       | 12.7727          | +            | chr1;motif_1 |
| chr1 | 31718 | 31737     | 12.7576             | +          | chr1;motif_2 |
| chr1 | 33517 | 33536     | 12.7727             | -            | chr1;motif_3 |
| chr1 | 43588 | 43607     | 13.1818             | +            | chr1;motif_4 |

   - The chr/start/stop are genomic locations for motif hits in the genome.
   - Score is the score attributed to how well the motif matches that genomic location assigned by the fimo scan.
   - Strand the motif hit was on (+/-).
   - The motif_id matches the 2nd column of the distance table for a given TF.

**5. plots**
   - *single_position_BaseDistribution as previously described shows the position specific mononucleotide base distributions from both experimental and simulated data. Both are plotted as a qc metric and should be roughly identical, with a sharp increase in GC composition around position 0.
   - sample_probability_givenX_BaseDistribution shows plots related to the position specific conditional probabilities for A/T/C/G.
   - rbg_significance, rbg_gc_content, rbg_elliptic_fit all plot the observed MD-score (y-axis) vs the expected MD-score (x-axis). The only difference in these three plots is the coloration. Significance marks enrichment/depletion. GC_content marks the range of GC content of the TF motifs. Elliptical fit shows which points (inliers) are fit to the linear regression, and thus which residuals are fit to the normal distribution to assess significance.
     
**6. scores**
   - There are two MD-score output files, simulated_traditional_md_score.txt and experimental_traditional_md_score.txt. These files contain intermediate information regarding the calculation of MD-scores.
   - **md_score_experimental_vs_simulated_significance.txt** is the main results file.
Example results file:

| tf      | percent_gc     | md_score_exp | md_score_sim | elliptic_outlier | pval |significance|
| -------------- | ------------ | -------- | ------------- | ------------ | -|-|
| NRF1_HUMAN.H11MO.0.A | 0.69 | 0.40       | 0.47             | -1       | 9.1e-78 | Depleted|
| ELK4_HUMAN.H11MO.0.A | 0.63 | 0.41       | 0.41             | -1       | 5.8e-68 | Depleted|
| FOXC2_HUMAN.H11MO.0.D | 0.37 | 0.09       | 0.08            | 1       | 0.13 | Not Significant|
| KAISO_HUMAN.H11MO.1.A | 0.61 | 0.38       | 0.16             | -1       | 2.3e-52 | Enriched|

   - tf is the TF name (from the MEME file)
   - percent_gc is the GC content of the motif derived from the MEME file.
   - md_score_exp/md_score_sim are the MD-scores calculated from the experimental/observed data, and from the simulated/expected data.
   - elliptical outlier is 1/-1 (True/False) where all values that are FALSE (ie are not outliers) are used to generate the linear regression. The residuals of the FALSE elliptical outlier points are fit to a normal distribution. It's this distribution that is contrasted with all residuals to attribute significance values as shown in column 6-pval.
   - The final column, significance, contains three possible values, Depleted (the experimental/observed MD-score is lower than expectation), Enriched (the experimental/observed MD-score is lower than expectation) or Not Significant. The user can specify the p-value cut-off with the -p flag. The default is 0.05.

**7. temp** - contains intermediate files. It is recommended to remove temp right away.

## Contact Information ##

[Taylor Jones](tcaron360@gmail.com)
