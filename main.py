######################################### Imports #########################################
import sys
from sequence_generator import run_sequence_generator
from fimo_scanner import run_fimo_scanner
from md_score import run_md_score

######################################### Run Main #########################################
def run(outdir, annotation, genome, sample, motifs, background_file=None, pre_scan=None, sequence_num=20000, window=1500, chrom_num=10, threshold_fimo=0.000001, cpus=1, seed=True, experimental_fimo=False, whole_genome_fimo=False, verbose=False):   
 #     #print("--------------Check needed modules---------------\npython/3.6.3\nbedtools/2.25.0\nsamtools/1.8\nmeme/5.0.3")
    #actually check here and if not present- exit
    if verbose == True:
        print('--------------Generating Sequences--------------')
    run_sequence_generator(verbose=verbose, outdir=outdir, sample=sample,
                           sequence_num=sequence_num, chrom_num=chrom_num, window=window,
                           genome=genome, annotation=annotation, seed=seed,
                           experimental_fimo=experimental_fimo, whole_genome_fimo=whole_genome_fimo)

    if verbose == True:
        print('--------------Running FIMO Scan--------------')
    run_fimo_scanner(verbose=verbose, outdir=outdir, sample=sample,
                           cpus=cpus, motifs=motifs, threshold_fimo=threshold_fimo,
                           background_file=background_file, genome=genome,
                           experimental_fimo=experimental_fimo, whole_genome_fimo=whole_genome_fimo)
    
    if verbose == True: 
        print('--------------Calculating MD-Scores--------------')
    run_md_score(verbose=verbose, outdir=outdir, sample=sample, window=window, cpus=cpus)
    
    print('done')
    sys.exit(0)
    

    
    
    
    
    