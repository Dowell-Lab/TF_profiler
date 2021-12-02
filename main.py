######################################### Imports #########################################
import sys
from os import path
import datetime
from sequence_generator import run_sequence_generator
from fimo_scanner import run_fimo_scanner
from distance_module import run_distance_calculation

######################################### Run Main #########################################
def run(verbose, outdir, sample, genome, annotation, 
        sequence_num, chrom_num, motifs, background_file, seed, cpus=1, window=1500,
        mononucleotide_generation=False, dinucleotide_generation=True, skip_simulated_fimo=False,
        experimental_fimo=False, pre_scan=None,
        threshold_fimo='1e-5'): 
#################################### Sequence Generation ####################################
    if verbose == True:
        print('--------------Generating Sequences--------------')
        print('Start time: %s' % str(datetime.datetime.now()))
    run_sequence_generator(verbose=verbose, 
                           outdir=outdir, sample=sample, genome=genome,
                           annotation=annotation, cpus=cpus, seed=seed,
                           sequence_num=sequence_num, chrom_num=chrom_num,
                           window=window, 
                           mononucleotide_generation=mononucleotide_generation,
                           dinucleotide_generation=dinucleotide_generation,
                           experimental_fimo=experimental_fimo, pre_scan=pre_scan)
    if verbose == True:
        print('--------------Sequence Generation Complete--------------')
        print('Stop time: %s' % str(datetime.datetime.now()))
######################################### FIMO Scan ######################################### 
    if verbose == True:
        print('--------------Running FIMO Scan--------------')
        print('Start time: %s' % str(datetime.datetime.now()))  
    
    ###Dinucleotide simulated scan
    if skip_simulated_fimo == True:
        if verbose == True:
            print('Not performing dinucleotide simulated sequence FIMO scan.')
    elif dinucleotide_generation==True:
        run_fimo_scanner(verbose=verbose, outdir=outdir, sample=sample, cpus=cpus, motifs=motifs,
                         threshold_fimo=threshold_fimo, background_file=background_file, 
                         seq_type='simulated')
    else:
        if verbose == True:
            print('Skipping dinucleotide simulated sequence FIMO scan.')
            
    ###Mononucleotide simulated scan
    if skip_simulated_fimo == True:
        if verbose == True:
            print('Not performing mononucleotide simulated sequence FIMO scan.')
    elif mononucleotide_generation == True:
        run_fimo_scanner(verbose=verbose, outdir=outdir, sample=sample, cpus=cpus, motifs=motifs,
                         threshold_fimo=threshold_fimo, background_file=background_file, 
                         seq_type='mononucleotide_simulated')
    else:
        if verbose == True:
            print('Skipping mononucleotide simulated sequence FIMO scan.')
    
    ###Experimental scan
    if pre_scan is not None:
            if verbose==True:
                print('Skipping experimental scan to use pre-scanned set of bedfiles.')
                print('Please verify that the pre-scanned experimental set (ie a whole-genome scan) uses the same fimo parameters (background/threshold) as the simulated sets and motif hits are in bedfile format.')           
    elif experimental_fimo == True:
        run_fimo_scanner(verbose=verbose, outdir=outdir, sample=sample, cpus=cpus, motifs=motifs,
                         threshold_fimo=threshold_fimo, background_file=background_file, 
                         seq_type='experimental')
    else:
        if verbose == True:
            print('Skipping experimental sequence FIMO scan.')
   
    if verbose == True:
        print('--------------FIMO Scan(s) Complete--------------')
        print('Stop time: %s' % str(datetime.datetime.now()))        
######################################## Distance Calculation ######################################### 
    if verbose == True: 
        print('--------------Calculating Distances --------------')
    if mononucleotide_generation == True:
        run_distance_calculation(verbose=verbose, outdir=outdir,annotation=annotation, sample=sample, window=window, 
                             cpus=cpus, seq_type='mononucleotide_simulated', pre_scan=pre_scan)                      
    if dinucleotide_generation == True:
        run_distance_calculation(verbose=verbose, outdir=outdir, annotation=annotation, sample=sample, window=window, 
                             cpus=cpus, seq_type='simulated', pre_scan=pre_scan)
    if experimental_fimo == True:
        run_distance_calculation(verbose=verbose, outdir=outdir, annotation=annotation, sample=sample, window=window, 
                             cpus=cpus, seq_type='experimental', pre_scan=pre_scan) 
                      
                      