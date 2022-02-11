######################################### Imports #########################################
import sys
from os import path
import datetime
from sequence_generator import run_sequence_generator
from fimo_scanner import run_fimo_scanner
from distance_module import run_distance_calculation
from scoring_module import run_scoring_module
from statistics_module import run_statistics_module

######################################### Run Main #########################################
def run(verbose, outdir, sample, genome, annotation, 
        sequence_num, chrom_num, motifs, background_file, seed, cpus=1, window=1500,
        mononucleotide_generation=False, dinucleotide_generation=True, simulated_pre_scan=None,
        experimental_fimo=False, pre_scan=None, rerun=False,
        threshold_fimo='1e-5', traditional_md=True): 
################################ Sequence Generation ####################################
#     if verbose == True:
#         print('--------------Generating Sequences--------------')
#         print('Start time: %s' % str(datetime.datetime.now()))
#     run_sequence_generator(verbose=verbose, 
#                            outdir=outdir, sample=sample, genome=genome,
#                            annotation=annotation, cpus=cpus, seed=seed,
#                            sequence_num=sequence_num, chrom_num=chrom_num,
#                            window=window, 
#                            mononucleotide_generation=mononucleotide_generation,
#                            dinucleotide_generation=dinucleotide_generation,
#                            experimental_fimo=experimental_fimo, pre_scan=pre_scan, rerun=rerun)
#     if verbose == True:
#         print('--------------Sequence Generation Complete--------------')
#         print('Stop time: %s' % str(datetime.datetime.now()))
# ####################################### FIMO Scan ######################################### 
#     if verbose == True:
#         print('--------------Running FIMO Scan--------------')
#         print('Start time: %s' % str(datetime.datetime.now()))  
    
#     ###Dinucleotide simulated scan
#     if simulated_pre_scan is not None:
#         if verbose == True:
#             print('Not performing dinucleotide simulated sequence FIMO scan.')
#     elif dinucleotide_generation==True:
#         if verbose == True:
#             print('--------------Simulated Sequence Scan--------------')
#         run_fimo_scanner(verbose=verbose, outdir=outdir, sample=sample, cpus=cpus, motifs=motifs,
#                          threshold_fimo=threshold_fimo, background_file=background_file, rerun=rerun, 
#                          seq_type='simulated')
#     else:
#         if verbose == True:
#             print('Skipping dinucleotide simulated sequence FIMO scan.')
            
#     ###Mononucleotide simulated scan
#     if simulated_pre_scan is not None:
#         if verbose == True:
#             print('Not performing mononucleotide simulated sequence FIMO scan.')
#     elif mononucleotide_generation == True:
#         if verbose == True:
#             print('--------------Mononucleotide Simulated Sequence Scan--------------')
#         run_fimo_scanner(verbose=verbose, outdir=outdir, sample=sample, cpus=cpus, motifs=motifs,
#                          threshold_fimo=threshold_fimo, background_file=background_file, rerun=rerun,
#                          seq_type='mononucleotide_simulated')
#     else:
#         if verbose == True:
#             print('Skipping mononucleotide simulated sequence FIMO scan.')
    
#     ###Experimental scan
#     if pre_scan is not None:
#         if verbose==True:
#             print('Skipping experimental scan to use pre-scanned set of bedfiles.')
#             print('Please verify that the pre-scanned experimental set (ie a whole-genome scan) uses the same fimo parameters (background/threshold) as the simulated sets and motif hits are in the correct bedfile format.')           
#     elif experimental_fimo == True:
#         if verbose == True:
#             print('--------------Experimental Sequence Scan--------------')
#         run_fimo_scanner(verbose=verbose, outdir=outdir, sample=sample, cpus=cpus, motifs=motifs,
#                          threshold_fimo=threshold_fimo, background_file=background_file, rerun=rerun,
#                          seq_type='experimental')
#     else:
#         if verbose == True:
#             print('Skipping experimental sequence FIMO scan.')
   
#     if verbose == True:
#         print('--------------FIMO Scan(s) Complete--------------')
#         print('Stop time: %s' % str(datetime.datetime.now()))        
# ######################################## Distance Calculation ######################################### 
#     if verbose == True: 
#         print('--------------Calculating Distances--------------')
#     if mononucleotide_generation == True:
#         if verbose == True:
#             print('--------------Mononucleotide Simulated Distances--------------') 
#         run_distance_calculation(verbose=verbose, outdir=outdir,annotation=annotation, sample=sample, window=window, 
#                              cpus=cpus, seq_type='mononucleotide_simulated', pre_scan=pre_scan)                      
#     if dinucleotide_generation == True:
#         if verbose == True:
#             print('--------------Simulated Distances--------------')  
#         run_distance_calculation(verbose=verbose, outdir=outdir, annotation=annotation, sample=sample, window=window, 
#                              cpus=cpus, seq_type='simulated', pre_scan=pre_scan)
#     if experimental_fimo == True:
#         if verbose == True:
#             print('--------------Experimental Distances--------------')        
#         run_distance_calculation(verbose=verbose, outdir=outdir, annotation=annotation, sample=sample, window=window, 
#                              cpus=cpus, seq_type='experimental', pre_scan=pre_scan) 

######################################## Scoring Module ######################################### 
#     if verbose == True: 
#         print('--------------Scoring Motif Displacement--------------')
#     if mononucleotide_generation == True:
#         run_scoring_module(verbose=verbose, outdir=outdir,sample=sample, window=window, cpus=cpus,
#                            seq_type='mononucleotide_simulated')                      
#     if dinucleotide_generation == True:
#         run_scoring_module(verbose=verbose, outdir=outdir,sample=sample, window=window,cpus=cpus,
#                            seq_type='simulated')  
#     if experimental_fimo == True:
#         run_scoring_module(verbose=verbose, outdir=outdir,sample=sample, window=window,cpus=cpus,
#                            seq_type='experimental')

######################################## Statistics Module ######################################### 
    if verbose == True:
        print('--------------Determining Significance and Plotting--------------')
    run_statistics_module(verbose=verbose, outdir=outdir,sample=sample, traditional_md=traditional_md)








                      