######################################### Imports #########################################
import sys
from os import path
import datetime
from sequence_generator import run_sequence_generator
from sequence_generator import run_dinucleotide_generator
from sequence_generator import run_mononucleotide_generator
from sequence_generator import run_experimental_fimo_formatter
from fimo_scanner import run_fimo_scanner
from distance_module import run_distance_calculation
from scoring_module import run_scoring_module
from statistics_module import run_statistics_module

######################################### Run Main #########################################
def run(verbose, outdir, sample, genome, annotation, 
        sequence_num, chrom_num, motifs, background_file, seed, cpus=1, window=1500,
        mononucleotide_generation=False, dinucleotide_generation=True, simulated_pre_scan=None,
        experimental_fimo=False, pre_scan=None, continue_run=False,
        threshold_fimo='1e-5', traditional_md=True): 
#################################### Sequence Generation ####################################
    if verbose == True:
        print('--------------Generating Sequences--------------')
        print('Start time: %s' % str(datetime.datetime.now()))
    if dinucleotide_generation == False and mononucleotide_generation == False and simulated_pre_scan is not None:
        if verbose == True:
            print('No sequences were simulated.\n'+
                  'By default, if a simulated pre-scan is listed then no sequences are generated.\n'+
                  'To simulate sequences set dinucleotide_generation or mononucleotide_generation to True.')
    else:
        if verbose == True:
            print('Performing sequence generation initiation functions.')
        ls, rs_list, sequence_num, chrom_num = run_sequence_generator(verbose=verbose, 
                               outdir=outdir, sample=sample, genome=genome,
                               annotation=annotation, cpus=cpus, seed=seed,
                               sequence_num=sequence_num, chrom_num=chrom_num,
                               window=window)
    ### Dinucleotide Simulation ###
    if simulated_pre_scan is not None:
        if verbose == True:
            print('Not generating dinucleotide sequences. Pre-scan was provided.')
    elif dinucleotide_generation==True:
        if verbose == True:
            print('--------------Dinucleotide Sequence Generation--------------')
        run_dinucleotide_generator(verbose=verbose, 
                               outdir=outdir, sample=sample, sequence_num=sequence_num, chrom_num=chrom_num,
                               cpus=cpus, window=window, ls=ls, rs_list=rs_list)
    ### Mononucleotide Simulation ###    
    if simulated_pre_scan is not None:
        if verbose == True:
            print('Not generating mononucleotide sequences. Pre-scan was provided.')
    elif mononucleotide_generation==True:
        if verbose == True:
            print('--------------Mononucleotide Sequence Generation--------------')
        run_mononucleotide_generator(verbose=verbose, 
                               outdir=outdir, sample=sample, sequence_num=sequence_num, chrom_num=chrom_num,
                               cpus=cpus, window=window, ls=ls, rs_list=rs_list)
    ### Experimental Pre-Formatting ###
    if pre_scan is not None:
        if verbose == True:
            print('Not re-formatting sequences for fimo scan. Pre-scan was provided.')
    elif experimental_fimo == True:
        if verbose == True:
            print('--------------Formatting Experimental Sequences for FIMO Scan--------------')
        run_experimental_fimo_formatter(verbose=verbose, 
                               outdir=outdir, sample=sample, window=window, ls=ls)    
    if verbose == True:
        print('--------------Sequence Generation Complete--------------')
        print('Stop time: %s' % str(datetime.datetime.now()))
###################################### FIMO Scan ######################################### 
    if verbose == True:
        print('--------------Running FIMO Scan(s)--------------')
        print('Start time: %s' % str(datetime.datetime.now()))  
    ### Dinucleotide Scan ###
    if simulated_pre_scan is not None:
        if verbose == True:
            print('Not performing dinucleotide simulated sequence FIMO scan. Pre-scan was provided.')
    elif dinucleotide_generation==True:
        if verbose == True:
            print('--------------Simulated Sequence Scan--------------')
        run_fimo_scanner(verbose=verbose, outdir=outdir, sample=sample, cpus=cpus, motifs=motifs,
                         threshold_fimo=threshold_fimo, background_file=background_file, continue_run=continue_run, 
                         seq_type='simulated')  
    ### Mononucleotide Scan ###
    if simulated_pre_scan is not None:
        if verbose == True:
            print('Not performing mononucleotide simulated sequence FIMO scan. Pre-scan was provided.')
    elif mononucleotide_generation == True:
        if verbose == True:
            print('--------------Mononucleotide Simulated Sequence Scan--------------')
        run_fimo_scanner(verbose=verbose, outdir=outdir, sample=sample, cpus=cpus, motifs=motifs,
                         threshold_fimo=threshold_fimo, background_file=background_file, continue_run=continue_run,
                         seq_type='mononucleotide_simulated')
    ### Experimental Scan ###
    if pre_scan is not None:
        if verbose==True:
            print('Skipping experimental scan to use pre-scanned set of bedfiles.')
            print('Please verify that the pre-scanned experimental set (ie a whole-genome scan) uses the same fimo parameters (background/threshold) as the simulated sets.')           
    elif experimental_fimo == True:
        if verbose == True:
            print('--------------Experimental Sequence Scan--------------')
        run_fimo_scanner(verbose=verbose, outdir=outdir, sample=sample, cpus=cpus, motifs=motifs,
                         threshold_fimo=threshold_fimo, background_file=background_file, continue_run=continue_run,
                         seq_type='experimental')   
    if verbose == True:
        print('--------------FIMO Scan(s) Complete--------------')
        print('Stop time: %s' % str(datetime.datetime.now()))        
######################################## Distance Calculation ######################################### 
    if verbose == True:
        print('--------------Calculating Distances--------------')
        print('Start time: %s' % str(datetime.datetime.now()))  
    ### Simulated Distances ###
    ### If simulated_pre_scan is specified then it will be used by default ###
    ### If simulated_pre_scan is none then generated dinucleotide sequences will be run ###
    if simulated_pre_scan is not None or dinucleotide_generation == True:
        if verbose == True:
            print('--------------Simulated Distances--------------')  
        run_distance_calculation(verbose=verbose, outdir=outdir, annotation=annotation, sample=sample, window=window, 
                             cpus=cpus, seq_type='simulated', 
                                 pre_scan=pre_scan, simulated_pre_scan=simulated_pre_scan)            
    ### Mononucleotide Distances ###
    if simulated_pre_scan is not None:
        print('Simulated pre-scan is complete.')
    elif mononucleotide_generation == True:
        if verbose == True:
            print('--------------Mononucleotide Simulated Distances--------------') 
        run_distance_calculation(verbose=verbose, outdir=outdir,annotation=annotation, sample=sample, window=window, 
                             cpus=cpus, seq_type='mononucleotide_simulated', 
                                 pre_scan=pre_scan, simulated_pre_scan=simulated_pre_scan)  
    ### Experimental Distances ###
    ### If pre_scan is specified then it will be used by default ###
    ### If pre_scan is none then the newly scanned experimental sequences will be run ###
    if pre_scan is not None or experimental_fimo == True:
        if verbose == True:
            print('--------------Experimental Distances--------------')        
        run_distance_calculation(verbose=verbose, outdir=outdir, annotation=annotation, sample=sample, window=window, 
                             cpus=cpus, seq_type='experimental', 
                                 pre_scan=pre_scan, simulated_pre_scan=simulated_pre_scan) 
    if verbose == True:
        print('--------------Calculating Distances Complete--------------')
        print('Stop time: %s' % str(datetime.datetime.now()))   
####################################### Scoring Module ######################################### 
    if verbose == True:
        print('--------------Scoring Motif Displacement--------------')    
        print('Start time: %s' % str(datetime.datetime.now()))                    
    ### Simulated Score ### 
    ### Calls simulated distances that are calculated from prescan in outdir + '/distances/simulated/' ###
    ### Otherwise the distances in outdir + '/distances/simulated/' are from dinucleotide_simulated data ###
    if simulated_pre_scan is not None or dinucleotide_generation == True:
        if verbose == True:
            print('--------------Simulated MD Score--------------')  
            print('Calling distances from ' + outdir + '/distances/simulated/')
        run_scoring_module(verbose=verbose, outdir=outdir,sample=sample, window=window,cpus=cpus,
                           seq_type='simulated')
    ### Mononucleotide Distances ###
    if simulated_pre_scan is not None:
        print('Simulated MD Score Calculation is complete.')
    elif mononucleotide_generation == True:
        if verbose == True:
            print('--------------Mononucleotide MD Score--------------') 
        run_scoring_module(verbose=verbose, outdir=outdir,sample=sample, window=window, cpus=cpus,
                           seq_type='mononucleotide_simulated') 
    ### Experimental Score ### 
    ### Calls experimental distances that are calculated from prescan in outdir + '/distances/experimental/' ###
    ### Otherwise the distances in outdir + '/distances/experimental/' are from the newly scanned experimental sequences ###     
    if pre_scan is not None or experimental_fimo == True:
        if verbose == True:
            print('--------------Experimental MD Score--------------')
            print('Calling distances from ' + outdir + '/distances/experimental/')
        run_scoring_module(verbose=verbose, outdir=outdir, sample=sample, window=window, cpus=cpus,
                           seq_type='experimental')
    if verbose == True:
        print('--------------Motif Displacement Scoring Complete--------------')
        print('Stop time: %s' % str(datetime.datetime.now())) 
####################################### Statistics Module ######################################### 
    if verbose == True:
        print('--------------Determining Significance and Plotting--------------')
        print('Start time: %s' % str(datetime.datetime.now()))  
    run_statistics_module(verbose=verbose, outdir=outdir,sample=sample, traditional_md=traditional_md)
 
    if verbose == True:
        print('--------------RBG Workflow Complete--------------')
        print('Stop time: %s' % str(datetime.datetime.now())) 








                      