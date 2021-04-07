# Outside imports
import sys
import time
import datetime
import os
import threading
import multiprocessing
from functools import partial
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Local imports
from qc_and_keys import *
from get_sequences import *
from counting import *
from base_plot import *
from file_organizer import *
from fimo_scanner import *
from md_score import *
from fasta_and_bed_processing import *
from sequence_generator import *

#maybe move output file to a log file to track the course of the run

def run(outdir, annotation, genome, sample, motifs, background_file=None, pre_scan=None, sequence_num=20000, window=1500, chrom_num=10, cutoff_fimo=0.000001, threads=1, seed=True, experimental_fimo=False, whole_genome_fimo=False):   
    print ("Log file for rbg work-flow...")    
    print("--------------Check needed modules---------------\npython/3.6.3\nbedtools/2.25.0\nsamtools/1.8\nmeme/5.0.3")
    #actually check here and if not present- exit
    print("--------------Making Directories---------------")
    mkdirectories(outdir, experimental_fimo, whole_genome_fimo)
    
    print("--------------Expanding Windows---------------")
    fasta_bed_process=FastaBedProcess(outdir, sample, sequence_num, chrom_num, window, genome, annotation)
    fasta_bed_process.get_windows()
    
    print("-------------Extracting Sequences-------------")
    seqs = ExtractSequences(genome, sample, outdir)
    seqs.get_sequences()

    print("-----------Counting Base Content and Plot Generation-------------")
    listseq = ListSequences(outdir, sample)
    window_seq = listseq.list_sequences()
    base_plot(window_seq, sample, outdir, sample, "All", int(window))

    print("-----------Start Simulated Sequence Generation-------------")
    start_time = int(time.time())
    print('Start time: %s' % str(datetime.datetime.now()))

    if seed == 0:
        np.random.seed()
        print("-----------Setting Seed for the Sequence Generator.-----------")
        print("Using random seed.")
        print("WARNING: Sequences generated using a random seed can not be reproduced.")
    elif seed is True:
        np.random.seed(global_seed())
    else:
        np.random.seed(seed)
        print("-----------Setting Seed for the Sequence Generator.-----------")
        print("User defined seed:", seed)

    print("-----------Reading Per-Base Sequence Frequency----------------")
    tsv = outdir + '/generated_sequences/' + sample + '_base_content.tsv'
    print('Base Content File: ' + tsv)
    position_prob = sequence_input(tsv)

    print("-----------Generating Sequences-------------------------------")
    sequences_generating = sequence_generator(bases=['A', 'T', 'G', 'C'], sequence_num=sequence_num, position_feq=position_prob)

    print("-----------Writing Sequences to Fasta File--------------------")
    write_fasta(sequences_generating, sample, outdir, int(window), sequence_num, chrom_num)
    print('Simulated Fasta File: ' + outdir + '/generated_sequences/' + str(sample) + "_simulated.fa")
    print("----------------Simulated Sequence Generation Complete------------------------------")
    stop_time = int(time.time())
    print('Stop time: %s' % str(datetime.datetime.now()))
    print("Total Run time :", stop_time-start_time, " seconds")
    
    print("-----------Indexing Fasta and Generating Chrm Sizes File--------------------")
    fasta_bed_process.index_and_chrm_sizes(seq_type='simulated')
    
    print("-----------Generating Bedfile for Simulated Sequences--------------------")
    fasta_bed_process.generate_bed(seq_type='simulated')
    
    print("---------FIMO Scan: Identifying Locations of Motif Hits in Simulated Genome----------")
    motif_list = fimo_motif_names(motifs)
    fimoscan= FIMOScan(motif_list, motifs, cutoff_fimo, background_file, genome, outdir, sample)

    print('Initializing ' + str(threads) + ' threads to run FIMO scan.')
    start_time = int(time.time())
    print('Start time: %s' % str(datetime.datetime.now()))
    pool = multiprocessing.Pool(threads)
    results = pool.map(fimoscan.fimo_scanner, motif_list) ##NEED TO PUT SEQ_TYPE HERE 'simulated'
    pool.close()
    stop_time = int(time.time())
    print('Stop time: %s' % str(datetime.datetime.now()))
    print("Total Run time :", (stop_time-start_time)/3600, " hrs")

    print("----------------Convert FIMO to Bed Format-----------------------")
    fimotobed(outdir, seq_type='simulated')
    print("Bed file columns: ['sequence_name', 'start', 'stop', 'identifier', 'motif_id', 'score', 'strand', 'matched_sequence']")
    
    if experimental_fimo == True:
        print("---------Formating Experimental Fasta for FIMO Scan----------")
        fasta_bed_process.reformat_expt_fasta()
        
        print("-----------Indexing Fasta and Generating Chrm Sizes File--------------------")
        fasta_bed_process.index_and_chrm_sizes(seq_type='experimental')

        print("-----------Generating Bedfile for Experimental Sequences--------------------")
        fasta_bed_process.generate_bed(seq_type='experimental')
        
        print("---------FIMO Scan: Identifying Locations of Motif Hits in Experimental Genome---------")
        print('Initializing ' + str(threads) + ' threads to run FIMO scan.')
        start_time = int(time.time())
        print('Start time: %s' % str(datetime.datetime.now()))
        pool = multiprocessing.Pool(threads)
        results = pool.map(fimoscan.fimo_scanner, motif_list) ##NEED TO PUT SEQ_TYPE HERE 'experimental'
        pool.close()
        pool.join()
        stop_time = int(time.time())
        print('Stop time: %s' % str(datetime.datetime.now()))
        print("Total Run time :", (stop_time-start_time)/3600, " hrs")
        
        print("----------------Convert FIMO to Bed Format-----------------------")
        fimotobed(outdir, seq_type='experimental')
        print("Bed file columns: ['sequence_name', 'start', 'stop', 'identifier', 'motif_id', 'score', 'strand', 'matched_sequence']")
        
    elif whole_genome_fimo == True:
        print("-----------Generating Chrm Sizes File--------------------")
        fasta_bed_process.chrm_sizes_whole_genome()
        
        print("---------FIMO Scan: Identifying Locations of Motif Hits in Whole Genome----------")
        print('Initializing ' + str(threads) + ' threads to run FIMO scan.')
        start_time = int(time.time())
        print('Start time: %s' % str(datetime.datetime.now()))
        pool = multiprocessing.Pool(threads)
        results = pool.map(fimoscan.whole_genome_scanner, motif_list)
        pool.close()
        pool.join()
        stop_time = int(time.time())
        print('Stop time: %s' % str(datetime.datetime.now()))
        print("Total Run time :", (stop_time-start_time)/3600, " hrs")
        
        print("----------------Convert FIMO to Bed Format-----------------------")
        fimotobed(outdir, seq_type='whole_genome')
        print("Bed file columns: ['sequence_name', 'start', 'stop', 'identifier', 'motif_id', 'score', 'strand', 'matched_sequence']")
        
    else:
        print('No experimental motifs scanned. Note: It is important that both the experimental motifs and simulated motifs are called with the same fimo parameters (meme file, cutoff and background file).')
        
#     print("---------Calculating the MD Scores for the Simulated Genome----------")
#     print('Initializing ' + str(threads) + ' threads to calculate MD Scores.')
#     start_time = int(time.time())
#     print('Start time: %s' % str(datetime.datetime.now()))
#     mdcalc(window, threads, outdir, sample)
#     stop_time = int(time.time())
#     print('Stop time: %s' % str(datetime.datetime.now()))
#     print("Total Run time :", (stop_time-start_time)/60, " minutes")

#     if experimental_fimo == True:   
#         print("---------Calculating the MD Scores for the Experimental Genome----------")
#         print('Initializing ' + str(threads) + ' threads to calculate MD Scores.')
#         start_time = int(time.time())
#         print('Start time: %s' % str(datetime.datetime.now()))       
#         mdcalc_expt(window, threads, outdir, sample)
#         stop_time = int(time.time())
#         print('Stop time: %s' % str(datetime.datetime.now()))
#         print("Total Run time :", (stop_time-start_time)/60, " minutes")           
#     elif whole_genome_fimo == True:
#         print("---------Calculating the MD Scores for the Whole Genome----------")
#         print('Initializing ' + str(threads) + ' threads to calculate MD Scores.')
#         start_time = int(time.time())
#         print('Start time: %s' % str(datetime.datetime.now()))         
#         mdcalc_whole_genome(annotation, window, threads, outdir, sample, genome)
#         stop_time = int(time.time())
#         print('Stop time: %s' % str(datetime.datetime.now()))
#         print("Total Run time :", (stop_time-start_time)/60, " minutes")     
#     elif pre_scan != None:
#         print("---------Calculating the MD Scores from Motif Database Provided----------")
#         print("The motifs generated here MUST be from the same genome provided by the -g flag")
#         print("-----------Generating Chrm Sizes File--------------------")
#         fasta_bed_process.chrm_sizes_whole_genome()
#         print('Initializing ' + str(threads) + ' threads to calculate MD Scores.')
#         start_time = int(time.time())
#         print('Start time: %s' % str(datetime.datetime.now()))         
#         mdcalc_prescan(annotation, window, threads, outdir, sample, genome, pre_scan)
#         stop_time = int(time.time())
#         print('Stop time: %s' % str(datetime.datetime.now()))
#         print("Total Run time :", (stop_time-start_time)/60, " minutes")        
#     else:
#         print('No experimental motifs provided. Further downstream analyses halted')
#         sys.exit(0)

#     sys.exit(0)
    
    
    
    
    
    
    