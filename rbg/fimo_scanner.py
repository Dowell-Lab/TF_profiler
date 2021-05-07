######################################### Imports #########################################
import os
from os import path
import sys
import time
import datetime
import subprocess
from functools import partial
import multiprocessing
import glob
import pandas as pd

######################################### Sequence Generator Main #########################################
def run_fimo_scanner(verbose, outdir, sample, cpus, motifs, 
                     threshold_fimo, background_file, genome, 
                     experimental_fimo, whole_genome_fimo):
    if verbose == True: 
        print("---------FIMO Scan: Identifying Locations of Motif Hits in Simulated Genome----------")
        print('Initializing ' + str(cpus) + ' cpus to run FIMO scan.')
        print('Start time: %s' % str(datetime.datetime.now()))
        start_time = int(time.time())
    motif_list = fimo_motif_names(verbose=verbose,motifs=motifs)
    fimo_dirs(verbose, outdir=outdir, seq_type='simulated')
    seq_type='simulated'
    pool = multiprocessing.Pool(cpus)
    results = pool.map(partial(scanner, 
                               inputs=[verbose, outdir, sample, motifs, threshold_fimo, background_file, genome, seq_type]), motif_list)
    if verbose == True: 
        stop_time = int(time.time())
        print('Stop time: %s' % str(datetime.datetime.now()))
        print("Total Run time :", (stop_time-start_time)/3600, " hrs")
        print("----------------Convert Simulated FIMO to Bed Format-----------------------")
    fimotobed(verbose, outdir, seq_type='simulated')
    if verbose == True: 
        print("Bed file columns: ['sequence_name','start', 'stop', 'motif_id','score', 'strand','identifier']")

    if experimental_fimo == True:
        if verbose == True: 
            print("---------FIMO Scan: Identifying Locations of Motif Hits in Experimental Genome----------")
            print('Start time: %s' % str(datetime.datetime.now()))
            start_time = int(time.time())
        fimo_dirs(verbose, outdir=outdir, seq_type='experimental')
        seq_type='experimental'
        pool.map(partial(scanner,
                         inputs=[verbose, outdir, sample, motifs, threshold_fimo, background_file, genome, seq_type]), motif_list)
        if verbose == True: 
            stop_time = int(time.time())
            print('Stop time: %s' % str(datetime.datetime.now()))
            print("Total Run time :", (stop_time-start_time)/3600, " hrs")
            print("----------------Convert Experimental FIMO to Bed Format-----------------------")
        fimotobed(verbose, outdir, seq_type='experimental')
    
    elif whole_genome_fimo == True:
        if verbose == True: 
            print("---------FIMO Scan: Identifying Locations of Motif Hits in Whole Genome----------")
            print('Start time: %s' % str(datetime.datetime.now()))
            start_time = int(time.time())
        fimo_dirs(verbose, outdir=outdir, seq_type='whole_genome')
        seq_type='whole_genome'
        pool.map(partial(scanner,
                         inputs=[verbose, outdir, sample, motifs, threshold_fimo, background_file, genome, seq_type]), motif_list)
        if verbose == True: 
            stop_time = int(time.time())
            print('Stop time: %s' % str(datetime.datetime.now()))
            print("Total Run time :", (stop_time-start_time)/3600, " hrs")
            print("----------------Convert Whole Genome FIMO to Bed Format-----------------------")
        fimotobed(verbose, outdir, seq_type='whole_genome')
    else:
        if verbose == True:
            print('No experimental motifs scanned.')
            print('Note: It is important that both the experimental motifs and simulated motifs are called with the same fimo parameters (meme file, cutoff and background file).')
        
######################################### FIMO Functions #########################################

def fimo_dirs(verbose, outdir, seq_type):
    if (path.exists(outdir + '/temp') == False):
        os.system('mkdir -p ' + outdir + '/temp') 
    try:
        os.system('mkdir -p ' + outdir + '/temp/' + seq_type + '_fimo_out')
    except OSError:
        print ('Creation of the directory %s failed' % outdir + '/temp/' + seq_type + '_fimo_out')
        sys.exit(1)
    else:
        if verbose == True:
            print('Successfully created the directory %s' % outdir + '/temp/' + seq_type + '_fimo_out')
    
    if (path.exists(outdir + '/motifs') == False):
        os.system('mkdir -p ' + outdir + '/motifs')        
    try:
        os.system('mkdir -p ' + outdir + '/motifs/' + seq_type)
    except OSError:
        print ('Creation of the directory %s failed' % outdir + '/motifs/' + seq_type)
        sys.exit(1)
    else:
        if verbose == True:
            print('Successfully created the directory %s' % outdir + '/motifs/' + seq_type)
            
def fimo_motif_names(verbose, motifs):
    '''Extracts motif names from a MEME formatted motif database'''
    motif_list = list()
    with open(motifs) as F:
        for line in F:
            if 'MOTIF' in line:
                line = line.strip('\n').split()
                motif_name = line[-1]
                motif_list.append(motif_name)
    first = [motif_list[0:1]]
    last =  [motif_list[-1:]]
    if verbose == True:
        print('There are ' + str(len(motif_list)) + " motifs in this meme file. " + str(first) + '...' + str(last))
    return motif_list

def scanner(motif_list, inputs):
    verbose, outdir, sample, motifs, threshold_fimo, background_file, genome, seq_type = inputs
    if seq_type == 'whole_genome':
        fasta = genome
    else:
        fasta = outdir + '/generated_sequences/' + sample + '_' + seq_type + '.fa'

    if verbose == True:
        fimo_verbosity = '--verbosity 2 '
    else:
        fimo_verbosity = '--verbosity 1 '

    if background_file is not None:
        os.system('fimo ' + fimo_verbosity + '--thresh ' + str(threshold_fimo) + ' --bfile ' + background_file + ' --motif ' + motif_list + ' --oc ' + outdir + '/temp/' + seq_type + '_fimo_out/' + motif_list + ' ' + motifs + ' ' + fasta)
    else:
        os.system('fimo ' + fimo_verbosity + '--thresh ' + str(threshold_fimo) + ' --motif ' + motif_list + ' --oc ' + outdir + '/temp/' + seq_type + '_fimo_out/' + motif_list + ' ' + motifs + ' ' + fasta)
        
def fimotobed(verbose, outdir, seq_type):
    motif_dirs = glob.glob(outdir + '/temp/' + seq_type + '_fimo_out/*')
    for motif_dir in motif_dirs:
        if path.isdir(motif_dir):
            motif_name = motif_dir.split('/')[-1]
            df = pd.read_csv('%s/fimo.tsv' % (motif_dir), sep ='\t')
            if verbose == True:
                print("Post-processing FIMO output for %s" % motif_name)
            df.drop(df.tail(3).index,inplace=True)
            if df.empty == True:
                if verbose == True:
                    print('Skipping ' + motif_name + ' -no motif hits.')
            else:
                df = df.sort_values('sequence_name').reset_index()
                df['identifier'] = df.apply(lambda row: identifier(row), axis=1)

                df = df[['sequence_name','start', 'stop', 'motif_id','score', 'strand','identifier']]
                df['start'] = df['start'].astype(int) 
                df['stop'] = df['stop'].astype(int)                 
                df.to_csv(outdir + '/motifs/' + seq_type + '/' + motif_name + '.sorted.bed', sep='\t', header=None, index=False)
                
def identifier(row):
    ident= str(row['sequence_name']) + ':' + str(row['start']) + '-' + str(row['stop']) + '_' + str(row['strand'] + ';' + str(row['motif_id']))
    return ident
