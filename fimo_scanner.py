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
import numpy as np

######################################### Sequence Generator Main #########################################
def run_fimo_scanner(verbose, outdir, sample, cpus, motifs, 
                     threshold_fimo, background_file, continue_run, seq_type):
    
    fimo_dirs(verbose, outdir=outdir, seq_type=seq_type)
    motif_list = fimo_motif_names(verbose=verbose,motifs=motifs)
    if continue_run == True:
        scanned_motif_list=check_scanned_motifs(outdir=outdir, sample=sample, seq_type=seq_type)
        mset = set(motif_list)
        sset = set(scanned_motif_list)
        motif_list=list(mset.difference(sset))
        if verbose==True:
            print(str(len(sset)) + ' motifs were previously scanned.')
            print('Continuing scan for ' + str(len(motif_list)) + ' motifs.')
    
    if len(motif_list) == 0:
        print('There are no motifs to be scanned in this meme file.')
    else:
        if verbose == True: 
            print("---------Calculating GC Content of Motifs from Motif Database (.meme file)----------")
        get_gc(verbose, outdir=outdir, sample=sample, motifs=motifs, motif_list=motif_list, alphabet=['A','C','G','T'])
        if verbose == True: 
            print("---------FIMO Scan: Identifying Locations of Motif Hits -------------")
            print("This scan is for the " + seq_type + " genome.")
            print('Initializing ' + str(cpus) + ' cpus to run FIMO scan.')
            print('Start time: %s' % str(datetime.datetime.now()))
            start_time = int(time.time())
        seq_type=seq_type
        pool = multiprocessing.Pool(cpus)
        results = pool.map(partial(scanner, 
                                   inputs=[verbose, outdir, sample, motifs, threshold_fimo, background_file, seq_type]), motif_list)
        pool.close()
        pool.join()
        if verbose == True:
            print('FIMO scan complete at: %s' % str(datetime.datetime.now()))
        if verbose == True: 
            print("Bed file columns: ['sequence_name','start', 'stop','score', 'strand', 'motif_region_name']")
            stop_time = int(time.time())
            print("The scan for the " + seq_type + " genome is complete.")
            print('Stop time: %s' % str(datetime.datetime.now()))
            print("Total Run time :", (stop_time-start_time)/3600, " hrs")

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
            print(outdir + '/temp/' + seq_type + '_fimo_out exists.')
    
    if (path.exists(outdir + '/motifs') == False):
        os.system('mkdir -p ' + outdir + '/motifs')        
    try:
        os.system('mkdir -p ' + outdir + '/motifs/' + seq_type)
    except OSError:
        print ('Creation of the directory %s failed' % outdir + '/motifs/' + seq_type)
        sys.exit(1)
    else:
        if verbose == True:
            print(outdir + '/motifs/' + seq_type + ' exists.')
            
def fimo_motif_names(verbose, motifs):
    '''Extracts motif names from a MEME formatted motif database'''
    motif_list = list()
    with open(motifs) as F:
        for line in F:
            if 'MOTIF' in line:
                line = line.strip('\n').split()
                motif_name = line[-1]
                motif_list.append(motif_name)
    if verbose == True:
        print('There are ' + str(len(motif_list)) + " motifs in this meme file.")
    return motif_list

def check_scanned_motifs(outdir, sample, seq_type):
    scanned_motif_list = []
    motif_path = outdir + '/motifs/' + seq_type + '/*'
    motif_filenames = glob.glob(motif_path)
    motif_count = len(motif_filenames)
    for filename in motif_filenames:
        filename_no_path = filename.split('/')[-1]
        filename_no_path=filename_no_path.replace('.sorted.bed','')
        scanned_motif_list.append(filename_no_path)
    return scanned_motif_list

def get_gc(verbose, outdir, sample, motifs, motif_list, alphabet=['A','C','G','T']):
    '''Obtain a pssm model from a meme formatted database file.'''
    full_path = os.path.abspath(__file__)
    parent_path = os.path.dirname(full_path)
    meme_name=str(motifs)
    meme_name = meme_name.split('/')[-1]
    meme_name = meme_name.replace('.meme', '')    
    if (path.exists(parent_path + '/assets/'+ meme_name +'_motif_gc_percentage.txt') == False):
        gc_out={}
        for motif in motif_list:
            motif_hit = False
            PSSM = []
            with open(motifs,'r') as F:
                for line in F:
                    if 'MOTIF' in line:
                        if motif in line:
                            motif_hit = True
                        else:
                            motif_hit = False
                    elif motif_hit and 'URL' not in line and 'letter-probability' not in line and line != '\n':
                        acgt_probabilities = [float(x) for x in line.strip('\n').split()]
                        total_prob = sum(acgt_probabilities)
                        acgt_probabilities = [x/total_prob for x in acgt_probabilities] #Convert to probabilities
                        PSSM.append(acgt_probabilities)

                gc = 0
                for base in PSSM:
                    gc += base[alphabet.index('C')]
                    gc += base[alphabet.index('G')]

                gc = gc/float(len(PSSM))
            gc_out[motif] = motif,gc,(len(PSSM))

            df_gc_out = pd.DataFrame.from_dict(gc_out)
            df_gc_out = df_gc_out.transpose()
            df_gc_out.columns =['motif', 'percent_gc','motif_length']
    
            df_gc_out.to_csv(parent_path + '/assets/'+ meme_name +'_motif_gc_percentage.txt', header=None, index=False, sep='\t')
        else:
            if verbose == True:
                print('Motif length and GC content already calculated.')
                print('File: ' + parent_path + '/assets/'+ meme_name +'_motif_gc_percentage.txt')

def scanner(motif_list, inputs):
    verbose, outdir, sample, motifs, threshold_fimo, background_file, seq_type = inputs
    
    fasta = outdir + '/generated_sequences/' + sample + '_' + seq_type + '.fa'

    if verbose == True:
        fimo_verbosity = '--verbosity 1 ' ##change back to 2 eventually
    else:
        fimo_verbosity = '--verbosity 1 '

    if background_file is not None:
        os.system('fimo ' + fimo_verbosity + '--thresh ' + str(threshold_fimo) + ' --bfile ' + background_file + ' --motif ' + motif_list + ' --oc ' + outdir + '/temp/' + seq_type + '_fimo_out/' + motif_list + ' ' + motifs + ' ' + fasta)
    else:
        os.system('fimo ' + fimo_verbosity + '--thresh ' + str(threshold_fimo) + ' --motif ' + motif_list + ' --oc ' + outdir + '/temp/' + seq_type + '_fimo_out/' + motif_list + ' ' + motifs + ' ' + fasta)
        
    df = pd.read_csv(outdir + '/temp/' + seq_type + '_fimo_out/'+motif_list+'/fimo.tsv', sep ='\t')
    df.drop(df.tail(3).index,inplace=True)
    if df.empty == True:
        pass
    else:
        df = df.sort_values(by=['sequence_name', 'start']).reset_index()
        df['start'] = df['start'].astype(int) 
        df['stop'] = df['stop'].astype(int)
        df['count'] = (np.arange(len(df)))
        df['count'] = (df['count']+1).apply(str)
        df['motif_region_name'] = df['sequence_name'] + ';motif_' + df['count']
        df.drop(['count'], axis=1, inplace=True)
        df = df[['sequence_name','start', 'stop','score', 'strand', 'motif_region_name']]
        df.to_csv(outdir + '/motifs/' + seq_type + '/' + motif_list + '.sorted.bed', sep='\t', header=None, index=False)
