######################################### Imports #########################################
import subprocess
from functools import partial
import multiprocessing
import pandas as pd
import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import time
import datetime
import os      
from os import path
import warnings
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

######################################### MD Score Main #########################################

def run_md_score(verbose, outdir, sample, window, cpus):
    if verbose == True: 
        print('--------------Beginning MD-Score Calculation- Simulated---------------')
        print('Initializing ' + str(cpus) + ' threads to calculate MD Scores.')
        start_time = int(time.time())
        print('Start time: %s' % str(datetime.datetime.now()))
    md_dirs(verbose, outdir=outdir, seq_type='simulated')
    
    if verbose==True:
        print('--------------Pulling in Annotation and Getting List of Motifs---------------')
    annotation_df = read_annotation(verbose=verbose, sample=sample, 
                                        outdir=outdir, seq_type='simulated')
    tf_list = get_tfs(verbose=verbose, outdir=outdir, seq_type='simulated') 
    
    for tf in tf_list:
        if verbose == True:
            print('Processing ' + tf + '.')
        motif_df = read_motif(verbose=verbose, outdir=outdir, seq_type='simulated', tf=tf)
        
        if verbose == True:
            print("---------Pulling Chromosomes From Annotation File and Motif File----------") 
        chrs = get_chrs(verbose=verbose, motif_df=motif_df, annotation_df=annotation_df, tf=tf)
        if verbose == True:
            print("---------Calculating Motif Distance Scores----------")
        pool = multiprocessing.Pool(cpus)
        motif_distance_dfs = pool.map(partial(get_distances, 
                               inputs=[verbose, motif_df, annotation_df]), chrs)
        motif_distance_df = pd.concat(motif_distance_dfs, axis=0)

        motif_score_df = calculate_motif_distance_score(verbose=verbose, outdir=outdir, sample=sample,
                                                        motif_distance_df=motif_distance_df, 
                                                        window=window, distance_weight=0.1,tf=tf, seq_type='simulated')

    if verbose == True:
        print("---------Simulated MD-Score Calculation Complete----------")
        stop_time = int(time.time())
        print('Stop time: %s' % str(datetime.datetime.now()))
        print('Total Run time :', (stop_time-start_time)/60, ' minutes')
        print('---------Compiling MD-scores---------')  
    pull_scores(verbose=verbose, outdir=outdir, seq_type='simulated', tf_list=tf_list, sample=sample)
      
    if verbose == True: 
        print('--------------Beginning MD-Score Calculation- Experimental---------------')
        print('Initializing ' + str(cpus) + ' threads to calculate MD Scores.')
        start_time = int(time.time())
        print('Start time: %s' % str(datetime.datetime.now()))
    md_dirs(verbose, outdir=outdir, seq_type='experimental')
    
    if verbose==True:
        print('--------------Pulling in Annotation and Getting List of Motifs---------------')
    annotation_df = read_annotation(verbose=verbose, sample=sample, 
                                        outdir=outdir, seq_type='experimental')
    tf_list = get_tfs(verbose=verbose, outdir=outdir, seq_type='experimental') 
    
    for tf in tf_list:
        if verbose == True:
            print('Processing ' + tf + '.')
        motif_df = read_motif(verbose=verbose, outdir=outdir, seq_type='experimental', tf=tf)
        
        if verbose == True:
            print("---------Pulling Chromosomes From Annotation File and Motif File----------") 
        chrs = get_chrs(verbose=verbose, motif_df=motif_df, annotation_df=annotation_df, tf=tf)
        if verbose == True:
            print("---------Calculating Motif Distance Scores----------")
        pool = multiprocessing.Pool(cpus)
        motif_distance_dfs = pool.map(partial(get_distances, 
                               inputs=[verbose, motif_df, annotation_df]), chrs)
        motif_distance_df = pd.concat(motif_distance_dfs, axis=0)

        motif_score_df = calculate_motif_distance_score(verbose=verbose, outdir=outdir, sample=sample,
                                                        motif_distance_df=motif_distance_df, 
                                                        window=window, distance_weight=0.1,tf=tf, seq_type='experimental')

    if verbose == True:
        print("---------Experimental MD-Score Calculation Complete----------")
        stop_time = int(time.time())
        print('Stop time: %s' % str(datetime.datetime.now()))
        print('Total Run time :', (stop_time-start_time)/60, ' minutes')
        print('---------Compiling MD-scores---------')
    pull_scores(verbose=verbose, outdir=outdir, seq_type='experimental', tf_list=tf_list, sample=sample)
            
######################################### MD Score Functions #########################################

def md_dirs(verbose, outdir, seq_type):
    if (path.exists(outdir + '/temp') == False):
        os.system('mkdir -p ' + outdir + '/temp') 
    try:
        os.system('mkdir -p ' + outdir + '/temp/' + seq_type + '_motif_scores')
    except OSError:
        print ('Creation of the directory %s failed' % outdir + '/temp/' + seq_type + '_motif_scores')
        sys.exit(1)
    else:
        if verbose == True:
            print('Successfully created the directory %s' % outdir + '/temp/' + seq_type + '_motif_scores')

def get_tfs(verbose, outdir, seq_type):
    tf_list = []
    motif_directory = os.fsencode(outdir + '/motifs/' + seq_type)
    for motif in os.listdir(motif_directory):
        if motif.endswith(b'.bed'):
            motif = motif.decode('utf-8')
            motif=motif.split('/')[-1].split('.')[-3]
            tf_list.append(motif)
        else:
            continue
    if verbose==True:
        print('There are ' + str(len(tf_list)) + ' motifs in the ' + seq_type + ' motif directory.')
    return tf_list

def read_motif(verbose, outdir, seq_type, tf):
    if verbose == True:
        print('Reading motif file ' + tf + ' and centering regions...')
    motif_file = outdir + '/motifs/' + seq_type + '/' + tf + '.sorted.bed'
    motif_df = pd.read_csv(motif_file, sep='\t', header=None)
    motif_df.columns = ['chr','start', 'stop', 'motif_id','score', 'strand','identifier', 'motif_region_name']
    motif_df['center'] = round((motif_df['stop'] + motif_df['start'])/2)
    motif_df['motif_length'] = motif_df['stop']-motif_df['start']
    return motif_df

def read_annotation(verbose, sample, outdir, seq_type):
    if verbose == True:
        print('Reading annotation file and centering regions...')
    annotation_file = outdir + '/annotations/' + sample + '_' + seq_type + '_window.bed'
    annotation_df = pd.read_csv(annotation_file, sep='\t', header=None)
    annotation_df.columns = ['chr', 'start', 'stop', 'region_name']
    annotation_df['center'] = round((annotation_df['stop'] + annotation_df['start'])/2)
    return annotation_df

def get_chrs(verbose, motif_df, annotation_df, tf):
    if verbose == True:
        print('Getting chromosome list...')
        print('Checking if chromosomes in ' + tf + ' match the annotation file.')
    all_chr_motif = list(motif_df['chr'].unique())
    all_chr_annotation = list(annotation_df['chr'].unique())
    chrs = list(set(all_chr_motif + all_chr_annotation))
    if verbose == True:
        print('There are ' + str(len(chrs)) + ' chromosomes in this dataset. They are:')
        print(chrs)
        #checking for odd ball chromosomes
        mtf = set(all_chr_motif)
        ann = set(all_chr_annotation)
        motif_not_annotation = mtf.difference(ann)
        annotation_not_motif = ann.difference(mtf)
        if len(motif_not_annotation) != 0:
            print(motif_not_annotation + ' chromosomes are unique to the motif file')
        if len(annotation_not_motif) != 0:
            print(annotation_not_motif + ' chromosomes are unique to the annotation file')
        else:
            print('The chromosomes match between the motif and annotation files.')
    return chrs

def findbedregions(row, single_chr_annotation):
    center = row['center']
    hits_df = single_chr_annotation[single_chr_annotation['start'] < center]
    hits_df = hits_df[hits_df['stop'] > center]
    distance = list(hits_df['center'] - row['center'])
    region_name = list(hits_df['region_name'])
    return [distance[0], region_name[0]]

def get_distances(chrs, inputs):
    verbose, motif_df, annotation_df = inputs
    if verbose == True:
        print('Calculating distances for chromosome ' + chrs + '.') 
    distances =[]
    single_chr_motif = motif_df[motif_df['chr'] == chrs]
    if verbose == True:
        before_drop = (len(single_chr_motif))
    single_chr_motif = single_chr_motif.sort_values(by='score', ascending=False)
    single_chr_motif = single_chr_motif.drop_duplicates(subset=['center'], keep='first')
    if verbose == True:
        after_drop = (len(single_chr_motif))
        dropped = before_drop - after_drop
        if dropped != 0:
            print(str(dropped) + ' hits were dropped on chromosome ' + chrs + ' due to duplicate hits.\n')  
        else:
            print('No hits were dropped on chromosome ' + chrs + '.\n') 
    single_chr_annotation = annotation_df[annotation_df['chr'] == chrs]

    single_chr_motif['dis_rn'] = single_chr_motif.apply(lambda row: findbedregions(row, single_chr_annotation), axis=1)
    single_chr_motif_distances = single_chr_motif[['dis_rn']]
    distances.append(single_chr_motif_distances)

    distance_df = pd.concat(distances, ignore_index=False)
    motif_distance_df = motif_df.merge(distance_df, left_index=True, right_index=True)
    motif_distance_df[['distance','region_name']] = pd.DataFrame(motif_distance_df.dis_rn.tolist(), index= motif_distance_df.index)
    motif_distance_df = motif_distance_df.drop(['dis_rn'], axis=1)
    return motif_distance_df

def no0s(row):
    if (row['distance'] == 0):
        return row['distance'] + 0.99999
    else:
        return row['distance']

def calculate_motif_distance_score(verbose, outdir, sample, motif_distance_df, window, distance_weight, tf, seq_type):
        rhf = motif_distance_df.groupby('region_name').count().reset_index()
        rhf = rhf[['region_name', 'chr']]
        rhf.columns = ['region_name', 'region_hit_frequency']
        motif_distance_df=motif_distance_df.merge(rhf, on='region_name')
        motif_distance_df['distance'] = motif_distance_df.apply(lambda row: no0s(row), axis=1)

        ##Calculating distance scores
        #this sets the weight of the distance score
        #this parameter adjusts the slope of the negative exponential function
        #ie small_window^exponent = distance weight
        #where small_window is 10% of the size of the large window
        small_window = window*0.1
        exponent = math.log(distance_weight, small_window)
        motif_distance_df['distance_score'] = (abs(motif_distance_df['distance']))**exponent
        motif_distance_df.to_csv(outdir + '/temp/' + seq_type + '_motif_scores/' + sample + '_' + tf + '.txt', sep='\t', index=False)
        if verbose ==  True:
            print('Successfully calculated distance scores for ' + tf + '.')
           
def pull_scores(verbose, outdir, seq_type, tf_list, sample):
    if (path.exists(outdir + '/results') == False):
        os.system('mkdir -p ' + outdir + '/results') 
    dd = {}
    if verbose == True:
        print('Reading scores for...')
    for tf in tf_list:
        if verbose == True:
            print(tf)
        df = pd.read_csv(outdir + '/temp/' + seq_type + '_motif_scores/' + sample + '_' + tf + '.txt', 
                         sep='\t')
        total_distance_score = float(df[['distance_score']].sum())
        motif_id = list(df[['motif_id']].loc[0])
        motif_id = str(motif_id[0])
        motif_hit_occurrences = int(len(df))
        dd[motif_id] = total_distance_score, motif_hit_occurrences
        out = pd.DataFrame.from_dict(dd, orient='index').reset_index()
        out.to_csv(outdir + '/results/' + sample + '_' + seq_type + '_md_scores.txt', 
                   sep="\t", header=['motif_id', 'total_distance_score', 'motif_hit_occurrences'], index=False)
    if verbose == True:
        print('Generation of '+ seq_type + '_md_scores.txt is complete.')