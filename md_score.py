######################################### Imports #########################################
import pandas as pd
import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import time
import datetime
import os
from os import path
######################################### MD Score Main #########################################
def run_md_score(verbose, outdir, sample, window):
    if verbose == True: 
        print('--------------Beginning MD-Score Calculation---------------')
#         print('Initializing ' + str(cpus) + ' threads to calculate MD Scores.')
        start_time = int(time.time())
        print('Start time: %s' % str(datetime.datetime.now()))
    md_dirs(verbose, outdir=outdir, seq_type='experimental')
    
    if verbose==True:
        print('--------------Normalizing Motif Quality Scores---------------')
    ##func to read in each motif file present in a dir
    motif_file = '/Users/tajo5912/rbg/squid/motifs/experimental/SNAI2_M09145_1_mini.sorted.bed'
    motif_df = calculate_motif_quality_score(motif_file=motif_file)
    
    if verbose == True:
        print("---------Calculating Motif Distance Scores----------")
    annotation_file = '/Users/tajo5912/rbg/squid/annotations/squid_experimental_window_mini.bed'
    distance_dict = get_distances(verbose=verbose, motif_df=motif_df, 
                                  annotation_file=annotation_file) #I want to thread this by chromosome
    motif_score_df = calculate_motif_distance_score(verbose=verbose, outdir=outdir, sample=sample,
                                                    distance_dict=distance_dict, 
                                                    window=window, distance_weight=0.1,seq_type='experimental')
    
    if verbose == True:
        print("---------Calculating Region Score----------")
    region_score_df = calculate_region_score(verbose=verbose, outdir=outdir, sample=sample, 
                                             motif_score_df=motif_score_df,seq_type='experimental')
    
    if verbose == True:
        print("---------MD-Score Calculation Complete----------")
        stop_time = int(time.time())
        print('Stop time: %s' % str(datetime.datetime.now()))
        print('Total Run time :', (stop_time-start_time)/60, ' minutes')

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
            
    try:
        os.system('mkdir -p ' + outdir + '/temp/' + seq_type + '_region_scores')
    except OSError:
        print ('Creation of the directory %s failed' % outdir + '/temp/' + seq_type + '_region_scores')
        sys.exit(1)
    else:
        if verbose == True:
            print('Successfully created the directory %s' % outdir + '/temp/' + seq_type + '_region_scores')

def calculate_motif_quality_score(motif_file):
    motif_df = pd.read_csv(motif_file, sep='\t', header=None)
    motif_df.columns = ['chr','start', 'stop', 'motif_id','score', 'strand','identifier', 'motif_region_name']
    motif_df['center'] = round((motif_df['stop'] + motif_df['start'])/2)
    motif_df['motif_length'] = motif_df['stop']-motif_df['start']
    motif_df['normalized_quality_score'] = motif_df['score']/((motif_df['score']).max())
    return motif_df

def get_distances(verbose, motif_df, annotation_file):
    annotation_df = pd.read_csv(annotation_file, sep='\t', header=None)
    annotation_df.columns = ['chr', 'start', 'stop', 'region_name']
    annotation_df['center'] = round((annotation_df['stop'] + annotation_df['start'])/2)
    
    distance_dict = {}
    for index, motif_row in motif_df.iterrows():
        for index, annotation_row in annotation_df.iterrows():
            if (motif_row['center'] >= annotation_row['start']) & (motif_row['center'] <= annotation_row['stop']):
                distance = (annotation_row['center'] - motif_row['center']) #motif is the center or annotation is the center?--doesn't matter for the score but does for the plots
                distance_dict[motif_row['motif_region_name']] = annotation_row['region_name'], distance, motif_row['normalized_quality_score']
    if verbose==True:
        print(distance_dict)
    return distance_dict

def calculate_motif_distance_score(verbose, outdir, sample, distance_dict, window, distance_weight, seq_type):
    motif_score_df = pd.DataFrame.from_dict(distance_dict)
    motif_score_df= motif_score_df.transpose().reset_index()
    cols=['motif_region_name', 'region_name', 'distance', 'normalized_quality_score']
    motif_score_df.columns = cols
    rhf=motif_score_df.groupby('region_name').count().reset_index()
    rhf=rhf[cols[0:2]]
    rhf.columns = ['region_hit_frequency', 'region_name']
    motif_score_df=motif_score_df.merge(rhf, on='region_name')
    
    ##Calculating distance scores
    #this sets the weight of the distance score
    #this parameter adjusts the slope of the negative exponential function
    #ie small_window^exponent = distance weight
    #where small_window is 10% of the size of the large window
    small_window = window*0.1
    exponent = math.log(distance_weight, small_window)
    motif_score_df['normalized_distance_score'] = (abs(motif_score_df['distance']))**exponent
    motif_score_df.to_csv(outdir + '/temp/' + seq_type + '_motif_scores/' + 'TFTEST.txt', sep='\t', index=False)
    return motif_score_df

def calculate_region_score(verbose, outdir, sample, motif_score_df, seq_type): 
    motif_score_df['motif_hit_score'] = motif_score_df['normalized_quality_score']*motif_score_df['normalized_distance_score']*(1/motif_score_df['region_hit_frequency'])
    region_score_df = motif_score_df[['region_name', 'motif_hit_score']]
    region_score_df = region_score_df.groupby('region_name').sum().reset_index()
    region_score_df.columns = ['region_name', 'region_score']
    region_score_df.to_csv(outdir + '/temp/' + seq_type + '_region_scores/' + 'TFTEST.txt', sep='\t', index=False)
    return region_score_df