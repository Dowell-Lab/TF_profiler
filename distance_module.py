######################################### Imports #########################################
import subprocess
from functools import partial
import multiprocessing
import pandas as pd
import math
import numpy as np
import glob
import time
import datetime
import os      
from os import path
import warnings
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
######################################### MD Score Main #########################################
def run_distance_calculation(verbose, outdir, sample, annotation, window, cpus, seq_type, pre_scan):
    if verbose == True: 
        print('--------------Beginning Distance Calculation---------------')
        print('Initializing ' + str(cpus) + ' threads to calculate MD Scores.')
        start_time = int(time.time())
        print('Start time: %s' % str(datetime.datetime.now()))
    md_dirs(verbose, outdir=outdir, seq_type=seq_type)
    
    if verbose==True:
        print('--------------Pulling in Annotation and Getting List of Motifs---------------')
    annotation_df = read_annotation(verbose=verbose, sample=sample, 
                                        outdir=outdir, pre_scan=pre_scan, annotation=annotation, seq_type=seq_type)
    tf_list = get_scanned_tfs(verbose=verbose, outdir=outdir, sample=sample, pre_scan=pre_scan, seq_type=seq_type) 
    
    for tf in tf_list:
        if verbose == True:
            print('Processing ' + tf + '.')
        motif_df = read_motif(verbose=verbose, outdir=outdir, seq_type=seq_type, pre_scan=pre_scan, tf=tf)
        
        if verbose == True:
            print("---------Pulling Chromosomes From Annotation File and Motif File----------") 
        chrs = get_chrs(verbose=verbose, motif_df=motif_df, annotation_df=annotation_df, tf=tf)
        if verbose == True:
            print("---------Calculating Motif Distance Scores----------")
        pool = multiprocessing.Pool(cpus)
        motif_distance_dfs = pool.map(partial(get_distances, 
                               inputs=[verbose, motif_df, annotation_df, window, tf]), chrs)
        pool.close()
        pool.join()
        motif_distance_df = pd.concat(motif_distance_dfs, axis=0)
        motif_score_df = calculate_rhf(verbose=verbose, outdir=outdir, sample=sample,
                                                        motif_distance_df=motif_distance_df, 
                                                        window=window,tf=tf, seq_type=seq_type)        

    if verbose == True:
        print("---------Simulated MD-Score Calculation Complete----------")
        stop_time = int(time.time())
        print('Stop time: %s' % str(datetime.datetime.now()))
        print('Total Run time :', (stop_time-start_time)/60, ' minutes')
    
######################################### MD Score Functions #########################################
def md_dirs(verbose, outdir, seq_type):
    if (path.exists(outdir + '/distances') == False):
        os.system('mkdir -p ' + outdir + '/distances') 
    try:
        os.system('mkdir -p ' + outdir + '/distances/' + seq_type + '_motif_scores')
    except OSError:
        print ('Creation of the directory %s failed' % outdir + '/distances/' + seq_type + '_motif_scores')
        sys.exit(1)
    else:
        if verbose == True:
            print(outdir + '/distances/' + seq_type + '_motif_scores exists.')
def get_scanned_tfs(verbose, outdir, sample, pre_scan, seq_type):
    tf_list = []
    if seq_type == 'experimental' and pre_scan is not None:
        tf_motif_path = pre_scan
    else:
        tf_motif_path = outdir + '/motifs/' + seq_type + '/*'
    motif_filenames = glob.glob(tf_motif_path)
    motif_count = len(motif_filenames)
    if verbose == True:
        print("Processing %s motif files in %s" % (motif_count, tf_motif_path))
    for filename in motif_filenames:
        filename_no_path = filename.split('/')[-1]
        filename_no_path=filename_no_path.replace('.sorted.bed','')
        tf_list.append(filename_no_path)
    if verbose == True:
        print('There are ' + str(len(tf_list)) + " motifs with hits in this dataset.")
    return tf_list

def read_annotation(verbose, sample, outdir, pre_scan, annotation, seq_type):
    if verbose == True:
        print('Reading annotation file and centering regions...')
    if seq_type == 'experimental' and pre_scan is not None:
        annotation_file = annotation
    else:
        annotation_file = outdir + '/annotations/' + sample + '_' + seq_type + '_window.bed'
        annotation_df = pd.read_csv(annotation_file, sep='\t', header=None)
        annotation_df.columns = ['chr', 'start', 'stop', 'region_name']
        annotation_df['center'] = round((annotation_df['stop'] + annotation_df['start'])/2)
        return annotation_df

def read_motif(verbose, outdir, seq_type, pre_scan, tf):
    if verbose == True:
        print('Reading motif file ' + tf + ' and centering regions...')
    if seq_type == 'experimental' and pre_scan is not None:
        motif_file = pre_scan + '/' + tf + '.sorted.bed'
        motif_df = pd.read_csv(motif_file, sep='\t', header=None)
        motif_df=motif_df[[0,1,2]] 
        motif_df.columns = ['chr','start', 'stop']
        motif_df['count'] = (np.arange(len(motif_df)))
        motif_df['count'] = (motif_df['count']+1).apply(str)
        motif_df['motif_region_name'] = motif_df['chr'] + ';motif_' + motif_df['count']
        motif_df.drop(['count'], axis=1, inplace=True)
    else:
        motif_file = outdir + '/motifs/' + seq_type + '/' + tf + '.sorted.bed'
        motif_df = pd.read_csv(motif_file, sep='\t', header=None)
        motif_df.columns = ['chr','start', 'stop', 'motif_id','score', 'strand','identifier', 'motif_region_name']
    motif_df['center'] = round((motif_df['stop'] + motif_df['start'])/2)
    motif_df = motif_df.sort_values(by=['chr', 'start'])
    return motif_df

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
            print(str(motif_not_annotation) + ' chromosomes are unique to the motif file')
        if len(annotation_not_motif) != 0:
            print(str(annotation_not_motif) + ' chromosomes are unique to the annotation file')
        else:
            print('The chromosomes match between the motif and annotation files.')
    return chrs

####TO BE EDITED
#some sort of intersect

###module to actually get distances

def find_bed_regions(row, single_chr_annotation):
    center = row['center']
    hits_df = single_chr_annotation[single_chr_annotation['start'] < center]
    hits_df = hits_df[hits_df['stop'] > center]
    distance = list(hits_df['center'] - row['center'])
    region_name = list(hits_df['region_name'])
    return [distance[0], region_name[0]]

def get_distances(chrs, inputs):
    verbose, motif_df, annotation_df, window, tf = inputs
    if verbose == True:
        print('Calculating ' + tf + ' distances for chromosome ' + chrs + '.') 
    
    #selecting a single chromosome for the annotation
    single_chr_annotation = annotation_df[annotation_df['chr'] == chrs]    
    
    if chrs in motif_df.values:
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
    
        distances =[]
        single_chr_motif['dis_rn'] = single_chr_motif.apply(lambda row: find_bed_regions(row, single_chr_annotation), axis=1)
        single_chr_motif_distances = single_chr_motif[['dis_rn']]
        distances.append(single_chr_motif_distances)

        distance_df = pd.concat(distances, ignore_index=False)
        motif_distance_df = motif_df.merge(distance_df, left_index=True, right_index=True)
        motif_distance_df[['distance','region_name']] = pd.DataFrame(motif_distance_df.dis_rn.tolist(), index= motif_distance_df.index)
        motif_distance_df = motif_distance_df.drop(['dis_rn'], axis=1)
        motif_distance_df=motif_distance_df[['motif_id', 'identifier','distance','region_name']]   
    else:
        if verbose == True:
            print('There are no motif hits on chromosome ' + chrs + '.')
    return motif_distance_df

def calculate_rhf(verbose, outdir, sample, motif_distance_df, window, tf, seq_type):
        rhf = motif_distance_df.groupby('region_name').count().reset_index()
        rhf = rhf[['region_name', 'motif_id']]
        rhf.columns = ['region_name', 'region_hit_frequency']
        motif_distance_df=motif_distance_df.merge(rhf, on='region_name')

        motif_distance_df.to_csv(outdir + '/distances/' + seq_type + '_motif_scores/' + tf + '.txt', 
                                 sep='\t', index=False)
        if verbose ==  True:
            print('Successfully calculated distance scores for ' + tf + '.')
