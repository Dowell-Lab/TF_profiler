######################################### Imports #########################################
import os
from os import path
import time
import datetime
import pandas as pd
import glob
import multiprocessing
from functools import partial

######################################### MD Score Main #########################################
def run_scoring_module(verbose, outdir, sample, window, cpus, seq_type):
    if verbose==True:
        print('--------------Pulling in Annotation and Getting List of Motifs---------------')
    scoring_dirs(verbose=verbose, outdir=outdir, seq_type=seq_type)
    tf_list = get_distance_tfs(verbose=verbose, outdir=outdir, sample=sample, seq_type=seq_type)
    #if traditional_md ==True: ##build in flag for different scoring options
    if verbose == True: 
        print('--------------Beginning Traditional MD Score Calculation---------------')
        print('Initializing ' + str(cpus) + ' threads to calculate MD score.')
        start_time = int(time.time())
        print('Start time: %s' % str(datetime.datetime.now()))
        pool = multiprocessing.Pool(cpus)
        tf_lhits = pool.map(partial(calculate_traditional_md_score, 
                               inputs=[verbose, outdir, sample, window, seq_type]), tf_list)
        pool.close()
        pool.join()
        traditional_md_score_df = pd.DataFrame.from_records(tf_lhits, columns=['tf', 'small_hits', 'large_hits', 'md_score'])
        traditional_md_score_df.to_csv(outdir+'/scores/' + seq_type + '_traditional_md_score.txt', sep='\t', index=False)
    if verbose == True:
        print("---------Traditional MD Score Calculation Complete----------")
        stop_time = int(time.time())
        print('Stop time: %s' % str(datetime.datetime.now()))
        print('Total Run time :', (stop_time-start_time)/60, ' minutes')

######################################### MD Score Functions #########################################
def scoring_dirs(verbose, outdir, seq_type):
    try:
        os.system('mkdir -p ' + outdir + '/scores')
    except OSError:
        print ('Creation of the directory %s failed' % outdir + '/scores')
        sys.exit(1)
    else:
        if verbose == True:
            print(outdir + '/scores exists.')
    
def get_distance_tfs(verbose, outdir, sample, seq_type):
    tf_list = []
    tf_distance_path = outdir + '/distances/'+seq_type+'/*'
    motif_filenames = glob.glob(tf_distance_path)
    motif_count = len(motif_filenames)
    if verbose == True:
        print("Processing %s motif files in %s" % (motif_count, tf_distance_path))
    for filename in motif_filenames:
        filename_no_path = filename.split('/')[-1]
        filename_no_path=filename_no_path.replace('_distances.txt','')
        tf_list.append(filename_no_path)
    if verbose == True:
        print('There are ' + str(len(tf_list)) + " motifs with hits in this dataset.")
    return tf_list

def calculate_traditional_md_score(tf_list, inputs):
    verbose, outdir, sample, window, seq_type = inputs
    distance_df = pd.read_csv(outdir+'/distances/'+seq_type+'/'+tf_list+'_distances.txt', sep='\t')
    hlarge=len(distance_df)+1
    hsmall= len(distance_df[(distance_df['distance'] <= window*0.1) & (distance_df['distance'] >= -window*0.1)])+1
    md_score = hsmall/hlarge
    return tf_list, hsmall, hlarge, md_score

