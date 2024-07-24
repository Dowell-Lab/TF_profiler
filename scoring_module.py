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
def run_scoring_module(verbose, outdir, sample, window, cpus, seq_type, simulated_pre_scan, annotation):
    if verbose==True:
        print('--------------Pulling in Annotation and Getting List of Motifs---------------')
    scoring_dirs(verbose=verbose, outdir=outdir, seq_type=seq_type)
    tf_list = get_distance_tfs(verbose=verbose, outdir=outdir, sample=sample, seq_type=seq_type, simulated_pre_scan=simulated_pre_scan)
    sequence_num, promoter_num = get_sequence_numbers(outdir=outdir, sample=sample, annotation=annotation)

    if verbose == True: 
        print('--------------Beginning Traditional MD Score Calculation---------------')
        print('Initializing ' + str(cpus) + ' threads to calculate MD score.')
        start_time = int(time.time())
        print('Start time: %s' % str(datetime.datetime.now()))
        pool = multiprocessing.Pool(cpus)
        tf_lhits = pool.map(partial(calculate_traditional_md_score, 
                               inputs=[verbose, outdir, sample, window, 
                                       seq_type, sequence_num, promoter_num, simulated_pre_scan]), tf_list)
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
    
def get_distance_tfs(verbose, outdir, sample, seq_type, simulated_pre_scan):
    tf_list = []
    if simulated_pre_scan is not None and seq_type=='simulated':
        tf_distance_path=simulated_pre_scan+'/*_distances.txt'
    else:
        tf_distance_path = outdir + '/distances/'+seq_type+'/*_distances.txt'
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

def get_sequence_numbers(outdir, sample, annotation):
    with open(annotation, 'r') as an:
        sequence_num = len(an.readlines())
        
    directory = os.path.abspath(__file__)
    directory = os.path.dirname(directory)
                                                                          
    assetbed = os.path.join(directory, '/assets/hg38_refseq_merge_1000bp_TSSs.bed') #mm10_refseq_unique_TSSs_1000bp_merge.sorted.bed

    os.system('bedtools intersect -wa -u -a '+ annotation + 
         ' -b ' + assetbed + ' > ' +
         outdir + '/annotations/'+sample+'_promoters.bed')
        
    with open(outdir + '/annotations/'+sample+'_promoters.bed', 'r') as pro:
        promoter_num = len(pro.readlines())
    return sequence_num, promoter_num

def calculate_traditional_md_score(tf_list, inputs):
    verbose, outdir, sample, window, seq_type, sequence_num, promoter_num, simulated_pre_scan = inputs
    if simulated_pre_scan is not None and seq_type=='simulated':
        distance_df = pd.read_csv(simulated_pre_scan+'/'+tf_list+'_distances.txt', sep='\t')
    else:
        distance_df = pd.read_csv(outdir+'/distances/'+seq_type+'/'+tf_list+'_distances.txt', sep='\t')
    
    distance_df = distance_df.sort_values(by=['region_id', 'distance_rank', 'quality_rank'])
    distance_df=distance_df.drop_duplicates(subset=['region_id'], keep='first')
    
    hlarge=len(distance_df)+1
    hsmall= len(distance_df[(distance_df['distance'] <= window*0.1) & (distance_df['distance'] >= -window*0.1)])+1
    md_score = hsmall/hlarge
    return tf_list, hsmall, hlarge, md_score


