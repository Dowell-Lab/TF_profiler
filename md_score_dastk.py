##Adapted from DAStk by Ignacio Tripodi
######################################### Imports #########################################
import sys
import os
from os import path
import pandas as pd
import csv
import glob
import multiprocessing
import numpy as np
import pandas as pd
from functools import partial
from operator import itemgetter
import time
import datetime

######################################### MD Score Main #########################################
def run_md_score_dastk(window, cpus, outdir, sample, experimental_fimo, verbose):
    if verbose == True:
        print("---------Calculating the DAStk MD Scores for the Simulated Genome----------")
        print('Initializing ' + str(cpus) + ' cpus to calculate MD Scores.')
        start_time = int(time.time())
        print('Start time: %s' % str(datetime.datetime.now()))
    mdcalc(window, cpus, outdir, sample, verbose, seq_type='simulated')
    norm_name(verbose, outdir, sample, seq_type='simulated')
    if verbose == True:
        stop_time = int(time.time())
        print('Stop time: %s' % str(datetime.datetime.now()))
        print("Total Run time :", (stop_time-start_time)/60, " minutes")
    if experimental_fimo == True: 
        if verbose == True:
            print("---------Calculating the DAStk MD Scores for the Experimental Genome----------")
            print('Initializing ' + str(cpus) + ' cpus to calculate MD Scores.')
            start_time = int(time.time())
            print('Start time: %s' % str(datetime.datetime.now()))       
        mdcalc(window, cpus, outdir, sample, verbose, seq_type='experimental')
        norm_name(verbose, outdir, sample, seq_type='experimental')
        if verbose == True:
            stop_time = int(time.time())
            print('Stop time: %s' % str(datetime.datetime.now()))
            print("Total Run time :", (stop_time-start_time)/60, " minutes")     
        else:
            if verbose == True:
                print('No experimental motifs provided. Further downstream analyses halted')
            sys.exit(0)
######################################### MD Score Functions #########################################
def mdcalc(window, cpus, outdir, sample, verbose, seq_type):
    bed_name = outdir + '/annotations/' + sample + '_' + seq_type + '_centered.bed' 
    input_file_statinfo = os.stat(bed_name)    
    if input_file_statinfo.st_size !=0:
        bed = open(bed_name)
        bed_csv_reader = csv.reader(bed, delimiter='\t')
        bed_line = next(bed_csv_reader)
        bed_peak_count = 0        
    else:
        raise ValueError("The annotation file (%s_%s_centered.bed) is empty!" % (sample, seq_type))

    chromosomes = outdir + '/generated_sequences/' + sample + '_' + seq_type + '.chrom.sizes'
    if chromosomes:
        chr_df = pd.read_csv(chromosomes, header=None, comment='#', sep="\t", usecols=[0,1], names=['chrom', 'size'], na_filter=False, dtype={'chrom':'str', 'size':'int'})
        chr_list = list(chr_df['chrom'])
        CHROMOSOMES = [word for word in chr_list if len(word) <= 6]
        if verbose == True:
            print('Counting over chromosomes: ' + str(CHROMOSOMES))
    else:
        raise ValueError("Chromosome size file missing from generated_sequences(sample_' + seq_type + '.chrom.sizes)!") 

    bed.seek(0)
    bed_csv_reader = csv.reader(bed, delimiter='\t')
    bed_line = next(bed_csv_reader)
    while(bed_line[0][0] == '#'):
        bed_line = next(bed_csv_reader)    

    motif_stats = []
    sorted_motif_stats = []
    tf_motif_path = outdir + '/motifs/' + seq_type + '/*'
    motif_filenames = glob.glob(tf_motif_path)
    motif_count = len(motif_filenames)
    if verbose == True:
        print("Processing %s motif files in %s" % (motif_count, tf_motif_path))
    for filename in motif_filenames:
        filename_no_path = filename.split('/')[-1]
        if os.path.getsize(filename) > 0 and os.path.basename(filename).endswith(tuple(['.bed'])):
            [md_score, small_window, large_window, motif_site_count, heat] = get_md_score(filename, int(cpus), bed_name, CHROMOSOMES, window, verbose)
            if verbose == True:
                print('The %s MD-score for %s is %.6f' % (seq_type, filename_no_path, md_score))
            motif_stats.append({ 'motif_file': filename_no_path,'md_score': md_score, 'small_window': small_window, 'large_window': large_window, 'motif_site_count': motif_site_count, 'heat': heat })
    # sort the stats dictionary by MD-score, descending order
    sorted_motif_stats = sorted(motif_stats, key=itemgetter('md_score'), reverse=True)

    md_score_fp = open("%s/results/%s_%s_dastk_md_scores.txt" % (outdir, sample, seq_type), 'w')
    for stat in sorted_motif_stats:
        md_score_fp.write("%s,%s,%s,%s,%s,%s\n" % (stat['motif_file'], stat['md_score'], stat['small_window'], stat['large_window'], stat['motif_site_count'], stat['heat']))
    md_score_fp.close()
        
def is_in_window(motif_interval, bed_median, window_size):
    start = int(motif_interval.start)
    end = int(motif_interval.end)
    if (end >= (bed_median - window_size) and end <= (bed_median + window_size)) \
        or (start >= (bed_median - window_size) and start <= (bed_median + window_size)):
        return True
    else:
        return False
    
def find_motifs_in_chrom(current_chrom, files):
    tf_motif_filename, bed_name, window, verbose = files
    H = window          # in bps, the MD-score parameter (large window)
    h = window * 0.1           # in bps, the MD-score parameter (small window)

    wdf = pd.read_csv(bed_name, header=None, comment='#', sep="\t", usecols=[0, 1, 2], \
                          names=['chrom', 'start', 'end'], na_filter=False, dtype={'chrom':'str', 'start':'int', 'end':'int'})
    w_iter = wdf[(wdf.chrom == current_chrom)].itertuples()
    motif_df = pd.read_csv(tf_motif_filename, header=None, comment='#', sep="\t", usecols=[0, 1, 2], \
                           names=['chrom', 'start', 'end'], na_filter=False, dtype={'chrom':'str', 'start':'int', 'end':'int'})
    if len(motif_df) == 0:
        return None
    motif_iter = motif_df[(motif_df.chrom == current_chrom)].itertuples()
    last_motif = None
    g_h = 0
    g_H = 0
    total_motif_sites = 0
    tf_distances = []
    try:
        motif_region = next(motif_iter)   # get the first motif in the list
    except StopIteration:
        if verbose == True:
            print('No motifs for chromosome ' + current_chrom + ' on file ' + tf_motif_filename)
        return None

    peak_count_overlapping_motif = 0
    for bed_peak in w_iter:
        motifs_within_region = True

        bed_median = bed_peak.start + (bed_peak.end - bed_peak.start)/2
        # check the last motif, too, in case any of them falls within the region of
        # interest of two sequential ATAC-Seq peaks
        if last_motif:
            if is_in_window(last_motif, bed_median, h):
                g_h = g_h + 1
            if is_in_window(last_motif, bed_median, H):
                g_H = g_H + 1
                tf_median = last_motif.start + (last_motif.end - last_motif.start)/2
                tf_distances.append(bed_median - tf_median)
                try:
                    motif_region = next(motif_iter)   # get the next motif in the list
                    total_motif_sites += 1
                except StopIteration:
                    pass
                last_motif = motif_region
                if motif_region.start > (bed_median + H):
                    motifs_within_region = False

        # Move to the next putative motif sites until we get one past our evaluation window
        while motifs_within_region:
            total_motif_sites += 1
            # account for those within the smaller window (h)
            if is_in_window(motif_region, bed_median, h):
                g_h = g_h + 1
            # account for those within the larger window (H)
            if is_in_window(motif_region, bed_median, H):
                g_H = g_H + 1
                tf_median = motif_region.start + (motif_region.end - motif_region.start)/2
                tf_distances.append(bed_median - tf_median)

            # if we still haven't shifted past this peak...
            if motif_region.start <= (bed_median + H):
                try:
                    motif_region = next(motif_iter)   # get the next motif in the list
                    total_motif_sites += 1
                except StopIteration:
                    # No more TF motifs for this chromosome
                    break
                last_motif = motif_region
            else:
                motifs_within_region = False

    # Count any remaining TF motif sites after the last peak
    while(len(motif_region) > 0):
        try:
            motif_region = next(motif_iter)   # this gets the next motif in the list
        except StopIteration:
            break

    return [tf_distances, g_h, g_H, total_motif_sites]

def get_md_score(tf_motif_filename, cpus, bed_name, CHROMOSOMES, window, verbose):

    HISTOGRAM_BINS = 150
    pool = multiprocessing.Pool(cpus)
    results = pool.map(partial(find_motifs_in_chrom, files=[tf_motif_filename, bed_name, window, verbose]), CHROMOSOMES)
    pool.close()
    pool.join()

    results_matching_motif = [x for x in results if x is not None]
    if len(results_matching_motif) > 0:
        sums = np.sum(results_matching_motif, axis=0)
        overall_g_h = sums[1]
        overall_g_H = sums[2]
        overall_motif_sites = sums[3]

        # Calculate the heatmap for this motif's barcode
        tf_distances = sums[0]
        heatmap, xedges = np.histogram(tf_distances, bins=HISTOGRAM_BINS)
        str_heatmap = np.char.mod('%d', heatmap)
        if overall_g_H >= 0:
            return [float(overall_g_h + 1)/(overall_g_H + 1), (overall_g_h + 1), (overall_g_H + 1), (overall_motif_sites + 1), ';'.join(str_heatmap)]
    else:
        return None    

def norm_name(verbose, outdir, sample, seq_type):
    df = pd.read_csv(outdir + '/results/' + sample + '_' + seq_type + '_dastk_md_scores.txt', sep=',', 
                      names= ['motif_id', 'md_score', 'small_window', 'large_window',
                               'motif_site_count', 'heat'])
    df=df.replace(to_replace='.sorted.bed', value='', regex=True)
    df.to_csv(outdir + '/results/' + sample + '_' + seq_type + '_dastk_md_scores.txt', sep='\t', index=False)
    if verbose == True:
        print('MD Score file reformated.')