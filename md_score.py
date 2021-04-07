##Adapted from DAStk by Ignacio Tripodi
import os
import sys
from os import path
import pandas as pd
import csv
import glob
from md_functions import *

def mdcalc(window, threads, outdir, sample):
    bed_name = outdir + '/annotations/' + sample + '_simulated_centered.bed' 
    input_file_statinfo = os.stat(bed_name)    
    if input_file_statinfo.st_size !=0:
        bed = open(bed_name)
        bed_csv_reader = csv.reader(bed, delimiter='\t')
        bed_line = next(bed_csv_reader)
        bed_peak_count = 0        
    else:
        raise ValueError("The annotation file (sample_simulated_centered.bed) is empty!")

    chromosomes = outdir + '/generated_sequences/' + sample + '_simulated.chrom.sizes'
    if chromosomes:
        chr_df = pd.read_csv(chromosomes, header=None, comment='#', sep="\t", usecols=[0,1], names=['chrom', 'size'], na_filter=False, dtype={'chrom':'str', 'size':'int'})
        chr_list = list(chr_df['chrom'])
        CHROMOSOMES = [word for word in chr_list if len(word) <= 6]
        print('Counting over chromosomes: ' + CHROMOSOMES)
    else:
        raise ValueError("Chromosome size file missing from generated_sequences(sample_simulated.chrom.sizes)!") 

    bed.seek(0)
    bed_csv_reader = csv.reader(bed, delimiter='\t')
    bed_line = next(bed_csv_reader)
    while(bed_line[0][0] == '#'):
        bed_line = next(bed_csv_reader)    

    motif_stats = []
    sorted_motif_stats = []
    tf_motif_path = outdir + '/motifs/simulated/*'
    motif_filenames = glob.glob(tf_motif_path)
    motif_count = len(motif_filenames)
    print("Processing %s motif files in %s" % (motif_count, tf_motif_path))
    for filename in motif_filenames:
        filename_no_path = filename.split('/')[-1]
        if os.path.getsize(filename) > 0 and os.path.basename(filename).endswith(tuple(['.bed'])):
            [md_score, small_window, large_window, motif_site_count, heat] = get_md_score(filename, int(threads), bed_name, CHROMOSOMES, window)
            print('The MD-score for %s is %.6f' % (filename_no_path, md_score))
            motif_stats.append({ 'motif_file': filename_no_path,'md_score': md_score, 'small_window': small_window, 'large_window': large_window, 'motif_site_count': motif_site_count, 'heat': heat })
    # sort the stats dictionary by MD-score, descending order
    sorted_motif_stats = sorted(motif_stats, key=itemgetter('md_score'), reverse=True)

    md_score_fp = open("%s/temp/md_scores/%s_simulated_md_scores.txt" % (outdir, sample), 'w')
    for stat in sorted_motif_stats:
        md_score_fp.write("%s,%s,%s,%s,%s,%s\n" % (stat['motif_file'], stat['md_score'], stat['small_window'], stat['large_window'], stat['motif_site_count'], stat['heat']))
    md_score_fp.close()
    
    
    
## Experimental
def mdcalc_expt(window, threads, outdir, sample):
    bed_name = outdir + '/annotations/' + sample + '_experimental_genome_centered.bed' 
    input_file_statinfo = os.stat(bed_name)    
    if input_file_statinfo.st_size !=0:
        bed = open(bed_name)
        bed_csv_reader = csv.reader(bed, delimiter='\t')
        bed_line = next(bed_csv_reader)
        bed_peak_count = 0        
    else:
        raise ValueError("The annotation file (sample_experimental_genome_centered) is empty!")

    chromosomes = outdir + '/generated_sequences/' + sample + '_experimental.chrom.sizes'
    if chromosomes:
        chr_df = pd.read_csv(chromosomes, header=None, comment='#', sep="\t", usecols=[0,1], names=['chrom', 'size'], na_filter=False, dtype={'chrom':'str', 'size':'int'})
        chr_list = list(chr_df['chrom'])
        CHROMOSOMES = [word for word in chr_list if len(word) <= 6]
    else:
        raise ValueError("Chromosome size file missing from generated_sequences(sample_experimental.chrom.sizes)!") 

    bed.seek(0)
    bed_csv_reader = csv.reader(bed, delimiter='\t')
    bed_line = next(bed_csv_reader)
    while(bed_line[0][0] == '#'):
        bed_line = next(bed_csv_reader)    

    motif_stats = []
    sorted_motif_stats = []
    tf_motif_path = outdir + '/motifs/experimental/*'
    motif_filenames = glob.glob(tf_motif_path)
    motif_count = len(motif_filenames)
    print("Processing %s motif files in %s" % (motif_count, tf_motif_path))
    for filename in motif_filenames:
        filename_no_path = filename.split('/')[-1]
        if os.path.getsize(filename) > 0 and os.path.basename(filename).endswith(tuple(['.bed'])):
            [md_score, small_window, large_window, motif_site_count, heat] = get_md_score(filename, int(threads), bed_name, CHROMOSOMES, window)
            print('The MD-score for %s is %.6f' % (filename_no_path, md_score))
            motif_stats.append({ 'motif_file': filename_no_path,'md_score': md_score, 'small_window': small_window, 'large_window': large_window, 'motif_site_count': motif_site_count, 'heat': heat })
    # sort the stats dictionary by MD-score, descending order
    sorted_motif_stats = sorted(motif_stats, key=itemgetter('md_score'), reverse=True)

    md_score_fp = open("%s/temp/md_scores/%s_experimental_md_scores.txt" % (outdir, sample), 'w')
    for stat in sorted_motif_stats:
        md_score_fp.write("%s,%s,%s,%s,%s,%s\n" % (stat['motif_file'], stat['md_score'], stat['small_window'], stat['large_window'], stat['motif_site_count'], stat['heat']))
    md_score_fp.close()


###Whole genome
def mdcalc_whole_genome(annotation, window, threads, outdir, sample, genome):
    bed_name = annotation
    input_file_statinfo = os.stat(bed_name)    
    if input_file_statinfo.st_size !=0:
        bed = open(bed_name)
        bed_csv_reader = csv.reader(bed, delimiter='\t')
        bed_line = next(bed_csv_reader)
        bed_peak_count = 0        
    else:
        raise ValueError("The annotation file is empty!")

        
    g = genome.split('/')[-1].split('.')[-2]    
    chromosomes = outdir + "/generated_sequences/" + g + ".chrom.sizes"
    if chromosomes:
        chr_df = pd.read_csv(chromosomes, header=None, comment='#', sep="\t", usecols=[0,1], names=['chrom', 'size'], na_filter=False, dtype={'chrom':'str', 'size':'int'})
        chr_list = list(chr_df['chrom'])
        CHROMOSOMES = [word for word in chr_list if len(word) <= 6]
    else:
        raise ValueError("Chromosome size file missing from generated_sequences(genome.chrom.sizes)!") 

    bed.seek(0)
    bed_csv_reader = csv.reader(bed, delimiter='\t')
    bed_line = next(bed_csv_reader)
    while(bed_line[0][0] == '#'):
        bed_line = next(bed_csv_reader)    

    motif_stats = []
    sorted_motif_stats = []
    tf_motif_path = outdir + '/motifs/whole_genome/*'
    motif_filenames = glob.glob(tf_motif_path)
    motif_count = len(motif_filenames)
    print("Processing %s motif files in %s" % (motif_count, tf_motif_path))
    for filename in motif_filenames:
        filename_no_path = filename.split('/')[-1]
        if os.path.getsize(filename) > 0 and os.path.basename(filename).endswith(tuple(['.bed'])):
            [md_score, small_window, large_window, motif_site_count, heat] = get_md_score(filename, int(threads), bed_name, CHROMOSOMES, window)
            print('The MD-score for %s is %.6f' % (filename_no_path, md_score))
            motif_stats.append({ 'motif_file': filename_no_path,'md_score': md_score, 'small_window': small_window, 'large_window': large_window, 'motif_site_count': motif_site_count, 'heat': heat })
    # sort the stats dictionary by MD-score, descending order
    sorted_motif_stats = sorted(motif_stats, key=itemgetter('md_score'), reverse=True)

    md_score_fp = open("%s/temp/md_scores/%s_experimental_md_scores.txt" % (outdir, sample), 'w')
    for stat in sorted_motif_stats:
        md_score_fp.write("%s,%s,%s,%s,%s,%s\n" % (stat['motif_file'], stat['md_score'], stat['small_window'], stat['large_window'], stat['motif_site_count'], stat['heat']))
    md_score_fp.close()


#Pre_scan
def mdcalc_prescan(annotation, window, threads, outdir, sample, genome, pre_scan):
    bed_name = annotation 
    input_file_statinfo = os.stat(bed_name)    
    if input_file_statinfo.st_size !=0:
        bed = open(bed_name)
        bed_csv_reader = csv.reader(bed, delimiter='\t')
        bed_line = next(bed_csv_reader)
        bed_peak_count = 0        
    else:
        raise ValueError("The annotation file is empty!")

        
    g = genome.split('/')[-1].split('.')[-2]    
    chromosomes = outdir + "/generated_sequences/" + g + ".chrom.sizes"
    if chromosomes:
        chr_df = pd.read_csv(chromosomes, header=None, comment='#', sep="\t", usecols=[0,1], names=['chrom', 'size'], na_filter=False, dtype={'chrom':'str', 'size':'int'})
        chr_list = list(chr_df['chrom'])
        CHROMOSOMES = [word for word in chr_list if len(word) <= 6]
    else:
        raise ValueError("Chromosome size file missing from generated_sequences(genome.chrom.sizes)!") 

    bed.seek(0)
    bed_csv_reader = csv.reader(bed, delimiter='\t')
    bed_line = next(bed_csv_reader)
    while(bed_line[0][0] == '#'):
        bed_line = next(bed_csv_reader)    

    motif_stats = []
    sorted_motif_stats = []
    tf_motif_path = pre_scan + '/*'
    motif_filenames = glob.glob(tf_motif_path)
    motif_count = len(motif_filenames)
    print("Processing %s motif files in %s" % (motif_count, tf_motif_path))
    for filename in motif_filenames:
        filename_no_path = filename.split('/')[-1]
        if os.path.getsize(filename) > 0 and os.path.basename(filename).endswith(tuple(['.bed'])):
            [md_score, small_window, large_window, motif_site_count, heat] = get_md_score(filename, int(threads), bed_name, CHROMOSOMES, window)
            print('The MD-score for %s is %.6f' % (filename_no_path, md_score))
            motif_stats.append({ 'motif_file': filename_no_path,'md_score': md_score, 'small_window': small_window, 'large_window': large_window, 'motif_site_count': motif_site_count, 'heat': heat })
    # sort the stats dictionary by MD-score, descending order
    sorted_motif_stats = sorted(motif_stats, key=itemgetter('md_score'), reverse=True)

    md_score_fp = open("%s/temp/md_scores/%s_experimental_md_scores.txt" % (outdir, sample), 'w')
    for stat in sorted_motif_stats:
        md_score_fp.write("%s,%s,%s,%s,%s,%s\n" % (stat['motif_file'], stat['md_score'], stat['small_window'], stat['large_window'], stat['motif_site_count'], stat['heat']))
    md_score_fp.close()