import os
from os import path
import glob
import sys
import pandas as pd

# Author fimo_motif_names: Jonathan Rubin    
def fimo_motif_names(motifs):
    '''Extracts motif names from a MEME formatted motif database
    Parameters'''
    motif_list = list()
    with open(motifs) as F:
        for line in F:
            if 'MOTIF' in line:
                line = line.strip('\n').split()
                motif_name = line[-1]
                motif_list.append(motif_name)
    first = [motif_list[0:1]]
    last =  [motif_list[-1:]]
    print('There are ' + str(len(motif_list)) + " motifs in this meme file. " + str(first) + '...' + str(last))
    return motif_list

def fimotobed(outdir):
    def identifier(row):
        return str(row['sequence_name']) + ':' + str(row['start']) + '-' + str(row['stop']) + '_' + str(row['strand'] + ';' + str(row['motif_id']))

    def percent_gc(row):
        r = pd.Series(row['matched_sequence'])
        z = r.str.count('G') + r.str.count('C')
        x = z/len(row['matched_sequence'])
        return x

    motif_dirs = glob.glob(outdir + '/temp/sim_fimo_out/*')
    for motif_dir in motif_dirs:
        if path.isdir(motif_dir):
            motif_name = motif_dir.split('/')[-1]
            df = pd.read_csv('%s/fimo.tsv' % (motif_dir), sep ='\t')
            print("Post-processing FIMO output for %s" % motif_name)
            df.drop(df.tail(3).index,inplace=True)
            if df.empty == True:
                print('Skipping ' + motif_name + ' -no motif hits.')
            else:
                df = df.sort_values('sequence_name').reset_index()
                df['identifier'] = df.apply(lambda row: identifier(row), axis=1)
                df['percent_gc'] = df.apply(lambda row: percent_gc(row), axis=1)
                df['mean_percent_gc'] = df['percent_gc'].mean()
                df = df[['sequence_name', 'start', 'stop', 'identifier', 'motif_id','score', 'strand', 
                         'matched_sequence', 'percent_gc', 'mean_percent_gc']]
                df['start'] = df['start'].astype(int) 
                df['stop'] = df['stop'].astype(int)                 
                df.to_csv(outdir + '/motifs/simulated/' + motif_name + '.sorted.bed', sep='\t', header=None, index=False)

def gc_content(outdir):
    directory = os.fsencode(outdir + '/motifs/simulated')
    output = pd.DataFrame({'motif_id': [], 'identifier' : [], 'mean_percent_gc' : []})
    for file in os.listdir(directory):
        if file.endswith(b'.bed'):
            file = file.decode('utf-8')
            df = pd.read_csv(outdir + '/motifs/simulated/' + file, sep="\t", header = None, 
                             names = ['sequence_name', 'start', 'stop', 'identifier', 'motif_id', 'score', 'strand', 'matched_sequence', 'percent_gc', 'mean_percent_gc'])
            df = df[['motif_id','identifier', 'mean_percent_gc']].loc[0]
            output=output.append(df)
            output.to_csv(outdir + '/annotations/simulated_mean_gc_percent.txt', sep="\t", header=None, index=False)
        else:
            continue
            
########### Experimental Genome
def fimotobed_expt(outdir):
    def identifier(row):
        return str(row['sequence_name']) + ':' + str(row['start']) + '-' + str(row['stop']) + '_' + str(row['strand'] + ';' + str(row['motif_id']))

    def percent_gc(row):
        r = pd.Series(row['matched_sequence'])
        z = r.str.count('G') + r.str.count('C')
        x = z/len(row['matched_sequence'])
        return x

    motif_dirs = glob.glob(outdir + '/temp/expt_fimo_out/*')
    for motif_dir in motif_dirs:
        if path.isdir(motif_dir):
            motif_name = motif_dir.split('/')[-1]
            print(motif_name)
            df = pd.read_csv('%s/fimo.tsv' % (motif_dir), sep ='\t')
            print("Post-processing FIMO output for %s" % motif_name)
            df.drop(df.tail(3).index,inplace=True)
            if df.empty == True:
                print('Skipping ' + motif_name + ' -no motif hits.')
            else:
                df = df.sort_values('sequence_name').reset_index()
                df['identifier'] = df.apply(lambda row: identifier(row), axis=1)
                df['percent_gc'] = df.apply(lambda row: percent_gc(row), axis=1)
                df['mean_percent_gc'] = df['percent_gc'].mean()
                df = df[['sequence_name', 'start', 'stop', 'identifier', 'motif_id', 'score', 'strand', 
                         'matched_sequence', 'percent_gc', 'mean_percent_gc']]
                df['start'] = df['start'].astype(int) 
                df['stop'] = df['stop'].astype(int) 
                df.to_csv(outdir + '/motifs/experimental/' + motif_name + '.sorted.bed', sep='\t', header=None, index=False)

def gc_content_expt(outdir):
    directory = os.fsencode(outdir + '/motifs/experimental')
    output = pd.DataFrame({'motif_id': [], 'identifier' : [], 'mean_percent_gc' : []})
    for file in os.listdir(directory):
        if file.endswith(b'.bed'):
            file = file.decode('utf-8')
            df = pd.read_csv(outdir + '/motifs/experimental/' + file, sep="\t", header = None, 
                             names = ['sequence_name', 'start', 'stop', 'identifier', 'motif_id', 'score', 'strand', 'matched_sequence', 'percent_gc', 'mean_percent_gc'])
            df = df[['motif_id','identifier', 'mean_percent_gc']].loc[0]
            output=output.append(df)
            output.to_csv(outdir + '/annotations/experimental_mean_gc_percent.txt', sep="\t", header=None, index=False)
        else:
            continue
########### Whole Genome
def fimotobed_whole_genome(outdir):
    def identifier(row):
        return str(row['sequence_name']) + ':' + str(row['start']) + '-' + str(row['stop']) + '_' + str(row['strand'] + ';' + str(row['motif_id']))

    def percent_gc(row):
        r = pd.Series(row['matched_sequence'])
        z = r.str.count('G') + r.str.count('C')
        x = z/len(row['matched_sequence'])
        return x

    motif_dirs = glob.glob(outdir + '/temp/whole_genome_fimo_out/*')
    for motif_dir in motif_dirs:
        if path.isdir(motif_dir):
            motif_name = motif_dir.split('/')[-1]
            print(motif_name)
            df = pd.read_csv('%s/fimo.tsv' % (motif_dir), sep ='\t')
            print("Post-processing FIMO output for %s" % motif_name)
            df.drop(df.tail(3).index,inplace=True)
            if df.empty == True:
                print('Skipping ' + motif_name + ' -no motif hits.')
            else:
                df = df.sort_values('sequence_name').reset_index()
                df['identifier'] = df.apply(lambda row: identifier(row), axis=1)
                df['percent_gc'] = df.apply(lambda row: percent_gc(row), axis=1)
                df['mean_percent_gc'] = df['percent_gc'].mean()
                df = df[['sequence_name', 'start', 'stop', 'identifier', 'motif_id', 'score', 'strand', 
                         'matched_sequence', 'percent_gc', 'mean_percent_gc']]
                df['start'] = df['start'].astype(int) 
                df['stop'] = df['stop'].astype(int) 
                df.to_csv(outdir + '/motifs/whole_genome/' + motif_name + '.sorted.bed', sep='\t', header=None, index=False)

def gc_content_whole_genome(outdir):
    directory = os.fsencode(outdir + '/motifs/whole_genome')
    output = pd.DataFrame({'motif_id': [], 'identifier' : [], 'mean_percent_gc' : []})
    for file in os.listdir(directory):
        if file.endswith(b'.bed'):
            file = file.decode('utf-8')
            df = pd.read_csv(outdir + '/motifs/whole_genome/' + file, sep="\t", header = None, 
                             names = ['sequence_name', 'start', 'stop', 'identifier', 'motif_id', 'score', 'strand', 'matched_sequence', 'percent_gc', 'mean_percent_gc'])
            df = df[['motif_id','identifier', 'mean_percent_gc']].loc[0]
            output=output.append(df)
            output.to_csv(outdir + '/annotations/whole_genome_mean_gc_percent.txt', sep="\t", header=None, index=False)
        else:
            continue