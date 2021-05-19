######################################### Imports #########################################
import os
from os import path
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import math
import datetime
import time

######################################### Sequence Generator Main #########################################
def run_sequence_generator(verbose, outdir, sample,
                           sequence_num, chrom_num, window, 
                           genome, annotation, seed,
                           experimental_fimo, whole_genome_fimo):
    if verbose == True: 
        print('--------------Expanding Windows---------------')
    window_annotation(verbose, annotation=annotation, outdir=outdir, sample=sample, window=window)
    
    if verbose == True:
        print('-------------Extracting Sequences-------------')
    get_sequences(verbose, genome=genome, outdir=outdir, sample=sample)

    if verbose == True: 
        print('-----------Counting Base Content and Plot Generation-------------')
    window_seq = list_sequences(outdir=outdir, sample=sample)
    base_plot(verbose, seq=window_seq, plottitle=sample, outdir=outdir,
              sample=sample,  window=window)
    tsv = (outdir + '/generated_sequences/' + sample +'_base_content.tsv')

    if verbose == True: 
        print('-----------Start Simulated Sequence Generation-------------')
        start_time = int(time.time())
    npseed = set_seed(verbose, seed=seed)
    
    if verbose == True: 
        print('-----------Reading Per-Base Sequence Frequency----------------')
        print('Base Content File: ' + tsv)
    position_freq = base_composition_input(base_composition=tsv)

    if verbose == True: 
        print('-----------Generating Sequences-------------------------------')
    generated_sequences = sequence_generator(verbose, bases=['A', 'T', 'G', 'C'], position_freq=position_freq, 
                                             sequence_num=sequence_num)

    if verbose == True: 
        print('-----------Writing Sequences to Fasta File--------------------')
    write_fasta(verbose, generated_sequences=generated_sequences, sample=sample, outdir=outdir, 
                window=window, sequence_num=sequence_num, chrom_num=chrom_num)
    
    if verbose == True: 
        print('Simulated Fasta File: ' + outdir + '/generated_sequences/' + str(sample) + '_simulated.fa')
        print('----------------Simulated Sequence Generation Complete------------------------------')
        stop_time = int(time.time())
        print('Stop time: %s' % str(datetime.datetime.now()))
        print('Total Run time :', stop_time-start_time, ' seconds')

    if verbose == True: 
        print('-----------Indexing Fasta, Generating Chrm Sizes and Creating Bedfiles--------------------')
    index_and_chrm_sizes(verbose, outdir=outdir, sample=sample, genome=genome, seq_type='simulated')
    generate_bed(verbose, outdir=outdir, sample=sample, sequence_num=sequence_num, chrom_num=chrom_num, 
                 window=window, annotation=annotation, seq_type='simulated')
    if experimental_fimo == True:
        reformat_expt_fasta(verbose, outdir=outdir, sample=sample, chrom_num=chrom_num, window=window)
        index_and_chrm_sizes(verbose, outdir=outdir, sample=sample, genome=genome, seq_type='experimental')
        generate_bed(verbose, outdir=outdir, sample=sample, sequence_num=sequence_num, chrom_num=chrom_num, 
             window=window, annotation=annotation, seq_type='experimental')
    if whole_genome_fimo == True:
        index_and_chrm_sizes(verbose, outdir=outdir, sample=sample, genome=genome, seq_type='whole_genome') 
        
######################################### Sequence Generator Functions #########################################

def window_annotation(verbose, annotation, outdir, sample, window):
    '''This function takes in the annotation file from Tfit and redefines mu and extends the window
    Calls windower- a function that windows and saves bedfiles.
    '''
    if verbose == True:    
        print('Annotation File: ' + annotation)
    
    try:
        os.system('mkdir -p ' + outdir + '/annotations')
    except OSError:
        print ('Creation of the directory %s failed' % outdir + '/annotations')
        sys.exit(1)
    else:
        if verbose == True:
            print('Successfully created the directory %s' % outdir + '/annotations')
    os.system('rsync ' + annotation + ' ' + outdir + '/annotations/')        
    windower(bed=annotation, outdir=outdir, sample=sample, window=int(window), seq_type='experimental', ident=False)
    
def get_sequences(verbose, genome, outdir, sample):
    '''This function pulls the sequences out of the windowed annotation file and outputs them in 
    generated_sequences for further use'''
    if verbose == True:
        print('Genome Used: ' + genome + '\nOutdir: ' + outdir + '\nSample: ' + sample)
    
    try:
        os.system('mkdir -p ' + outdir + '/temp')
    except OSError:
        print ('Creation of the directory %s failed' % outdir + '/temp')
        sys.exit(1)
    else:
        if verbose == True:
            print('Successfully created the directory %s' % outdir + '/temp')
            
    os.system('bedtools getfasta -fi ' + genome + 
              ' -bed ' + outdir + '/annotations/' + sample + '_experimental_window.bed' + 
              ' -fo ' + outdir + '/temp/' + sample + '_experimental_window_sequences.fa')
    if (path.exists(outdir + '/temp/' + sample + '_experimental_window_sequences.fa') == False):
        print('Extracted experimental sequences failed. Make sure bedtools is installed.')
        sys.exit(1)
    else:
        if verbose == True:
            print('Extracted Experimental Sequences Output: ' + outdir + '/temp/' + sample + '_experimental_window_sequences.fa')

def list_sequences(outdir, sample):
    '''This function strings experimental sequences together in a list'''
    fasta_sequences = outdir + '/temp/' + sample + '_experimental_window_sequences.fa'
    sequences = []
    sequence_names = []

    #add sequences and sequence names to lists                                                                                                                     
    with open(fasta_sequences) as fa:
        for line in fa:
            line = line.strip('\n')
            if '>' not in line:
                sequences.append(line)
            else:
                sequence_names.append(line)

    #return list with sequence         
    return sequences

def base_plot(verbose, seq, plottitle, outdir, sample, window):
    '''plot line plot for each bases content across sequences for all sequences on the same grid.
    This script calls count_bases.
    Parameters:
    plottitle : str - title for the plot represented as a string'''
    
    try:
        os.system('mkdir -p ' + outdir + '/generated_sequences')
    except OSError:
        if verbose == True:
            print ('Creation of the directory %s failed' % outdir + '/generated_sequences')
            sys.exit(1)
    else:
        if verbose == True:
            print('Successfully created the directory %s' % outdir + '/generated_sequences')
    
    counts = count_bases(verbose, sequences=seq, outdir=outdir, sample=sample, window=window)

    ##get positions                                                   
    starting=int(window)
    stopping=int(window)+1
    positions = np.arange(-starting,stopping,1)
                
    ### Line plot for each base on one grid ###
    plt.figure(figsize=(12,10))
    gs = plt.GridSpec(1, 1)
    ax0 = plt.subplot(gs[0])
    color_list = ['blue', 'purple', 'orange', 'red']
    base_list=['A','T','C','G']
    for i in range(0,4):
        ax0.plot(positions,counts[i],color=color_list[i], alpha=0.75, label=base_list[i])
    ax0.set_xlabel('Distance (bp)',fontsize=30,fontweight='bold')
    ax0.set_ylabel('Base Content',fontsize=30,fontweight='bold')
    plt.xticks(fontsize = 25)
    plt.yticks(fontsize = 25)
    plt.legend(bbox_to_anchor=(1.05, 1),fontsize = 25, loc=2, borderaxespad=0.)
    plt.suptitle(plottitle, fontsize=40, fontweight='bold')
    plt.savefig(outdir + '/generated_sequences/' + sample + '_BaseDistribution.png',bbox_inches='tight')
    plt.cla()
                                                                                                         
    ### Smooth the frequencies for better visualization ### 
    plt.figure(figsize=(12,10))
    gs1 = plt.GridSpec(1, 1)
    ax1 = plt.subplot(gs1[0])

    for i in range(0,4):
        smoothed_count = savgol_filter(tuple(np.array(counts)), 61, 3) # window size 61, polynomial order 3     
        ax1.plot(positions,smoothed_count[i],color=color_list[i], alpha=0.75, label=base_list[i])
    ax1.set_xlabel('Distance (bp)',fontsize=30,fontweight='bold')
    ax1.set_ylabel('Base Content',fontsize=30,fontweight='bold')
    plt.xticks(fontsize = 25)
    plt.yticks(fontsize = 25)
    plt.legend(bbox_to_anchor=(1.05, 1),fontsize = 25, loc=2, borderaxespad=0.)
    plt.suptitle(plottitle, fontsize=40, fontweight='bold')
    plt.savefig(outdir + '/generated_sequences/' + sample + '_SmoothedBaseDistribution.png',bbox_inches='tight')
    plt.cla()
    
    if verbose == True:
        print('Generated Base Composition Plots: ' + outdir + '/generated_sequences/' + sample + '_SmoothedBaseDistribution.png and ' + sample + '_BaseDistribution.png')


def count_bases(verbose, sequences, outdir, sample, window):
    '''Called by base_plot. This script calculates per position base composition across multiple sequences of even length.
    anew,tnew,cnew,gnew,nnew : list of lists normalized base counts for every position across all sequences'''
    if verbose==True:
        print('Sequence List: Contains ' + str(len(sequences)) + ' elements. ')

    #since the window is flanking a single base we accout for that by 2(w)+1
    sequence_length = int(2*int(window)+1)
    if verbose==True:
        print('Full Window Length: +/- ' + str(window) + ' = ' + str(sequence_length))

    ##initialize lists with length of sequences                                                                                                                                  
    a = [0]*sequence_length
    t = [0]*sequence_length
    c = [0]*sequence_length
    g = [0]*sequence_length
    n = [0]*sequence_length

    ##for each positions count the occurance of each base                                                                                                                        
    ##across all sequences in the the input list                                                                                                                                 
    for i in range(sequence_length):
        ##initialize counters                                                                                                                                                    
        count_a = 0
        count_t = 0
        count_c = 0
        count_g = 0
        count_n = 0
        for j in sequences:
            if j[i] == 'a' or j[i] == 'A':
                count_a = count_a + 1
                a[i] = count_a
            elif j[i] == 't' or j[i] == 'T':
                count_t = count_t + 1
                t[i] = count_t
            elif j[i] == 'g' or j[i] == 'G':
                count_g = count_g + 1
                g[i] = count_g
            elif j[i] == 'c' or j[i] == 'C':
                count_c = count_c + 1
                c[i] = count_c
            elif j[i] == 'n' or j[i] == 'N':
                count_n = count_n + 1
                n[i] = count_n

    ##evenly distribute Ns across all bases                                                                                                                                      
    nnew = [x / 4 for x in n]

    anew = [ai + bi for ai,bi in zip(a,nnew)]
    tnew = [ai + bi for ai,bi in zip(t,nnew)]
    gnew = [ai + bi for ai,bi in zip(g,nnew)]
    cnew = [ai + bi for ai,bi in zip(c,nnew)]

    ##get the base frequencies of all bases                                                                                                                                      
    anew = [x / len(sequences) for x in anew]
    tnew = [x / len(sequences) for x in tnew]
    cnew = [x / len(sequences) for x in cnew]
    gnew = [x / len(sequences) for x in gnew]

    base_df = pd.DataFrame({'A': anew,
                            'T': tnew,
                            'G': cnew,
                            'C': gnew})

    base_df.to_csv(outdir + '/generated_sequences/' + sample +'_base_content.tsv', sep='\t')
    if verbose == True:
        print('Generated Positional Base Frequencies: ' + outdir + '/generated_sequences/' + sample +'_base_content.tsv')
    return anew, tnew, cnew, gnew, nnew

def set_seed(verbose, seed):
    '''Sets numpy seed for the sequence generation'''
    if verbose == True: 
        print('Start time: %s' % str(datetime.datetime.now()))
        if seed == 0:
            npseed = np.random.seed()
            if verbose == True: 
                print('-----------Setting Seed for the Sequence Generator.-----------')
                print('Using random seed.')
                print('WARNING: Sequences generated using a random seed can not be reproduced.')
        elif seed is True:
            npseed = np.random.seed(global_seed(verbose=verbose))
        else:
            npseed = np.random.seed(seed)
            if verbose == True:
                print('-----------Setting Seed for the Sequence Generator.-----------')
                print('User defined seed:', seed)
            return npseed
                
def global_seed(verbose):
    '''function used by set seed to set seed to the time'''
    if verbose==True:
        print('-----------Setting Seed for the Sequence Generator------------')
        print('Time in fractional seconds of:', time.time())
        print('Time as an interger:',int(time.time()))
        print('Setting seed on the clock')
    print('np.random.random(',int(time.time()),') => ',np.random.random())
    return(int(time.time()))

def base_composition_input(base_composition):
    '''takes in base_compositions that are in .tsv format and 
    frequencies per position to use in simulating base_compositions using a 
    1st order Markov Model.
    '''

    position_feq = []

    with open(base_composition) as bc:

        ##remove the header line
                ##remove the header line
        lines = bc.readlines()[1:]
        for i in range(len(lines)):
            pos_prob = []
            try:
                pos_prob.append(lines[i].strip('\n').split(',')[1])
                pos_prob.append(lines[i].strip('\n').split(',')[2])
                pos_prob.append(lines[i].strip('\n').split(',')[3])
                pos_prob.append(lines[i].strip('\n').split(',')[4])
                position_feq.append(pos_prob)
            except IndexError:
                
                pos_prob.append(lines[i].strip('\n').split('\t')[1])
                pos_prob.append(lines[i].strip('\n').split('\t')[2])
                pos_prob.append(lines[i].strip('\n').split('\t')[3])
                pos_prob.append(lines[i].strip('\n').split('\t')[4])
                position_feq.append(pos_prob)
            else:
                print("Check the input file. \n It should be tab or comma separated")


    return(position_feq)

def sequence_generator(verbose, bases, position_freq, sequence_num):
    '''takes in frequencies per position and simulates sequences using a 
    1st order Markov Model.
    '''
    
    first = [position_freq[0:1]]
    last =  [position_freq[-1:]]
    if verbose == True:
        print('Position Frequency: Contains ' + str(len(position_freq)) + ' elements. ' + str(first) + '...' + str(last))
    
    sequences = np.empty([sequence_num, len(position_freq)], dtype=str)
    
    for i in range(len(position_freq)):
        column = np.random.choice(bases, sequence_num, p=position_freq[i])
        sequences[:,i] = column
        
    joined_sequences = [''.join(row) for row in sequences]
    
    return joined_sequences

def write_fasta(verbose, generated_sequences, sample, outdir, window, sequence_num, chrom_num):
    ''' writes sequences generated into fasta file format and outputs them in generated sequences'''
    seq_per_chrom, seq_last_chrom = define_seq_structure(verbose, sequence_num=sequence_num, chrom_num=chrom_num, window=window)
    
    file=[]
    for i in range(len(generated_sequences)):
        file.append(str(generated_sequences[i]) + str('N')*20)      

    dd = {}
    for i in range(0,chrom_num-1):
        dd['>sim_' + str(i + 1)] = str(file[(seq_per_chrom)*i:(seq_per_chrom)*(i+1)]).replace("', '","").replace("['","").replace("']","")
    dd['>sim_' + str(chrom_num)] = str(file[(seq_per_chrom)*(chrom_num-1):]).replace("', '","").replace("['","").replace("']","")

    df = pd.DataFrame.from_dict(dd, orient='index')

    df = df.reset_index()
    df.columns = range(df.shape[1])
    df = df.stack()
    df = pd.DataFrame(df).reset_index(drop=True)

    return df.to_csv(outdir + '/generated_sequences/' + str(sample) + '_simulated.fa', header=None, index=False, sep='\t')

def index_and_chrm_sizes(verbose, outdir, sample, genome, seq_type):
    '''generates index and chromosome size files for the genome of interest'''
    if seq_type == 'whole_genome':
        g = genome.split('/')[-1].split('.')[-2]
        os.system('cut -f1,2 ' + genome + '.fai > ' + outdir + '/generated_sequences/' + g + '.chrom.sizes')
        if (path.exists(outdir + '/generated_sequences/' + g + '.chrom.sizes') == False):
            print('Creation of chromosome size file failed.')
            sys.exit(1)
        else:
            if verbose == True:
                print('Successfully created chromosome size file %s' % outdir + '/generated_sequences/' + g + '.chrom.sizes')
    
    else:
        os.system('samtools faidx ' + outdir + '/generated_sequences/' + sample + '_' + seq_type + '.fa')
        if (path.exists(outdir + '/generated_sequences/' + sample + '_' + seq_type + '.fa.fai') == False):
            print('Creation of index file failed. Make sure samtools is installed.')
            sys.exit(1)
        else:
            if verbose == True:
                print('Successfully created index file %s' % outdir + '/generated_sequences/' + sample + '_' + seq_type + '.fa.fai')

        os.system('cut -f1,2 ' + outdir + '/generated_sequences/' + sample + '_' + seq_type + '.fa.fai > ' + outdir + '/generated_sequences/' + sample + '_' + seq_type + '.chrom.sizes')

        if (path.exists(outdir + '/generated_sequences/' + sample + '_' + seq_type + '.chrom.sizes') == False):
            print('Creation of chromosome size file failed.')
            sys.exit(1)
        else:
            if verbose == True:
                print('Successfully created chromosome size file %s' % outdir + '/generated_sequences/' + sample + '_' + seq_type + '.chrom.sizes')      
    
           
def reformat_expt_fasta(verbose, outdir, sample, chrom_num, window):
    '''reformats the experimental sequences pulled from the whole genome
    by get_sequences() into a more 'genome-like' orientation on large chromosomes'''
    df = pd.read_csv(outdir + '/temp/' + sample + '_experimental_window_sequences.fa', header=None)
    df = df[df.index % 2 != 0]
    l_expt_seq = df.values.tolist()
    flat_list = [item for sublist in l_expt_seq for item in sublist]
    strN = str('N')*20
    file = ['{}{}'.format(i,strN) for i in flat_list]
    
    
    sequence_num = len(file)
    if verbose == True:
        print('Verifying the reformating of the experimental fasta...')
    seq_per_chrom, seq_last_chrom = define_seq_structure(verbose, sequence_num=sequence_num, chrom_num=chrom_num, window=window) 
    
    dd = {}
    for i in range(0,chrom_num-1):
        dd['>exp_' + str(i + 1)] = str(file[(seq_per_chrom)*i:(seq_per_chrom)*(i+1)]).replace("', '","").replace("['","").replace("']","")
    dd['>exp_' + str(chrom_num)] = str(file[(seq_per_chrom)*(chrom_num-1):]).replace("', '","").replace("['","").replace("']","")
    
    df = pd.DataFrame.from_dict(dd, orient='index')
    df.columns = range(df.shape[1])
    df = df.reset_index()
    df = df.stack()
    df = pd.DataFrame(df).reset_index(drop=True)
    df.to_csv(outdir + '/generated_sequences/' + str(sample) + '_experimental.fa', header=None, index=False, sep='\t')

def generate_bed(verbose, outdir, sample, sequence_num, chrom_num, window, annotation, seq_type):
    '''generates bed files defining mu for the simulated and experimental 'genomes'
    Calls windower'''
    if seq_type == 'experimental': 
        #getting the sequence number
        annotation_file = pd.read_csv(annotation, header=None)
        sequence_num = len(annotation_file)
    elif seq_type == 'simulated':
        sequence_num = sequence_num
    if verbose == True:
        print('Verifying the structure of the bed files...')
    seq_per_chrom, seq_last_chrom = define_seq_structure(verbose, sequence_num, chrom_num, window)
    
    #creating the locations
    dbed = {}
    for i in range(0,(seq_per_chrom)):
        dbed['r_' + str(i)] = [int(window-1)+i*(2*window+21), int(window+1)+i*(2*window+21)]

    dfbed = pd.DataFrame.from_dict(dbed, orient='index')
    dfbed.columns = ['start', 'stop']
    #duplicating the dataframe for all chromosomes except the last
    dfbed = pd.concat([dfbed]*(chrom_num-1))

    #creating the locations - last chrom
    dlast = {}
    for i in range(0,(seq_last_chrom)):
        dlast['z_' + str(i)] = [int(window-1)+i*(2*window+21), int(window+1)+i*(2*window+21)]
    dflast = pd.DataFrame.from_dict(dlast, orient='index')
    dflast.columns = ['start', 'stop']    

    #adding the chromosome name to the genomic location
    if seq_type == 'experimental':
        l=[]
        for i in range(int(seq_per_chrom)):
            for c in range(int(chrom_num-1)):
                l.append(str('exp_' + str(c+1)))
        l = sorted(l)
        dfbed['chr'] = l
    elif seq_type == 'simulated':
        l=[]
        for i in range(int(seq_per_chrom)):
            for c in range(int(chrom_num-1)):
                l.append(str('sim_' + str(c+1)))
        l = sorted(l)
        dfbed['chr'] = l
    
    #adding the chromosome name to the genomic location- last chrom
    if seq_type == 'experimental':
        l=[]
        for i in range(int(seq_last_chrom)):
            l.append(str('exp_' + str(chrom_num)))
        l = sorted(l)
        dflast['chr'] = l
    elif seq_type == 'simulated':
        l=[]
        for i in range(int(seq_last_chrom)):
            l.append(str('sim_' + str(chrom_num)))
        l = sorted(l)
        dflast['chr'] = l
    
    #combining the bed locations
    df = pd.concat([dfbed, dflast])
    df = df[['chr', 'start', 'stop']]
    df['count'] = (np.arange(len(df)))
    df['count'] = (df['count']+1).apply(str)
    df['region_name'] = df['chr'] + ';region_' + df['count']
    df.drop(['count'], axis=1, inplace=True)
    #saving the new annotation
    df.to_csv(outdir + '/annotations/' + str(sample) + '_' + seq_type + '_centered.bed', header=None, index=False, sep='\t')
    pull_bed= outdir + '/annotations/' + str(sample) + '_' + seq_type + '_centered.bed'
    windower(bed=pull_bed, outdir=outdir, sample=sample, window=window, seq_type=seq_type, ident=True)

######################################### Functions Called by Other Functions #########################################    
def windower(bed, outdir, sample, window, seq_type, ident):
    '''this function windows bedfiles'''
    bed_df = pd.read_csv(bed, sep ='\t',header=None)
    if ident == True:
        bed_df = bed_df.loc[:, 0:3]
        bed_df.columns = ['chr', 'start', 'stop', 'region_name']    
    elif ident==False:
        bed_df = bed_df.loc[:, 0:2]
        bed_df.columns = ['chr', 'start', 'stop']
        

    ##redefine mu to get new start and stop coordinates
    bed_df['start_new'] = bed_df.apply(lambda x: round((x['start'] + x['stop'])/2), axis=1)

    bed_df['stop_new'] = bed_df.apply(lambda x: x['start_new'] + 1, axis = 1)

    ##the -1500 position from 'origin'
    bed_df['start'] = bed_df.apply(lambda x: x['start_new'] - int(window), axis=1)

    ##the 1500 position from the 'origin'
    bed_df['stop'] = bed_df.apply(lambda x: x['stop_new'] + int(window), axis=1)

    ##saving the new annotation
    if ident == True:
        bed_df.to_csv(outdir + '/annotations/' + sample + '_' + seq_type + '_window.bed', sep='\t',
                columns=['chr','start','stop', 'region_name'],
                header = False, index = False)
    elif ident == False:
        bed_df.to_csv(outdir + '/annotations/' + sample + '_' + seq_type + '_window.bed', sep='\t',
                    columns=['chr','start','stop'],
                    header = False, index = False)

def define_seq_structure(verbose, sequence_num, chrom_num, window):
    '''defines chromosome structuring for generated 'genomes' and bed files'''
    seq_length = ((2*int(window))+21)
    tot_bases = int(sequence_num*seq_length)
    seq_per_chrom = math.floor(sequence_num/chrom_num)
    seq_last_chrom = seq_per_chrom + sequence_num%chrom_num
    #this is to check if there is the correct amount of sequences per chromosome
    scount = seq_last_chrom + (seq_per_chrom)*(chrom_num-1)

    base_per_chrom = seq_length*seq_per_chrom 
    base_last_chrom = seq_length*seq_last_chrom
    #this is to check if there is the correct amount of bases per chromosome
    bscount = base_last_chrom + (base_per_chrom)*(chrom_num-1)

    if verbose == True:
        print('Individual Sequence Length: ' + str(seq_length))
        print('Number of Sequences: ' + str(sequence_num)) 
        print('Number of Chromosomes: ' + str(chrom_num) + '\n')
        print('Sequences per Chromosome: ' + str(seq_per_chrom))
        print('Sequences on Last Chromosome: ' + str(seq_last_chrom) + '\n')
        print('Bases per Chromosome: ' + str(base_per_chrom))
        print('Bases on Last Chromosome: ' + str(base_last_chrom))
        print('Total Bases: ' + str(tot_bases))
        
    return seq_per_chrom, seq_last_chrom