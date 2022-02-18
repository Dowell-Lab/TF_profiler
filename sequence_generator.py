######################################### Imports #########################################
import os
from os import stat
from os import path
import sys
import pandas as pd
import numpy as np
import math
import time
from functools import partial
import multiprocessing

import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

import warnings
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

######################################### Sequence Generator Main #########################################
def run_sequence_generator(verbose, outdir, sample, genome, annotation,
                           sequence_num, chrom_num, cpus, seed, window):    
    if verbose == True:
        start_time = float(time.time())
        print('--------------Expanding Windows and Extracting Sequences---------------')
    window_annotation(verbose, annotation=annotation, 
                      outdir=outdir, sample=sample, window=window)
    get_sequences(verbose, genome=genome, outdir=outdir, sample=sample, plot_dinucleotide=False)
    ls = list_sequences(outdir=outdir, sample=sample)

    if sequence_num is not None:
        sequence_num = int(sequence_num)
    else:
        sequence_num = int(len(ls))
        if verbose == True:
            print('There are ' + str(sequence_num) + ' bidirectionals in the annotation file provided.')
            print(str(sequence_num) + ' sequences will be simulated.')
    if chrom_num is not None:
        chrom_num = int(chrom_num)
    else:
        chrom_num = int(math.ceil(sequence_num*(window*2+21)/1150000)) #11500000
        if verbose == True:
            print('The 10% of the average size of a human chromosome is ~11,500,000 bases.')
            print('Therefore, we are generating ' + str(chrom_num) + ' artificial chromosomes.')    
    
    if verbose == True:
        print('We are generating ' + str(sequence_num) + ' ' + str(window*2) + ' base long sequences on ' + str(chrom_num) + ' artificial chromosomes.')
    
    rs_list=set_seed(verbose=verbose, seed=seed, cpus=cpus, sequence_num=sequence_num)
    return ls, rs_list, sequence_num, chrom_num
    
def run_dinucleotide_generator(verbose, outdir, sample,
                           sequence_num, chrom_num, cpus, window, ls, rs_list):
    if verbose == True:
        start_prob = float(time.time())
        print('-----------Calculating Dinucleotide Dependant Base Position Probabilities-------------')
    first_position = first_dinuc(ls=ls, outdir=outdir, sample=sample)
    position_givenA=conditional_probability_counter(nuc='A', ls=ls, window=window, outdir=outdir, sample=sample)
    position_givenC=conditional_probability_counter(nuc='C', ls=ls, window=window, outdir=outdir, sample=sample)
    position_givenG=conditional_probability_counter(nuc='G', ls=ls, window=window, outdir=outdir, sample=sample)
    position_givenT=conditional_probability_counter(nuc='T', ls=ls, window=window, outdir=outdir, sample=sample)
    if verbose == True:
        stop_prob = float(time.time())
        prob_time = (stop_prob-start_prob)/60
        print('Probability calculations took ' + str(prob_time) + ' min.') 
        print('--------------------Generating Dinucleotide Sequences----------------------')
        print('Initiating ' + str(cpus) + ' cpus.')
    pool = multiprocessing.Pool(cpus)
    seqs = pool.map(partial(sequence_generator, 
                           inputs=[verbose, window, first_position, position_givenA, position_givenC, position_givenG, position_givenT]), rs_list)
    pool.close()
    pool.join()
    seqs=[seq for seqs_per_cpu in seqs for seq in seqs_per_cpu]
    write_fasta(verbose, generated_sequences=seqs, sample=sample, outdir=outdir, 
                window=window, sequence_num=sequence_num, chrom_num=chrom_num, seq_type='simulated')
    if verbose == True:
        stop_time = float(time.time())
        seq_time = (stop_time-stop_prob)/60
        print('Dinucleotide sequence generation took ' + str(seq_time) + ' min.') 
        print('-----------Indexing Dinucleotide Sequences and Creating Bedfiles--------------------')
    index_and_chrm_sizes(verbose, outdir=outdir, sample=sample, seq_type='simulated')
    generate_bed(verbose, outdir=outdir, sample=sample, 
                 sequence_num=sequence_num, chrom_num=chrom_num, 
                 window=window, ls=ls, seq_type='simulated')

    if verbose == True:
        print('-----------Plotting positional biases-----------')
    ###plotting mononucleotide probabilites from the experimental data
    mono_probabilities=mono_probability_counter(ls=ls, window=window, outdir=outdir, sample=sample)
    plot_positional_bias(outdir=outdir, sample=sample, window=window, 
                         base='not_conditional', probabilities=mono_probabilities)

    ###plotting conditional probabilities from the experimental data
    bases=['A','C','G','T']
    conditional_probabilities_list= [position_givenA,position_givenC,position_givenG,position_givenT]
    for i,conditional_probabilities in enumerate(conditional_probabilities_list):
        plot_positional_bias(outdir=outdir, sample=sample, window=window,
                     base=bases[i], probabilities=conditional_probabilities)

    ###plotting mononucleotide probabilities from the new dinucleotide simulated data
    get_sequences(verbose=verbose, genome=(outdir + '/generated_sequences/' + str(sample) + '_simulated.fa'), 
                  outdir=outdir, sample=sample, plot_dinucleotide=True)
    dbls = list_sequences(outdir=outdir, sample='dinucleotide')
    mono_probabilities_from_dinucleotide_sequences = mono_probability_counter(ls=dbls, window=window, outdir=outdir,
                                                 sample='dinucleotide')
    plot_positional_bias(outdir=outdir, sample='dinucleotide', window=window,
                         base='not_conditional', probabilities=mono_probabilities_from_dinucleotide_sequences)
        
def run_mononucleotide_generator(verbose, outdir, sample,
                           sequence_num, chrom_num, cpus, window, ls, rs_list):
    if verbose == True:
        start_prob = float(time.time())
        print('-----------Calculating Mononucleotide Dependant Base Position Probabilities-------------')
    mono_probabilities=mono_probability_counter(ls=ls, window=window, outdir=outdir, sample=sample)
    plot_positional_bias(outdir=outdir, sample=sample, window=window, 
                         base='not_conditional', probabilities=mono_probabilities)
    if verbose == True:
        stop_prob = float(time.time())
        prob_time = (stop_prob-start_prob)/60
        print('Probability calculations took ' + str(prob_time) + ' min.') 
        print('--------------------Generating Mononucleotide Sequences----------------------')
        print('Initiating ' + str(cpus) + ' cpus.')
    pool = multiprocessing.Pool(cpus)
    seqs = pool.map(partial(mono_sequence_generator, 
                           inputs=[window, mono_probabilities]), rs_list)
    pool.close()
    pool.join()        
    seqs=[seq for seqs_per_cpu in seqs for seq in seqs_per_cpu]
    write_fasta(verbose, generated_sequences=seqs, sample=sample, outdir=outdir, 
                window=window, sequence_num=sequence_num, chrom_num=chrom_num, seq_type='mononucleotide_simulated')
    if verbose == True:
        stop_time = float(time.time())
        seq_time = (stop_time-stop_prob)/60
        print('Mononucleotide sequence generation took ' + str(seq_time) + ' min.') 
        print('-----------Indexing Mononucleotide Sequences and Creating Bedfiles--------------------')
    index_and_chrm_sizes(verbose, outdir=outdir, sample=sample, seq_type='mononucleotide_simulated')
    generate_bed(verbose, outdir=outdir, sample=sample, sequence_num=sequence_num, chrom_num=chrom_num, 
             window=window, ls=ls, seq_type='mononucleotide_simulated')

def run_experimental_fimo_formatter(verbose, outdir, sample, window, ls):
    if verbose == True:
        print('-----------Formating Experimental Sequences for FIMO Scan--------------------')
    chrom_num_exp = int(math.ceil(len(ls)*(window*2+21)/1150000))
    write_fasta(verbose, generated_sequences=ls, sample=sample, outdir=outdir, 
                window=window, sequence_num=len(ls), chrom_num=chrom_num_exp, seq_type='experimental')
    if verbose == True:
        print('-----------Indexing Experimental Sequences and Creating Bedfiles--------------------')
    index_and_chrm_sizes(verbose, outdir=outdir, sample=sample, seq_type='experimental')
    generate_bed(verbose, outdir=outdir, sample=sample, sequence_num=len(ls), chrom_num=chrom_num_exp, 
             window=window, ls=ls, seq_type='experimental')

######################################### Sequence Generator Functions #########################################
################################################## Setting up ##################################################
def set_seed(verbose, seed, cpus, sequence_num):
    '''Sets numpy seed for the sequence generation'''
    if seed is not None:
        if verbose == True:
            print('-----------Setting Seed for the Sequence Generator.-----------')
            print('User defined seed: '+ str(seed))
        try:
            seed=int(seed)
        except:
            if verbose==True:
                print('User defined seed is not a number! Setting seed to time.')
            seed=global_seed(verbose=verbose)
    else:
        seed=global_seed(verbose=verbose)

    if seed==0:
        if verbose==True:
            print('Seed is set to 0, this cannot be reproduced. Resetting seed to time.')
        seed=global_seed(verbose=verbose)
    
    rng_list=[]
    for c in range(cpus):
        rng = np.random.default_rng(seed+c) 
        rng_list.append(rng)
        
    ### set how many sequences will be generated per cpu
    ### this number of sequences will be paired with a rng 
    seq_per_cpu = math.floor(sequence_num/cpus)
    seq_last_cpu = math.floor(sequence_num/cpus) + sequence_num%cpus
    if cpus==1:
        spc_list=[seq_per_cpu]
        print('Only one cpu allocated to generate ' + str(sequence_num) + ' sequences.')
    elif cpus > 1:
        if verbose == True:
            print('We are generating ' + str(seq_per_cpu) + ' sequences per cpu.')
            print('With ' + str(sequence_num%cpus) + ' additional sequence(s) on the last cpu.')
        spc_list=[seq_per_cpu]*(cpus-1)
        spc_list.append(seq_last_cpu)

    ### combining the lists to loop through during sequence generation
    rs_list=[]
    for i,rng in enumerate(rng_list):
        rs=[rng,spc_list[i]]
        rs_list.append(rs)
    return rs_list
                
def global_seed(verbose):
    '''function used by set seed to set seed to the time'''
    if verbose==True:
        print('-----------Setting Seed for the Sequence Generator------------')
        t= time.time()
        print('Time in fractional seconds of:', str(t))
        print('Setting seed on the clock: ' + str(int(t)))
    return(int(t))

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
            print(outdir + '/annotations exists.')     
    try:
        os.system('mkdir -p ' + outdir + '/temp')
    except OSError:
        print ('Creation of the directory %s failed' % outdir + '/temp')
        sys.exit(1)
    else:
        if verbose == True:
            print(outdir + '/temp exists')
    try:
        os.system('mkdir -p ' + outdir + '/generated_sequences')
    except OSError:
        if verbose == True:
            print ('Creation of the directory %s failed' % outdir + '/generated_sequences')
            sys.exit(1)
    else:
        if verbose == True:
            print(outdir + '/generated_sequences exists')
    try:
        os.system('mkdir -p ' + outdir + '/plots')
    except OSError:
        print ('Creation of the directory %s failed' % outdir + '/plots')
        sys.exit(1)
    else:
        if verbose == True:
            print(outdir + '/plots exists.')  
    windower(bed=annotation, outdir=outdir, sample=sample, window=int(window), spef_dir='/temp/', seq_type='experimental', ident=False)

def get_sequences(verbose, genome, outdir, sample, plot_dinucleotide):
    '''This function pulls the sequences out of the windowed annotation file and outputs them in 
    generated_sequences for further use'''
    if plot_dinucleotide==True:
        annotation_file=outdir + '/annotations/' + sample + '_simulated_window.bed'
        os.system('bedtools getfasta -fi ' + genome + 
          ' -bed ' + annotation_file + 
          ' -fo ' + outdir + '/temp/dinucleotide_window_sequences.fa')
        if os.stat(outdir + '/temp/dinucleotide_window_sequences.fa').st_size == 0:
            print('Extraction for plotting failed.')
    else:
        if (path.exists(outdir + '/temp/' + sample + '_window_sequences.fa') == False) or (os.stat(outdir + '/temp/' + sample + '_window_sequences.fa').st_size == 0):
            annotation_file=outdir + '/temp/' + sample + '_experimental_window.bed'   
            os.system('bedtools getfasta -fi ' + genome + 
                      ' -bed ' + annotation_file + 
                      ' -fo ' + outdir + '/temp/' + sample + '_window_sequences.fa')
            if (path.exists(outdir + '/temp/' + sample + '_window_sequences.fa') == False):
                print('Extracted experimental sequences failed. Make sure bedtools/2.25.0 is installed.')
                sys.exit(1)
            elif os.stat(outdir + '/temp/' + sample + '_window_sequences.fa').st_size == 0:
                print('Extracted experimental sequences failed. The output file is empty!')
                sys.exit(1)
            else:
                if verbose == True:
                    print('Extracted Experimental Sequences Output: ' + outdir + '/temp/' + sample + '_window_sequences.fa')
        else:
            print('Experimental Sequences were previously extracted.')


def list_sequences(outdir, sample):
    '''This function strings experimental sequences together in a list'''
    fasta_sequences = outdir + '/temp/' + sample + '_window_sequences.fa'
    sequences = []
    sequence_names = []                                                                                                        
    with open(fasta_sequences) as fa:
        for line in fa:
            line = line.strip('\n')
            if '>' not in line:
                sequences.append(line)
            else:
                sequence_names.append(line)      
    return sequences

######################################### Dinucleotide Sequences #########################################
def first_dinuc(ls, outdir, sample):
    '''For dinucleotide dependant sequence generation. This function calculates the probability of the first dinucleotide pair in the window'''        
    dinuc_list = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
                  'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    first_position=[1]*len(dinuc_list)
    for j in ls:
        for n in range(len(dinuc_list)):
            if j[0]+j[1] == dinuc_list[n]:
                first_position[n]=first_position[n]+1
            else:
                continue

    first_position = [x / (sum(first_position)) for x in first_position]

    df = pd.DataFrame(first_position, index=dinuc_list)
    df.columns=['dinucleotide_probabilities']
    df.to_csv(outdir + '/generated_sequences/' + sample + '_position1_dinucleotide_probabilities.tsv', sep='\t')
    return first_position

def conditional_probability_counter(nuc, ls, window, outdir, sample):
    '''For dinucleotide dependant sequence generation. This function calculates the probability of the next nucleotide in the sequence given information on the previous nucleotide'''
    ll=[]
    seq_length = window*2+1
    pos=[x+1 for x in range(1,seq_length)]

    for i in range(1,seq_length-1):
        ll.append([1,1,1,1])
        for j in ls:
            if j[i] == nuc:
                seq_counter = j
            else:
                continue

            if seq_counter[i+1] == 'A':
                ll[i-1][0]=ll[i-1][0]+1
            if seq_counter[i+1] == 'C':
                ll[i-1][1]=ll[i-1][1]+1
            if seq_counter[i+1] == 'G':
                ll[i-1][2]=ll[i-1][2]+1
            if seq_counter[i+1] == 'T':
                ll[i-1][3]=ll[i-1][3]+1
            else:
                continue

        ll[i-1] = [x / (sum(ll[i-1])) for x in ll[i-1]]
    base_df = pd.DataFrame.from_records(ll, columns=['A','C','G','T'])
    base_df.to_csv(outdir + '/generated_sequences/' + sample +'_conditional_probabilites_given'+ nuc + '.tsv', sep='\t', index=False)
    return ll

def sequence_generator(rs_list, inputs):
    '''For dinucleotide dependant sequence generation. This function generates the simulated sequences'''
    verbose, window, first_position, position_givenA, position_givenC, position_givenG, position_givenT = inputs
    rng=rs_list[0]
    sequence_num=rs_list[1]
    
    print(rng)
    print(sequence_num)
              
    dinuc_list = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
                  'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    bases=['A','C','G','T']
                  
    #checking CONDITIONAL sequence lengths
    #should be total sequence length minus 2 since the first 2 are generated from first position
    if (len(position_givenA) == window*2-1) == False:
        sys.exit()
    elif (len(position_givenA) == len(position_givenC)) == False:
        sys.exit()
    elif (len(position_givenA) == len(position_givenG)) == False:
        sys.exit()
    elif (len(position_givenA) == len(position_givenT)) == False:
        sys.exit()
    else:
        sequence_length= (window*2+1)

    sequences = np.empty([sequence_num, sequence_length], dtype=str)
    column = rng.choice(dinuc_list, sequence_num, p=first_position)    

    for i in range(sequence_num):
        sequences[:,0][i] = column[i][0]
        sequences[:,1][i] = column[i][1]
        for j in range(2,sequence_length):
            if str(sequences[:,j-1][i]) == 'A':
                choice=(rng.choice(bases, 1, p=position_givenA[j-2]))
                sequences[:,j][i] = choice[0]
            elif str(sequences[:,j-1][i]) == 'C':
                choice=(rng.choice(bases, 1, p=position_givenC[j-2]))
                sequences[:,j][i] = choice[0]
            elif str(sequences[:,j-1][i]) == 'G':
                choice=(rng.choice(bases, 1, p=position_givenG[j-2]))
                sequences[:,j][i] = choice[0]
            elif str(sequences[:,j-1][i]) == 'T':
                choice=(rng.choice(bases, 1, p=position_givenT[j-2]))
                sequences[:,j][i] = choice[0]
            else:
                print('ERROR')
                sys.exit()

    joined_sequences = [''.join(row) for row in sequences]
    return joined_sequences

######################################### Mononucleotide Sequences #########################################
def mono_probability_counter(ls, window, outdir, sample):
    '''For mononucleotide dependant sequence generation. This function calculates the probability of each nucleotide at a given position'''
    ll=[]
    seq_length = window*2+1
    ##for each positions count the occurance of each base                                                                                                                        
    ##across all sequences in the the input list                                                                                                                                 
    for i in range(seq_length):
        ll.append([1,1,1,1])
        for j in ls:
            if j[i] == 'A':
                ll[i][0]=ll[i][0]+1
            if j[i] == 'C':
                ll[i][1]=ll[i][1]+1
            if j[i] == 'G':
                ll[i][2]=ll[i][2]+1
            if j[i] == 'T':
                ll[i][3]=ll[i][3]+1
            else:
                continue
        ll[i] = [x / (sum(ll[i])) for x in ll[i]]
    if sample == 'dinucleotide':
        print('Calculating positional base composition based on dinucleotide generated sequences.')
    base_df = pd.DataFrame.from_records(ll, columns=['A','C','G','T'])
    base_df.to_csv(outdir + '/generated_sequences/' + sample +'_mononucleotide_probabilites.tsv', sep='\t', index=False)
    return ll

def mono_sequence_generator(rs_list, inputs):
    '''For mononucleotide dependant sequence generation. This function generates the simulated sequences'''
    window, mono_probabilities = inputs
    
    rng=rs_list[0]
    sequence_num=rs_list[1]
    
    print(rng)
    print(sequence_num)
    
    bases=['A','C','G','T']
                              
    #checking CONDITIONAL sequence lengths
    #should be total sequence length minus 2 since the first 2 are generated from first position

    if (len(mono_probabilities) == window*2+1) == False:
        sys.exit()
    else:
        sequence_length= (window*2+1)

    sequences = np.empty([sequence_num, sequence_length], dtype=str)

    for i in range(sequence_length):
        column = rng.choice(bases, sequence_num, p=mono_probabilities[i])
        sequences[:,i] = column

    joined_sequences = [''.join(row) for row in sequences]
    return joined_sequences    

######################################### Saving and Formating Sequences #########################################
def write_fasta(verbose, generated_sequences, sample, outdir, window, sequence_num, chrom_num, seq_type):
    ''' writes sequences generated into fasta file format and outputs them in generated sequences'''
    seq_per_chrom, seq_last_chrom = define_seq_structure(verbose, sequence_num=sequence_num, chrom_num=chrom_num, window=window)
    if seq_type =='simulated' or seq_type =='mononucleotide_simulated':
        chrom = '>sim'
    elif seq_type == 'experimental':
        chrom= '>exp'
    else:
        chrom='>chr'

    file=[]
    for i in range(len(generated_sequences)):
        file.append(str(generated_sequences[i]) + str('N')*20)      

    dd = {}
    for i in range(0,chrom_num-1):
        dd[chrom + str(i + 1)] = str(file[(seq_per_chrom)*i:(seq_per_chrom)*(i+1)]).replace("', '","").replace("['","").replace("']","")
    dd[chrom + str(chrom_num)] = str(file[(seq_per_chrom)*(chrom_num-1):]).replace("', '","").replace("['","").replace("']","")

    df = pd.DataFrame.from_dict(dd, orient='index')

    df = df.reset_index()
    df.columns = range(df.shape[1])
    df = df.stack()
    df = pd.DataFrame(df).reset_index(drop=True)
    df.to_csv(outdir + '/generated_sequences/' + str(sample) + '_' + seq_type + '.fa', header=None, index=False, sep='\t')

def index_and_chrm_sizes(verbose, outdir, sample, seq_type):
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

def generate_bed(verbose, outdir, sample, sequence_num, chrom_num, window, ls, seq_type):
    '''generates bed files defining mu for the simulated and experimental 'genomes'
    Calls windower'''
    if seq_type == 'simulated' or seq_type == 'mononucleotide_simulated':
        sequence_num = sequence_num
        chrom='sim'
    elif seq_type == 'experimental': 
        #getting the sequence number
        sequence_num = len(ls)
        chrom='exp'
    else:
        sequence_num = sequence_num
        chrom='chr'
    
    if verbose == True:
        print('Verifying the structure of the bed files...')
    seq_per_chrom, seq_last_chrom = define_seq_structure(verbose, sequence_num, chrom_num, window)
        
    #creating the locations
    dbed = {}
    for i in range(0,(seq_per_chrom)):
        dbed['r_' + str(i)] = [int(window-1)+i*(2*window+21), int(window+1)+i*(2*window+21)]

    dfbed = pd.DataFrame.from_dict(dbed, orient='index')
    dfbed.columns = ['start', 'stop']
    #list for chromosome names
    l=[]
    if chrom_num == 1:
        df=dfbed
        #defining chromosom names
        for i in range(int(seq_per_chrom)):
            l.append(str(chrom + str(chrom_num)))
    else:
        #duplicating the dataframe for all chromosomes except the last
        dfbed = pd.concat([dfbed]*(chrom_num-1))
        #creating the locations - last chrom
        dlast = {}
        for i in range(0,(seq_last_chrom)):
            dlast['z_' + str(i)] = [int(window-1)+i*(2*window+21), int(window+1)+i*(2*window+21)]
        dflast = pd.DataFrame.from_dict(dlast, orient='index')
        dflast.columns = ['start', 'stop'] 

        #combining the bed locations
        df = pd.concat([dfbed, dflast])
    
        #adding the chromosome name to the genomic location
        l=[]
        for i in range(int(seq_per_chrom)):
            for c in range(int(chrom_num-1)):
                l.append(str(chrom + str(c+1)))
        l = sorted(l)
        #adding the chromosome name to the genomic location- last chrom
        for i in range(int(seq_last_chrom)):
            l.append(str(chrom + str(chrom_num)))
        
    df['chr'] = l
    df = df[['chr', 'start', 'stop']]
    df['count'] = (np.arange(len(df)))
    df['count'] = (df['count']+1).apply(str)
    df['region_name'] = df['chr'] + ';region_' + df['count']
    df.drop(['count'], axis=1, inplace=True)
    df = df.sort_values(by=['chr', 'start'])
    #saving the new annotation
    df.to_csv(outdir + '/annotations/' + str(sample) + '_' + seq_type + '_centered.bed', header=None, index=False, sep='\t')
    pull_bed= outdir + '/annotations/' + str(sample) + '_' + seq_type + '_centered.bed'
    windower(bed=pull_bed, outdir=outdir, sample=sample, window=window, spef_dir='/annotations/', seq_type=seq_type, ident=True)

######################################### Functions Called by Other Functions #########################################    
def windower(bed, outdir, sample, window, seq_type, spef_dir, ident):
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
    bed_df['stop'] = bed_df.apply(lambda x: x['stop_new'] + int(window), axis=1)
    
    bed_df = bed_df.sort_values(by=['chr', 'start'])
    bed_df = bed_df[bed_df['start'] >= 0]
    ##saving the new annotation
    
    if ident == True:
        bed_df.to_csv(outdir + spef_dir + sample + '_' + seq_type + '_window.bed', sep='\t',
                columns=['chr','start','stop', 'region_name'],
                header = False, index = False)
    elif ident == False:
        bed_df.to_csv(outdir + spef_dir + sample + '_' + seq_type + '_window.bed', sep='\t',
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

#     if verbose == True:
#         print('Individual Sequence Length: ' + str(seq_length))
#         print('Number of Sequences: ' + str(sequence_num)) 
#         print('Number of Chromosomes: ' + str(chrom_num) + '\n')
#         print('Sequences per Chromosome: ' + str(seq_per_chrom))
#         print('Sequences on Last Chromosome: ' + str(seq_last_chrom) + '\n')
#         print('Bases per Chromosome: ' + str(base_per_chrom))
#         print('Bases on Last Chromosome: ' + str(base_last_chrom))
#         print('Total Bases: ' + str(tot_bases))

    return seq_per_chrom, seq_last_chrom

###################################### Plotting positional bias ######################################
def plot_positional_bias(outdir, sample, window, base, probabilities):
    ##get positions                                                   
    starting=int(window)
    
    if base == 'not_conditional':
        stopping=int(window)+1
        out_name='single_position'
        
    else:
        stopping=int(window)-1
        out_name='probability_given'+base
    
    positions = np.arange(-starting,stopping,1)

    ### Line plot for each base on one grid ###
    plt.figure(figsize=(12,10))
    gs = plt.GridSpec(1, 1)
    ax = plt.subplot(gs[0])
    color_list = ['blue', 'red', 'orange', 'purple']
    bases=['A','C','G','T']

    aprobs=[]
    cprobs=[]
    gprobs=[]
    tprobs=[]
    for i in range(len(probabilities)):
        aprobs.append(probabilities[i][0])
        cprobs.append(probabilities[i][1])
        gprobs.append(probabilities[i][2])
        tprobs.append(probabilities[i][3])
    probs=[aprobs,cprobs,gprobs,tprobs]

    for i,base in enumerate(bases):
        ax.plot(positions,probs[i],color=color_list[i], alpha=0.75, label=bases[i])

    ax.set_xlabel('Distance (bp)',fontsize=30,fontweight='bold')
    ax.set_ylabel('Base Content',fontsize=30,fontweight='bold')
    plt.xticks(fontsize = 25)
    plt.yticks(fontsize = 25)
    plt.legend(bbox_to_anchor=(1.05, 1),fontsize = 25, loc=2, borderaxespad=0.)
    plt.suptitle(sample+'_'+out_name, fontsize=40, fontweight='bold')
    plt.savefig(outdir + '/plots/' + sample + '_' + out_name + '_BaseDistribution.png',bbox_inches='tight')

    # ## Smooth the frequencies for better visualization ### 
    plt.figure(figsize=(12,10))
    gs1 = plt.GridSpec(1, 1)
    ax1 = plt.subplot(gs1[0])

    for i in range(0,4):
        smoothed_probs = savgol_filter(tuple(np.array(probs)), 61, 3) # window size 61, polynomial order 3     
        ax1.plot(positions,smoothed_probs[i],color=color_list[i], alpha=0.75, label=bases[i])
    ax1.set_xlabel('Distance (bp)',fontsize=30,fontweight='bold')
    ax1.set_ylabel('Base Content',fontsize=30,fontweight='bold')
    plt.xticks(fontsize = 25)
    plt.yticks(fontsize = 25)
    plt.legend(bbox_to_anchor=(1.05, 1),fontsize = 25, loc=2, borderaxespad=0.)
    plt.suptitle(sample+'_'+out_name, fontsize=40, fontweight='bold')
    plt.savefig(outdir + '/plots/' + sample + '_' + out_name + '_SmoothedBaseDistribution.png',bbox_inches='tight')
