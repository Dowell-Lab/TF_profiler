#Global seed, sequence input and sequence generator by Rutendo F. Sigauke 
import os
import numpy as np
import pandas as pd
import time
import math
from qc_and_keys import *

def global_seed():
    print("-----------Setting Seed for the Sequence Generator------------")
    print("Time in fractional seconds of:", time.time())
    print("Time as an interger:",int(time.time()))
    print("Setting seed on the clock")
    return(int(time.time()))
    print("np.random.random(",int(time.time()),") => ",np.random.random())


def sequence_input(sequence):
    '''takes in sequences from base_content that are in .csv format and 
    frequencies per position to use in simulating sequences using a 
    1st order Markov Model.
    '''

    position_feq = []

    with open(sequence) as seq:

        ##remove the header line
        lines = seq.readlines()[1:]
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
        
def sequence_generator(bases,position_feq, sequence_num):
    '''takes in frequencies per position and simulates sequences using a 
    1st order Markov Model.
    '''
    
    first = [position_feq[0:1]]
    last =  [position_feq[-1:]]
    print('Position Frequency: Contains ' + str(len(position_feq)) + " elements. " + str(first) + '...' + str(last))
    
    sequences = np.empty([sequence_num, len(position_feq)], dtype=str)
    
    for i in range(len(position_feq)):
        column = np.random.choice(bases, sequence_num, p=position_feq[i])
        sequences[:,i] = column
        
    joined_sequences = [''.join(row) for row in sequences]
    
    return joined_sequences

def write_fasta(generated_sequences, sample, outdir, window, sequence_num, chrom_num):
    ''' writes sequences generated into fasta file
    '''
    seq_per_chrom, seq_last_chrom = check_seq_structure(sequence_num=sequence_num, chrom_num=chrom_num, window=window)
    
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

    return df.to_csv(outdir + "/generated_sequences/" + str(sample) + "_simulated.fa", header=None, index=False, sep='\t')