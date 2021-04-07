import math
import sys
import os

#dict_seq_type = {'simulated' : 'sim', 'experimental: exp'}

def check_seq_structure(sequence_num, chrom_num, window):
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

    print("Individual Sequence Length: " + str(seq_length))
    print("Number of Experimental Sequences: " + str(sequence_num)) 
    print("Number of Chromosomes: " + str(chrom_num) + "\n")
    print("Sequences per Chromosome: " + str(seq_per_chrom))
    print("Sequences on Last Chromosome: " + str(seq_last_chrom))
    if sequence_num-scount != 0:
        print('Error! The total number of sequences does not equal the sum of the sequences divided into chromosomes!')
        sys.exit(1)
    else:
        print('Number of sequences per chromosome are properly structured.\n')
    print("Bases per Chromosome: " + str(base_per_chrom))
    print("Bases on Last Chromosome: " + str(base_last_chrom))
    print("Total Bases in Experiment: " + str(tot_bases))
    if tot_bases-bscount != 0:
        print('Error! The total number of bases does not equal the sum of the bases divided into chromosomes!')
        sys.exit(1)
    else:
        print('Number of bases per chromosome are properly structured.\n')
    return seq_per_chrom, seq_last_chrom