import os
import pandas as pd
import math
from qc_and_keys import *

class FastaBedProcess:
    def __init__(self, outdir, sample, sequence_num, chrom_num, window, genome, annotation):
        self.outdir = outdir
        self.sample = sample
        self.sequence_num = sequence_num
        self.chrom_num = chrom_num
        self.window = window
        self.genome = genome
        self.annotation = annotation

    ##Author of get_windows: Rutendo F. Sigauke
    def get_windows(self):
        '''This function takes in bed files from Tfit and redefines mu and extends the window
        '''
        ##takes in Tfit regions/calls or any other bed file with chr, start, stop columns
        print("Annotation File: " + self.annotation)
        os.system("rsync " + self.annotation + ' ' + self.outdir + '/annotations/')

        bed = pd.read_csv(self.annotation, sep ='\t',header=None)

        ##select the coordinate colmuns only    
        bed_df = bed.loc[:, 0:2]
        bed_df.columns = ["chr", "start", "stop"]

        ##redefine mu to get new start and stop coordinates
        bed_df["start_new"] = bed_df.apply(lambda x: round((x["start"] + x["stop"])/2), axis=1)

        bed_df["stop_new"] = bed_df.apply(lambda x: x["start_new"] + 1, axis = 1)

        ##the -1500 position from "origin"
        bed_df["start"] = bed_df.apply(lambda x: x["start_new"] - int(self.window), axis=1)

        ##the 1500 position from the "origin"
        bed_df["stop"] = bed_df.apply(lambda x: x["stop_new"] + int(self.window), axis=1)

        ##saving the new annotation
        bed_df.to_csv(self.outdir + '/annotations/' + self.sample + '_experimental_window.bed', sep='\t',
                        columns=["chr","start","stop"],
                        header = False, index = False)  
        
#####Creating indexes and chromosome size files for downstream processing
    def index_and_chrm_sizes(self, seq_type):
        try:
            os.system("samtools faidx " + self.outdir + "/generated_sequences/" + self.sample + "_" + seq_type + ".fa")
        except OSError:
            print ("Creation of index file %s failed" % self.outdir + "/generated_sequences/" + self.sample + "_" + seq_type + ".fa.fai")
        else:
            print ("Successfully created index file %s" % self.outdir + "/generated_sequences/" + self.sample + "_" + seq_type + ".fa.fai")
    
        try:
            os.system("cut -f1,2 " + self.outdir + "/generated_sequences/" + self.sample + "_simulated.fa.fai > " + self.outdir + "/generated_sequences/" + self.sample + "/" + seq_type + ".fa.fai")
        except OSError:
            print ("Creation of chromosome size file %s failed" % self.outdir + "/generated_sequences/" + self.sample + "/" + seq_type + ".chrom.sizes")
        else:
            print ("Successfully created chromosome size file %s" % self.outdir + "/generated_sequences/" + self.sample + "/" + seq_type + ".chrom.sizes")   
            
    def chrm_sizes_whole_genome(self):
        g = self.genome.split('/')[-1].split('.')[-2]
        try:
            os.system("cut -f1,2 " + self.genome + ".fai > " + self.outdir + "/generated_sequences/" + g + ".chrom.sizes")
        except OSError:
            print ("Creation of chromosome size file %s failed" % self.outdir + "/generated_sequences/" + g + ".chrom.sizes")
        else:
            print ("Successfully created chromosome size file %s" % self.outdir + "/generated_sequences/" + g + ".chrom.sizes")

######For experimental genome if experimental_fimo == True, reformat the experimental fasta into a genome of only windowed bidirectionals                
    def reformat_expt_fasta(self):
        df = pd.read_csv(self.outdir + "/generated_sequences/" + self.sample + "_experimental_window_sequences.fa", header=None)
        df = df[df.index % 2 != 0]
        l_expt_seq = df.values.tolist()
        flat_list = [item for sublist in l_expt_seq for item in sublist]
        strN = str('N')*20
        file = ["{}{}".format(i,strN) for i in flat_list]

        sequence_num = len(file)
        seq_per_chrom, seq_last_chrom = check_seq_structure(sequence_num=sequence_num, chrom_num=self.chrom_num, window=self.window)

        dd = {}
        for i in range(0,self.chrom_num-1):
            dd['>exp_' + str(i + 1)] = str(file[(seq_per_chrom)*i:(seq_per_chrom)*(i+1)]).replace("', '","").replace("['","").replace("']","")
        dd['>exp_' + str(self.chrom_num)] = str(file[(seq_per_chrom)*(self.chrom_num-1):]).replace("', '","").replace("['","").replace("']","")
        df = pd.DataFrame.from_dict(dd, orient='index')
        df.columns = range(df.shape[1])
        df = df.reset_index()
        df = df.stack()
        df = pd.DataFrame(df).reset_index(drop=True)
        return df.to_csv(self.outdir + "/generated_sequences/" + str(self.sample) + "_experimental.fa", header=None, index=False, sep='\t')

######generate bedfiles for "new" genomes. This will produce the center and window
    def generate_bed(self, seq_type):
        if seq_type == "experimental": 
            #getting the sequence number
            annotation_file = pd.read_csv(self.annotation, header=None)
            sequence_num = len(annotation_file)
        elif seq_type == "simulated":
            sequence_num = self.sequence_num
        seq_per_chrom, seq_last_chrom = check_seq_structure(sequence_num, self.chrom_num, self.window)
        #creating the locations for the all chromosomes except the last
        dfbed=annotation_creater(seq_per_chrom=seq_per_chrom, chrom_num=self.chrom_num, window=self.window, seq_type=seq_type)
        #creating the locations for the last chromosome
        dflast=annotation_creater(seq_per_chrom=seq_last_chrom, chrom_num=self.chrom_num, window=self.window, seq_type=seq_type)
        #combining the bed locations
        df = pd.concat([dfbed, dflast])
        df = df[['chr', 'start', 'stop']]
        #saving the new annotation
        df.to_csv(self.outdir + "/annotations/" + str(self.sample) + "_" + seq_type + "_centered.bed", header=None, index=False, sep='\t')
        
        #creating the windowed bed
        df["start_new"] = df.apply(lambda x: round((x["start"] + x["stop"])/2), axis=1)
        df["stop_new"] = df.apply(lambda x: x["start_new"] + 1, axis = 1)

        ##the -1500 position from "origin"
        df["start"] = df.apply(lambda x: x["start_new"] - int(self.window), axis=1)

        ##the 1500 position from the "origin"
        df["stop"] = df.apply(lambda x: x["stop_new"] + int(self.window), axis=1)

        ##saving the new annotation
        df.to_csv(self.outdir + '/annotations/' + self.sample + '_' + seq_type + '_window.bed', sep='\t', columns=["chr","start","stop"], header = False, index = False)
        return print('Bed File Genereated: ' + self.outdir + "/annotations/" + str(self.sample) + "_" + seq_type + "_centered.bed and " + str(self.sample) + "_" + seq_type + "_window.bed")  
       
###function used to make bed files  -- called by generate bed  
def annotation_creater(seq_per_chrom, chrom_num, window, seq_type):
    dbed = {}
    for i in range(0,(seq_per_chrom)):
        dbed['r_' + str(i)] = [int(window-1)+i*(2*window+21), int(window+1)+i*(2*window+21)]

    dfbed = pd.DataFrame.from_dict(dbed, orient='index')
    dfbed.columns = ['start', 'stop']
    #duplicating the dataframe for all chromosomes except the last
    dfbed = pd.concat([dfbed]*(chrom_num-1))
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
    return dfbed