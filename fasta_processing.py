import os
import pandas as pd
import math

class FastaProcess:
    def __init__(self, outdir, sample, sequence_num, chrom_num, window, genome, annotation):
        self.outdir = outdir
        self.sample = sample
        self.sequence_num = sequence_num
        self.chrom_num = chrom_num
        self.window = window
        self.genome = genome
        self.annotation = annotation
        
    def index_and_chrm_sizes(self):
        try:
            os.system("samtools faidx " + self.outdir + "/generated_sequences/" + self.sample + "_simulated.fa")
        except OSError:
            print ("Creation of index file %s failed" % self.outdir + "/generated_sequences/" + self.sample + "_simulated.fa.fai")
        else:
            print ("Successfully created index file %s" % self.outdir + "/generated_sequences/" + self.sample + "_simulated.fa.fai")

        try:
            os.system("cut -f1,2 " + self.outdir + "/generated_sequences/" + self.sample + "_simulated.fa.fai > " + self.outdir + "/generated_sequences/" + self.sample + "_simulated.chrom.sizes")
        except OSError:
            print ("Creation of chromosome size file %s failed" % self.outdir + "/generated_sequences/" + self.sample + "_simulated.chrom.sizes")
        else:
            print ("Successfully created chromosome size file %s" % self.outdir + "/generated_sequences/" + self.sample + "_simulated.chrom.sizes")
              
    def sim_bed(self):
        #calculating the number of sequences per chromosome
        seq_per_chrom = math.floor(self.sequence_num/self.chrom_num)
        seq_last_chrom = seq_per_chrom + self.sequence_num%self.chrom_num
        scount = seq_last_chrom + (seq_per_chrom)*(self.chrom_num-1)
        #checking fasta parameters
        seq_length = ((2*int(self.window))+21)
        tot_bases = int(self.sequence_num*seq_length)
        #this is to check if there is the correct amount of sequences per chromosome
        scount = seq_last_chrom + (seq_per_chrom)*(self.chrom_num-1)
        
        print("----------Checking bed parameters------------")
        print("Individual Sequence Length: " + str(seq_length))
        print("Number of Experimental Sequences: " + str(self.sequence_num)) 
        print("Number of Chromosomes: " + str(self.chrom_num) +'\n')
        print("Sequences per Chromosome: " + str(seq_per_chrom))
        print("Sequences on Last Chromosome: " + str(seq_last_chrom))
        if self.sequence_num-scount != 0:
            print('Error! The total number of sequences does not equal the sum of the sequences divided into chromosomes!')
            sys.exit(1)
        else:
            print('Number of sequences per chromosome are properly structured.\n')        
        
        #creating the locations for the all chromosomes except the last
        dbed = {}
        for i in range(0,(seq_per_chrom)):
            dbed['r_' + str(i)] = [int(self.window-1)+i*(2*self.window+21), int(self.window+1)+i*(2*self.window+21)]

        dfbed = pd.DataFrame.from_dict(dbed, orient='index')
        dfbed.columns = ['start', 'stop']
        #duplicating the dataframe for all chromosomes except the last
        dfbed = pd.concat([dfbed]*(self.chrom_num-1))
        #adding the chromosome name to the genomic location
        l_sim=[]
        for i in range(int(seq_per_chrom)):
            for c in range(int(self.chrom_num-1)):
                l_sim.append(str('sim_' + str(c+1)))
        l_sim = sorted(l_sim)
        dfbed['sim'] = l_sim

        #creating the locations for the last chromosome
        dlast = {}
        for i in range(0,(seq_last_chrom)):
            dlast['l_' + str(i)] = [int(self.window-1)+i*(2*self.window+21), int(self.window+1)+i*(2*self.window+21)]
        dflast = pd.DataFrame.from_dict(dlast, orient='index')
        dflast.columns = ['start', 'stop']
        #adding the chromosome name to the genomic location
        l_last = []
        for i in range(int(seq_last_chrom)):
            l_last.append(str('sim_' + str(self.chrom_num)))
        dflast['sim'] = l_last

        #combining the bed locations
        df = pd.concat([dfbed, dflast])
        df = df[['sim', 'start', 'stop']]
        #saving the new annotation
        df.to_csv(self.outdir + "/annotations/" + str(self.sample) + "_simulated_centered.bed", header=None, index=False, sep='\t')

        #creating the windowed bed
        df["start_new"] = df.apply(lambda x: round((x["start"] + x["stop"])/2), axis=1)
        df["stop_new"] = df.apply(lambda x: x["start_new"] + 1, axis = 1)

        ##the -1500 position from "origin"
        df["start"] = df.apply(lambda x: x["start_new"] - int(self.window), axis=1)

        ##the 1500 position from the "origin"
        df["stop"] = df.apply(lambda x: x["stop_new"] + int(self.window), axis=1)

        ##saving the new annotation
        df.to_csv(self.outdir + '/annotations/' + self.sample + '_simulated_window.bed', sep='\t', columns=["sim","start","stop"], header = False, index = False)
        return print('Simulated Bed Files Genereated: ' + self.outdir + "/annotations/" + str(self.sample) + "_simulated_centered.bed and " + str(self.sample) + "_simulated_window.bed")

###### For experimental genome if experimental_fimo == True                
    def reformat_expt_fasta(self):
        df = pd.read_csv(self.outdir + "/generated_sequences/" + self.sample + "_experimental_window_sequences.fa", header=None)
        df = df[df.index % 2 != 0]
        l_expt_seq = df.values.tolist()
        flat_list = [item for sublist in l_expt_seq for item in sublist]
        strN = str('N')*20
        file = ["{}{}".format(i,strN) for i in flat_list]

        sequence_num = len(file)
        seq_length = ((2*int(self.window))+21)
        tot_bases = int(sequence_num*seq_length)
        seq_per_chrom = math.floor(sequence_num/self.chrom_num)
        seq_last_chrom = seq_per_chrom + sequence_num%self.chrom_num
        #this is to check if there is the correct amount of sequences per chromosome
        scount = seq_last_chrom + (seq_per_chrom)*(self.chrom_num-1)

        base_per_chrom = seq_length*seq_per_chrom 
        base_last_chrom = seq_length*seq_last_chrom
        #this is to check if there is the correct amount of bases per chromosome
        bscount = base_last_chrom + (base_per_chrom)*(self.chrom_num-1)

        print("Individual Sequence Length: " + str(seq_length))
        print("Number of Experimental Sequences: " + str(sequence_num)) 
        print("Number of Chromosomes: " + str(self.chrom_num) +'\n')
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
            print('Number of bases per chromosome are properly structured.')

        dd = {}
        for i in range(0,self.chrom_num-1):
            dd['>exp_' + str(i + 1)] = str(file[(seq_per_chrom)*i:(seq_per_chrom)*(i+1)]).replace("', '","").replace("['","").replace("']","")
        dd['>exp_' + str(self.chrom_num)] = str(file[(seq_per_chrom)*(self.chrom_num-1):]).replace("', '","").replace("['","").replace("']","")
        df = pd.DataFrame.from_dict(dd, orient='index')
        df.columns = range(df.shape[1])
        df=df.reset_index()
        df = df.stack()
        df = pd.DataFrame(df).reset_index(drop=True)
        return df.to_csv(self.outdir + "/generated_sequences/" + str(self.sample) + "_experimental.fa", header=None, index=False, sep='\t')


    def index_and_chrm_sizes_expt(self):
            try:
                os.system("samtools faidx " + self.outdir + "/generated_sequences/" + self.sample + "_experimental.fa")
            except OSError:
                print ("Creation of index file %s failed" % self.outdir + "/generated_sequences/" + self.sample + "_experimental.fa.fai")
            else:
                print ("Successfully created index file %s" % self.outdir + "/generated_sequences/" + self.sample + "_experimental.fa.fai")

            try:
                os.system("cut -f1,2 " + self.outdir + "/generated_sequences/" + self.sample + "_experimental.fa.fai > " + self.outdir + "/generated_sequences/" + self.sample + "_experimental.chrom.sizes")
            except OSError:
                print ("Creation of chromosome size file %s failed" % self.outdir + "/generated_sequences/" + self.sample + "_experimental.chrom.sizes")
            else:
                print ("Successfully created chromosome size file %s" % self.outdir + "/generated_sequences/" + self.sample + "_experimental.chrom.sizes")

    def expt_bed(self):
        #getting the sequence number
        annotation_file = pd.read_csv(self.annotation, header=None)
        sequence_num = len(annotation_file)
        #calculating the number of sequences per chromosome
        seq_per_chrom = math.floor(sequence_num/self.chrom_num)
        seq_last_chrom = seq_per_chrom + sequence_num%self.chrom_num
        #checking fasta parameters
        seq_length = ((2*int(self.window))+21)
        #this is to check if there is the correct amount of sequences per chromosome
        scount = seq_last_chrom + (seq_per_chrom)*(self.chrom_num-1)
        
        print("----------Checking bed parameters------------")
        print("Individual Sequence Length: " + str(seq_length))
        print("Number of Experimental Sequences: " + str(sequence_num)) 
        print("Number of Chromosomes: " + str(self.chrom_num) +'\n')
        print("Sequences per Chromosome: " + str(seq_per_chrom))
        print("Sequences on Last Chromosome: " + str(seq_last_chrom))
        if sequence_num-scount != 0:
            print('Error! The total number of sequences does not equal the sum of the sequences divided into chromosomes!')
            sys.exit(1)
        else:
            print('Number of sequences per chromosome are properly structured.\n')


        #creating the locations for the all chromosomes except the last
        dbed = {}
        for i in range(0,(seq_per_chrom)):
            dbed['r_' + str(i)] = [int(self.window-1)+i*(2*self.window+21), int(self.window+1)+i*(2*self.window+21)]

        dfbed = pd.DataFrame.from_dict(dbed, orient='index')
        dfbed.columns = ['start', 'stop']
        #duplicating the dataframe for all chromosomes except the last
        dfbed = pd.concat([dfbed]*(self.chrom_num-1))
        #adding the chromosome name to the genomic location
        l_expt=[]
        for i in range(int(seq_per_chrom)):
            for c in range(int(self.chrom_num-1)):
                l_expt.append(str('exp_' + str(c+1)))
        l_expt = sorted(l_expt)
        dfbed['expt'] = l_expt

        #creating the locations for the last chromosome
        dlast = {}
        for i in range(0,(seq_last_chrom)):
            dlast['l_' + str(i)] = [int(self.window-1)+i*(2*self.window+21), int(self.window+1)+i*(2*self.window+21)]
        dflast = pd.DataFrame.from_dict(dlast, orient='index')
        dflast.columns = ['start', 'stop']
        #adding the chromosome name to the genomic location
        l_last = []
        for i in range(int(seq_last_chrom)):
            l_last.append(str('exp_' + str(self.chrom_num)))
        dflast['expt'] = l_last

        #combining the bed locations
        df = pd.concat([dfbed, dflast])
        df = df[['expt', 'start', 'stop']]
        #saving the new annotation
        df.to_csv(self.outdir + "/annotations/" + str(self.sample) + "_experimental_genome_centered.bed", header=None, index=False, sep='\t')

        #creating the self.windowed bed
        df["start_new"] = df.apply(lambda x: round((x["start"] + x["stop"])/2), axis=1)
        df["stop_new"] = df.apply(lambda x: x["start_new"] + 1, axis = 1)

        ##the -1500 position from "origin"
        df["start"] = df.apply(lambda x: x["start_new"] - int(self.window), axis=1)

        ##the 1500 position from the "origin"
        df["stop"] = df.apply(lambda x: x["stop_new"] + int(self.window), axis=1)

        ##saving the new annotation
        df.to_csv(self.outdir + '/annotations/' + self.sample + '_experimental_genome_window.bed', sep='\t', columns=["expt","start","stop"], header = False, index = False)    
        return print('Experimental Bed Files Genereated: ' + self.outdir + "/annotations/" + str(self.sample) + "_experimental_genome_centered.bed and " + str(self.sample) + "_experimental_genome_window.bed")


###### For whole genome if whole_genome_fimo == True
    def chrm_sizes_whole_genome(self):
        g = self.genome.split('/')[-1].split('.')[-2]
        try:
            os.system("cut -f1,2 " + self.genome + ".fai > " + self.outdir + "/generated_sequences/" + g + ".chrom.sizes")
        except OSError:
            print ("Creation of chromosome size file %s failed" % self.outdir + "/generated_sequences/" + g + ".chrom.sizes")
        else:
            print ("Successfully created chromosome size file %s" % self.outdir + "/generated_sequences/" + g + ".chrom.sizes")
