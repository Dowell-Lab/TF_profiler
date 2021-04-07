# Author: Rutendo F. Sigauke 
import os
class ExtractSequences:
    '''Extract sequences from fasta based on bed coordinates using bedtools'''
    # Initializer / Instance Attributes
    def __init__(self, genome, sample, outdir):
        self.genome = genome
        self.sample = sample
        self.outdir = outdir

    # instance method
    def get_sequences(self):
        print('Genome Used: ' + self.genome + '\nOutdir: ' + self.outdir + '\nSample: ' + self.sample)
        os.system("bedtools getfasta -fi " + self.genome + " -bed " + self.outdir + '/annotations/' + self.sample + "_experimental_window.bed" + " -fo " + self.outdir + "/generated_sequences/" + self.sample + "_experimental_window_sequences.fa")
        print('Extracted Experimental Sequences Output: ' + self.outdir + "/generated_sequences/" + self.sample + "_experimental_window_sequences.fa")