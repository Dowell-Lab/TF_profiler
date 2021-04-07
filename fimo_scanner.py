import os
import sys
import time
import subprocess
import multiprocessing
from fimo_functions import *


class FIMOScan:
    '''Setting up for parallel FIMO scan'''
    def __init__(self, motif_list, motifs, cutoff_fimo, background_file, genome, outdir, sample):
        self.motif_list = motif_list
        self.motifs = motifs
        self.cutoff_fimo = cutoff_fimo
        self.background_file = background_file
        self.genome = genome
        self.outdir = outdir
        self.sample = sample   

    def sim_scanner(self, mls):
        sim_fasta = self.outdir + "/generated_sequences/" + self.sample + "_simulated.fa"
        os.system('mkdir -p ' + self.outdir + '/temp/sim_fimo_out')
#         print("fimo --verbosity 1 --thresh " + str(self.cutoff_fimo) + " --bfile " + self.background_file + " --motif " + mls + " --oc " + self.outdir + '/temp/sim_fimo_out/' + mls + " " + self.motifs + " " + sim_fasta)
        if self.background_file is not None:
            os.system("fimo --verbosity 1 --thresh " + str(self.cutoff_fimo) + " --bfile " + self.background_file + " --motif " + mls + " --oc " + self.outdir + '/temp/sim_fimo_out/' + mls + " " + self.motifs + " " + sim_fasta)
        else:
            os.system("fimo --verbosity 1 --thresh " + str(self.cutoff_fimo) + " --motif " + mls + " --oc " + self.outdir + '/temp/sim_fimo_out/' + mls + " " + self.motifs + " " + sim_fasta)
    
    def expt_scanner(self, mls):
        expt_fasta = self.outdir + "/generated_sequences/" + self.sample + "_experimental.fa"
        os.system('mkdir -p ' + self.outdir + '/temp/expt_fimo_out')
#         print("fimo --verbosity 1 --thresh " + str(self.cutoff_fimo) + " --bfile " + self.background_file + " --motif " + mls + " --oc " + self.outdir + '/temp/expt_fimo_out/' + mls + " " + self.motifs + " " + expt_fasta)
        if self.background_file is not None:
            os.system("fimo --verbosity 1 --thresh " + str(self.cutoff_fimo) + " --bfile " + self.background_file + " --motif " + mls + " --oc " + self.outdir + '/temp/expt_fimo_out/' + mls + " " + self.motifs + " " + expt_fasta)
        else:
            os.system("fimo --verbosity 1 --thresh " + str(self.cutoff_fimo) + " --motif " + mls + " --oc " + self.outdir + '/temp/expt_fimo_out/' + mls + " " + self.motifs + " " + expt_fasta)

  
    def whole_genome_scanner(self, mls):
        os.system('mkdir -p ' + self.outdir + '/temp/whole_genome_fimo_out')
#         print("fimo --verbosity 1 --thresh " + str(self.cutoff_fimo) + " --bfile " + self.background_file + " --motif " + mls + " --oc " + self.outdir + '/temp/whole_genome_fimo_out/' + mls + " " + self.motifs + " " + self.genome)
        if self.background_file is not None:
            os.system("fimo --verbosity 1 --thresh " + str(self.cutoff_fimo) + " --bfile " + self.background_file + " --motif " + mls + " --oc " + self.outdir + '/temp/whole_genome_fimo_out/' + mls + " " + self.motifs + " " + self.genome)
        else:
            os.system("fimo --verbosity 1 --thresh " + str(self.cutoff_fimo) + " --motif " + mls + " --oc " + self.outdir + '/temp/whole_genome_fimo_out/' + mls + " " + self.motifs + " " + self.genome)
    
   
