import argparse
from argparse import RawTextHelpFormatter
from multiprocessing import Pool
import main
import sys

if __name__ == "__main__":
    
    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')    

    p = argparse.ArgumentParser(description = 'This work flow functions to assess active transcription factors in a ATAC-seq or PRO-seq data set. First, it calculates base content per position over a 3000bp window (given by a bed file). Subsequently, it generates sequences using a first or second order MM to obtain a simulated background sequence. Next, it can run FIMO to generate motif calls within the generated sequence (and the original sequence provided) or take a set of pre-scanned sequences. Next, it generates an MD score file for both the original sequence and the simulated sequence. Lastly, it outputs a text file with TF activation information along with various plots displaying this information.', formatter_class=RawTextHelpFormatter)    
    
    general = p.add_argument_group('General Arguments')
    sequence_gen = p.add_argument_group('Sequence Generation Arguments')
    fimo = p.add_argument_group('FIMO Arguments')
    
    #### General options
    general.add_argument('-v', '--verbose',dest="verbose", type=str2bool, nargs='?', const=True, default=False, help = 'will generate verbose output file', metavar="True/False", required=False)
 
    general.add_argument('-a', '--annotation', dest="annotation", help = 'input bed file of bidirectionals or ATAC peaks, ends with .bed or .sorted.bed', metavar="annotation.bed", required=True)
    general.add_argument('-o', '--outdir', dest="outdir", help = 'directory for output', metavar="/full/path/to/output", required=True)
    general.add_argument('-s', '--sample', dest="sample", help = 'name of the sample to be run (str)', metavar="name_of_sample", required=True)
    general.add_argument('-c', '--cpus', type=int, dest='cpus', metavar='int', help='number of CPUs for multiprocessing. Default=1', default=1, required=False)
        
    ### Sequence generation options
    sequence_gen.add_argument('-l', '--mononucleotide_generation',dest="mononucleotide_generation", type=str2bool, nargs='?', const=True, default=False, help = 'If False mononucleotide simulated sequences will not be generated. By default these sequences will also be scanned in FIMO unless the skip_simulated_fimo flag is used (Default: False)', metavar="True/False", required=False)
    sequence_gen.add_argument('-d', '--dinucleotide_generation',dest="dinucleotide_generation", type=str2bool, nargs='?', const=True, default=True, help = 'If False dinucleotide simulated sequences will not be generated. By default these sequences will also be scanned in FIMO unless the skip_simulated_fimo flag is used (Default: True)', metavar="True/False", required=False) 
    sequence_gen.add_argument('-g', '--genome', dest="genome", help = 'reference genome in fasta format. Genome index must be in the same directory (genome.fa.fai).', metavar="genome.fa", required=True)  
    
    sequence_gen.add_argument('-e', '--seed',dest="seed", default=None, help = 'seed for initializing the random number generator. To set a specific seed -e int. The default is "None" so the seed is based on the clock.', metavar="int", required=False)
    sequence_gen.add_argument('-w', '--window', dest="window",type=int, default=1500, help = 'window to extract sequences. Default=1500', metavar="int", required=False)
    sequence_gen.add_argument('-n', '--sequence_num', dest="sequence_num",type=int, help = 'number of simulated sequences to be generated. Default: set equal to number of peaks called in the annotation file provided', metavar="int", required=False)
    sequence_gen.add_argument('-i', '--chrom_num', dest="chrom_num",type=int, help = 'number of simulated chromosomes. Default: approximately ~11.5 million bases per chromosome (10 percent of the average size of a human chromosome).', metavar="int", required=False)

    
#     ###FIMO Scanning options
    fimo.add_argument('-p', '--pre_scan', dest="pre_scan", default=None, help = 'directory containing pre-scanned motif hits in bed format over the whole genome obtained using the same fimo parameters (background and threshold) as the simulated dataset. If path is set experimental genome scan will be skipped and pre-scanned motif hits will be used instead.', metavar="/full/path/to/pre-scanned/motifs", required=False)                   
    fimo.add_argument('-x', '--experimental_fimo', dest="experimental_fimo", type=str2bool, nargs='?', const=True, default=False, help = 'will run fimo over only the annotated regions from the experimental dataset. True will increase run time. Recommended if you are only looking at a single dataset. If False, provide destination of the pre-scanned genome using the "--pre_scan" flag. Default: False.', metavar="True/False", required=False)
    fimo.add_argument('-k', '--skip_simulated_fimo', dest="skip_simulated_fimo", type=str2bool, nargs='?', const=True, default=False, help = 'if True then the fimo scan for the simulated sequence(s) will be skipped', metavar="True/False", required=False)
    fimo.add_argument('-m', '--motifs', dest="motifs", help = 'meme file for TFs', metavar="motif_database.meme", required=True)
 
    
    fimo.add_argument('-t', '--threshold_fimo',dest="threshold_fimo",type=float, default='1e-5', help = 'threshold for motifs called by fimo. Default: 1e-5', metavar="float", required=False)
    fimo.add_argument('-b', '--background_file', dest="background_file", help = 'background base composition of a given genome. This flag is HIGHLY recommended. See background options in useful files.', metavar="background.csv", default=None, required=False)    
        ###Can I default in a useful files folder for a certain background?
        
#     ###Distance
#     p.add_argument('-q', '--distance', dest="distance", type=str2bool, nargs='?', const=True, default=False, help = 'will calculate distances between a given motif hit and mu. Default: False.', metavar="True/False", required=False)    
    
#     ###Scoring
#     p.add_argument('-j', '--score', dest="score", type=str2bool, nargs='?', const=True, default=False, help = 'will run negative exponential scoring. Default: False.', metavar="True/False", required=False)
#     p.add_argument('-k', '--dastk', dest="dastk", type=str2bool, nargs='?', const=True, default=False, help = 'will run dastk scoring as well as new scoring method. Default: False.', metavar="True/False", required=False)
#     #p.add_argument('-q', '--score_distances', dest="score_distances", help='score based on a pre-generated distance file')  

    args = p.parse_args()

main.run(args.verbose, args.outdir, args.sample, args.genome, args.annotation, 
        args.sequence_num, args.chrom_num, args.motifs, args.background_file, args.seed, args.cpus, args.window,
        args.mononucleotide_generation, args.dinucleotide_generation, args.skip_simulated_fimo,
        args.experimental_fimo, args.pre_scan,
        args.threshold_fimo)

