import argparse
from argparse import RawTextHelpFormatter
from multiprocessing import Pool
import main

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

    parser = argparse.ArgumentParser(description = 'This work flow functions to assess active transcription factors in a ATAC-seq or PRO-seq data set. First, it calculates base content per position over a window (given by a bed file). Subsequently, it generates sequences using a first order MM to obtain a simulated background sequence. Next, it runs FIMO to generate motif calls within the generated sequence (and the original sequence provided if needed). Next, it generates an MD score file for both the original sequence and the simulated sequence. Lastly, it outputs a text file with TF activation information along with various plots displaying this information.', formatter_class=RawTextHelpFormatter)
    
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    
    optional.add_argument('-v', '--verbose',dest="verbose", type=str2bool, nargs='?', const=True, default=False, help = 'will generate verbose output file', metavar="True/False", required=False)

#seq generation    Add in dinucleotide switch?
    required.add_argument('-g', '--genome', dest="genome", help = 'reference genome in fasta format. Genome index must be in the same directory (genome.fa.fai).', metavar="genome.fa")
    required.add_argument('-a', '--annotation', dest="annotation", help = 'input bed file of bidirectionals or ATAC peaks, ends with .bed or .sorted.bed', metavar="annotation.bed") ##Add in option to give multiple annotation files and use mumerge
    required.add_argument('-o', '--outdir', dest="outdir", help = 'directory for output', metavar="/full/path/to/output")
    required.add_argument('-s', '--sample', dest="sample", help = 'name of the sample to be run (str)', metavar="name_of_sample")

    optional.add_argument('-d', '--seed',dest="seed",type=int, default=True, help = 'seed for initializing the random number generator. To set a specific seed -d int. The dafault is "True" so the seed is based on the clock. If "0" a random seed will be selected.', metavar="int", required=False)
    optional.add_argument('-w', '--window', dest="window",type=int, default=1500, help = 'window to extract sequences (Default=1500)', metavar="int", required=False)
    optional.add_argument('-n', '--sequence_num', dest="sequence_num",type=int, default=20000, help = 'number of simulated sequences to be generated. For best results, set equal to number of peaks called in the annotation file provided', metavar="int", required=False) #change default to be equal to the number of experimental sequences (ie len annotation)
    optional.add_argument('-i', '--chrom_num', dest="chrom_num",type=int, default=10, help = 'number of simulated chromosomes. It is ideal to have approximately ~15 million bases per chromosome', metavar="int", required=False) #change default to hit this 15mill approximation
    
#fimo    
    required.add_argument('-m', '--motifs', dest="motifs", help = 'meme file for TFs', metavar="motif_database.meme")
    optional.add_argument('-t', '--threshold_fimo',dest="threshold_fimo",type=float, default=0.000001, help = 'threshold for motifs called by fimo. Default: 1e-6', metavar="float", required=False)
    optional.add_argument('-b', '--background_file', dest="background_file", help = 'background base composition of a given genome. This flag is HIGHLY recommended otherwise a uniform background distribution is assumed ie A/T/G/C = 25/25/25/25', metavar="background.csv", default=None, required=False)    
    optional.add_argument('-e', '--experimental_fimo', dest="experimental_fimo", type=str2bool, nargs='?', const=True, default=False, help = 'will run fimo over only the annotated regions from the experimental dataset. True will increase run time. Recommended if you are only looking at a single dataset. If False, provide destination of the pre-scanned genome using the "--pre_scan" flag or scan the whole genome. The -e, -x and -p flags are mutually exclusive. Default: False.', metavar="True/False", required=False)
    optional.add_argument('-x', '--whole_genome_fimo', dest="whole_genome_fimo", type=str2bool, nargs='?', const=True, default=False, help = 'will run fimo over the full genome provided by the "--genome" flag to generate the directory that can be called by the "--pre_scan" flag. Given the same fimo threshold, meme file, genome and background distribution this will not change. True will greatly increase run time. It is recommended to run this only once given the same fimo parameters as simulated and recall the data using the "--pre_scan" flag. The -e, -x and -p flags are mutually exclusive. Default: False.', metavar="True/False", required=False)
    optional.add_argument('-p', '--pre_scan', dest="pre_scan", default=None, help = 'directory containing pre-scanned motifs over the whole genome obtained using the same fimo parameters as the simulated dataset. The -e, -x and -p flags are mutually exclusive.', metavar="/full/path/to/pre-scanned/motifs", required=False)
    optional.add_argument('-c', '--cpus', type=int, dest='cpus', metavar='int', help='number of CPUs for multiprocessing. Default=1', default=1, required=False)

#scoring options    
    optional.add_argument('-k', '--dastk', dest="dastk", type=str2bool, nargs='?', const=True, default=False, help = 'will run dastk scoring as well as new scoring method. Default: False.', metavar="True/False", required=False)
        
#add histogram bin size change for md heat map parameter

    args = parser.parse_args()

#     ##I'm not really sure why this doesn't work........
#     if (args.experimental_fimo==True) and (args.whole_genome_fimo==True):
#         raise ValueError("The --experimental_fimo and --whole_genome_fimo arguments are mutually exclusive. Please provide only one of them.")    
    
#     if (args.experimental_fimo==True) and (args.pre_scan!=''):
#         raise ValueError("The --experimental_fimo and --pre_scan arguments are mutually exclusive. Please provide only one of them.") 
        
#     if (args.whole_genome_fimo==True) and (args.pre_scan!=None):
#         raise ValueError("The --whole_genome_fimo and --pre_scan arguments are mutually exclusive. Please provide only one of them.")
    
main.run(args.outdir, args.annotation, args.genome, args.sample, args.motifs, args.background_file, args.pre_scan, args.sequence_num, args.window, args.chrom_num, args.threshold_fimo, args.cpus, args.seed, args.experimental_fimo, args.whole_genome_fimo, args.dastk, args.verbose)