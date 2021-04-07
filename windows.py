# Author: Rutendo F. Sigauke   
import pandas as pd
import os

class BedWindows:

    # Initializer / Instance Attributes
    def __init__(self, annotation, outdir, sample, window = None):        
        self.annotation = annotation
        self.outdir = outdir
        self.sample = sample
        self.window = window

    def get_windows(self):
        '''This function takes in bed files from Tfit and  redefines mu and extends the window
        '''
        ##takes in Tfit regions/calls or any other bed file with chr, start, stop columns
        print("Annotation File: " + self.annotation)
        os.system("rsync " + self.annotation + ' ' + self.outdir + '/annotations/')
#         os.system("mv " self.outdir + '/annotations/' + self.annotation + " " + self.outdir + '/annotations/' + self.sample + ".bed")
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
