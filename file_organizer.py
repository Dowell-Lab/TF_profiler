import os

def mkdirectories(outdir, experimental_fimo, whole_genome_fimo):
    try:
        os.system("mkdir -p " + outdir + "/temp")
    except OSError:
        print ("Creation of the directory %s failed" % outdir + "/temp")
    else:
        print ("Successfully created the directory %s" % outdir + "/temp")

    try:
        os.system("mkdir -p " + outdir + "/temp/md_scores")
    except OSError:
        print ("Creation of the directory %s failed" % outdir + "/temp/md_scores")
    else:
        print ("Successfully created the directory %s" % outdir + "/temp/md_scores")

    try:
        os.system("mkdir -p " + outdir + "/annotations")
    except OSError:
        print ("Creation of the directory %s failed" % outdir + "/annotations")
    else:
        print ("Successfully created the directory %s" % outdir + "/annotations")

    try:
        os.system("mkdir -p " + outdir + "/final_output")
    except OSError:
        print ("Creation of the directory %s failed" % outdir + "/final_output")
    else:
        print ("Successfully created the directory %s" % outdir + "/final_output")

    try:
        os.system("mkdir -p " + outdir + "/generated_sequences")
    except OSError:
        print ("Creation of the directory %s failed" % outdir + "/generated_sequences")
    else:
        print ("Successfully created the directory %s" % outdir + "/generated_sequences")

    try:
        os.system("mkdir -p " + outdir + "/motifs")
    except OSError:
        print ("Creation of the directory %s failed" % outdir + "/motifs")
    else:
        print ("Successfully created the directory %s" % outdir + "/motifs")

    try:
        os.system("mkdir -p " + outdir + "/motifs/simulated")
    except OSError:
        print ("Creation of the directory %s failed" % outdir + "/motifs/simulated")
    else:
        print ("Successfully created the directory %s" % outdir + "/motifs/simulated")

    if experimental_fimo == True:
        try:
            os.system("mkdir -p " + outdir + "/motifs/experimental")
        except OSError:
            print ("Creation of the directory %s failed" % outdir + "/motifs/experimental")
        else:
            print ("Successfully created the directory %s" % outdir + "/motifs/experimental")
    elif whole_genome_fimo == True:
        try:
            os.system("mkdir -p " + outdir + "/motifs/whole_genome")
        except OSError:
            print ("Creation of the directory %s failed" % outdir + "/motifs/whole_genome")
        else:
            print ("Successfully created the directory %s" % outdir + "/motifs/whole_genome")