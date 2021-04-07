# Author of get_gc: Jonathan Rubin  
def get_gc(motif=None, motif_database=None, alphabet=['A','C','G','T']):
    '''
    Obtain a pssm model from a meme formatted database file. Warning: If there are multiple
    motif matches, this function will return the last match in the database.

    '''
    motif_hit = False
    PSSM = list()
    with open(motif_database,'r') as F:
        for line in F:
            if 'MOTIF' in line:
                if motif in line:
                    motif_hit = True
                else:
                    motif_hit = False
            elif motif_hit and 'URL' not in line and 'letter-probability' not in line and line != '\n':
                acgt_probabilities = [float(x) for x in line.strip('\n').split()]
                total_prob = sum(acgt_probabilities)
                acgt_probabilities = [x/total_prob for x in acgt_probabilities] #Convert to probabilities
                PSSM.append(acgt_probabilities)

    gc = 0
    for base in PSSM:
        gc += base[alphabet.index('C')]
        gc += base[alphabet.index('G')]
    
    gc = gc/float(len(PSSM))
    
    return gc