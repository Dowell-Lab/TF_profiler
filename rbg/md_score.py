def fimo_parse(fimo_file=None, largewindow=None, retain='distance', 
    '''Parses a fimo output file and writes into a new file that is formatted
        in a way that can be parsed within existing TFEA functions
    Parameters
    ----------
    largewindow : float
        the size of the larger window to perform TFEA. Specified by user in
        config file
    tempdir : string
        full path to temp directory in output directory (created by TFEA)
    fimo_file : string
        full path to the fimo output file to be parsed by this function
    motif_file : string
        the name of the motif being parsed, this function will create a file
        using this motif_file string
        
    Returns
    -------
    outname : string
        the full path to the file to be used by other TFEA functions
    '''
    d = dict()
    with open(fimo_file) as F:
        header = F.readline().strip('\n').split('\t')
        if len(header) > 1:
            start_index = header.index('start')
            stop_index = header.index('stop')
            name_index = header.index('sequence_name')
            score_index = header.index('score')
            for line in F:
                line = line.strip('\n').split('\t')
                rank = line[name_index].split(',')[-1]
                rank = int(rank)
                start = line[start_index]
                stop = line[stop_index]
                distance = ((int(start)+int(stop))/2)-int(largewindow)
                score = line[score_index]
                if rank not in d:
                    d[rank] = [rank, score, distance]
                elif retain == 'score':
                    prev_score = float(d[rank][-2])
                    if prev_score < float(score):
                        d[rank] = [rank, score, distance]
                elif retain == 'distance':
                    prev_distance = float(d[rank][-1])
                    if prev_distance < float(distance):
                        d[rank] = [rank, score, distance]
    distances = list()
    for rank in range(1, linecount+1):
        if rank in d:
            distances.append(d[rank][-1])
        else:
            distances.append('.')

    return distances