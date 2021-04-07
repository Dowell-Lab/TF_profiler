##From DAStk by Ignacio Tripodi
import multiprocessing
import numpy as np
import pandas as pd
from functools import partial
from operator import itemgetter


def is_in_window(motif_interval, bed_median, window_size):
    start = int(motif_interval.start)
    end = int(motif_interval.end)
    if (end >= (bed_median - window_size) and end <= (bed_median + window_size)) \
        or (start >= (bed_median - window_size) and start <= (bed_median + window_size)):
        return True
    else:
        return False
    
def find_motifs_in_chrom(current_chrom, files):
    tf_motif_filename, bed_name, window = files
    H = window          # in bps, the MD-score parameter (large window)
    h = window * 0.1           # in bps, the MD-score parameter (small window)

    wdf = pd.read_csv(bed_name, header=None, comment='#', sep="\t", usecols=[0, 1, 2], \
                          names=['chrom', 'start', 'end'], na_filter=False, dtype={'chrom':'str', 'start':'int', 'end':'int'})
    w_iter = wdf[(wdf.chrom == current_chrom)].itertuples()
    motif_df = pd.read_csv(tf_motif_filename, header=None, comment='#', sep="\t", usecols=[0, 1, 2], \
                           names=['chrom', 'start', 'end'], na_filter=False, dtype={'chrom':'str', 'start':'int', 'end':'int'})
    if len(motif_df) == 0:
        return None
    motif_iter = motif_df[(motif_df.chrom == current_chrom)].itertuples()
    last_motif = None
    g_h = 0
    g_H = 0
    total_motif_sites = 0
    tf_distances = []
    try:
        motif_region = next(motif_iter)   # get the first motif in the list
    except StopIteration:
        print('No motifs for chromosome ' + current_chrom + ' on file ' + tf_motif_filename)
        return None

    peak_count_overlapping_motif = 0
    for bed_peak in w_iter:
        motifs_within_region = True

        bed_median = bed_peak.start + (bed_peak.end - bed_peak.start)/2
        # check the last motif, too, in case any of them falls within the region of
        # interest of two sequential ATAC-Seq peaks
        if last_motif:
            if is_in_window(last_motif, bed_median, h):
                g_h = g_h + 1
            if is_in_window(last_motif, bed_median, H):
                g_H = g_H + 1
                tf_median = last_motif.start + (last_motif.end - last_motif.start)/2
                tf_distances.append(bed_median - tf_median)
                try:
                    motif_region = next(motif_iter)   # get the next motif in the list
                    total_motif_sites += 1
                except StopIteration:
                    pass
                last_motif = motif_region
                if motif_region.start > (bed_median + H):
                    motifs_within_region = False

        # Move to the next putative motif sites until we get one past our evaluation window
        while motifs_within_region:
            total_motif_sites += 1
            # account for those within the smaller window (h)
            if is_in_window(motif_region, bed_median, h):
                g_h = g_h + 1
            # account for those within the larger window (H)
            if is_in_window(motif_region, bed_median, H):
                g_H = g_H + 1
                tf_median = motif_region.start + (motif_region.end - motif_region.start)/2
                tf_distances.append(bed_median - tf_median)

            # if we still haven't shifted past this peak...
            if motif_region.start <= (bed_median + H):
                try:
                    motif_region = next(motif_iter)   # get the next motif in the list
                    total_motif_sites += 1
                except StopIteration:
                    # No more TF motifs for this chromosome
                    break
                last_motif = motif_region
            else:
                motifs_within_region = False

    # Count any remaining TF motif sites after the last peak
    while(len(motif_region) > 0):
        try:
            motif_region = next(motif_iter)   # this gets the next motif in the list
        except StopIteration:
            break

    return [tf_distances, g_h, g_H, total_motif_sites]

def get_md_score(tf_motif_filename, threads, bed_name, CHROMOSOMES, window):

    HISTOGRAM_BINS = 150
    pool = multiprocessing.Pool(threads)
    results = pool.map(partial(find_motifs_in_chrom, files=[tf_motif_filename, bed_name, window]), CHROMOSOMES)
    pool.close()
    pool.join()

    results_matching_motif = [x for x in results if x is not None]
    if len(results_matching_motif) > 0:
        sums = np.sum(results_matching_motif, axis=0)
        overall_g_h = sums[1]
        overall_g_H = sums[2]
        overall_motif_sites = sums[3]

        # Calculate the heatmap for this motif's barcode
        tf_distances = sums[0]
        heatmap, xedges = np.histogram(tf_distances, bins=HISTOGRAM_BINS)
        str_heatmap = np.char.mod('%d', heatmap)
        if overall_g_H >= 0:
            return [float(overall_g_h + 1)/(overall_g_H + 1), (overall_g_h + 1), (overall_g_H + 1), (overall_motif_sites + 1), ';'.join(str_heatmap)]
    else:
        return None