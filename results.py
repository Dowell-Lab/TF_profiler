######################################### Imports #########################################
import pandas as pd

######################################### Results Main #########################################
def run_results(verbose, outdir, sample, window):
    if verbose == True:
        print('Combining MD Scores.')
    combine_md_scores(outdir=outdir, sample=sample) 
#     if verbose == True:
#         print('Some quick plotting.')
#     quick_plots(outdir=outdir, sample=sample)
    
######################################### Results Functions #########################################
def combine_md_scores(outdir, sample):
    dfed = pd.read_csv(outdir + '/results/' + sample + '_experimental_dastk_md_scores.txt', sep='\t')
    dfsd = pd.read_csv(outdir + '/results/' + sample + '_simulated_dastk_md_scores.txt', sep='\t')
    dfe = pd.read_csv(outdir + '/results/' + sample + '_experimental_md_scores.txt', sep='\t')
    dfs = pd.read_csv(outdir + '/results/' + sample + '_simulated_md_scores.txt', sep='\t')

    dfdas = dfed.merge(dfsd, on='motif_id', suffixes = ('_dastk_experimental', '_dastk_simulated'))
    dfnew = dfe.merge(dfs, on='motif_id', suffixes = ('_experimental', '_simulated'))
    df = dfnew.merge(dfdas, on='motif_id')
    gc = pd.read_csv(outdir + '/generated_sequences/' + sample + '_motif_gc_percentage.txt', 
                 sep='\t', names=['motif_id', 'gc_percentage'])
    df = df.merge(gc, on='motif_id')
    df.to_csv(outdir + '/results/' + sample + '_md_scores.txt', sep='\t', index=False)
    
# def quick_plots(outdir, sample):
#     df= pd.read_csv(outdir + '/results/' + sample + '_md_scores.txt', sep='\t')
    
#     #scores
#     fig = px.scatter(df, x="total_distance_score_simulated", y="total_distance_score_experimental", 
#                      template="simple_white", color='gc_percentage')
#     fig.write_html(outdir + '/results/' + sample + '_total_distance_score.html')

#     fig = px.scatter(df, x="md_score_dastk_simulated", y="md_score_dastk_experimental", 
#                      template="simple_white", color='gc_percentage')
#     fig.write_html(outdir + '/results/' + sample + '_md_score_dastk.html')
    
#     #hit occurances
#     fig = px.scatter(df, x="motif_hit_occurrences_simulated", y="motif_hit_occurrences_experimental", 
#                      template="simple_white", color='gc_percentage')
#     fig.write_html(outdir + '/results/' + sample + '_hit_occurances.html')
    
#     fig = px.scatter(df, x="motif_hit_occurrences_simulated", y="large_window_dastk_simulated", 
#                      template="simple_white", color='gc_percentage')
#     fig.write_html(outdir + '/results/' + sample + '_simulated_hits.html')

#     fig = px.scatter(df, x="motif_hit_occurrences_experimental", y="large_window_dastk_experimental", 
#                      template="simple_white", color='gc_percentage')
#     fig.write_html(outdir + '/results/' + sample + '_experimental_hits.html')
    
    
    
