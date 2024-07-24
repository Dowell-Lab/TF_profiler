######################################### Imports #########################################
import os
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from sklearn.linear_model import LinearRegression
from sklearn.covariance import EllipticEnvelope
from scipy.stats import norm

######################################### MD Score Main #########################################
def run_statistics_module(verbose, outdir, sample, window, motifs, traditional_md, pval_cutoff, outliers_fraction, plot_barcode, seq_type):
    md_score_df = import_md_scores(outdir=outdir,sample=sample, motifs=motifs, traditional_md=traditional_md, seq_type=seq_type)
    md_score_df = define_unchanging_md_scores(verbose=verbose, md_score_df=md_score_df, 
                                                  outliers_fraction=outliers_fraction)
    md_score_df, slope, intercept, sigma_cutoff = determine_significance(verbose=verbose, outdir=outdir, 
                                                                         md_score_df=md_score_df, pval_cutoff=pval_cutoff)

    color_types=['significance','gc_content', 'elliptic_fit']
    for color_type in color_types:
        plot_rbg(verbose=verbose, outdir=outdir, sample=sample, 
                 color_type=color_type, md_score_df=md_score_df, 
                 slope=slope, intercept=intercept, 
                 sigma_cutoff=sigma_cutoff)
    if plot_barcode != 'none':
        plot_barcodes(verbose=verbose, outdir=outdir, sample=sample,
                      window=window, pval_cutoff=pval_cutoff, md_score_df=md_score_df, 
                      plot_barcode=plot_barcode)

######################################### MD Score Functions #########################################
def import_md_scores(outdir,sample,motifs,traditional_md,seq_type):
    
    exp_md_score_df = pd.read_csv(outdir+'/scores/experimental_traditional_md_score.txt', sep='\t')

    if seq_type=='simulated':
        sim_md_score_df = pd.read_csv(outdir+'/scores/simulated_traditional_md_score.txt', sep='\t')
    elif seq_type=='mononucleotide_simulated':
        sim_md_score_df = pd.read_csv(outdir+'/scores/mononucleotide_simulated_traditional_md_score.txt', sep='\t')
    
    md_score_df = exp_md_score_df.merge(sim_md_score_df, on='tf', suffixes=('_exp', '_sim'))

    full_path = os.path.abspath(__file__)
    parent_path = os.path.dirname(full_path)    
    meme_name=str(motifs)
    meme_name = meme_name.split('/')[-1]
    meme_name = meme_name.replace('.meme', '')
    gc = pd.read_csv(parent_path + '/assets/'+ meme_name +'_motif_gc_percentage.txt', sep='\t', 
                    names=['tf','percent_gc','length'])
    
    md_score_df=md_score_df.merge(gc, on = 'tf')
    md_score_df=md_score_df[['tf','percent_gc','md_score_exp','md_score_sim']]
    return md_score_df


def define_unchanging_md_scores(verbose, md_score_df, outliers_fraction):
    X = np.array(md_score_df[[('md_score_sim'), ('md_score_exp')]])
    n_samples = len(X)
    n_outliers = int(outliers_fraction * n_samples)
    n_inliers = n_samples - n_outliers
    if verbose == True:
        print('At ' + str(outliers_fraction*100) + ' percent outliers there are '+ str(n_outliers) + ' total.')
        print(str(n_inliers) + ' MD-scores remain for subsequent linear regression analysis.')

    # Anomaly detection method
    algorithm = EllipticEnvelope(contamination=outliers_fraction)
    
    # Fit
    algorithm.fit(X)
    y_pred_elliptic = algorithm.fit(X).predict(X)
    md_score_df['elliptic_outlier'] = y_pred_elliptic    
    return md_score_df

def determine_significance(verbose, outdir, md_score_df, pval_cutoff):
    fit_unchanging_md_score = md_score_df[md_score_df['elliptic_outlier'] == 1]
    X = fit_unchanging_md_score['md_score_sim'].values.reshape(-1, 1)  # values converts it into a numpy array
    Y = fit_unchanging_md_score['md_score_exp'].values.reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column
    linear_regressor = LinearRegression(fit_intercept=True)  # create object for the class
    linear_regressor.fit(X, Y)  # perform linear regression
    y_pred = linear_regressor.predict(X)  # make predictions
    slope=float(linear_regressor.coef_)
    intercept=float(linear_regressor.intercept_)

    if verbose==True:
        print('The slope of the linear regression is ' + str(slope) +'.')
        print('The intercept is '+ str(intercept) + '.')
    
    fit_unchanging_md_score = md_score_df[md_score_df['elliptic_outlier'] == 1]
    residuals = (fit_unchanging_md_score['md_score_exp'] - fit_unchanging_md_score['md_score_sim']*slope)

    (mu, sigma) = norm.fit(residuals)

    if verbose == True:
        print('Normal distribution information for the residuals:')
        print('mu = ' + str(mu))
        print('sigma = '+ str(sigma))
        
    z_cutoff = norm.ppf(pval_cutoff)
    sigma_cutoff = z_cutoff*sigma+mu

    all_residuals= (md_score_df['md_score_exp'] - md_score_df['md_score_sim']*slope)
    z_scores = (all_residuals - mu)/sigma
    p_values = norm.sf(abs(z_scores))
    md_score_df['pval'] = p_values
    md_score_df=md_score_df.sort_values(by='pval')
    
    md_score_df['significance'] = md_score_df.apply(lambda row: significance_caller(row, pval_cutoff), axis=1)
    md_score_df.to_csv(outdir+'/scores/md_score_experimental_vs_simulated_significance.txt',sep='\t',index=False)
    return md_score_df, slope, intercept, sigma_cutoff

def significance_caller(row, pval_cutoff):
    if (row['pval'] <= pval_cutoff and (row['md_score_exp'] - row['md_score_sim'] > 0)):
        return 'Enriched'
    if (row['pval'] <= pval_cutoff and (row['md_score_exp'] - row['md_score_sim'] < 0)):
        return 'Depleted'
    else:
        return 'Not Significant'
       
def plot_rbg(verbose, outdir, sample, color_type, md_score_df, slope, intercept, sigma_cutoff):
    os.system('mkdir -p ' + outdir + '/plots')
    if verbose=='True':
        print('Plotting rbg. Coloring is based on '+ color_type +'.')
    plt.clf()
    plt.figure(figsize=(12,10))
    gs = plt.GridSpec(1, 1)
    ax = plt.subplot(gs[0])
    max_x=max(md_score_df['md_score_sim'])
    max_y=max(md_score_df['md_score_exp'])

    plt.suptitle(sample+'_rbg_'+color_type, fontsize=40, fontweight='bold')
    plt.xticks(fontsize = 25)
    plt.yticks(fontsize = 25)
    plt.xlim([0, max_x*1.05])
    plt.ylim([0, max_y*1.05])
    ax.set_xlabel('Nucleotide Background MD-Score',fontsize=30,fontweight='bold')
    ax.set_ylabel('Observed MD-Score',fontsize=30,fontweight='bold')
    
    if color_type == 'elliptic_fit':
        md_score_df_group_outlier = md_score_df.groupby('elliptic_outlier')
        for outlier,md_score_df_outlier in md_score_df_group_outlier:
            if outlier == -1:
                color='grey'
            elif outlier == 1:
                color='blue'
            ax.scatter(md_score_df_outlier['md_score_sim'], 
               md_score_df_outlier['md_score_exp'], color=color, alpha=0.7)
        ax.plot([0, max_x*1.1], [intercept, max_x*slope*1.1+intercept], 
                color='black', linestyle='--', alpha = 0.7)
    
    if color_type=='significance':
        md_score_df_group_significance = md_score_df.groupby('significance')
        for significance,md_score_df_significance in md_score_df_group_significance:
            if significance == 'Not Significant':
                color='blue'
            elif significance =='Enriched':
                color='green'
            elif significance == 'Depleted':
                color='red'
            else:
                color='black'
            ax.scatter(md_score_df_significance['md_score_sim'], 
                       md_score_df_significance['md_score_exp'], 
                       label=significance, color=color, alpha=0.7)
        #here I am plotting the trend line and trend line +/-pval_cutoff, max_x*slope=expected_md_score_exp 
        #*1.1 is simply to make the line trail off the plot slightly for the ~aesthetic~
#         ax.plot([0, max_x], [intercept, max_x*slope+intercept], 
#                 color='black', linestyle='--', alpha = 0.7)
#         ax.plot([0, max_x], [-sigma_cutoff+intercept, -sigma_cutoff+intercept+max_x*slope], 
#                 color='black', linestyle='--', alpha = 0.3)
#         ax.plot([0, max_x], [sigma_cutoff+intercept, sigma_cutoff+intercept+max_x*slope], 
#                 color='black', linestyle='--', alpha = 0.3)
        plt.legend(bbox_to_anchor=(1.05, 1),fontsize = 25, loc=2, borderaxespad=0.)

    if color_type == 'gc_content':
        gc=ax.scatter(md_score_df['md_score_sim'], 
                   md_score_df['md_score_exp'], 
                   c=md_score_df['percent_gc'], cmap='plasma')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        plt.colorbar(gc, cax=cax, orientation='vertical')
    
    plt.savefig(outdir + '/plots/rbg_'+ color_type + '.png', bbox_inches='tight')

def plot_barcodes(verbose, outdir, sample, window, pval_cutoff, md_score_df, plot_barcode):
    try:
        os.system('mkdir -p ' + outdir + '/plots/barcodes')
    except OSError:
        print('Creation of the directory %s failed' % outdir + '/plots/barcodes')
        sys.exit(1)
    else:
        if verbose == True:
            print(outdir + '/plots/barcodes exists.') 
    
    bins = int(window*0.1)
    
    ### plot barcode flag determines which barcodes are plotted
    if plot_barcode == 'significance':
        significant_md_score_df = md_score_df[md_score_df['pval'] <= pval_cutoff]
        tf_list=list(significant_md_score_df['tf'])
    elif plot_barcode == 'all':
        tf_list=list(md_score_df['tf'])
    else:
        all_tf_list = list(md_score_df['tf'])
        tf_list=[]
        tf_upper_case = str.upper(plot_barcode)
        tf_list.extend([tf for tf in all_tf_list if tf_upper_case in tf])

    if len(tf_list) == 0:
        if verbose == True:
            print('No TFs with an MD-score in both simulated and experimental sets match the user defined TF.')
    else:
        if verbose == True:
            print('Plotting barcodes for ' + str(len(tf_list)) + ' TFs.')
        for tf in tf_list:
            ###setting up experimental data for a given tf
            exp_distance_df = pd.read_csv(outdir+'/distances/experimental/'+tf+'_distances.txt', sep='\t')
            exp_heat, xedges = np.histogram(exp_distance_df['distance'],bins=bins)
            exp_n_distance=len(exp_distance_df)
            max_exp=max(exp_heat)
            exp_md_score=round(float(md_score_df['md_score_exp'][md_score_df['tf'] == tf]), 4)

            ###setting up simulated data for a given tf
            sim_distance_df = pd.read_csv(outdir+'/distances/simulated/'+tf+'_distances.txt', sep='\t')
            sim_heat, xedges = np.histogram(sim_distance_df['distance'],bins=bins)
            sim_n_distance=len(sim_distance_df)
            max_sim=max(sim_heat)
            sim_md_score=round(float(md_score_df['md_score_sim'][md_score_df['tf'] == tf]), 4)

            pval=round(float(md_score_df['pval'][md_score_df['tf'] == tf]), 4)

            plt.clf()
            plt.figure(figsize=(12,10))

            ###plotting experimental data
            ax = plt.subplot(1, 2, 1)
            heat_m = np.nan * np.empty(shape=(int(bins/4), bins))
            for row in range(int(bins/4)):
                heat_m[row] = exp_heat
                ax.matshow(heat_m, cmap='YlOrRd', vmin=0, vmax=max_exp)
            ax.axis('off')
            ax.text(bins/2, 43, 'Experimental MD-score = ' + str(exp_md_score), 
                     ha='center', size=14, zorder=0)
            ax.text(bins/2, 50, 'N = ' + str(exp_n_distance), 
                     ha='center', size=14, zorder=0)

            ###plotting simulated data
            ax1 = plt.subplot(1, 2, 2)
            heat_m = np.nan * np.empty(shape=(int(bins/4), bins))
            for row in range(int(bins/4)):
                heat_m[row] = sim_heat
                ax1.matshow(heat_m, cmap='YlOrRd', vmin=0, vmax=max_sim)
            ax1.axis('off')
            ax1.text(bins/2, 43, 'Simulated MD-score = ' + str(sim_md_score), 
                     ha='center', size=14, zorder=0)
            ax1.text(bins/2, 50, 'N = ' + str(sim_n_distance), 
                     ha='center', size=14, zorder=0)

            ###plot title including the significance value
            plt.suptitle(tf + '; p-value = ' + str(pval), 
                     y = 0.59, x = 0.5, size=20, zorder=0)
            plt.tight_layout()

            ###plot export
            plt.savefig(outdir+'/plots/barcodes/'+tf+'.png', bbox_inches='tight')
