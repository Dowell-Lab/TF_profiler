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
def run_statistics_module(verbose, outdir, sample, traditional_md, pval_cutoff, outliers_fraction):
    md_score_df = import_md_scores(outdir=outdir,sample=sample, traditional_md=traditional_md)
    md_score_df, xx, yy, Z = define_unchanging_md_scores(verbose=verbose, md_score_df=md_score_df, 
                                                  outliers_fraction=outliers_fraction)
    md_score_df, slope, intercept, sigma_cutoff = determine_significance(verbose=verbose, outdir=outdir, 
                                                                         md_score_df=md_score_df, pval_cutoff=pval_cutoff)

    color_types=['significance','gc_content', 'elliptic_fit']
    for color_type in color_types:
        plot_rbg(verbose=verbose, outdir=outdir, sample=sample, 
                 color_type=color_type, md_score_df=md_score_df, 
                 slope=slope, intercept=intercept, 
                 sigma_cutoff=sigma_cutoff,
                 xx=xx, yy=yy, Z=Z)

######################################### MD Score Functions #########################################
def import_md_scores(outdir,sample,traditional_md):
    if traditional_md == True:
        exp_md_score_df = pd.read_csv(outdir+'/scores/experimental_traditional_md_score.txt', sep='\t')
        sim_md_score_df = pd.read_csv(outdir+'/scores/simulated_traditional_md_score.txt', sep='\t')
    #     sim_md_score_df = pd.read_csv(outdir+'/scores/mononucleotide_simulated_traditional_md_score.txt', sep='\t')
    
    md_score_df = exp_md_score_df.merge(sim_md_score_df, on='tf', suffixes=('_exp', '_sim'))
    
    gc = pd.read_csv(outdir+'/generated_sequences/'+sample+'_motif_gc_percentage.txt', sep='\t', 
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
    
    # Get ellipse shape
    xx, yy = np.meshgrid(np.linspace(0.05, 0.15),
                 np.linspace(0.05, 0.15))
    Z = algorithm.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    
    return md_score_df, xx, yy, Z

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
       
def plot_rbg(verbose, outdir, sample, color_type, md_score_df, slope, intercept, sigma_cutoff, xx, yy, Z):
    if verbose=='True':
        print('Plotting rbg. Coloring is based on '+ color_type +'.')
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
        plt.contour(xx, yy, Z, levels=[0], linewidths=2, colors='black', alpha = 0.5)
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
        ax.plot([0, max_x*1.1], [intercept, max_x*slope*1.1+intercept], 
                color='black', linestyle='--', alpha = 0.7)
        ax.plot([0, max_x*1.1], [-sigma_cutoff+intercept, -sigma_cutoff+intercept+max_x*slope*1.1], 
                color='black', linestyle='--', alpha = 0.3)
        ax.plot([0, max_x*1.1], [sigma_cutoff+intercept, sigma_cutoff+intercept+max_x*slope*1.1], 
                color='black', linestyle='--', alpha = 0.3)
        plt.legend(bbox_to_anchor=(1.05, 1),fontsize = 25, loc=2, borderaxespad=0.)

    if color_type == 'gc_content':
        gc=ax.scatter(md_score_df['md_score_sim'], 
                   md_score_df['md_score_exp'], 
                   c=md_score_df['percent_gc'], cmap='plasma')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        plt.colorbar(gc, cax=cax, orientation='vertical')
    
    plt.savefig(outdir + '/plots/rbg_'+ color_type + '.png', bbox_inches='tight')

