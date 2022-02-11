######################################### Imports #########################################
import os
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from sklearn.linear_model import LinearRegression
from scipy.stats import norm

######################################### MD Score Main #########################################
def run_statistics_module(verbose, outdir, sample, traditional_md):
    if traditional_md == True:
        md_score_df=import_md_scores(outdir=outdir,sample=sample)
        md_score_df, slope, sigma = get_linear_regression_and_normal_distribution(verbose=verbose, outdir=outdir, md_score_df=md_score_df)
        color_types=['significance','gc_content']
        for color_type in color_types:
            plot_rbg(verbose=verbose, outdir=outdir, sample=sample, 
                     color_type=color_type, md_score_df=md_score_df, slope=slope, sigma=sigma)

######################################### MD Score Functions #########################################

def import_md_scores(outdir,sample):
    exp_md_score_df = pd.read_csv(outdir+'/scores/experimental_traditional_md_score.txt', sep='\t')
    sim_md_score_df = pd.read_csv(outdir+'/scores/simulated_traditional_md_score.txt', sep='\t')
#     sim_md_score_df = pd.read_csv(outdir+'/scores/mononucleotide_simulated_traditional_md_score.txt', sep='\t')
    md_score_df = exp_md_score_df.merge(sim_md_score_df, on='tf', suffixes=('_exp', '_sim'))
    
    gc = pd.read_csv(outdir+'/generated_sequences/'+sample+'_motif_gc_percentage.txt', sep='\t', 
                    names=['tf','percent_gc','length'])
    md_score_df=md_score_df.merge(gc, on = 'tf')
    return md_score_df


def get_linear_regression_and_normal_distribution(verbose, outdir, md_score_df):
    X = md_score_df['md_score_sim'].values.reshape(-1, 1)  # values converts it into a numpy array
    Y = md_score_df['md_score_exp'].values.reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column
    linear_regressor = LinearRegression(fit_intercept=0)  # create object for the class
    linear_regressor.fit(X, Y)  # perform linear regression
    Y_pred = linear_regressor.predict(X)  # make predictions
    slope=float(linear_regressor.coef_)
    
    if verbose==True:
        print('The slope of the linear regression is ' + str(slope) +'.')
        print('The intercept is set to '+ str(linear_regressor.intercept_) + '.')
    
    sim_md = list(md_score_df['md_score_sim'])
    (mu, sigma) = norm.fit(sim_md)
    
    if verbose == True:
        print('The normal distribution for the simulated data is as follows:')
        print('Mu = ' + str(mu))
        print('Sigma = '+ str(sigma))
    md_score_df['expected_md_score_exp'] = md_score_df['md_score_sim']*slope
    md_score_df['significance'] = md_score_df.apply(lambda row: significance_caller(row, sigma), axis=1)
    md_score_df=md_score_df.sort_values(by='significance', ascending=False)
    md_score_df.to_csv(outdir+'/scores/traditional_md_score_experimental_vs_simulated_significance.txt',sep='\t',index=False)
    return md_score_df, slope, sigma

def significance_caller(row, sigma):
    if (row['md_score_exp'] >= row['expected_md_score_exp']+sigma*2):
        return 'Enriched'
    if (row['md_score_exp']<= row['expected_md_score_exp']-sigma*2):
        return 'Depleted'
    else:
        return ''
       
def plot_rbg(verbose, outdir, sample, color_type, md_score_df, slope, sigma):
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

    if color_type=='significance':
        md_score_df_group_significance = md_score_df.groupby("significance")
        for significance,md_score_df_significance in md_score_df_group_significance:
            if significance == '':
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
        #here I am plotting the trend line and trend line +/-2sigma, max_x*slope=expected_md_score_exp 
        #*1.1 is simply to make the line trail off the plot slightly for the ~aesthetic~
        ax.plot([0, max_x*1.1], [0, max_x*slope*1.1], 
                color='black', linestyle='--', alpha = 0.7)
        ax.plot([0, max_x*1.1], [-2*sigma, -2*sigma+max_x*slope*1.1], 
                color='black', linestyle='--', alpha = 0.3)
        ax.plot([0, max_x*1.1], [2*sigma, 2*sigma+max_x*slope*1.1], 
                color='black', linestyle='--', alpha = 0.3)
        plt.legend(bbox_to_anchor=(1.05, 1),fontsize = 25, loc=2, borderaxespad=0.)

    if color_type == 'gc_content':
        gc=ax.scatter(md_score_df['md_score_sim'], 
                   md_score_df['md_score_exp'], 
                   c=md_score_df['percent_gc'], cmap='plasma')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        plt.colorbar(gc, cax=cax, orientation='vertical')
    
    plt.savefig(outdir + '/plots/' + sample + '_rbg_'+ color_type + '.png', bbox_inches='tight')


