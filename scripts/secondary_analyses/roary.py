#!/usr/bin/env python3

# Plotting imports
# %matplotlib inline

import sys
import logging
import matplotlib.pyplot as plt
import seaborn as sns

# Other imports
# import os
import pandas as pd
import numpy as np
import prince


working_dir = sys.argv[1]

try:
    imported_gene_matrix = pd.read_csv(
        f"{working_dir}/secondary_analyses/roary/roary_output/gene_presence_absence.csv", low_memory=False
    )
except BaseException:
    logging.exception("Roary file processing failed")
    sys.exit()

# imported_gene_matrix.info(verbose=True)
roary_igm = imported_gene_matrix.copy()
# Set index (group name)
roary_igm.set_index('Gene', inplace=True)
# Drop the other info columns
roary_igm.drop(list(roary_igm.columns[:13]), axis=1, inplace=True)
# Transform it in a presence/absence matrix (1/0)
roary_igm.replace('.{2,100}', 1, regex=True, inplace=True)
roary_igm.replace(np.nan, 0, regex=True, inplace=True)
roary_igm = roary_igm.drop('Inference',axis=1)
roary_igm.rename(columns=lambda x: x.split('_')[0], inplace=True)
roary_igm

################### seaborn plotting global parameters
sns.set_style('white',{'axes.linewidth': 20.0})
sns.set_context("paper", 
                font_scale = 1.8
               )
# plt.rcParams['xtick.major.size'] = 20
# plt.rcParams['xtick.major.width'] = 4
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True

################### making boxplot
def box_plot(data,x_label,y_label):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    # box_fig = plt.figure(figsize=(7, 10))
    sns.set_theme(style="white", rc=custom_params)
    sns.set_context("paper", 
                font_scale = 1.8
               )
    ax = sns.boxplot(y=data)
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True
    # plt.ylabel(y_label,fontsize=20)
    # plt.xlabel(x_label,fontsize=20)
    plt.savefig(f"{working_dir}/secondary_analyses/roary/Pseudogene_number_boxplot.png", bbox_inches="tight", dpi=100)

box_plot(roary_igm.sum(axis=0),'Isolates','Number of pseudogenes')

##############################
collect_vs_pseuds = roary_igm.sum().to_frame()
collect_vs_pseuds['Genomes'] = collect_vs_pseuds.index
collect_vs_pseuds.rename(columns={0:'Number of pseudogenes'}, inplace=True)

################### plotting histogram for pseudogenes
import numpy as np
from numpy import mean, size, zeros, where, transpose
from numpy.random import normal
from matplotlib.pyplot import hist
# from scipy import linspace
import array
    
from matplotlib.pyplot import figure,  plot, xlabel, ylabel,hist

def calculate_FD(pd_column):
    x = pd_column.to_list()
    # print(len(x))

    # my metric/tweaks to help estimate the range of bins that should be looked at using Freedman-Diaconis’s Rule of 
    # estimating bins to divide the amount of data giving the final metric for defining the range
    q3, q1 = np.percentile(x, [75 ,25])
    IQR = q3 - q1
    FD = (max(x) - min(x))/(2 * IQR / (len(x)**(1/3)))
    # print(FD)

    normalized_len_data = int(len(x)/FD)
    # print(normalized_len_data)

    return normalized_len_data, FD

def remove_outlier(df_in, col_name):
    q1 = df_in[col_name].quantile(0.25)
    q3 = df_in[col_name].quantile(0.75)
    iqr = q3-q1 #Interquartile range
    # print(f"IQR is {iqr}")
    fence_low  = q1-1.5*iqr
    fence_high = q3+1.5*iqr
    # print(f"upper boundary = {fence_high}", f"lower boundary = {fence_low}")
    df_out = df_in.loc[(df_in[col_name] > fence_low) & (df_in[col_name] < fence_high)]
    return df_out

def optimal_bin_width(pd_column, normalized_len_data, FD):
    plt.rcParams.update({'font.size': 15})

    x = pd_column.to_list()
    # print(len(x))

    x_max = max(x)
    x_min = min(x)
    N_MIN = int(FD-normalized_len_data)  #Minimum number of bins (integer)
            #N_MIN must be more than 1 (N_MIN > 1).
    N_MAX = int(FD+normalized_len_data)  #Maximum number of bins (integer)
    N = range(N_MIN,N_MAX) # #of Bins
    N = np.array(N)
    D = (x_max-x_min)/N    #Bin size vector
    C = zeros(shape=(size(D),1))
    # print(x_max,"\n",x_min, "\n",C,"\n",N,"\n",D)

#Computation of the cost function
    # fig = figure(figsize=(15,10))
    # ylabel("Number of events per bin")
    # xlabel("node_degree")
    N[N < 0] = 0
    for i in range(size(N)):
        edges = np.linspace(x_min,x_max,N[i]+1) # Bin edges
        ki = hist(x,edges) # Count # of events in bins
        ki = ki[0]
        k = mean(ki) #Mean of event count
        v = sum((ki-k)**2)/N[i] #Variance of event count
        C[i] = (2*k-v)/((D[i])**2) #The cost Function

#Optimal Bin Size Selection
    C = np.nan_to_num(C, nan=10.0)
    # print("*************",C,"*************")
    cmin = min(C)
    # print(cmin)
    idx  = where(C==cmin)
    idx = int(idx[0])
    optD = D[idx]
    # print(optD)
    # print('\n',optD,'\n')

    # edges = np.linspace(x_min,x_max,N[idx]+1)
    # fig1 = figure(figsize=(15,10))
    # ax = fig1.add_subplot(111)
    # ax.hist(x,edges)
# title(u"Histogram")
#     ylabel("Frequency")
#     xlabel("node_degree")
# # savefig('Hist.png')
#     fig2 = figure(figsize=(15,10))
#     plot(D,C,'.b',optD,cmin,'*r')
#     xlabel("Bin size")
#     ylabel("Estimated cost")

    return round(optD)

# pd.to_numeric(collect_vs_pseuds['Number of pseudogenes'])
# print(collect_vs_pseuds['Number of pseudogenes'])

def draw_pseudogenome_histogram(collect_vs_pseuds):

    normalized_len_data, FD = calculate_FD(collect_vs_pseuds['Number of pseudogenes'])

    def draw_hist(optimal_bin):
        plt.figure(figsize=(15,12))
        sns.set_style("ticks",{'axes.linewidth': 20.0})
        sns.set_context('talk', 
                    font_scale = 1.5)
        ax = sns.histplot(collect_vs_pseuds['Number of pseudogenes'],bins=optimal_bin,)
        sns.despine()
        plt.savefig(f"{working_dir}/secondary_analyses/roary/Pseudogene_number_histogram_{optimal_bin}.png", bbox_inches="tight", dpi=100)

    if normalized_len_data == 0:

        #normal seaborn plot
        plt.figure(figsize=(15,12))
        sns.set_style("ticks",{'axes.linewidth': 20.0})
        sns.set_context('talk', 
                    font_scale = 1.5)
        ax = sns.histplot(collect_vs_pseuds['Number of pseudogenes'])
        sns.despine()
        plt.savefig(f"{working_dir}/secondary_analyses/roary/Pseudogene_number_histogram_none.png", bbox_inches="tight", dpi=100)

        ###removing outlier
        outlier_removed = remove_outlier(collect_vs_pseuds,'Number of pseudogenes')
        normalized_len_data_out, FD_out = calculate_FD(outlier_removed['Number of pseudogenes'])
        if normalized_len_data_out > 0:
            optimal_bin_out = optimal_bin_width(outlier_removed['Number of pseudogenes'], normalized_len_data_out, FD_out)
            draw_hist(optimal_bin_out)
        else:
            pass
    elif normalized_len_data > 0:
        optimal_bin = optimal_bin_width(collect_vs_pseuds['Number of pseudogenes'], normalized_len_data, FD)
        draw_hist(optimal_bin)
        
draw_pseudogenome_histogram(collect_vs_pseuds)

################### plotting histogram for core and assessory pseudogenes
def plot_roary_histogram(roary_igm):
    plt.figure(figsize=(9, 6))
    plt.hist(roary_igm.sum(axis=1), roary_igm.shape[1],
         histtype="stepfilled", alpha=.7)

    plt.xlabel('Genomes',
    #fontdict={'fontsize':15}
    )
    plt.ylabel('Number of genes',
    #fontdict={'fontsize':15}
    )

    # sns.despine(top=True,
    #             right=True)
    sns.despine(right=True)
    plt.savefig(f"{working_dir}/secondary_analyses/roary/Pseudogene_roary_histogram.png", bbox_inches="tight", dpi=100)

plot_roary_histogram(roary_igm)

################### Pie chart plot
def draw_pie_chart(roary):
    # Plot the pangenome pie chart
    plt.figure(figsize=(10, 10))

    core = roary[roary.sum(axis=1) == roary.shape[1]].shape[0]
    softcore = roary[(roary.sum(axis=1) < roary.shape[1]) &
                     (roary.sum(axis=1) >= roary.shape[1]*0.95)].shape[0]
    shell = roary[(roary.sum(axis=1) < roary.shape[1]*0.95) &
                     (roary.sum(axis=1) >= roary.shape[1]*0.15)].shape[0]
    cloud = roary[roary.sum(axis=1) < roary.shape[1]*0.15].shape[0]

    total = roary.shape[0]

    def my_autopct(pct):
        val=int(pct*total/100.0)
        return '{v:d}'.format(v=val)

    a=plt.pie([core, softcore, shell, cloud],
          labels=['core\n(%d strains)'%roary.shape[1],
                  'soft-core\n(%d <= strains < %d)'%(roary.shape[1]*.95,
                                                     roary.shape[1]),
                  'shell\n(%d <= strains < %d)'%(roary.shape[1]*.15,
                                                 roary.shape[1]*.95),
                  'cloud\n(strains < %d)'%(roary.shape[1]*.15)],
          explode=[0.1, 0.05, 0.02, 0], radius=0.9,
          colors=[(0, 0, 1, float(x)/total) for x in (core, softcore, shell, cloud)],
          autopct=my_autopct)
    plt.savefig(f'{working_dir}/secondary_analyses/roary/Pseudogene_roary_pie.png', bbox_inches="tight", dpi=100)
draw_pie_chart(roary_igm)

####################
# sns.set(color_codes=True)
def draw_clustermap():
    
# sns.set(font_scale=1.5)
    gene_absence_presence_clustermap = sns.clustermap(roary_igm.T,cmap=plt.cm.Blues,
              figsize=(20, 15),col_cluster=False,yticklabels=False,xticklabels=False,
              dendrogram_ratio=(.2, .01), cbar_pos=(0, 0, .02, .1),
               cbar_kws = dict(ticks=[0,1]),
               robust=True,
            #   row_colors= [colors_loc,colors_yr,],
                  tree_kws={'linewidths':1.5,#'colors':[colmap[s] for s in countries_column]
                           })
    ax = gene_absence_presence_clustermap.ax_heatmap
# ax.tick_params(right=False)
    ax.set_xlabel("Pseudogenes", fontsize=20, labelpad=10.0)
    ax.set_ylabel("Isolates", fontsize=20, labelpad=10.0)

    plt.savefig(f'{working_dir}/secondary_analyses/roary/Pseudogene_roary_clustermap.png', bbox_inches="tight", dpi=100)

draw_clustermap()

############## DImensionality reduction ############################
############## MCA analysis
#### TODO Replace with a package conda can import from or find a way to get conda to recognize package
def mca_analysis (show_row_points=True,show_column_points=False,row_labels=True,show_column_labels=False):

    X=roary_igm.where(roary_igm == 0.0, 'present')
    X=X.where(X != 0.0, 'absent')

    # sns.set_style('white')
    mca = prince.MCA(n_components=2,n_iter=10, copy=True,check_input=True,engine='auto',random_state=42)
    mca = mca.fit(X.T)
    print(mca)

    ax = mca.plot_coordinates(X=X.T,ax=None,
                              figsize=(12, 10),show_row_points=show_row_points,
                              row_points_size=10,show_row_labels=row_labels,
                              show_column_points=show_column_points,column_points_size=10,
                              show_column_labels=show_column_labels)
    # print(row_coords,'\n',col_coords)
    plt.legend('', frameon=False)
    plt.savefig(f"{working_dir}/secondary_analyses/roary/Pseudogene_mca.png", bbox_inches="tight", dpi=100)
    return ax

mca_analysis(True,False,False,False)

########### UMAP
import umap

# roary_T_reset_index = roary_igm.T.reset_index()
# print(roary_T_reset_index)
save_path = f'{working_dir}/secondary_analyses/roary/Pseudogene_umap.png'
def umap_plotter(embedding, save_path):
    plt.figure(figsize=(12, 10))
#     sns.set(style='ticks')
    sns.set_style('white',{'axes.linewidth': 20.0})
    sns.set_context("paper", 
                    font_scale = 1.5
                   )
    # plt.rcParams['xtick.major.size'] = 20
    # plt.rcParams['xtick.major.width'] = 4
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True
    
    ax = sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1],
                    # hue=needed_metabolic_summary_cleaned.geo_loc_name,
                         #palette=country_color_dict,
#                          s=needed_metabolic_summary_cleaned.CD_rescaled,
                         edgecolor="black")
    # texts = []
    # for index,row in outsiders.iterrows():
    #     texts.append(plt.text(embedding[:, 0][index], y=embedding[:, 1][index],
    #                         s=outsiders.loc[index,'name'], 
    #                         fontsize=8,))
    sns.despine()
    # legend1 = plt.legend(loc=(1,0.3), title="Country", ncol=1)
    # produce a legend with a cross section of sizes from the scatter
#     handles, labels = ax.legend_elements(prop="sizes", alpha=0.6)
#     legend2 = ax.legend(handles, labels, loc="upper right", title="Sizes")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    # adjust_text(texts)
    plt.savefig(save_path, bbox_inches="tight", dpi=100)

def draw_umap(data, n_neighbors=15, min_dist=0.1, n_components=2, metric='euclidean', title='',save_path=save_path):
    fit = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        metric=metric,
        random_state=42
    )
    u = fit.fit_transform(data);
    fig = plt.figure()
    if n_components == 1:
        ax = fig.add_subplot(111)
        ax.scatter(u[:,0], range(len(u)),)# c=data.geo_loc_name)
    if n_components == 2:
        ax = fig.add_subplot(111)
        ax.scatter(u[:,0], u[:,1], )#c=data.geo_loc_name)
        umap_plotter(u,save_path)
    if n_components == 3:
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(u[:,0], u[:,1], u[:,2],)# c=data.geo_loc_name, s=100)
    plt.title(title, fontsize=18)

n_neighbours = 0.6 * roary_igm.T.shape[0]
n_neighbours = round(n_neighbours)
# print(n_neighbours)
draw_umap(roary_igm.T,n_neighbors=n_neighbours, 
    #title='n_neighbors = {}'.format(n),
    metric='correlation')

plt.close('all')