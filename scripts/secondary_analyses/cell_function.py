#!/usr/bin/env python3

# Plotting imports
import logging
import matplotlib.pyplot as plt
# import matplotlib
import seaborn as sns
import sys

# Other imports
import os
import pandas as pd
import numpy as np
# import prince
# from sklearn.preprocessing import StandardScaler
# from roary import draw_umap

sns.set_style('white')
############# KEGG ###################
working_dir = sys.argv[1]
app_dir = sys.argv[2]
cell_met_dir = f"{working_dir}/secondary_analyses/cell_and_metabolism"
os.makedirs(cell_met_dir,exist_ok=True)


################# COG ##################
def build_cog_matrix(path,element):
    root_folder = path
    full_cog_table = pd.DataFrame()
    count = 0
    try:
        for root, _, files in os.walk(root_folder):
            for file in files:
                if file.endswith("count.tsv"):
                    cog_table = pd.read_csv(f'{root}/{file}',sep='\t')
                    rename = root.split('/')[-1]
    #                 print(rename)
                    cog_table.rename(columns={'COUNT':rename},inplace=True)
                    cog_table[rename] = ((cog_table[rename]/sum(cog_table[rename]))*100).round(2)
                    if count == 0:
                        full_cog_table = pd.concat([full_cog_table,cog_table])
                        cols = full_cog_table.columns.tolist()
                        cols = [cols[2],cols[0],cols[3],cols[1]]
                        full_cog_table = full_cog_table[cols]
                        count +=1
                    else:
                        cog_table = cog_table[['DESCRIPTION',rename]]
                        full_cog_table = full_cog_table.merge(cog_table, on='DESCRIPTION',how='outer')
    #                     full_cog_table = full_cog_table.apply(lambda x: (x/sum(x)*100).round(2), axis=1)
    except BaseException:
        logging.exception(f"{element.capitalize()} COG post-analysis failed")
        # sys.exit()
        pass

    if 'pseudo' in path:
        full_cog_table.columns = full_cog_table.columns.str.split('_',expand=True).get_level_values(0)
    return full_cog_table

# cog_dir = f"{working_dir}/pan_files/pseudogenome/COG"
# full_cog_table = build_cog_matrix(cog_dir)

def plot_cog_heatmap(full_cog_table,save_path=cell_met_dir,element='Pseudogenome',x_ticks='F'):
    drop_filter = full_cog_table.filter(['COLOR','LETTER','DESCRIPTION'])
    heatmap = sns.clustermap(full_cog_table.drop(drop_filter, axis=1), 
                             xticklabels=x_ticks,
                             yticklabels=True, 
                                figsize=(30,25), 
                             row_cluster=False, 
    #                          col_cluster=True, 
                                 cbar_pos=(1, 1, 0.2, 0.04),
                                cbar_kws= {'orientation': 'horizontal',
                                           'label': 'Number of genes'}, 
                                 dendrogram_ratio=0.1,
                              tree_kws={'linewidths':2,})
    heatmap.ax_heatmap.set_yticklabels(full_cog_table.DESCRIPTION.to_list(), fontsize = 30,
                                      fontdict={'bbox':{'ec':'black','alpha':0.01}, 
                                                'fontweight':'bold',
                                               'rotation':'horizontal'})
#     heatmap.ax_heatmap.set_xlabel(x_label, fontsize=35, fontweight='bold')
    heatmap.ax_heatmap.set_ylabel('Clusters of Orthologous Genes (COGs)', fontsize=35, fontweight='bold')
    heatmap.ax_heatmap.set_xlabel(element, fontsize=35, fontweight='bold')
    cbar = heatmap.ax_heatmap.collections[0].colorbar
    # here set the labelsize by 20
    cbar.ax.tick_params(labelsize=20)
    cbar.ax.set_xlabel(f'Percentage per {element.lower()}',fontsize=25)
    #     plt.setp(heatmap.ax_heatmap.get_yticklabels(), bbox={'boxstyle':'round','ec':'black','alpha':0.3})
    # plt.setp(heatmap.ax_heatmap.get_yticklabels(), bbox={'ec':'black','alpha':0.01}, fontweight='bold')
#     plt.setp(heatmap.ax_heatmap.get_xticklabels(), fontsize=10, fontweight='bold')
    plt.savefig(save_path + "/" + element + "_heatmap.png", bbox_inches="tight", dpi=100)

# plot_cog_heatmap(full_cog_table=full_cog_table,element='Pseudogenome')


for element in ['pseudogenome','whole_genome','functional_genome']:
    try:
        cog_dir = f"{working_dir}/pan_files/{element}/COG"
        full_cog_table = build_cog_matrix(cog_dir,element)
        full_cog_table.to_csv(f"{cell_met_dir}/{element.capitalize()}_full_cog_table.csv")
        plot_cog_heatmap(full_cog_table=full_cog_table,element=element.capitalize())
        plt.clf()
    except BaseException:
        logging.exception(f"{element.capitalize()} COG post-analysis failed")
        # sys.exit()
        pass

###Functional genome
# try:
#     cog_dir = f"{working_dir}/pan_files/functional_genome/COG"
#     full_cog_table = build_cog_matrix(cog_dir)
#     plot_cog_heatmap(full_cog_table=full_cog_table,element='functional_genome')
# except BaseException:
#         logging.exception("functional_genome COG post-analysis failed")
#         # sys.exit()
#         pass
