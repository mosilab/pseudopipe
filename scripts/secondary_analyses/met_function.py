#!/usr/bin/env python3

# Plotting imports
import logging
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import sys

# Other imports
import os
import pandas as pd
import numpy as np
import prince
from sklearn.preprocessing import StandardScaler
from roary import draw_umap

sns.set_style('white')
############# KEGG ###################
working_dir = sys.argv[1]
app_dir = sys.argv[2]
cell_met_dir = f"{working_dir}/secondary_analyses/cell_and_metabolism"
os.makedirs(cell_met_dir,exist_ok=True)

def process_metabolic_summary():
    try:
        metabolic_summary_df = pd.read_csv(f"{working_dir}/pan_files/pseudogenome/KEGG/metabolic_summary__module_completeness.tab",sep='\t',index_col=1)
        metabolic_summary_df_dropped = metabolic_summary_df.drop(['module','pathway group'],axis=1)
        metabolic_summary_df_dropped.rename(columns=lambda x: x.split('_')[0], inplace=True)
        metabolic_summary_prepped = metabolic_summary_df_dropped.T
        return metabolic_summary_prepped
    except BaseException:
        logging.exception("Pseudogenome KEGG metabolic file processing failed")
        # sys.exit()
        pass

metabolic_summary_prepped = process_metabolic_summary()

################## PCA analysis
# X = metabolic_summary_prepped.dropna().set_index('name')
try:
    X = metabolic_summary_prepped.dropna()

    pca = prince.PCA(n_components=2,n_iter=10,
                    rescale_with_mean=True,rescale_with_std=True,copy=True,
                    check_input=True,engine='auto',random_state=42)
    pca = pca.fit(X)
    ax = pca.plot_row_coordinates(X,
                                ax=None,figsize=(12, 10),x_component=0,
                                y_component=1,labels=None,
                                #   color_labels=needed_metabolic_summary_cleaned.geo_loc_name,
                                ellipse_outline=False,ellipse_fill=True,show_points=True)
    plt.savefig(f"{cell_met_dir}/PCA.png", bbox_inches="tight", dpi=100)
    plt.clf()
except BaseException:
        logging.exception("Pseudogenome KEGG metabolic PCA failed")
        # sys.exit()
        pass

################## UMAP

scaled_data = StandardScaler().fit_transform(X.values)
draw_umap(scaled_data, n_neighbors=100, 
    #title='n_neighbors = {}'.format(n),
    metric='correlation',save_path=f"{cell_met_dir}/UMAP.png")
plt.clf()

########### KEGG heatmaps

def ma_table_cleaning(tab_file_path):
    proteome = pd.read_csv(tab_file_path,
                            sep="\t")
    proteome_cols = [idx.split('.')[0] for idx in list(proteome.columns)]
    proteome_counter=0
    prot_col_dict={}
    for i in range (len(list(proteome.columns))):
        prot_col_dict[list(proteome.columns)[proteome_counter]]=proteome_cols[proteome_counter]
        proteome_counter+=1
    proteome = proteome.rename(columns=prot_col_dict)
    return proteome

def pseuds_ma_table_cleaning(tab_file_path):
    pseudogenome = pd.read_csv(tab_file_path,sep="\t")
    pseudogenome.columns = pseudogenome.columns.str.split('_',expand=True).get_level_values(0)
    
    return pseudogenome

# final in-use
def plot_heatmap(metabolism_matrix, module_information, module_names, cluster, prefix, x_lab, save_path):
    matplotlib.rcParams['pdf.fonttype'] = 42
#     matplotlib.use('Agg')
    print("Creating output_file matrix and heatmap... ")
    # Check clustering asked
    if cluster is None:
        row_cluster = False
        col_cluster = False
    elif cluster == 'rows':
        row_cluster = True
        col_cluster = False
    elif cluster == 'cols':
        row_cluster = False
        col_cluster = True
    elif cluster == 'both':
        row_cluster = True
        col_cluster = True
    else:
        exit("Wrong clustering option, select between 'both', 'cols', 'rows' or don't pass this option.")
#     # Extract data from metabolism matrix and color dictionary
#     for genome, modules in metabolic_annotation.items():
#         # Fill metabolism matrix with completeness from the metabolic annotation dictionary
#         for module in modules:
#             metabolism_matrix.loc[module[0],genome] = module[1]
    # Get module id and its corresponding module name and the module and its corresponding color
    module_id_name = {}
    module_colors = {}
    for module, information in module_information.items():
        module_id_name[module] = information[0]
        module_colors[information[0]] = (information[1],information[2])
#     print(module_colors)
    # Get modules that are above 50% complete in at least one genome
    ylabel_text = 'Modules (at least 50% complete in at least one genome)'
    metabolism_matrix_retained = metabolism_matrix.loc[(metabolism_matrix >= 50).any(1)]
    metabolism_matrix_retained[metabolism_matrix_retained < 0] = 0
    table_empty = (metabolism_matrix_retained == 0).all()
    if table_empty.all() == True:
        print("There are no modules above 50% completeness. Plotting all modules regardless of completeness.")
        metabolism_matrix_retained = metabolism_matrix.loc[(metabolism_matrix > 0).any(1)]
        ylabel_text = 'Modules'
    colors_for_ticks = []
    colors_for_legend = {}
#     metabolism_matrix_retained_relabel = metabolism_matrix_retained.rename(index=module_id_name)
    if 'name' not in metabolism_matrix_retained.columns:
        metabolism_matrix_retained.insert(0,'name',module_names)
    metabolism_matrix_retained_relabel = metabolism_matrix_retained.set_index('name')

    # Plot heatmap with clustering
    # Check if clustering is possible with the number of genomes and modules
    table_dim = metabolism_matrix_retained_relabel.shape
    if row_cluster == True and table_dim[0] < 3:
        row_cluster = False
    if col_cluster == True and table_dim[1] < 3:
        col_cluster = False
    sns.set(font_scale=1.6)
    heatmap = sns.clustermap(metabolism_matrix_retained_relabel, xticklabels=True, yticklabels=True, 
                            figsize=(30,25), row_cluster=row_cluster, col_cluster=col_cluster, 
                             cbar_pos=(0.4, 1, 0.2, 0.04),
                            cbar_kws= {'orientation': 'horizontal', 'label': 'Module Completeness (%)'}, 
                             dendrogram_ratio=0.1,
                          tree_kws={'linewidths':1.5,})
    # cbar = heatmap.ax_cbar
# here set the labelsize by 20
    # cbar.set_yticklabels(labels=cbar.get_ylabel, fontdict={"fontsize":20})
    for module in heatmap.ax_heatmap.yaxis.get_ticklabels():
#         print (module.get_text())
        colors_for_ticks.append(module_colors[module.get_text()][1])
        if module_colors[module.get_text()][0] not in colors_for_legend:
            colors_for_legend[module_colors[module.get_text()][0]] = module_colors[module.get_text()][1]
    [mod.set_color(col) for (col,mod) in zip(colors_for_ticks,heatmap.ax_heatmap.yaxis.get_ticklabels())]
    for pathway, color in colors_for_legend.items():
        heatmap.ax_heatmap.scatter(1,1, color=color, label=pathway, zorder=0, s=45)
#     heatmap.ax_heatmap.set_xticklabels(heatmap.ax_heatmap.get_xmajorticklabels(), fontsize = 20)
    heatmap.ax_heatmap.set_yticklabels(heatmap.ax_heatmap.get_ymajorticklabels(), fontsize = 18)
    heatmap.ax_heatmap.set_xlabel(x_lab, fontsize=25, fontweight='bold')
    heatmap.ax_heatmap.set_ylabel(ylabel_text, fontsize=25, fontweight='bold')
    plt.figlegend(loc=(0.8,0.5), ncol=1, fontsize=18, markerscale=3)
#     plt.setp(heatmap.ax_heatmap.get_yticklabels(), bbox={'boxstyle':'round','ec':'black','alpha':0.3})
    plt.setp(heatmap.ax_heatmap.get_yticklabels(), bbox={'ec':'black','alpha':0.01}, fontweight='bold')
    plt.setp(heatmap.ax_heatmap.get_xticklabels(), fontsize=5, fontweight='bold')
    plt.savefig(f"{save_path}/{prefix}" + "_KEGG_heatmap.png", bbox_inches="tight", dpi=100)
    print("################Finished###############")
    # plt.show()
    return metabolism_matrix_retained_relabel, module_colors

def plot_heatmap_80(module_colors, module_group_matrix, metabolism_matrix_dropped_relabel, prefix, x_lab, save_path):
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.use('Agg')
    print("Grouping by metabolism type and plotting... ")
#     print(metabolism_matrix_dropped_relabel)
    for module in list(metabolism_matrix_dropped_relabel.index):
        for genome in list(metabolism_matrix_dropped_relabel.columns):
            if metabolism_matrix_dropped_relabel.loc[module,genome] >= 80:
                module_group_matrix.loc[module_colors[module][0],genome] += 1
#     print(module_group_matrix)
    module_group_matrix_transp = module_group_matrix.T
    emptyness = (module_group_matrix_transp == 0).all()
    if emptyness.all() == True:
        print("There are no modules above 80% completeness. No barplot will be generated.")
    else:
        module_group_matrix_transp = module_group_matrix_transp.loc[(module_group_matrix_transp >= 1).any(1)]
        module_group_matrix_transp = module_group_matrix_transp.loc[:, (module_group_matrix_transp != 0).any(0)]
        color_dict = {}
        color_list = []
        for pair in module_colors.values():
            if pair[0] not in color_dict:
                color_dict[pair[0]] = pair[1]
#         print(list(module_group_matrix_transp.columns))
        for pathway in list(module_group_matrix_transp.columns):
            color_list.append(color_dict[pathway])
#         print(len(color_list))

        Figure, Axis = plt.subplots(1,1, figsize=(35,20))
        sns.heatmap(module_group_matrix_transp.T,ax=Axis,xticklabels=1, yticklabels=1
#                     cmap=color_list,
#                     linewidths=2, 
                    #figsize=(25,15), legend=False, width=0.8
                   )
#         Axis.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize='medium', markerscale=0.3)
        Axis.set_xlabel(x_lab, fontsize=30, fontweight='bold')
        Axis.set_ylabel('Number of Modules (>=80% complete)', fontsize=30, fontweight='bold')
        Axis.tick_params(axis='x', labelsize=8)
        Axis.tick_params(axis='y', labelsize=20)
#         Axis.set_facecolor('white')
#         Figure.suptitle('Metabolism Module Category per Genome', fontsize=40)
        Figure.subplots_adjust()
        Figure.savefig(f"{save_path}/{prefix}" + "_KEGG_barplot.png", bbox_inches="tight", #transparent=True, 
                       dpi=100)
    print("Finished")
    # plt.show()
    return module_group_matrix_transp

def making_metabolic_plots(proteome,cluster,file_prefix,axis_lable,save_path):
    proteome_subset = proteome.drop(['module','name','pathway group'],axis=1)
    metabolism_matrix=proteome_subset.copy()
    
    module_information = {}
    module_ids = []
    group_names = []
    with open(f"{app_dir}/scripts/secondary_analyses/00.Module_Names.txt") as correspondence:
        for line in correspondence:
            line = line.strip().split("\t")
            module_information[line[0]] = [line[1], line[2], line[3]]
            # module_id_name[line[0]] = line[1]
            module_ids.append(line[0])
            # module_groups[line[0]] = line[2]
            if line[2] not in group_names:
                group_names.append(line[2])
    
    metabolism_matrix_dropped_relabel, module_colors = plot_heatmap(metabolism_matrix,
                                                                             module_information,
                                                                             proteome.name,
                                                                       cluster,
                                                                       file_prefix, 
                                                                       axis_lable,save_path)
    
    genome_names = list(proteome.columns[:3]) + list(sorted(proteome.columns[3:]))
    # Create module pathway group matrix
    module_group_matrix_80 = pd.DataFrame(0, index=group_names, columns=genome_names)
    # Create module matrix
    metabolism_matrix_80 = pd.DataFrame(0, index=module_ids, columns=genome_names)
    plot_heatmap_80(module_colors, module_group_matrix_80, 
                       metabolism_matrix_dropped_relabel, file_prefix, axis_lable,save_path)


element_dict = {'pseudogenome':'Pseudogenome','whole_genome':'Whole genome','functional_genome':'Functional genome'}
for element in ['pseudogenome','whole_genome','functional_genome']:
    met_summary_file_path = f"{working_dir}/pan_files/{element}/KEGG/metabolic_summary__module_completeness.tab"
    if os.path.isfile(met_summary_file_path):
        if element == 'pseudogenome':
            element_KEGG = pseuds_ma_table_cleaning()
        else:
            element_KEGG = ma_table_cleaning(met_summary_file_path)
        # capitalized_element = element.capitalize()
        making_metabolic_plots(element_KEGG, 'cols', element, element_dict[element],cell_met_dir)
        plt.clf()
    else:
        # raise FileNotFoundError("KEGG Metabolic Summary file not found. This analysis will be skipped!")
        pass

# whole_genome_KEGG = ma_table_cleaning(f"{working_dir}/pan_files/whole_genome/KEGG/metabolic_summary__module_completeness.tab")
# making_metabolic_plots(pseudogenome_KEGG,'cols','whole_genome','Whole genome')

# functional_genome_KEGG = ma_table_cleaning(f"{working_dir}/pan_files/functional_genome/KEGG/metabolic_summary__module_completeness.tab")
# making_metabolic_plots(pseudogenome_KEGG,'cols','functional_genome','Functional genome')
