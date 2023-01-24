#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import sys
import logging

working_dir = sys.argv[1]


try:
  data = pd.read_csv(f"{working_dir}/secondary_analyses/dnds/dn_ds_analysis_summary.csv")
except BaseException:
    logging.exception("dn/ds file processing failed")
    sys.exit()

fig,ax = plt.subplots(figsize=(12, 10))
fig.tight_layout()
#     sns.set(style='ticks')
sns.set_style('white',{'axes.linewidth': 20.0})
sns.set_context("paper", 
                font_scale = 2
               )
# plt.rcParams['xtick.major.size'] = 20
# plt.rcParams['xtick.major.width'] = 4
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True

fg = sns.regplot(y='genes_average',x='all_average',data=data,lowess=True,
               ax=ax,#scatter=False,
                 fit_reg=False,
                 scatter_kws={#"color":'green',
                   'clip_on': False,
               },truncate=False,
                 #robust=True,
                 #line_kws={"color": "red"}
                )
# ax1 = sns.scatterplot(y='genes_average',x='all_average',data=data,
#                      hue='geo_loc_name',palette=country_color_dict,s=25#edgecolor="black"
#                      )
fg1 = sns.regplot(y='pseudogenes_average',x='all_average',data=data,lowess=True,
               ax=ax,#scatter=False,
                  fit_reg=False,
                  scatter_kws={#"color":'green',
                   'clip_on': False,
               },truncate=False,
                 #robust=True,
                 #line_kws={"color": "red"}
                )
# ax2 = sns.scatterplot(y='pseudogenes_average',x='all_average',data=data,
#                      hue='geo_loc_name',palette=country_color_dict,s=25#edgecolor="black"
#                      )
# ax2 = sns.lineplot(y='pseudogenes_total_num',x='phylo_order',data=phylo_pseudo_df_dropped)
sns.despine()
fg.set(xlim=(data.all_average.min()-0.02, data.all_average.max()+0.02))

# handles, labels = ax.get_legend_handles_labels()
# unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
# country_legend = plt.legend(*zip(*unique),loc=(1,0.6), title="Country", ncol=1)
# plt.gca().add_artist(country_legend)

#     handles2, labels2 = ax2.get_legend_handles_labels()
orange_patch = mpatches.Patch(color='tab:orange', label='Pseudogenes')
blue_patch = mpatches.Patch(color='tab:blue', label='Functional genes')
plt.legend(handles=[blue_patch,orange_patch],
           #handles2, labels2,
           loc=(0.8,0.8), title="Genomic element", ncol=1)

ax.set_xlabel('Whole genome average dn/ds')
ax.set_ylabel('Average dn/ds')
ax.figure.savefig(f"{working_dir}/secondary_analyses/dnds/dnds.png", bbox_inches="tight", dpi=100)

plt.close('all')