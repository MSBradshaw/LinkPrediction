import pandas as pd
import matplotlib.pyplot as plt
import os

# add Scripts to path
import sys
sys.path.append('../Scripts/')
from plot_ppi_percent_vs_mnr import get_mnrs



file2color = {'Female':'#b66dff','Male':'#006ddb','Cancer':'#920000','PedCancer':'#ff6db6','RareDisease':'#004949', 'UltraRareDisease':'#009999','African':'#db6d00','EastAsian':'#22cf22','European':'#8f4e00','Latino':'#490092','RandomGenes':'#676767','RandomDiseases':'#616161'}
gene_groups = ['Female','Male','Cancer','PedCancer','RandomGenes']

# for RotatE get the HuRI filtered, string 100, string 50, string 25 MNR for each group
res = {'ppi':[],'mnr':[],'group':[],'model':[]}
ppis = ['original_monarch','monarch_filtered']
models = ['ComplEx','RotatE','TransE']
# check if Data/ppi_percent_vs_mnr.pkl exists, if so load it and plot

for model in models:
    for i,p in enumerate(ppis):
        tmp_median, tmp_names = get_mnrs(f'RankResults/{p}/{model}/') 
        for j in range(len(tmp_median)):
            res['ppi'].append(p)
            res['mnr'].append(tmp_median[j])
            res['group'].append(tmp_names[j])
            res['model'].append(model)

df = pd.DataFrame(res)

file_prefixes = ['Female','Male','Cancer','PedCancer','RareDisease', 'UltraRareDisease','African','EastAsian','European','Latino','RandomGenes','RandomDiseases']
file_names = ['Female','Male','Cancer','Pediatric Cancer','Rare Disease', 'Ultra Rare Disease','African','East Asian','European','Latino','Random Genes','Random Diseases']

prefix2name = {prefix:name for prefix,name in zip(file_prefixes,file_names)}

fig, ax = plt.subplots(1,4,figsize=(20,5))

lines = [] # for legend
for i, model in enumerate(df['model'].unique()):
    for group in df['group'].unique():
        tmp = df[(df['model']==model) & (df['group']==group)]
        # assert there are two rows left
        assert tmp.shape[0] == 2
        # x = ppi == original_monarch
        # y = ppi == monarch_filtered
        x = tmp[tmp['ppi']=='original_monarch']['mnr'].values[0]
        y = tmp[tmp['ppi']=='monarch_filtered']['mnr'].values[0]
        ax[i].scatter(x,y, color=file2color[group], label=group)
        lines.append(plt.Line2D([0], [0], color=file2color[group], lw=4))
        # set x range to 0.5 to 1.0
        ax[i].set_xlim(0.5,1.0)
        # set y range to 0.5 to 1.0
        ax[i].set_ylim(0.5,1.0)
        # plot a horizontal line of slope 1
        ax[i].plot([0.5,1.0],[0.5,1.0], color='black', linestyle='dashed', linewidth=1)
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
    ax[i].set_title(model)
    ax[i].set_xlabel('MNR Full KG')
    ax[i].set_ylabel('MNR Filtered KG')

# sort legend so it goes gene_groups first then the rest, with a black space in between
gene_lines = [lines[i] for i in range(len(df['group'].unique())) if df['group'].unique()[i] in gene_groups]
not_gene_lines = [lines[i] for i in range(len(df['group'].unique())) if df['group'].unique()[i] not in gene_groups]
not_gene_names = [df['group'].unique()[i] for i in range(len(df['group'].unique())) if df['group'].unique()[i] not in gene_groups]

names = gene_groups + ['   '] + not_gene_names
lines = gene_lines + [plt.Line2D([0], [0], color='white', lw=4)] + not_gene_lines
names = [ prefix2name[x] if x in prefix2name else '' for x in names ]
ax[-1].legend(lines, names, loc='center', frameon=False)
ax[-1].axis('off')

plt.tight_layout()
plt.savefig('Figures/full_vs_filtered_mnr_scatter.png',dpi=300)
plt.show()

