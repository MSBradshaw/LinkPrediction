import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os
import pickle

def get_mnrs(dpath):
    medians = []
    file_prefixes = ['Female','Male','Cancer','PedCancer','RareDisease', 'UltraRareDisease','African','EastAsian','European','Latino','RandomGenes','RandomDiseases']
    file_names = ['Female','Male','Cancer','Pediatric Cancer','Rare Disease', 'Ultra Rare Disease','African','East Asian','European','Latino','Random Genes','Random Diseases']
    gene_groups = ['Female','Male','Cancer','PedCancer','Random']
    # ensure dpath ends with '/'
    if dpath[-1] != '/':
        dpath += '/'
    # 'RankResults/original_monarch/TransE/'
    # list all files in dpath
    files = [dpath+f for f in os.listdir(dpath) if f.endswith('.tsv')]
    # bins 0.0 to 1.0 by 0.05
    bins = [i/20 for i in range(21)]
    for i,prefix in enumerate(file_prefixes):
        f = f"{dpath}{prefix}_ranking_data.tsv"
        try:
            df = pd.read_csv(f, sep='\t')
        except:
            print(f'Error reading {f}')
            continue

        median = df['rank'].median()
        medians.append(median)
    return medians, file_prefixes

def count_hgnc_edges(el_path):
    # read in the edge list
    el = pd.read_csv(el_path, sep='\t')
    el.columns = ['subject','predicate','object']
    # if HGNC prefix in subject and object keep
    el = el[el['subject'].str.contains('HGNC') & el['object'].str.contains('HGNC')]
    # count the number of edges
    return el.shape[0]

def main():
    file2color = {'Female':'#b66dff','Male':'#006ddb','Cancer':'#920000','PedCancer':'#ff6db6','RareDisease':'#004949', 'UltraRareDisease':'#009999','African':'#db6d00','EastAsian':'#22cf22','European':'#8f4e00','Latino':'#490092','RandomGenes':'#676767','RandomDiseases':'#616161'}
    gene_groups = ['Female','Male','Cancer','PedCancer','RandomGenes']
    # get the number of HGNC to HGNC edges in Monarch KG Filtered and each of the ppis
    mkg_edges = count_hgnc_edges('ELs_for_Rotate/Monarch_KG_Filtered/train.txt')
    huri_edges = count_hgnc_edges('ELs_for_Rotate/Monarch_HuRI_Filtered/train.txt')
    s100_edges = count_hgnc_edges('ELs_for_Rotate/Monarch_STRING_t100_new/train.txt')
    s50_edges = count_hgnc_edges('ELs_for_Rotate/Monarch_STRING_t50_new/train.txt')
    s25_edges = count_hgnc_edges('ELs_for_Rotate/Monarch_STRING_t25_new/train.txt')
    ppi_counts = [mkg_edges,huri_edges,s100_edges,s50_edges,s25_edges,]
    # for RotatE get the HuRI filtered, string 100, string 50, string 25 MNR for each group
    res = {'ppi':[],'mnr':[],'group':[],'model':[],'ppi_percent':[]}
    ppis = ['original_monarch','HuRI_filtered','string_filtered_t100', 'string_filtered_t50', 'string_filtered_t25']
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
                res['ppi_percent'].append(ppi_counts[i]/mkg_edges*100)
    df = pd.DataFrame(res)
    # pickle the df
    # df.to_pickle('Data/ppi_percent_vs_mnr.pkl')

    fig, ax = plt.subplots(2,4)
    # set fig size
    fig.set_size_inches(18,5)
    for i,model in enumerate(models):
        # for each group, plot x = ppi_percent y = mnr, color using file2color
        for j,group in enumerate(df['group'].unique()):
            tmp = df[(df['model'] == model) & (df['group'] == group)]
            # sort by ppi_percent
            tmp = tmp.sort_values(by='ppi_percent')
            index = 1
            if group in gene_groups:
                index = 0
            ax[index,i].plot(tmp['ppi_percent'],tmp['mnr'],color=file2color[group],label=group)
        for j in range(2): 
            ax[j,i].set_title(model)
            ax[j,i].set_xlabel('PPI Percent')
            ax[j,i].set_ylabel('Median Normalized Rank')
            # remove top and right spines
            ax[j,i].spines['top'].set_visible(False)
            ax[j,i].spines['right'].set_visible(False)
            # set the y range to be 0.0 to 1.0
            ax[j,i].set_ylim(0.0,1.0)
    
    # make a custom legend and place it in 1,0
    # make a set of lines for the gene group and one for the non gene group

    g_lines = [plt.Line2D([0], [0], color=file2color[group], linewidth=3, linestyle='-') for group in df['group'].unique() if group in gene_groups]
    g_labels = [group for group in df['group'].unique() if group in gene_groups]
    ax[0, -1].legend(g_lines, g_labels, loc='center', title_fontsize='medium',frameon=False)
    ax[0, -1].axis('off')

    ng_lines = [plt.Line2D([0], [0], color=file2color[group], linewidth=3, linestyle='-') for group in df['group'].unique() if group not in gene_groups]
    ng_labels = [group for group in df['group'].unique() if group not in gene_groups]
    ax[1, -1].legend(ng_lines, ng_labels, loc='center', title_fontsize='medium',frameon=False)
    ax[1, -1].axis('off')

    plt.tight_layout()
    plt.savefig('Figures/ppi_percent_vs_mnr.png',dpi=300)


if __name__ == '__main__':
    main()