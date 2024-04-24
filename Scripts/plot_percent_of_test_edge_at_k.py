import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
import concurrent.futures
import random

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='Input dirrectory')
    parser.add_argument('-t', type=str, help='title prefix')
    parser.add_argument('-o', type=str, help='Output image')
    return parser.parse_args()

def plot_kurve_of_tsv_in_dir(dpath,output,title_prefix):
    file2color = {'Female':'#b66dff','Male':'#006ddb','Cancer':'#920000','PedCancer':'#ff6db6','RareDisease':'#004949', 'UltraRareDisease':'#009999','African':'#db6d00','EastAsian':'#22cf22','European':'#8f4e00','Latino':'#490092','RandomGenes':'#676767','RandomDiseases':'#616161'}
    file_prefixes = ['Female','Male','Cancer','PedCancer','RareDisease', 'UltraRareDisease','African','EastAsian','European','Latino','RandomGenes','RandomDiseases']
    file_names = ['Female','Male','Cancer','Pediatric Cancer','Rare Disease', 'Ultra Rare Disease','African','East Asian','European','Latino','Random Genes','Random Diseases']
    gene_groups = ['Female','Male','Cancer','PedCancer','RandomGenes']
    # ensure dpath ends with '/'
    if dpath[-1] != '/':
        dpath += '/'
    # 'RankResults/original_monarch/TransE/'
    # list all files in dpath
    files = [dpath+f for f in os.listdir(dpath) if f.endswith('.tsv')]
    # bins 0.0 to 1.0 by 0.05
    bins = [i/20 for i in range(21)]
    fig, axes = plt.subplots( ncols=3, nrows=1, figsize=(13,5), gridspec_kw={'width_ratios': [5, 5, 3]})
    for i,prefix in enumerate(file_prefixes):
        f = f"{dpath}{prefix}_ranking_data.tsv"
        try:
            df = pd.read_csv(f, sep='\t')
        except:
            print(f'Error reading {f}')
            continue
        
        # create plotting data
        percent_hits = []
        ks = [1,5,10,15,20,25,50,100] + list(range(100,1001,100))
        for k in ks:
            hits = df[df['rank_integer'] <= k].shape[0]
            if df.shape[0] == 0:
                percent_hits.append(0)
            else:
                percent_hits.append(hits/df.shape[0] * 100)

        # plot
        if prefix in gene_groups:
            axes[0].plot(ks, percent_hits,  color=file2color[prefix], label=file_names[i])
        else:
            axes[1].plot(ks, percent_hits,  color=file2color[prefix], label=file_names[i])
    axes[0].set_title(title_prefix)
    axes[0].set_xlabel('k')
    axes[0].set_ylabel('Percent of test set edges')
    axes[1].set_xlabel('k')
    axes[1].set_ylabel('Percent of test set edges')
    # add median line
    # remove top and right spines
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    
    # create a custom legend that will go in axes[2]
    custom_lines = [plt.Line2D([0], [0], color=file2color[prefix], lw=4) for prefix in gene_groups]
    custom_lines += [plt.Line2D([0], [0], color='white', lw=4, linestyle='dashed')]
    custom_lines += [plt.Line2D([0], [0], color=file2color[prefix], lw=4) for prefix in file_prefixes if prefix not in gene_groups]
    axes[2].legend(custom_lines, gene_groups + ['   '] + [prefix for prefix in file_prefixes if prefix not in gene_groups], loc='center', frameon=False)
    axes[2].axis('off')

    axes[0].set_ylim(0,50)
    axes[1].set_ylim(0,50)
    plt.tight_layout()
    plt.savefig(output,dpi=300)
    plt.show()

def main():
    args = get_args()
    if args.i[-1] != '/':
        args.i += '/'
    plot_kurve_of_tsv_in_dir(args.i, args.o, args.t)


if __name__ == '__main__':
    main()