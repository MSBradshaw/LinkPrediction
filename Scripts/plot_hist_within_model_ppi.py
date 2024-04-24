import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import networkx as nx

def get_args():
    parser = argparse.ArgumentParser(description='Plot histogram of groupings within one model and KG')
    parser.add_argument('-i', type=str, help='input directory path')
    parser.add_argument('-o', type=str, help='output plot path')
    parser.add_argument('-t', type=str, help='title prefix')
    # ensure that input dir ends with /
    return parser.parse_args()

def plot_tsv_in_dir(dpath,output,by_median, title_prefix):
    file2color = {'Female':'#b66dff','Male':'#006ddb','Cancer':'#920000','PedCancer':'#ff6db6','RareDisease':'#004949', 'UltraRareDisease':'#009999','African':'#db6d00','EastAsian':'#22cf22','European':'#8f4e00','Latino':'#490092','RandomGenes':'#676767','RandomDiseases':'#616161'}
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
    fig, axes = plt.subplots(nrows=1, ncols=len(file_prefixes), figsize=(5*len(files),5))
    for i,prefix in enumerate(file_prefixes):
        f = f"{dpath}{prefix}_ranking_data.tsv"
        try:
            df = pd.read_csv(f, sep='\t')
        except:
            print(f'Error reading {f}')
            continue
        if by_median:
            # aggregate df according to
            if prefix in gene_groups:
                # group by head_label
                df = df.groupby('head_label').median('rank').reset_index()
            else:
                # group by tail_label
                df = df.groupby('tail_label').median('rank').reset_index()
        # plot
        try:
            axes[i].hist(df['rank'], bins=bins, color=file2color[prefix])
        except KeyError:
            print(f'Error plotting {f}')
            continue
        median = df['rank'].median()
        axes[i].set_title(title_prefix + '\n' + file_names[i] + f'\nmedian={median:.2f}')
        axes[i].set_xlabel('rank')
        axes[i].set_ylabel('count')
        # add median line
        axes[i].axvline(median, color='r', linestyle='dashed', linewidth=1)
        # remove top and right spines
        axes[i].spines['top'].set_visible(False)
        axes[i].spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(output,dpi=300)
    plt.show()


def main():    
    args = get_args()
    if not args.i.endswith('/'):
        args.i += '/'
    plot_tsv_in_dir(args.i, args.o, False, args.t)
    
if __name__ == '__main__':
    main()