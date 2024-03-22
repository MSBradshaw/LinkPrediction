import matplotlib.pyplot as plt
import argparse
import pandas as pd

def get_args():
    # input dir, title, outfile
    parser = argparse.ArgumentParser(description='Plot hits at k curve')
    parser.add_argument('--ppi_model_dir', type=str, help='Directory containing the ppi model results')
    parser.add_argument('--title', type=str, help='Title of the plot')
    parser.add_argument('--outfile', type=str, help='Output file')
    return parser.parse_args()

def plot_kurve(ppi_model_dir:str,title:str,outfile:str):
    if ppi_model_dir[-1] != '/':
        ppi_model_dir += '/'
    
    dirs2color = {'Female':'#b66dff','Male':'#006ddb','Cancer':'#920000','PedCancer':'#ff6db6','RareDisease':'#004949', 'UltraRareDisease':'#009999','African':'#db6d00','EastAsian':'#22cf22','European':'#8f4e00','Latino':'#490092','Random':'#676767','RandomDiseases':'#e6e6e6'}
    dirs = ['Female','Male','Cancer','PedCancer','RareDisease', 'UltraRareDisease','African','EastAsian','European','Latino','Random','RandomDiseases']
    dirs_names = ['Female','Male','Cancer','Pediatric Cancer','Rare Disease', 'Ultra Rare Disease','African','East Asian','European','Latino','Random Genes','Random Diseases']
    
    gene_dir_names = ['Female','Male','Cancer','Pediatric Cancer','Random']
    gene_groups = ['Female','Male','Cancer','PedCancer','Random']
    non_gene_dir_names = ['Rare Disease', 'Ultra Rare Disease','African','East Asian','European','Latino']
    
    fig, axes = plt.subplots(2,2,figsize=(10,10), gridspec_kw={'width_ratios': [3, 1]})
    for i,d in enumerate(dirs):
        try:
            df = pd.read_csv(ppi_model_dir + d + '/avg_hits_at_k.tsv',sep='\t')
        except FileNotFoundError:
            print(f'No file found for {d}')
            continue
        if d in gene_groups:
            axes[0,0].plot(df['k'], df['percent_hits@k'], label=dirs_names[i],color=dirs2color[d])
        else:
            axes[1,0].plot(df['k'], df['percent_hits@k'], label=dirs_names[i],color=dirs2color[d])
    axes[0,0].set_xlabel('k')
    axes[0,0].set_ylabel('mean % hits@k')
    
    axes[1,0].set_xlabel('k')
    axes[1,0].set_ylabel('mean % hits@k')

    axes[0,0].set_title(title)
    axes[1,0].set_title(title)
    
    axes[0,0].spines['top'].set_visible(False)
    axes[0,0].spines['right'].set_visible(False)
    axes[1,0].spines['top'].set_visible(False)
    axes[1,0].spines['right'].set_visible(False)
    # axes[0].legend()

    # create a custom legend to put in axes[1] of the dirs_names and colors
    custom_lines_genes = [plt.Line2D([0], [0], color=dirs2color[d], lw=2) for d in dirs if d in gene_groups]
    axes[0,1].legend(custom_lines_genes, gene_dir_names, loc='center', title='Group', frameon=False)

    custom_lines_non_genes = [plt.Line2D([0], [0], color=dirs2color[d], lw=2) for d in dirs if d not in gene_groups]
    axes[1,1].legend(custom_lines_non_genes, non_gene_dir_names, loc='center', title='Group', frameon=False)
    axes[0,1].axis('off')
    axes[1,1].axis('off')

    plt.tight_layout()
    plt.savefig(outfile,dpi=300)
    plt.clf()

if __name__ == '__main__':
    args = get_args()
    plot_kurve(args.ppi_model_dir,args.title,args.outfile)