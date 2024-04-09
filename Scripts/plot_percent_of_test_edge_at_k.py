import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
import concurrent.futures
import random

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='Input dirrectory')
    parser.add_argument('-t', type=str, help='Input dirrectory')
    parser.add_argument('-o', type=str, help='Output image')
    return parser.parse_args()

def calc_hits(f, input_dir):
    result = {'k':[], 'percent_of_test_edge_at_k':[],'group':[]}
    df = pd.read_csv(input_dir + f,sep='\t')
    group = f.split('.')[0]

    df['global_rank'] = df['score'].rank(pct=True)
    n_test_edges = df['in_test_set'].sum()
    if n_test_edges == 0:
        for k in [1,5,10,15,20,25,50,100] + list(range(100,1001,100)):
            result['k'].append(k)
            result['percent_of_test_edge_at_k'].append(0)
            result['group'].append(group)
        print('Warning no test edges in', f)
        return pd.DataFrame(result)

    for k in [1,5,10,15,20,25,50,100] + list(range(100,1001,100)):
        hits_across_qs_at_k = 0
        for q in df['query_term'].unique():
            df_q = df[df['query_term'] == q]
            df_q = df_q.sort_values(by='score', ascending=False)
            hits_across_qs_at_k += df_q.iloc[:k]['in_test_set'].sum()
        percent =  hits_across_qs_at_k / n_test_edges
        result['k'].append(k)
        
        result['percent_of_test_edge_at_k'].append(int(percent * 100))
        result['group'].append(group)

    return pd.DataFrame(result)


def main():
    args = get_args()
    input_dir = args.i
    output_image = args.o
    if input_dir[-1] != '/':
        input_dir += '/'
    
    dirs2color = {'Female':'#b66dff','Male':'#006ddb','Cancer':'#920000','PedCancer':'#ff6db6','RareDisease':'#004949', 'UltraRareDisease':'#009999','African':'#db6d00','EastAsian':'#22cf22','European':'#8f4e00','Latino':'#490092','Random':'#676767','RandomDiseases':'#616161'}
    dirs = ['Female','Male','Cancer','PedCancer','RareDisease', 'UltraRareDisease','African','EastAsian','European','Latino','Random','RandomDiseases']
    dirs_names = ['Female','Male','Cancer','Pediatric Cancer','Rare Disease', 'Ultra Rare Disease','African','East Asian','European','Latino','Random Genes','Random Diseases']

    gene_dir_names = ['Female','Male','Cancer','Pediatric Cancer','Random']
    gene_groups = ['Female','Male','Cancer','PedCancer','Random']
    non_gene_dir_names = ['Rare Disease', 'Ultra Rare Disease','African','East Asian','European','Latino']

    files = [f for f in os.listdir(input_dir) if f.endswith('.tsv')]

    result = {'k':[], 'percent_of_test_edge_at_k':[],'group':[]}
    
    for i,f in enumerate(files):
        pass

    # parallize calc_hits for every file in files
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        # Submit the tasks to the executor
        # Each task corresponds to calling my_function with a value from values_list
        future_results = {executor.submit(calc_hits, value, input_dir): value for value in files}
        
        res_df = None
        # Retrieve the results as they become available
        for future in concurrent.futures.as_completed(future_results):
            value = future_results[future]
            result = future.result()
            print('Got one!', value, result.shape)
            if res_df is None:
                res_df = result
            else:
                res_df = pd.concat([res_df, result], axis=0)

    fig, ax = plt.subplots(2)
    # set fig size
    fig.set_size_inches(6, 10)
    for group in res_df['group'].unique():
        group_df = res_df[res_df['group'] == group]
        if group in gene_groups:
            index = 0
        else:
            index = 1
        ax[index].plot(group_df['k'], group_df['percent_of_test_edge_at_k'], label=group, color=dirs2color[group])

    for i in range(0,2):
        ax[i].set_xlabel('k')
        ax[i].set_ylabel('Percent of test edges at k')
        # remove top right borders
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        # legend with no border
        ax[i].legend(frameon=False)
    ax[0].set_title(args.t)
    
    plt.tight_layout()
    plt.savefig(output_image,dpi=300)

# def main_og():
#     args = get_args()
#     input_dir = args.i
#     output_image = args.o
#     if input_dir[-1] != '/':
#         input_dir += '/'
    
#     dirs2color = {'Female':'#b66dff','Male':'#006ddb','Cancer':'#920000','PedCancer':'#ff6db6','RareDisease':'#004949', 'UltraRareDisease':'#009999','African':'#db6d00','EastAsian':'#22cf22','European':'#8f4e00','Latino':'#490092','Random':'#676767','RandomDiseases':'#e6e6e6'}
#     dirs = ['Female','Male','Cancer','PedCancer','RareDisease', 'UltraRareDisease','African','EastAsian','European','Latino','Random','RandomDiseases']
#     dirs_names = ['Female','Male','Cancer','Pediatric Cancer','Rare Disease', 'Ultra Rare Disease','African','East Asian','European','Latino','Random Genes','Random Diseases']

#     gene_dir_names = ['Female','Male','Cancer','Pediatric Cancer','Random']
#     gene_groups = ['Female','Male','Cancer','PedCancer','Random']
#     non_gene_dir_names = ['Rare Disease', 'Ultra Rare Disease','African','East Asian','European','Latino']

#     files = [f for f in os.listdir(input_dir) if f.endswith('.tsv')]

#     result = {'k':[], 'percent_of_test_edge_at_k':[],'group':[]}
#     for i,f in enumerate(files):
#         print(f)
#         df = pd.read_csv(input_dir + f,sep='\t')
#         group = f.split('.')[0]

#         df['global_rank'] = df['score'].rank(pct=True)
#         n_test_edges = df['in_test_set'].sum()
#         if n_test_edges == 0:
#             for k in [1,5,10,15,20,25,50,100] + list(range(100,1001,100)):
#                 result['k'].append(k)
#                 result['percent_of_test_edge_at_k'].append(0)
#                 result['group'].append(group)
#             print('Warning no test edges in', f)
#             continue

#         for k in [1,5,10,15,20,25,50,100] + list(range(100,1001,100)):
#             hits = df.iloc[:k]['in_test_set'].sum()
#             hits_across_qs_at_k = 0
#             for q in df['query_term'].unique():
#                 df_q = df[df['query_term'] == q]
#                 df_q = df_q.sort_values(by='score', ascending=False)
#                 hits_across_qs_at_k += df_q.iloc[:k]['in_test_set'].sum()
#             percent =  hits_across_qs_at_k / n_test_edges
#             result['k'].append(k)
            
#             result['percent_of_test_edge_at_k'].append(int(percent * 100))
#             result['group'].append(group)

#     res_df = pd.DataFrame(result)
#     res_df.head()

#     fig, ax = plt.subplots(2)
#     # set fig size
#     fig.set_size_inches(6, 10)
#     for group in res_df['group'].unique():
#         group_df = res_df[res_df['group'] == group]
#         if group in gene_groups:
#             index = 0
#         else:
#             index = 1
#         ax[index].plot(group_df['k'], group_df['percent_of_test_edge_at_k'], label=group, color=dirs2color[group])

#     for i in range(0,2):
#         ax[i].set_xlabel('k')
#         ax[i].set_ylabel('Percent of test edges at k')
#         # remove top right borders
#         ax[i].spines['top'].set_visible(False)
#         ax[i].spines['right'].set_visible(False)
#         # legend with no border
#         ax[i].legend(frameon=False)
#     ax[0].set_title(args.t)
    
#     plt.tight_layout()
#     plt.savefig(output_image,dpi=300)

if __name__ == '__main__':
    main()