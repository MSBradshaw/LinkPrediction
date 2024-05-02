import argparse
import pandas as pd
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description='Plot time trial stats')
    # input files, take many
    parser.add_argument('-i',dest='input_files', type=str, nargs='+', help='input files')
    parser.add_argument('--input_names',dest='input_names', type=str, nargs='+', help='input names')
    # output file
    parser.add_argument('--output',dest='output', type=str, help='output file')
    return parser.parse_args()

def time_string_to_seconds(time_str):
    # check that time_str is a string is not skip
    if not isinstance(time_str, str):
        return None 
    days = int(time_str.split(" days ")[0]) * 24 * 60 * 60
    time_str = time_str.split(" days ")[1].split(":")
    hours = int(time_str[0]) * 60 * 60
    minutes = int(time_str[1]) * 60
    seconds = float(time_str[2].split(".")[0])

    total_seconds = days + hours + minutes + seconds
    
    return total_seconds

def load(input_file):
    df = pd.read_csv(input_file,sep='\t')
    # convert time stamp like this: 0 days 00:14:46.993231 to seconds
    df['duration2'] = df['duration'].apply(lambda x: time_string_to_seconds(x))
    # convert duration 2 to hours
    df['duration2'] = df['duration2'] / 3600
    return df

def main():
    args = get_args()
    # replace the first space in each label with a new line character
    for i, name in enumerate(args.input_names):
        args.input_names[i] = name.replace(" ", "\n")
    dfs = []
    print(args.input_files)
    print(args.input_names)
    for i, input_file in enumerate(args.input_files):
        dfs.append(load(input_file))
        print("Loaded", args.input_names[i] , input_file, dfs[-1].shape, 'Median=', dfs[-1]['duration2'].median())
        # drop NaNs
        dfs[-1] = dfs[-1].dropna(subset=['duration2'])
        

    # plot a box plot for each input file all on the same axis, with the names as the labels
    fig, ax = plt.subplots()
    ax.boxplot([df['duration2'] for df in dfs], labels=args.input_names)
    ax.set_ylabel('Trial duration (hours)')
    # remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xticks(rotation=45)
    plt.tight_layout()
    # log y
    # plt.yscale('log')
    plt.savefig(args.output)


if __name__ == '__main__':
    main()

"""
python Scripts/plot_time_trial_stats.py  \
-i PyKeenOut/ComplEx_monarch_kg_1/trials.tsv PyKeenOut/rotate_monarch_kg_1/trials.tsv PyKeenOut/TransE_monarch_kg_1/trials.tsv PyKeenOut/complex_monarch_kg_filtered_1/trials.tsv PyKeenOut/rotate_monarch_kg_filtered_1/trials.tsv PyKeenOut/transe_monarch_kg_filtered_1/trials.tsv \
--input_names "ComplEx Monarch" "Rotate Monarch" "TransE Monarch" "ComplEx Monarch Filtered" "Rotate Monarch Filtered" "TransE Monarch Filtered" \
--output Figures/time_trial_stats.png
"""