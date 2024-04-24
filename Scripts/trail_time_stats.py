import pandas as pd
import sys
from datetime import datetime

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

# read in sys.argv[1] as a DataFrame
df = pd.read_csv(sys.argv[1],sep='\t')

# print(df['duration'])

# convert time stamp like this: 0 days 00:14:46.993231 to seconds
df['duration2'] = df['duration'].apply(lambda x: time_string_to_seconds(x))
# print(df['duration2'])

# calc the mean and median and print them
mean = df['duration2'].mean()
median = df['duration2'].median()
print("Mean:", mean)
print("Median:", median)


"""
python del.py PyKeenOut/complex_monarch_kg_filtered_1/trials.tsv
python del.py PyKeenOut/rotate_monarch_kg_filtered_1/trials.tsv
python del.py PyKeenOut/transe_monarch_kg_filtered_1/trials.tsv

python del.py PyKeenOut/ComplEx_monarch_kg_1/trials.tsv
python del.py PyKeenOut/rotate_monarch_kg_1/trials.tsv
python del.py PyKeenOut/TransE_monarch_kg_1/trials.tsv
"""