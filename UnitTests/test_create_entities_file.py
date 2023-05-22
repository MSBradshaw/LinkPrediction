import os
import pandas as pd

test_el_1 = 'TestFiles/STRING_HPO_2019_small_test.edgelist.txt'
test_output = 'TestFilesOutput/example_output_1.tsv'


def ss(line):
    return line.strip().split('\t')

def test_this():
    command="""
    python Scripts/create_entities_file.py --train {} --validate {} --test {} --mapping {} --output {}
    """.format('TestFiles/STRING_HPO_2019_small_test.triples.txt',
            'TestFiles/STRING_HPO_2019_small_test.triples.1.txt',
            'TestFiles/STRING_HPO_2019_small_test.triples.2.txt',
            'TestFiles/STRING_HPO_2019_small_test.triples.mapping.tsv',
            'TestFilesOutput/mapping.tsv')

    # run the command
    print('Running command: {}'.format(command))
    e = os.system(command)
    if e != 0:
        raise Exception('Command failed')

    # load all input files as dfs
    train = pd.read_csv('TestFiles/STRING_HPO_2019_small_test.triples.txt', sep='\t', header=None)
    validate = pd.read_csv('TestFiles/STRING_HPO_2019_small_test.triples.1.txt', sep='\t', header=None)
    test = pd.read_csv('TestFiles/STRING_HPO_2019_small_test.triples.2.txt', sep='\t', header=None)

    # load the mapping file as a dictionary
    mapping = { ss(line)[1]:int(ss(line)[0]) for line in open('TestFilesOutput/mapping.tsv', 'r')}
    
    # check that every node in training, valid, test is in the mapping
    assert train[0].apply(lambda x: x in mapping).all()
    assert train[2].apply(lambda x: x in mapping).all()
    assert validate[0].apply(lambda x: x in mapping).all()
    assert validate[2].apply(lambda x: x in mapping).all()
    assert test[0].apply(lambda x: x in mapping).all()
    assert test[2].apply(lambda x: x in mapping).all()

    # assert that every value in mapping is unique
    assert len(mapping.values()) == len(set(mapping.values()))

    print('Tests passed')

    # remove the example output files
    os.remove('TestFilesOutput/mapping.tsv')
    
if __name__ == '__main__':
    test_this()