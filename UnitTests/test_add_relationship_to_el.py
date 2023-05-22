import os
import pandas as pd

test_el_1 = 'TestFiles/STRING_HPO_2019_small_test.edgelist.txt'
test_output = 'TestFilesOutput/example_output_1.tsv'


def test_this():
    command="""
    python Scripts/add_relationship_to_el.py --input {} --output {}
    """.format(test_el_1,test_output)

    # run the command
    print('Running command: {}'.format(command))
    e = os.system(command)
    if e != 0:
        raise Exception('Command failed')

    # check that the output file has the same number of lines as the input one
    f1 = pd.read_csv(test_el_1, sep='\t', header=None)
    f2 = pd.read_csv(test_output, sep='\t', header=None)
    assert f1.shape[0] == f2.shape[0]
    assert f2.shape[1] == 3

    # assert that there are 10 instances of 'STRING2STRING' and 10 instances of 'HPO2HPO' and 10 instances of 'STRING2HPO' in f2
    assert f2[1].value_counts()['STRING2STRING'] == 10
    assert f2[1].value_counts()['HPO2HPO'] == 10
    assert f2[1].value_counts()['STRING2HPO'] == 10

    # assert there are no a column 1 value that is not 'STRING2STRING' or 'HPO2HPO' or 'STRING2HPO'
    assert len(f2[1].unique()) == 3

    # assert that every node has a prefix 'STRING:' or 'HP:'
    assert f2[0].apply(lambda x: x.startswith('STRING:') or x.startswith('HP:')).all()
    assert f2[2].apply(lambda x: x.startswith('STRING:') or x.startswith('HP:')).all()

    print('Tests passed')

    # remove the example output files
    os.remove(test_output)
    
if __name__ == '__main__':
    test_this()