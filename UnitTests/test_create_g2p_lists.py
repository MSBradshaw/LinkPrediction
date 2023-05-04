import os

test_el_1 = 'TestFiles/example_edgelist_1.tsv'
test_el_2 = 'TestFiles/example_edgelist_2.tsv'
test_el_1_map = 'TestFiles/example_edgelist_1.mapping.tsv'

def test_this():
    command="""
    python Scripts/create_g2p_lists.py --edgelist1 {} --edgelist2 {} --node_map1 {} --output {} --output_numbered {}
    """.format(test_el_1, test_el_2, test_el_1_map, 'TestFilesOutput/example_output.tsv', 'TestFiles/example_output_numbered.tsv')

    # run the command
    print('Running command: {}'.format(command))
    e = os.system(command)
    if e != 0:
        raise Exception('Command failed')


    # check the output
    # 'TestFilesOutput/example_output.tsv' should have 1 line, 'Mel\tHP:Scott'
    # this test that only new edges are added that are between prexisting nodes and that contain exactly one node with 'HP:'
    with open('TestFilesOutput/example_output.tsv') as f:
        for line in f:
            assert line.strip() == 'Mel\tHP:Scott'

    # check the numbered out put, it should have one line '2\t4'
    with open('TestFiles/example_output_numbered.tsv') as f:
        for line in f:
            assert line.strip() == '2\t4'

    print('Tests passed')

    # remove the example output files
    os.remove('TestFilesOutput/example_output.tsv')
    os.remove('TestFiles/example_output_numbered.tsv')

if __name__ == '__main__':
    test_this()