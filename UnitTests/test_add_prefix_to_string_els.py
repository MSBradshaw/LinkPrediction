import os

test_el_1 = 'TestFiles/example_edgelist_1.tsv'
test_el_2 = 'TestFiles/example_edgelist_2.tsv'
test_el_1_map = 'TestFiles/example_edgelist_1.mapping.tsv'

def test_this():
    command="""
    python Scripts/add_prefix_to_string_els.py --edgelist1 {} --edgelist2 {} --node_map1 {} --node_map2 {} --output_el1 {} --output_el2 {} --output_mapping1 {} --output_mapping2 {}
    """.format(test_el_1,test_el_2,test_el_1_map,test_el_1_map,
            'TestFilesOutput/example_output_1.tsv',
            'TestFilesOutput/example_output_2.tsv',
            'TestFilesOutput/example_output_1.mapping.tsv',
            'TestFilesOutput/example_output_2.mapping.tsv')

    # run the command
    print('Running command: {}'.format(command))
    e = os.system(command)
    if e != 0:
        raise Exception('Command failed')

    # check that the input files and output files have the same number of lines
    def same_num_file(f1,f2):
        with open(f1) as f1:
            with open(f2) as f2:
                return len(f1.readlines()) == len(f2.readlines())

    assert same_num_file(test_el_1,'TestFilesOutput/example_output_1.tsv')
    assert same_num_file(test_el_2,'TestFilesOutput/example_output_2.tsv')
    assert same_num_file(test_el_1_map,'TestFilesOutput/example_output_1.mapping.tsv')

    # check that each line has 2 HP: nodes or 2 STRING: nodes or one of each
    def check_file(f):
        with open(f) as f:
            for line in f:
                line = line.strip()
                assert line.count('HP:') == 2 or line.count('STRING:') == 2 or (line.count('HP:') == 1 and line.count('STRING:') == 1)

    check_file('TestFilesOutput/example_output_1.tsv')
    check_file('TestFilesOutput/example_output_2.tsv')

    print('Tests passed')

    # remove the example output files
    os.remove('TestFilesOutput/example_output_1.tsv')
    os.remove('TestFilesOutput/example_output_2.tsv')
    os.remove('TestFilesOutput/example_output_1.mapping.tsv')

if __name__ == '__main__':
    test_this()