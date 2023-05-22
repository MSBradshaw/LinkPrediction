import argparse

# get params function
# params: input file, output file
def get_params():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, help='input edgelist file')
    parser.add_argument('--output', '-o', type=str, help='output triple list file')
    args = parser.parse_args()
    return args

def main():
    args = get_params()
    with open(args.output, 'w') as out:
        for line in open(args.input, 'r'):
            if line.startswith('#'):
                continue
            # if the line is empty
            if line == '\n':
                continue
            row = line.strip().split('\t')
            # count the number of occurances of the 'STRING:' substring in the line
            string_count = line.count('STRING:')
            hpo_count = line.count('HP:')
            if string_count == 2:
                out.write('{}\t{}\t{}\n'.format(row[0], 'STRING2STRING', row[1]))
            elif hpo_count == 2:
                out.write('{}\t{}\t{}\n'.format(row[0], 'HPO2HPO', row[1]))
            elif string_count == 1 and hpo_count == 1:
                out.write('{}\t{}\t{}\n'.format(row[0], 'STRING2HPO', row[1]))
            else:
                print('ERROR: line does not match any of the patterns')
                print(line)
                exit(1)
            


if __name__ == '__main__':
    main()