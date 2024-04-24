import argparse

def get_args():
        parser = argparse.ArgumentParser(description='Count duplicate edges in a network')
        parser.add_argument('--el1', type=str, help='first edge list')
        parser.add_argument('--el2', type=str, help='second edge list')
        parser.add_argument('--triples', action='store_true', help='flag for triples', default=False)
        parser.add_argument('--non_directed', action='store_true', help='flag for directed', default=False)
        return parser.parse_args()

def read_el(filepath,triples,non_directed):
        el = []
        non_double_count = 0
        with open(filepath, 'r') as f:
                for line in f:
                        non_double_count+=1
                        if triples:
                                h,p,t= line.strip().split()
                                el.append(str((h,p,t)))
                                if not non_directed:
                                      el.append(str((t,p,h)))  
                        else:
                                h,t = line.strip().split()
                                el.append(str((h,t)))
                                if not non_directed:
                                      el.append(str((t,h)))
        return el, non_double_count

def main():
        # read in the edge lists
        args = get_args()
        el1, c1 = read_el(args.el1,args.triples,args.non_directed)
        el2, c2 = read_el(args.el2,args.triples,args.non_directed)
        # print number of edges in each list, report if double counting because of non-dirrected
        if args.non_directed:
                print('Double counting edges because of non-directed edges')
                print('\tnumber of non-twice counted edges in el1:',c1)
                print('\tnumber of non-twice counted edges in el2:',c2)
        print('number of edges in el1:',len(el1))
        print('number of edges in el2:',len(el2))
        # count uniq edges in each EL
        uniq1 = set(el1)
        uniq2 = set(el2)
        print('number of unique edges in el1:',len(uniq1))
        print('number of unique edges in el2:',len(uniq2))
        # count duplicates
        duplicates = len(uniq1.intersection(uniq2))
        print('Intersection of el1 and el2:', duplicates)

if __name__ == '__main__':
        main()
"""
python Scripts/count_duplicate_edges.py --el1 ELs_for_Rotate/Monarch_HuRI/test.txt --el2 ELs_for_Rotate/Monarch_HuRI/train.txt --non_directed --triples
python Scripts/count_duplicate_edges.py --el1 ELs_for_Rotate/Monarch_HuRI_Filtered/test.txt --el2 ELs_for_Rotate/Monarch_HuRI_Filtered/train.txt --non_directed --triples
"""