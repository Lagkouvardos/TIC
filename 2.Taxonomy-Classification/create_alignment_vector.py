import argparse


# return list with each line an element, clear of whitespace
def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


def find_differences(s):
    alphabet = ['A', 'G', 'T', 'C', 'a', 'g', 't', 'c']
    positions = list()
    for i, letter in enumerate(s):
        if letter in alphabet:
            positions.append(i)
    return positions


def create_vector(input_file, out_dir):
    position_vector_file = open(out_dir + "/alignment_vector.csv", "w+")
    position_vector = [0] * 50000
    filepath = input_file
    with open(filepath) as fp:
        line = fp.readline()
        while line:
            if(">" not in line):
                align_pos = find_differences(line)
                for j in align_pos:
                    position_vector[j] = position_vector[j] + 1
            line = fp.readline()
    ready_vec = ','.join(str(e) for e in position_vector)
    position_vector_file.write(ready_vec + '\n')
    position_vector_file.close()


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", required=True, help="Input FASTA file", type=str)
parser.add_argument("-o", "--out_dir", required=True, help="Output DIR", type=str)
args = parser.parse_args()

input_file = str(args.input_file)
out_dir = str(args.out_dir)

create_vector(input_file, out_dir)
