from tqdm import tqdm
from sys import argv
from os import system


def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


OUTPUT_FASTA_EXTRACTION = argv[1]
OUTPUT_CSV_ALI_CLASS = argv[2]

sina_dict = dict()
sina_csv = read_file(OUTPUT_CSV_ALI_CLASS)[1:]
print('Incorporating Taxonomy Information')
for i in tqdm(range(len(sina_csv))):
    curr_line = sina_csv[i]
    line_tokens = curr_line.split(',')
    original_header = line_tokens[0]
    new_taxonomy = line_tokens[6]
    sina_dict[original_header] = new_taxonomy


filepath = OUTPUT_FASTA_EXTRACTION
out_fasta = open('dataset_extracted_new_taxonomy.fasta', 'w+')
with open(filepath) as fp:
    line = fp.readline()
    while line:
        if line[0] == '>':
            new_header = line[:-1] + ';tax=' + sina_dict[line[1:-1]].replace(' ', '') + '\n'
            out_fasta.write(new_header)
        else:
            out_fasta.write(line)
        line = fp.readline()

out_fasta.close()

system('mv dataset_extracted_new_taxonomy.fasta ' + OUTPUT_FASTA_EXTRACTION)
