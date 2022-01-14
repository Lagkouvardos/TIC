import argparse
from tqdm import tqdm
from os import mkdir


def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


def check_valid_names(input_taxonomy):
    taxo_tokens = input_taxonomy.split(';')
    valid_names_list = list()
    phylum = ''
    class_s = ''
    order_o = ''
    family = ''
    genus = ''
    if len(taxo_tokens) >= 1:
        try:
            domain = taxo_tokens[0].split()[0]
        except BaseException:
            domain = ''
    if len(taxo_tokens) >= 2:
        try:
            phylum = taxo_tokens[1].split()[0]
        except BaseException:
            phylum = ''
    if len(taxo_tokens) >= 3:
        try:
            class_s = taxo_tokens[2].split()[0]
        except BaseException:
            class_s = ''
    if len(taxo_tokens) >= 4:
        try:
            order_o = taxo_tokens[3].split()[0]
        except BaseException:
            order_o = ''
    if len(taxo_tokens) >= 5:
        try:
            family = taxo_tokens[4].split()[0]
        except BaseException:
            family = ''
    if len(taxo_tokens) >= 6:
        try:
            genus = taxo_tokens[5].split()[0]
        except BaseException:
            genus = ''
    if domain:
        valid_names_list.append(domain)
    else:
        return('other.fasta')
    if phylum and not (phylum == 'uncultured' or phylum == 'unknown' or 'Incertae' in phylum):
        valid_names_list.append(phylum)
    else:
        return('_'.join(valid_names_list) + '.fasta')
    if class_s and not (class_s == 'uncultured' or class_s == 'unknown' or 'Incertae' in class_s):
        valid_names_list.append(class_s)
    else:
        return('_'.join(valid_names_list) + '.fasta')
    if order_o and not (order_o == 'uncultured' or order_o == 'unknown' or 'Incertae' in order_o):
        valid_names_list.append(order_o)
    else:
        return('_'.join(valid_names_list) + '.fasta')
    if family and not (family == 'uncultured' or family == 'unknown' or 'Incertae' in family):
        valid_names_list.append(family)
    else:
        return('_'.join(valid_names_list) + '.fasta')
    if genus and not (genus == 'uncultured' or genus == 'unknown' or 'Incertae' in genus):
        valid_names_list.append(genus)
    return('_'.join(valid_names_list) + '.fasta')


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data_dir", required=True, help="Data Directory")
parser.add_argument("-i", "--input", required=True, help="Input File")

args = parser.parse_args()
MAIN_DIR = args.data_dir + '/'
input_file = args.input

try:
    mkdir(MAIN_DIR)
except BaseException:
    print(MAIN_DIR + ' already present please select other directory or move the present directory to another location')
    exit()

input_contents = read_file(input_file)
for i in tqdm(range(0, len(input_contents), 2)):
    header = input_contents[i]
    sequence = input_contents[i+1]
    curr_taxonomy = header.split('tax=')[1]
    out_file_name = check_valid_names(curr_taxonomy).replace('(', '').replace(')', '')
    tokens = out_file_name.split('_')
    counter = 1
    new_out_file_name_tokens = list()
    for curr_token in tokens:
        if counter < 7:
            new_out_file_name_tokens.append(curr_token)
        counter += 1
    new_out_file_name = '_'.join(new_out_file_name_tokens)
    if '.fasta' not in new_out_file_name:
        new_out_file_name += '.fasta'
    out_file = open(MAIN_DIR + new_out_file_name, 'a+')
    new_header = header.split('tax=')[0] + 'tax=' + new_out_file_name.split('.')[0].replace('_', ';') + ';'
    out_file.write(new_header + '\n')
    out_file.write(sequence + '\n')
    out_file.close()
