from sys import argv


# return list with each line an element, clear of whitespace
def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


def count_seq():
    out_file = open(OUTPUT_FASTA_EXTRACTION, 'w+')
    with open(INPUT_FASTA_EXTRACTION) as file_one:
        for line in file_one:
            if line[0] == '>':
                curr_header = line
            else:
                trimmed_sequence = line[EXTRACTION_REGION_START: EXTRACTION_REGION_END]
                collapsed_sequence = trimmed_sequence.replace('-', '').replace('.', '')
                if len(collapsed_sequence) >= BASES_LOW_LIMIT:
                    out_file.write(curr_header + '\n')
                    out_file.write(collapsed_sequence + '\n')
    out_file.close()


print('>>> Extracting regions based on user input...')
INPUT_FASTA_EXTRACTION = argv[1]
EXTRACTION_REGION_START = int(argv[2])
EXTRACTION_REGION_END = int(argv[3])
EXTRACTION_REGION_LIMIT = argv[4]

if 'AUTOMATIC' == EXTRACTION_REGION_LIMIT:
    ecoli = read_file('3.Extraction-Denoising/ecoli_sina.ali')[1]
    cut_region = ecoli[EXTRACTION_REGION_START:EXTRACTION_REGION_END].replace('.', '').replace('-', '')
    EXTRACTION_REGION_LIMIT = len(cut_region)
else:
    EXTRACTION_REGION_LIMIT = int(EXTRACTION_REGION_LIMIT)

BASES_LOW_LIMIT = EXTRACTION_REGION_LIMIT * 0.8
OUTPUT_FASTA_EXTRACTION = argv[5]

count_seq()

print('\tDone')
