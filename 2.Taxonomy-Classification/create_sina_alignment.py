from os import system
from sys import argv
from glob import glob

USER_FASTQ_FOLDER = argv[1]
CLUSTERING_TOOL = argv[2]
THREADS = argv[3]

SILVA_ARB = argv[4]
SINA_EXECUTABLE = argv[5]
PDF_REGION_OUTPUT = argv[6]
OUTPUT_FASTA_ALI_CLASS = argv[7]


def merging_all_samples():
    all_uniques = glob(USER_FASTQ_FOLDER + '/*_unique.fasta')
    output_file = open(USER_FASTQ_FOLDER + '/merged.fasta', 'w+')
    for curr_sample in all_uniques:
        sample_name = curr_sample.split('/')[-1].split('_')[0]
        with open(curr_sample) as fp:
            while True:
                line = fp.readline()
                if not line:
                    break
                if line[0] == '>':
                    output_file.write(line[:-1] + ';barcodelabel=' + sample_name + '\n')
                else:
                    output_file.write(line)
    output_file.close()


merging_all_samples()

cmd = SINA_EXECUTABLE + ' --quiet --in=' + USER_FASTQ_FOLDER + '/merged.fasta --out=' + OUTPUT_FASTA_ALI_CLASS + ' --db='
cmd += SILVA_ARB + ' --turn all --search --meta-fmt csv --lca-fields=tax_slv  --fasta-write-dna --threads ' + str(THREADS)
system(cmd)
print('\tDone')

cmd = 'python3 2.Taxonomy-Classification/create_alignment_vector.py '
cmd += '-i ' + OUTPUT_FASTA_ALI_CLASS + ' -o 2.Taxonomy-Classification'
system(cmd)

print('>>> Plotting PDF of 50K aligned positions...')
cmd = 'Rscript 2.Taxonomy-Classification/ggplot_alignment_vector.R 2.Taxonomy-Classification ' + PDF_REGION_OUTPUT
system(cmd)
print('\tDone')
