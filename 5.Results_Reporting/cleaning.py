from sys import argv
from os import system

OUTPUT_FASTA_EXTRACTION = argv[1]
CLUSTERING_DIRECTORY = argv[2]
OUTPUT_ZOTUS_EXTRACTION = argv[3]
USER_FASTQ_FOLDER = argv[4]

file_list = ['3.Extraction-Denoising/good_ZOTUS.log',
             '3.Extraction-Denoising/other_ZOTUS.fa',
             '3.Extraction-Denoising/ZOTUs-Table.tab',
             '2.Taxonomy-Classification/alignment_vector.csv',
             OUTPUT_FASTA_EXTRACTION, CLUSTERING_DIRECTORY, OUTPUT_ZOTUS_EXTRACTION,
             USER_FASTQ_FOLDER + '/merged.fasta',
             'test_s.fasta', 'test_z.fasta', 'sina_output.csv',
             'sina_output.fasta']

for i in file_list:
    try:
        system('rm -r ' + i)
    except BaseException:
        continue
