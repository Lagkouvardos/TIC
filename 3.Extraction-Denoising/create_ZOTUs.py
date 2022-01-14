from os import system, remove, chdir, getcwd
from sys import argv

USER_FASTQ_FOLDER = argv[1]
CLUSTERING_TOOL = argv[2]
THREADS = argv[3]
MIN_ZOTU_SIZE = argv[4]
SORT_ME_RNA_DB1 = argv[5]
SORT_ME_RNA_DB2 = argv[6]
SORT_ME_RNA_TOOL = argv[7]
RAPID_NJ = argv[8]
OUTPUT_FASTA_EXTRACTION = argv[9]
OUTPUT_ZOTUS_EXTRACTION = argv[10]


def dereplication_merged():
    if 'usearch' in CLUSTERING_TOOL:
        cmd = CLUSTERING_TOOL + " -fastx_uniques " + OUTPUT_FASTA_EXTRACTION + ' -sizeout '
        cmd += '-sizein -threads ' + THREADS + ' -fastaout ' + USER_FASTQ_FOLDER + '/dereped.fasta '
        cmd += '2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    else:
        cmd = CLUSTERING_TOOL + " --derep_fulllength " + OUTPUT_FASTA_EXTRACTION + ' -sizeout '
        cmd += '-sizein --output ' + USER_FASTQ_FOLDER + '/dereped.fasta -threads ' + THREADS
        cmd += ' 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    system(cmd)


def sort_merged():
    if 'usearch' in CLUSTERING_TOOL:
        cmd = CLUSTERING_TOOL + " -sortbysize " + USER_FASTQ_FOLDER + '/dereped.fasta'
        cmd += ' -fastaout ' + USER_FASTQ_FOLDER + '/sorted.fasta -threads ' + THREADS
        cmd += ' 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    else:
        cmd = CLUSTERING_TOOL + " -sortbysize " + USER_FASTQ_FOLDER + '/dereped.fasta'
        cmd += ' -output ' + USER_FASTQ_FOLDER + '/sorted.fasta -threads ' + THREADS
        cmd += ' 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    system(cmd)


def unoise():
    if 'usearch' in CLUSTERING_TOOL:
        cmd = CLUSTERING_TOOL + " -unoise3 " + USER_FASTQ_FOLDER + '/sorted.fasta -minsize ' + MIN_ZOTU_SIZE
        cmd += ' -zotus 3.Extraction-Denoising/zotus.fasta -threads ' + THREADS
        cmd += ' 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    else:
        cmd = CLUSTERING_TOOL + " -cluster_unoise " + USER_FASTQ_FOLDER + '/sorted.fasta -minsize ' + MIN_ZOTU_SIZE
        cmd += ' -centroids 3.Extraction-Denoising/zotus.fasta -threads ' + THREADS
        cmd += ' 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    system(cmd)


def remove_chimeras():
    top_directory = getcwd()
    chdir('3.Extraction-Denoising/')
    print(">>> Filtering out non 16S ZOTUs... ")
    cmd = SORT_ME_RNA_TOOL + ' --ref ' + SORT_ME_RNA_DB1 + ' --ref ' + SORT_ME_RNA_DB2
    cmd += ' --reads zotus.fasta --fastx --aligned good_ZOTUS --other other_ZOTUS --threads ' + THREADS
    cmd += ' --workdir sortme -e 0.1 >/dev/null 2>&1'
    try:
        system(cmd)
    except BaseException:
        print("The command for removal of non 16S seqs failed\n")
    else:
        system("rm -rf sortme/kvdb")
        remove('zotus.fasta')
        print("\tDone")
    chdir(top_directory)


def create_zotu_table():
    print('>>> Creating ZOTU table...')
    top_directory = getcwd()
    chdir('3.Extraction-Denoising/')
    cmd = top_directory + "/0.Setup_and_Testing/usearch -otutab " + USER_FASTQ_FOLDER + '/merged.fasta -zotus good_ZOTUS.fa'
    cmd += " -otutabout ZOTUs-Table.tab -id 0.97  -threads " + THREADS
    cmd += ' 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    try:
        system(cmd)
    except BaseException:
        print("ZOTU table formation command failed\n")
    else:
        print("\tDone")
    chdir(top_directory)


def create_ASVs():
    dereplication_merged()
    sort_merged()
    unoise()
    remove_chimeras()
    create_zotu_table()
    remove(USER_FASTQ_FOLDER + '/dereped.fasta')
    remove(USER_FASTQ_FOLDER + '/sorted.fasta')
    system('mv 3.Extraction-Denoising/good_ZOTUS.fa ' + OUTPUT_ZOTUS_EXTRACTION)


create_ASVs()
