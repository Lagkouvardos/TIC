from sys import argv
from os import system, remove, chdir, getcwd
from statistics import stdev, mean


CLUSTERING_TOOL = argv[1]
MAXDIFF = argv[2]
USER_FASTQ_FOLDER = argv[3]
TRIM_SCORE = argv[4]
MINMERGELEN = argv[5]
MAXMERGELEN = argv[6]
FORWARD_TRIM = argv[7]
REVERSE_TRIM = argv[8]
EXPECTED_ERROR_RATE = argv[9]
THREADS = argv[10]
MINPCTID = argv[11]


# return list with each line an element, clear of whitespace
def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


def fastq_merge_pairs(forward_file, reverse_file):
    if 'usearch' in CLUSTERING_TOOL:
        cmd = CLUSTERING_TOOL + ' -fastq_mergepairs ' + forward_file + ' -reverse ' + reverse_file + ' '
        cmd += '-fastq_maxdiffs ' + MAXDIFF + ' -fastq_pctid ' + MINPCTID + ' -fastqout ' + USER_FASTQ_FOLDER
        cmd += '/merged.fasta -fastq_trunctail ' + TRIM_SCORE + ' -fastq_minmergelen ' + MINMERGELEN
        cmd += ' -fastq_maxmergelen ' + MAXMERGELEN + ' -report ' + USER_FASTQ_FOLDER + '/report.txt -threads ' + THREADS
        cmd += '2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    else:
        cmd = CLUSTERING_TOOL + ' --fastq_mergepairs ' + forward_file + ' --reverse ' + reverse_file + ' '
        cmd += '-fastq_maxdiffs ' + MAXDIFF + ' -fastq_maxdiffpct ' + MINPCTID + ' -fastqout ' + USER_FASTQ_FOLDER
        cmd += '/merged.fasta -fastq_truncqual ' + TRIM_SCORE + ' -fastq_minmergelen ' + MINMERGELEN
        cmd += ' -fastq_maxmergelen ' + MAXMERGELEN + ' -threads ' + THREADS
        cmd += ' 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    system(cmd)


def trimming():
    if 'usearch' in CLUSTERING_TOOL:
        cmd = CLUSTERING_TOOL + " -fastx_truncate " + USER_FASTQ_FOLDER + '/merged.fasta -stripleft '
        cmd += FORWARD_TRIM + ' -stripright ' + REVERSE_TRIM + ' -fastqout ' + USER_FASTQ_FOLDER + '/trimmed.fastq -threads ' + THREADS
        cmd += ' 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    else:
        cmd = CLUSTERING_TOOL + " -fastx_filter " + USER_FASTQ_FOLDER + '/merged.fasta -fastq_stripleft '
        cmd += FORWARD_TRIM + ' -fastq_stripright ' + REVERSE_TRIM + ' -fastqout ' + USER_FASTQ_FOLDER + '/trimmed.fastq -threads ' + THREADS
        cmd += ' 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    system(cmd)


def trim_one_side(forward_file):
    if 'usearch' in CLUSTERING_TOOL:
        cmd = CLUSTERING_TOOL + " -fastx_truncate " + forward_file + " -stripleft " + FORWARD_TRIM
        cmd += " -fastqout " + USER_FASTQ_FOLDER + "/trimmed.fastq -threads " + THREADS
        cmd += ' 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    else:
        cmd = CLUSTERING_TOOL + " -fastx_filter " + forward_file + " -fastq_stripleft " + FORWARD_TRIM
        cmd += " -fastqout " + USER_FASTQ_FOLDER + "/trimmed.fastq -threads " + THREADS
        cmd += ' 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    system(cmd)


def filtering():
    cmd = CLUSTERING_TOOL + " -fastq_filter " + USER_FASTQ_FOLDER + '/trimmed.fastq --fastq_maxee_rate '
    cmd += EXPECTED_ERROR_RATE + ' -fastaout ' + USER_FASTQ_FOLDER + '/filtered.fasta -threads ' + THREADS
    cmd += ' 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    system(cmd)


def filter_merged_one_side(forward_file):
    curr_mean, curr_sd = seqFileStats(forward_file)
    minLength = curr_mean - int(0.1*curr_mean) - int(FORWARD_TRIM)  # remove the primer triming size plus 10% of the mean size
    cmd = CLUSTERING_TOOL + ' -fastq_filter trimmed.fastq -fastq_truncqual 20 -threads ' + THREADS
    cmd += ' -fastq_maxee_rate ' + EXPECTED_ERROR_RATE + ' -fastq_trunclen ' + str(minLength)
    cmd += ' -fastaout filtered.fasta '
    cmd += '2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    system(cmd)


def dereplication():
    if 'usearch' in CLUSTERING_TOOL:
        cmd = CLUSTERING_TOOL + " -fastx_uniques " + USER_FASTQ_FOLDER + '/filtered.fasta -sizeout '
        cmd += '-minuniquesize 1 -threads ' + THREADS + ' -fastaout ' + USER_FASTQ_FOLDER + '/unique.fasta '
        cmd += '2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    else:
        cmd = CLUSTERING_TOOL + " --derep_fulllength " + USER_FASTQ_FOLDER + '/filtered.fasta -sizeout '
        cmd += '-minuniquesize 1 -threads ' + THREADS + ' --output ' + USER_FASTQ_FOLDER + '/unique.fasta '
        cmd += '2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt'
    system(cmd)


def seqFileStats(seqFileName):
    try:
        seqfile_fh = read_file(seqFileName)
    except BaseException:
        print("Cannot open " + seqFileName + " to read from.")
        exit()
    lengths = list()
    for i in range(len(seqfile_fh)):
        if not seqfile_fh[i]:
            continue
        if seqfile_fh[i][0] == '@':
            curr_len = len(seqfile_fh[i+1])
            lengths.append(curr_len)
    mean_length = int(mean(lengths))
    st_dev = int(stdev(lengths))
    return(mean_length, st_dev)


def process_samples():
    top_directory = getcwd()
    chdir(USER_FASTQ_FOLDER)
    reverse_file = ''
    tab_file_contents = read_file('mapping_file.ssv')
    for line in tab_file_contents:
        tokens = line.split()
        sample_name = tokens[0]
        forward_file = tokens[2]
        if tokens[1] == '2':
            reverse_file = tokens[3]
            print('>>> Processing Sample:' + sample_name)
            out_file = open(top_directory + '/log_file.txt', 'a+')
            out_file.write('\n')
            out_file.write('>>> Processing Sample:' + sample_name + '\n')
            out_file.write('-----------------------------------\n')
            out_file.close()
            print('\tMerging...')
            try:
                fastq_merge_pairs(forward_file, reverse_file)
            except BaseException:
                print('Merging command failed with a critical failure. Exiting...\n')
                exit(1)
            print("\tMerging pairs " + forward_file + ' and ' + reverse_file + ' from sample ' + sample_name + '... Done')
            if 'usearch' in CLUSTERING_TOOL:
                system('cat ' + USER_FASTQ_FOLDER + '/report.txt >> ' + top_directory + '/log_file.txt')
                remove(USER_FASTQ_FOLDER + '/report.txt')
            print('\tTrimming...')
            try:
                trimming()
            except BaseException:
                print('Trimming command failed with a critical failure. Exiting...\n')
                exit(1)
            else:
                print('\tDone')
            print('\tFiltering merged reads...')
            try:
                filtering()
            except BaseException:
                print('Filtering command failed. Exiting...\n')
                exit(1)
            else:
                print('\tDone')
        elif tokens[1] == '1':
            print('>>> Processing Sample:' + sample_name)
            out_file = open(top_directory + '/log_file.txt', 'a+')
            out_file.write('\n')
            out_file.write('>>> Processing Sample:' + sample_name + '\n')
            out_file.write('-----------------------------------\n')
            out_file.close()
            print('\tTrimming...')
            try:
                trim_one_side(forward_file)
            except BaseException:
                print('Trimming command failed with a critical failure. Exiting...\n')
                exit(1)
            else:
                print('\tDone')
            print('\tFiltering...')
            try:
                filter_merged_one_side(forward_file)
            except BaseException:
                print('Filtering command failed. Exiting...\n')
                exit(1)
            else:
                print('\tDone')
        print('\tDereplicating sequences...')
        try:
            dereplication()
        except BaseException:
            print('Dereplication command failed. Exiting...\n')
            exit(1)
        else:
            print('\tDone')
        system('mv unique.fasta ' + sample_name + '_unique.fasta')
    try:
        remove('merged.fasta')
    except BaseException:
        pass
    try:
        remove('trimmed.fastq')
    except BaseException:
        pass
    try:
        remove('filtered.fasta')
    except BaseException:
        pass
    chdir(top_directory)


process_samples()
