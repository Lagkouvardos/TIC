from sys import argv
from os import mkdir, system, getcwd
from glob import glob


def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


def fasta2dict(fil):
    dic = {}
    cur_scaf = ''
    cur_seq = []
    for line in open(fil):
        if line.startswith(">") and cur_scaf == '':
            cur_scaf = line[1:].rstrip()
        elif line.startswith(">") and cur_scaf != '':
            dic[cur_scaf] = ''.join(cur_seq)
            cur_scaf = line[1:].rstrip()
            cur_seq = []
        else:
            cur_seq.append(line.rstrip())
    dic[cur_scaf] = ''.join(cur_seq)
    new_dict = dict()
    for key, value in dic.items():
        new_key = key.split(';')[0]
        new_dict[new_key] = value
    return new_dict


OUTPUT_FOLDER = argv[1]
OUTPUT_ASV_FASTA_WITH_TAXONOMY = argv[2]
OUTPUT_ASV_TABLE = argv[3]
CLUSTERING_DIRECTORY = argv[4]
INPUT_FASTA_CLUSTERING = argv[5]
KRONA_TOOL = argv[6]
OUTPUT_SOTU_FASTA_WITH_TAXONOMY = argv[7]
SILVA_ARB = argv[8]
SINA_EXECUTABLE = argv[9]
USER_FASTQ_FOLDER = argv[10]

mkdir(OUTPUT_FOLDER)
zotus_seqs_dict = fasta2dict(INPUT_FASTA_CLUSTERING)
all_stats = glob(CLUSTERING_DIRECTORY + '/species_stats/**/*.stats', recursive=True)
output_fasta = open(OUTPUT_FOLDER + '/' + OUTPUT_ASV_FASTA_WITH_TAXONOMY, 'w+')
output_fasta_sotu_centroids = open(OUTPUT_FOLDER + '/' + OUTPUT_SOTU_FASTA_WITH_TAXONOMY, 'w+')
for curr_stats in all_stats:
    stats_contents = read_file(curr_stats)
    for line in stats_contents:
        if not line:
            continue
        zotu_with_taxonomy = line.split('\t')[1]
        clean_zotu_name = zotu_with_taxonomy.split(';')[0]
        output_fasta.write('>' + zotu_with_taxonomy + '\n')
        output_fasta.write(zotus_seqs_dict[clean_zotu_name] + '\n')
        status = line.split('\t')[0]
        if status == 'S':
            output_fasta_sotu_centroids.write('>' + zotu_with_taxonomy + '\n')
            output_fasta_sotu_centroids.write(zotus_seqs_dict[clean_zotu_name] + '\n')
output_fasta_sotu_centroids.close()
output_fasta.close()

zotus_with_taxonomy_contents = read_file(OUTPUT_FOLDER + '/' + OUTPUT_ASV_FASTA_WITH_TAXONOMY)
zotus_taxonomy_dict = dict()
for line in zotus_with_taxonomy_contents:
    if line[0] == '>':
        tokens = line.split('tax=')
        zotu_name = tokens[0].split(';')[0][1:]
        taxonomy = tokens[1]
        zotus_taxonomy_dict[zotu_name] = taxonomy


out_tab = open(OUTPUT_FOLDER + '/' + OUTPUT_ASV_TABLE, 'w+')
tab_contents = read_file('3.Extraction-Denoising/ZOTUs-Table.tab')
for line in tab_contents:
    if line[0] == '#':
        out_tab.write(line + '\tTaxonomy\n')
    clean_zotu_name = line.split('\t')[0]
    if clean_zotu_name in zotus_taxonomy_dict.keys():
        new_line = line + '\t' + zotus_taxonomy_dict[clean_zotu_name] + '\n'
        out_tab.write(new_line)
out_tab.close()

taxonomy_counters_dict = dict()
for value in zotus_taxonomy_dict.values():
    if value not in taxonomy_counters_dict.keys():
        taxonomy_counters_dict[value] = 1
    else:
        taxonomy_counters_dict[value] += 1

out_file = open(OUTPUT_FOLDER + '/for_krona.tab', 'w+')
for key, value in taxonomy_counters_dict.items():
    out_file.write(str(value) + '\t' + key.replace(';', '\t') + '\n')
out_file.close()

system(KRONA_TOOL + ' ' + OUTPUT_FOLDER + 'for_krona.tab 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt')
system('mv text.krona.html ' + OUTPUT_FOLDER + '/krona_plot.html')
system('rm ' + OUTPUT_FOLDER + 'for_krona.tab')


def create_annotation():
    annot_1 = open(OUTPUT_FOLDER + '/annot.txt', 'w+')
    guide = open(OUTPUT_FOLDER + '/guide.txt', 'w+')
    phyla = list()
    classes = list()
    orders = list()
    families = list()
    for i in taxonomy_counters_dict.keys():
        tokens = i.split(';')
        cut_taxonomy = '.'.join([tokens[0], tokens[1], tokens[2], tokens[3], tokens[4]])
        guide.write(cut_taxonomy + '\n')
        curr_phylum = tokens[0] + '.' + tokens[1]
        if curr_phylum not in phyla:
            phyla.append(curr_phylum)
        curr_class = curr_phylum + '.' + tokens[2]
        if curr_class not in classes:
            classes.append(curr_class)
        curr_order = curr_class + '.' + tokens[3]
        if curr_order not in orders:
            orders.append(curr_order)
        curr_family = curr_order + '.' + tokens[4]
        if curr_family not in families:
            families.append(curr_family)

    annot_1.write('title_font_size\t33\n')
    annot_1.write('total_plotted_degrees\t340\n')
    annot_1.write('annotation_background_alpha\t0.1\n')
    annot_1.write('start_rotation\t270\n')
    annot_1.write('internal_label\t1\tDomain\n')
    annot_1.write('internal_label\t2\tPhyla\n')
    annot_1.write('internal_label\t3\tClasses\n')
    annot_1.write('internal_label\t4\tOrders\n')
    annot_1.write('internal_label\t5\tFamilies\n')
    annot_1.write('internal_labels_rotation\t270\n')
    for i in phyla:
        phrase_1 = i + '\tclade_marker_shape\th\n'
        annot_1.write(phrase_1)
        if 'UNK' in i:
            phrase_1 = i + '\tclade_marker_color\tred\n'
            annot_1.write(phrase_1)
    for i in classes:
        phrase_1 = i + '\tclade_marker_shape\tp\n'
        annot_1.write(phrase_1)
        if 'UNK' in i:
            phrase_1 = i + '\tclade_marker_color\tred\n'
            annot_1.write(phrase_1)
    for i in orders:
        phrase_1 = i + '\tclade_marker_shape\td\n'
        annot_1.write(phrase_1)
        if 'UNK' in i:
            phrase_1 = i + '\tclade_marker_color\tred\n'
            annot_1.write(phrase_1)
    for i in families:
        phrase_1 = i + '\tclade_marker_shape\ts\n'
        annot_1.write(phrase_1)
        if 'FOTU' in i:
            phrase_1 = i + '\tclade_marker_color\tred\n'
            annot_1.write(phrase_1)
    annot_1.close()
    guide.close()


create_annotation()
cmd = "graphlan_annotate --annot " + OUTPUT_FOLDER + "/annot.txt " + OUTPUT_FOLDER + '/guide.txt '
cmd += OUTPUT_FOLDER + "/guide.xml"
system(cmd)
system("graphlan " + OUTPUT_FOLDER + "/guide.xml " + OUTPUT_FOLDER + "/step.png --dpi 300 --size 6.5")
system("rm " + OUTPUT_FOLDER + "/annot.txt " + OUTPUT_FOLDER + '/guide.txt ' + OUTPUT_FOLDER + "/guide.xml")

sotu_zotu_map_dict = dict()
zotus_taxonomy = read_file(OUTPUT_FOLDER + '/zotus_with_taxonomy.tab')[1:]
for line in zotus_taxonomy:
    tokens = line.split('\t')
    taxonomy = tokens[-1]
    sotu = taxonomy.split(';')[-2]
    zotu = tokens[0]
    if sotu in sotu_zotu_map_dict.keys():
        sotu_zotu_map_dict[sotu] = sotu_zotu_map_dict[sotu] + ',' + zotu
    else:
        sotu_zotu_map_dict[sotu] = zotu

taxonomies = list()
for line in zotus_taxonomy:
    tokens = line.split('\t')
    taxonomy = tokens[-1]
    if taxonomy not in taxonomies:
        taxonomies.append(taxonomy)

out_file = open(OUTPUT_FOLDER + '/sotus_with_taxonomy.tab', 'w+')
header = read_file(OUTPUT_FOLDER + '/zotus_with_taxonomy.tab')[0]
out_file.write(header + '\n')
for taxonomy in taxonomies:
    out_line = taxonomy.split(';')[-2]
    curr_taxonomy_sample_sizes = list()
    for line in zotus_taxonomy:
        tokens = line.split('\t')
        curr_taxonomy = tokens[-1]
        if curr_taxonomy == taxonomy:
            samples_reads = '\t'.join(tokens[1:-1])
            curr_taxonomy_sample_sizes.append(samples_reads)
    samples_num = curr_taxonomy_sample_sizes[0].count('\t') + 1
    for i in range(samples_num):
        athroisma = 0
        for line in curr_taxonomy_sample_sizes:
            token = line.split('\t')[i]
            athroisma += int(token)
        out_line += '\t' + str(athroisma)
    out_line += '\t' + taxonomy + '\n'
    out_file.write(out_line)
out_file.close()


out_file = open(OUTPUT_FOLDER + '/zotu_sotu_map.tab', 'w+')
out_file.write('ZOTUs\tSOTU\n')
contents = read_file(OUTPUT_FOLDER + '/zotus_with_taxonomy.tab')[1:]
for line in contents:
    taxonomy = line.split('\t')[-1]
    zotu = line.split('\t')[0]
    sotu = taxonomy.split(';')[-2]
    out_file.write(zotu + '\t' + sotu + '\n')
out_file.close()

out_file = open(OUTPUT_FOLDER + '/sotu_gotu_map.tab', 'w+')
out_file.write('SOTUs\tGOTU\n')
contents = read_file(OUTPUT_FOLDER + '/zotus_with_taxonomy.tab')[1:]
passed_sotus = list()
for line in contents:
    taxonomy = line.split('\t')[-1]
    sotu = taxonomy.split(';')[-2]
    gotu = taxonomy.split(';')[-3]
    if sotu not in passed_sotus:
        out_file.write(sotu + '\t' + gotu + '\n')
        passed_sotus.append(sotu)
out_file.close()

out_file = open(OUTPUT_FOLDER + '/gotu_fotu_map.tab', 'w+')
out_file.write('GOTUs\tFOTU\n')
contents = read_file(OUTPUT_FOLDER + '/zotus_with_taxonomy.tab')[1:]
passed_gotus = list()
for line in contents:
    taxonomy = line.split('\t')[-1]
    gotu = taxonomy.split(';')[-3]
    fotu = taxonomy.split(';')[-4]
    if gotu not in passed_gotus:
        out_file.write(gotu + '\t' + fotu + '\n')
        passed_gotus.append(gotu)
out_file.close()


def sina_alignment():
    cmd = SINA_EXECUTABLE + ' --in=' + OUTPUT_FOLDER + OUTPUT_ASV_FASTA_WITH_TAXONOMY + ' --out=test_z.fasta --db='
    cmd += SILVA_ARB + ' --turn all --fasta-write-dna >/dev/null 2>/dev/null'
    system(cmd)
    #
    cmd = SINA_EXECUTABLE + ' --in=' + OUTPUT_FOLDER + OUTPUT_SOTU_FASTA_WITH_TAXONOMY + ' --out=test_s.fasta --db='
    cmd += SILVA_ARB + ' --turn all --fasta-write-dna >/dev/null 2>/dev/null'
    system(cmd)


# from the alignment get every column that has 1 base at least aligned
# so tree creation is correct and goes faster
def extract_columns(input_file_name):
    good_points = list()
    filepath = input_file_name
    with open(filepath) as fp:
        line = fp.readline()
        while line:
            if not line[0] == '>':
                counter = 0
                for curr_char in line:
                    if curr_char in ['A', 'C', 'G', 'T']:
                        good_points.append(counter)
                    counter += 1
            line = fp.readline()
    #
    good_points = list(set(good_points))
    #
    out_file = open('for_tree.fasta', 'w+')
    with open(filepath) as fp:
        line = fp.readline()
        while line:
            if line[0] == '>':
                header = line.split(';')[0] + '\n'
            else:
                curr_seq = ''
                for number in good_points:
                    curr_seq += line[number]
                out_file.write(header)
                out_file.write(curr_seq + '\n')
            line = fp.readline()
    out_file.close()
    system('mv for_tree.fasta ' + input_file_name)


def create_trees():
    top_dir = getcwd()
    extract_columns('test_s.fasta')
    extract_columns('test_z.fasta')
    rapidnj = top_dir + "/0.Setup_and_Testing/rapidNJ/bin/rapidnj "
    cmd1 = rapidnj + "test_s.fasta -n -o t -x " + OUTPUT_FOLDER + "/sotu_nj.tre"
    cmd2 = rapidnj + "test_z.fasta -n -o t -x " + OUTPUT_FOLDER + "/zotu_nj.tre"
    # cmd3 = "/crc/crc/binaries/FastTree -gtr -gamma -quiet -nt < test_s.fasta > sotu_aml.tre"
    # cmd4 = "/crc/crc/binaries/FastTree -gtr -gamma -quiet -nt < test_z.fasta > zotu_aml.tre"
    system(cmd1)
    system(cmd2)
    # print(cmd3)
    # print(cmd4)


sina_alignment()
create_trees()
