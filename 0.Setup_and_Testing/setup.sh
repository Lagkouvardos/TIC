echo 'Downloading usearch version 11, 32-bit version'
wget http://129.187.44.105:8082/usearch
chmod +x usearch
echo '    Done'

echo 'Downloading vsearch version 2.17, 64-bit version'
wget http://129.187.44.105:8082/vsearch
chmod +x vsearch
echo '    Done'

echo 'Downloading SORT_ME_RNA_DATABASES'
wget http://129.187.44.105:8082/sort_me_rna_dbs/silva-bac-16s-id90.fasta
wget http://129.187.44.105:8082/sort_me_rna_dbs/silva-arc-16s-id95.fasta
echo '    Done'

echo 'Downloading SORT_ME_RNA tool'
wget http://129.187.44.105:8082/sortmerna-4.3.3-Linux.tar.gz
echo '    Done'

echo 'Installing SORT_ME_RNA tool'
echo 'when asked choose y and y'
tar -xvzf sortmerna-4.3.3-Linux.tar.gz
./sortmerna-4.3.3-Linux/sortmerna-4.3.3-Linux.sh
rm sortmerna-4.3.3-Linux.tar.gz
echo '    Done'

echo 'Downloading Silva ARB file'
wget http://129.187.44.105:8082/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.gz
gunzip -d SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.gz
echo '    Done'

echo 'Downloading sina aligner tool'
wget http://129.187.44.105:8082/sina-1.7.2-linux.tar.gz
tar -xvzf sina-1.7.2-linux.tar.gz
rm sina-1.7.2-linux.tar.gz
echo '    Done'

echo 'Downloading rapidNJ tool'
wget http://129.187.44.105:8082/rapidNJ.tar.gz
tar -xvzf rapidNJ.tar.gz
rm rapidNJ.tar.gz
cd rapidNJ
make
cd ..
echo '    Done'

echo 'Downloading Krona tool'
wget http://129.187.44.105:8082/KronaTools-2.8.tar
tar -xf KronaTools-2.8.tar
rm KronaTools-2.8.tar
sudo ./KronaTools-2.8/install.pl
echo '    Done'


echo 'Downloading Template Data'
wget http://129.187.44.105:8082/template_data.zip
unzip template_data.zip
rm template_data.zip
echo '    Done'


CWD=$(pwd)
echo 'SILVA_ARB:'$CWD'/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb' >> ../tool_and_db_options.ini
echo 'SINA_EXECUTABLE:'$CWD'/sina-1.7.2-linux/sina' >> ../tool_and_db_options.ini
echo 'SORT_ME_RNA_DB1:'$CWD'/silva-bac-16s-id90.fasta' >> ../tool_and_db_options.ini
echo 'SORT_ME_RNA_DB2:'$CWD'/silva-arc-16s-id95.fasta' >> ../tool_and_db_options.ini
echo 'SORT_ME_RNA_TOOL:'$CWD'/sortmerna-4.3.3-Linux/bin/sortmerna' >> ../tool_and_db_options.ini
echo 'CLUSTERING_TOOL:'$CWD'/vsearch' >> ../tool_and_db_options.ini
echo 'KRONA_TOOL:'$CWD'/KronaTools-2.8/scripts/ImportText.pl' >> ../tool_and_db_options.ini
echo 'RAPID_NJ:'$CWD'/rapidNJ/bin/rapidnj' >> ../tool_and_db_options.ini
echo 'AL_CLUSTERING_TOOL:'$CWD'/usearch' >> ../tool_and_db_options.ini

echo 'Installing python packages'
pip3 install -r python_requirements.txt

echo 'Installing R packages'
sudo Rscript install_packages.R
