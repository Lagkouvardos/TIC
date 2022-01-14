![TIC-LOGO-v2](https://user-images.githubusercontent.com/8244618/139100265-7fb3319d-93a9-47f4-a617-ffe82a33b77d.png)

# A list of Python scripts for the analysis of microbial profiles and illumination of the diversity contained within

TIC-Pipeline is a suite that processes FASTQ files to create denoised sequences based on the latest golden standards,
and then performs a diversity analysis of the ZOTUs.

## Introduction:
In 16S rRNA sequence-based diversity analysis, a common practice is the clustering of the sequences based on similarity cutoffs to obtain groups reflecting molecular species, genera or families. Due to the size of the available data, greedy algorithms are preferred for their time efficiency. Such algorithms rely only on pairwise similarities and tend to cluster together sequences with diverse phylogenetic background. Taxonomic classifiers use position specific taxonomic information in assigning probable taxonomy to a given sequence. We developed a tool that uses classifier assigned taxonomy to restrict the clustering among sequences that share the same taxonomic path. We tested “Taxonomy Informed Clustering” with the sequences from SILVA release 138* and show that TIC outperforms greedy algorithms like Usearch and Vsearch in terms of clusters purity and entropy**. We offer a complete and automated pipeline for use in diversity analysis context, based on this concept, and a pipeline for 16S rRNA amplicon dataset processing. This implementation first process raw reads down to denoised amplicons, taxonomically classify them and apply TIC to cluster them further to sOTUs, gOTUs  and fOTUs.  The resulting tables offer more accurate insights at different evolutionary levels views that will be useful in microbiome research.


## Description:

### Organization
TIC-Pipeline is composed of 5 steps that can be run independently or as a set.

    1.Sample_Processing
    2.Taxonomy_Classification
    3.Extraction_Denoising
    4.Taxonomy_Informed_Clustering
    5.Results_Reporting

Running them in the given order simplifies the process as the outputs of each step are often the inputs of the next. Before running any script, please make sure you have read and fully understood the corresponding section of the documentation for each step.

## Requirements:
TIC-Pipeline is a mix of bash, Python and R scripts. Due to relience on bash tools we propose the usage of a Linux system for simplifying running the scripts.

## Installation:
In the folder: 0.Setup_and_Testing we provide a script called setup.sh that downloads and installs all required tools
for the pipeline.
To run the script execute the following commands:
```
cd 0.Setup_and_Testing
chmod +x setup.sh
./setup.sh
```

## Testing:
Before running the pipeline with your data, we propose to check if all programs are correctly installed.
To do this change the TESTING_MODE in config_options.txt to YES
and run:
```
python3 TIC-Pipeline.py
```
Based on the output of this command, you can verify that the pipeline will find the necessary tools.


## Citation:
