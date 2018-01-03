# NGS_Scripts
Various scripts that I made to analyze data for the Stupar lab. (Software versions are not listed, but will be added in the near future). Almost all of these scripts were intended to be run in parallel across mutiple nodes using GNU parallel.

## For new users
To download these scripts, you need to first log into MSI.
```
ssh <UMNx500>@login.msi.umn.edu
<x500 password>UMN
```
You will then need to change to the directory when you store your software. If you don't have one yet, make one.
```
#to make a directory
mkdir Software
#to enter a directory
cd Software
```
Once in your desired directory, you can then clone this repository

`git clone https://github.com/MeeshCompBio/Soybean_Scripts.git`

You can also export the path to your ~/.bashrc if you don't want to write out the full file path to use the software
Something like this at the end of the file (use "pwd -P" if you don't know the full path)

```export PATH=<path/to/directory/>Soybean_Scripts```

### Castseq.sh
This is a bash script to merge old sequencing files on MSI that were originally split based on their sequencing lane.

### Fastq2FPKM
These are bash scripts to automate the handling of RNA-seq data. The script can either use single or paried-end reads and will take trim off illumina truseq adapters.
* To get a full list of options and parameters
```
<path/to/script/>Fastq2FPKMV3.sh -h
```
* Fastq2FPKMV3.sh was made to use trimmomatic, (STAR/HISAT2) and cufflinks. 
* Fastq2FPKMGmaxV2.sh was made to run cutadapt, tophat2 and cufflinks. This is the older version of the two scripts and tophat has since been depreciated.

### Fastq2Readmap
These are bash scripts to automate the handling of whole-genome sequencing data.

* To get a full list of options and parameters
```
<path/to/script/>Fastq2ReadmapGmaxV2.sh -h
```
* Fastq2FPKMV3.sh was made to use cutadapt, bwa/bowtie2 and cufflinks
* Fastq2ReadmapGmaxV1.sh was made to run cutadapt, and bwa. This is the older version of the two scripts and some of the files were hardcoded so please avoid using it unless you are on the Stupar Labs MSI.

### MSIVariantCaller
This is just and old varaint caller scipt for new memebers of the lab to see how to submit a job on a PBS queue system. GATK Unified Genotyper is now depreciated.

### Orphan.pl
Used to identify orphaned reads from a .bam file. To be used in conjuntion with Transgenemap.sh

### SRA_Rename
Used to rename SRA given file names to the original samples names uploaded to the database
* To get a full list of options and parameters
```
<path/to/script/>SRA_Rename.sh -h
```

### TransGeneMap
Theis a bash script to identify where transgenes are integrated into the genome using whole-genome sequencing and vector sequence.
* To get a full list of options and parameters
```
<path/to/script/>TransGenMap.sh -h
```

### VCFquery.sh
This is a bash script for Stupar lab members to access variant calls and indels for specific regions of the genome. Use the -h option to get more information on the options and how to run the script. Currently I have indels and variants called for the 106 Genomes and the NAM parent population.

`bash <path/to/script/>VCFquery.sh -h`



