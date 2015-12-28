#!/bin/bash
"""This is an example script I used to call variants in soy on the MSI servers""

#PBS -l walltime=8:00:00,nodes=1:ppn=24,mem=62gb
#PBS -V
#PBS -N <logfilename>
#PBS -M <youremail>
#PBS -m abe


set -euo pipefail

module load java

cd /home/stuparr/michnoj0/Scratch/Austin/


java -jar /home/stuparr/shared/Software/GATK-3.3/GenomeAnalysisTK.jar \
        -T UnifiedGenotyper \
        -dcov 200 \
        #number of threads
        -nt 24 \
        #where your reference genome is located
        -R /home/stuparr/shared/References/Gmax.a2.v1/assembly/Gmax_275_v2.0.fa \
        #use -I option for each one of your .bams you want to call SNP's for
        -I /home/stuparr/michnoj0/Scratch/<File1.sorted.bam> \
        -I /home/stuparr/michnoj0/Scratch/<File2.sorted.bam> \
        #--alleles is for what alleles you want SNP's to be called for
        --alleles /home/stuparr/michnoj0/Scratch/Genotypes_50K/GMAXV2.recode.vcf \
        --genotyping_mode GENOTYPE_GIVEN_ALLELES \
        # this line will return genotype calls for all sites even if there is no SNP
        --output_mode EMIT_ALL_SITES \
        #output file
        -o HiLowSucrose_FNLines_Parents_50KGmaxV2.vcf
