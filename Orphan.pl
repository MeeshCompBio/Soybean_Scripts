#! /usr/bin/perl

#THIS IS A SLIGHTLY MODIFIED VERSION OF THE SCRIPT ORIGINALLY MADE BY Anuj Srivastava
#I modified it so it could take STDIN as input so I could use samtools to stream in a .bam
#PLEASE CITE ORIGINAL PAPER FOR THIER WORK
#Srivastava, A., Philip, V. M., Greenstein, I., Rowe, L. B., Barter, M., Lutz, C., &amp; Reinholdt, L. G. (2014).
#Discovery of transgene insertion sites by high throughput sequencing of mate pair libraries. BMC
#Genomics, 15(1), 367. http://doi.org/10.1186/1471-2164- 15-367

use strict;
use warnings;
use Getopt::Long;

################Written By Anuj Srivastava#######################################
my $usage = <<'USAGE';
#####################extract_unmapped_mates.pl###################################

 usage: perl extract_mapped_mates.pl  bowtie2_out.sam  gene_name_to_search\|gene2  orphaned_mates_out_file_name 

                -gs= gene string to search [required; could be multiple gene separated by \|]
                -o = outname attached to each created file [required]

##################################################################################

USAGE

my (  $geneString, $out_append);

my $result = GetOptions(
    "gs=s"  => \$geneString,
    "o=s"   => \$out_append,
)
;

die $usage unless ($geneString); ##### Mandatory arguments
die "$usage must set -o  must provide an name attached to each file" unless ($out_append);               ##### Mandatory arguments

#system ("grep . $infile   >$out_append.results"); #### picking the hits mapped to transgene
#system ("awk '{print $3}' $infile >$out_append.results"); #### picking the hits mapped to transgene

open(UnmappedMates, ">$out_append.fastq") || die "cannot open the file";                   #### writing the file to store mates whose one end mapped to transgenes

my ($mateMapped1, $mateMapped2);
my $flag = 0;
my @samFieldsMate1 =(); my @samFieldsMate2 = ();

#modified this file to take in stardard in
while(my $readTGHIT = <STDIN>) #### reading the sam file with transgene hits only; this file is already sorted by name
{
  if(($readTGHIT =~ /^\s*(.*?\s+){11,}/) && ($flag == 0))
  {
    chomp $readTGHIT;
    @samFieldsMate1 = split("\t", $readTGHIT);
    $mateMapped1 = $samFieldsMate1[5];
    $flag = 1;
    next;
  }

  if(($readTGHIT =~ /^\s*(.*?\s+){11,}/) && ($flag == 1))
  {
    chomp $readTGHIT;
    @samFieldsMate2 = split("\t", $readTGHIT);
    $mateMapped2 = $samFieldsMate2[5];
    $flag = 0;

     if(($mateMapped1 =~ /\*/) && ($mateMapped2 !~ /\*/))
     {

       print UnmappedMates  "\@$samFieldsMate1[0]\n$samFieldsMate1[9]\n+\n$samFieldsMate1[10]\n"; #### writing orpahned mate1

     }

    if(($mateMapped1 !~ /\*/) && ($mateMapped2 =~ /\*/))
     {

      print  UnmappedMates "\@$samFieldsMate2[0]\n$samFieldsMate2[9]\n+\n$samFieldsMate2[10]\n"; #### writing orphaned mate2

     }

    @samFieldsMate1 = (); 
    @samFieldsMate2 = ();
    next;

  }
  

}
