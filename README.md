# Soybean_Scripts
Various scripts that I made to analyze data for the Stupar lab

##For new users
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

###VCFquery.sh
This is a bash script for Stupar lab members to access variant calls and indels for specific regions of the genome. Use the -h option to get more information on the options and how to run the script. Currently I have indels and variants called for the 106 Genomes and the NAM parent population.

`bash <path/to/script/>VCFquery.sh`



