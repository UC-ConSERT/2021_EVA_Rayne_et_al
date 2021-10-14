#!/bin/bash
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
NCPU=$4
mkdir log_files

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

## A very basic (not aesthetic) script for filtering GBS data with vcftools
## Files were demultiplexed and barcodes/adaptors trimmed with https://github.com/Lanilen/GBS-PreProcess; and SNPs were genotyped using Stacks2 (refmap.pl wrapper)


echo 'Initial filter for minimum mean depth of 10, maximum missingness and remove indels 0.5'
vcftools --vcf 02_wildunfiltered.vcf --max-missing 0.5 --remove-indels --min-meanDP 10 --recode --recode-INFO-all --stdout > 03_geno0.5minmeanDP10removedindels.vcf

echo 'Keeping individuals with <90% missingness'
mkdir missing-out
vcftools --vcf 03_geno0.5minmeanDP10removedindels.vcf --missing-indv --out ./missing-out/90
mawk '$5 > 0.9' ./missing-out/90.imiss | cut -f1 > remove90.indv
vcftools --vcf 03_geno0.5minmeanDP10removedindels.vcf --remove remove90.txt  --recode --recode-INFO-all --out 04_geno0.5miss90

echo 'Removing loci with fewer than 60% SNPs'
vcftools --vcf 04_geno0.5miss90.recode.vcf --max-missing 0.6 --recode --recode-INFO-all --stdout > 05_geno0.6miss90.vcf

echo 'Keeping individuals with <70% missingness'
vcftools --vcf 05_geno0.6miss90.vcf --missing-indv --out ./missing-out/70
mawk '$5 > 0.7' ./missing-out/70.imiss | cut -f1 > remove70.txt#vcftools --vcf 05_geno0.6miss90.vcf --remove remove70.txt  --recode --recode-INFO-all --out 06_geno0.6miss70
echo 'Removing loci with fewer than 70% SNPs'
vcftools --vcf 06_geno0.6miss70.recode.vcf --max-missing 0.7 --recode --recode-INFO-all --stdout > 07_geno0.7miss70.vcf

echo 'Keeping individuals with <50% missingness'
vcftools --vcf 07_geno0.7miss70.vcf --missing-indv --out ./missing-out/50
mawk '$5 > 0.5' ./missing-out/50.imiss | cut -f1 > remove50.txt
vcftools --vcf 07_geno0.7miss70.vcf --remove remove50.txt  --recode --recode-INFO-all --out 08_geno0.7miss50

echo 'Filtering for min depth 3, minor allele count 3 and max mean depth 150'
vcftools --vcf 08_geno0.7miss50.recode.vcf --minDP 3 --mac 3 --max-meanDP 150 --recode --recode-INFO-all --out 09_geno0.7miss50mac

echo 'Removing loci with fewer than 90% SNPs'
vcftools --vcf 09_geno0.7miss50mac.recode.vcf --max-missing 0.9 --recode --recode-INFO-all --stdout > 10_geno0.9miss50mac.vcf

echo 'Keeping individuals with <40% missingness'
vcftools --vcf 10_geno0.9miss50mac.vcf --missing-indv --out ./missing-out/40
mawk '$5 > 0.4' ./missing-out/40.imiss | cut -f1 > remove40.txt
vcftools --vcf 10_geno0.9miss50mac.vcf --remove remove40.txt  --recode --recode-INFO-all --out 11_geno0.9miss40mac

echo 'Removing loci with fewer than 95% SNPs'
vcftools --vcf 11_geno0.9miss40mac.recode.vcf --max-missing 0.95 --recode --recode-INFO-all --stdout > 12_geno0.95miss40mac.vcf

echo 'Keeping individuals with <25% missingness'
vcftools --vcf 12_geno0.95miss40mac.vcf --missing-indv --out ./missing-out/25
mawk '$5 > 0.25' ./missing-out/25.imiss | cut -f1 > remove25.txt
vcftools --vcf 12_geno0.95miss40mac.vcf --remove remove25.txt --recode --recode-INFO-all --out 13_geno0.95miss25mac

echo 'Discard loci with linkage disequilibrium where r2 > 0.6 in 1000bp window'
bgzip -c 13_geno0.95miss25mac.recode.vcf > 13_geno0.95miss25mac.recode.vcf.gz
tabix -fp vcf 13_geno0.95miss25mac.recode.vcf.gz
bcftools +prune -m 0.6 -w 1000 -Ov 13_geno0.95miss25mac.recode.vcf.gz -o 14_geno0.95miss25macLD.vcf
vcftools --vcf 14_geno0.95miss25macLD.vcf --out geno0.95miss25macLD

echo 'Keeping individuals with <10% missingness'
vcftools --vcf 14_geno0.95miss25macLD.vcf --missing-indv --out ./missing-out/10
mawk '$5 > 0.10' ./missing-out/10.imiss | cut -f1 > remove10.txt
vcftools --vcf 14_geno0.95miss25macLD.vcf --remove remove10.txt --recode --recode-INFO-all --out 15_geno0.95miss10mac

echo 'Retaining biallelic SNPs only'
bcftools view -c 1 -v snps 15_geno0.95miss10mac.vcf -o filtered.vcf
vcftools --vcf filtered.vcf --recode --recode-INFO-all --out filtered
rm filtered.vcf
vcftools --vcf filtered.recode.vcf --missing-indv --out ./missing-out/final
cut -f1 ./missing-out/final.imiss > popmap.txt
populations -V filtered.recode.vcf -M popmap.txt -p 19 -r 0.01  -O ./ --vcf # to get populations summary output from Stacks

echo 'Heck yeah, we are finished.'
