Notes from Jitender at \\jic-hpc-data\HPC-Home\jitender.md
# /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia
```bash
color-combi.py #unknown function
```
```bash
VCF2Dis #VCF2Dis: A new simple and efficient software to calculate p-distance matrix based Variant Call Format
210s.M_persicae.onlySNPs.vcf.gz #input to VCF2Dis - presumably
p_dis.mat #output from VCF2Dis
```
distmat calculates the evolutionary distance between every pair of sequences in a multiple sequence alignment. A variety of methods to estimate distance may be selected, and differ in how they correct the observed substitution rates to more accurately reflect the true evolutionary distance. An output file containing a distance matrix for the set of sequences is written. The distances are expressed in terms of the number of substitutions per 100 bases or amino acids.
```bash
generate-mperc-pdis-csv.py #Generate an ID based p-distmat in for R (p-distance matrix)
p_dis.mat #input for p_dis_mperc.csv
p_dis_mperc.csv #output of p_dis_mperc.csv
```

## /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/reference
```bash
find-tss-tes-from-exons-gen-bed-gene1.py #
mperc-gff3.db #
MYZPE13164_O_EIv2.1.annotation.gff3 #
#Myzus_persicae_O_v2.0.scaffolds.braker2.gff3
#Myzus_persicae_O_v2.0.scaffolds.dict
#Myzus_persicae_O_v2.0.scaffolds.fa.amb
#Myzus_persicae_O_v2.0.scaffolds.fa.ann
#Myzus_persicae_O_v2.0.scaffolds.fa.bwt
#Myzus_persicae_O_v2.0.scaffolds.fa.fai
#Myzus_persicae_O_v2.0.scaffolds.fa.gz
#Myzus_persicae_O_v2.0.scaffolds.fa.pac
#Myzus_persicae_O_v2.0.scaffolds.fa.sa
tair10-gff3.db #
tss-test-from-exons-gene.bed #
```

## step-size-one and steps-faster
outputs from something - unknown

## /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/position-wise-snp-diff

Folder contains work to differentiate genic and non-genic SNPs?

```bash
cat README.md
#Positional SNP diff
#-------------------
#
#- We find positional SNP diff
#- 11870914
#
#[cheemaj@NBI-HPC interactive TABLE]$ wc -l mperc-210-genic-regions.vcf
#11,870,914
# 7,216,966 mperc-210-genic-regions.vcf
#[cheemaj@NBI-HPC interactive TABLE]$

cat vcf_filenames.txt
#Prunus_group1
#O_containing_group
#Prunus_group1
#Prunus_group2
#G006_group

cat Prunus_group2.txt
#S55
#S66
#S56
#S67
#S69
#S54
#S57
#S62
#S53
#S52
#S65
#S58
#S68
#I1
#S24

cat Prunus_group2.txt | wc -l
#15
cat Prunus_group1.txt | wc -l
#32
cat O_group_sample_names.txt | wc -l
#66
cat G006_group.txt | wc -l
#7
# look to be unique (names of samples) but dont add to 209 or 210?

#Nic.txt looks to be similar to the group names files but isn't in vcf_filenames.txt
cat Nic.txt | wc -l
#40

prunis.txt #is in the same format as those above however is not unique, and does not equal prunus group1+group2?

nano PCA_file_host.csv
#contains information on each of 210 samples:
#,sample.id,EV1,EV2,EV3,EV4,host,Collection,Country,Continent,simple_host,Collected_by,dummy
#what are EV1,2,3,4?
```

```bash
jic-bass-samples-extract.py #purpose unknown
PCA_file_host.csv #input for jic-bass-samples-extract.py
```
How to view these, cannot download/open files from mobaxterm?:
```bash
igv_snapshot.png
igv_snapshot-region2.png
```

Why?
```bash
add-group-info-Roland.py
#This script defines several lists of strings containing names or IDs of different groups, such as group_o, group_G006, group_Nic, group_prunus, group_Prunus_group1, and group_Prunus_group2.

#Then it opens a file named Prunus_group2.txt, reads the lines in the file and appends each line (after removing leading and trailing white spaces) to the list group_o. It then prints the length of group_o and the list itself.

#Finally, it opens a CSV file named PCA_file_host.csv, reads the lines in the file, and for each line, creates a dictionary d mapping the header values to the corresponding field values. It then stores the Collected_by value of each dictionary in a dictionary named geoD with the sample.id value as the key. It prints the keys of the geoD dictionary.
```

### /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/position-wise-snp-diff/TABLE
```bash
cat protocol.md

#[cheemaj@NBI-HPC interactive TABLE]$ pwd
#/jic/scratch/groups/Matthew-Hartley/cheemaj/scratch/scratch-work/jitender/saskia/position-wise-snp-diff/TABLE
#[cheemaj@NBI-HPC interactive TABLE]$ jit subset-genes-gff3-scaffolds1-6.py
#Done
#[cheemaj@NBI-HPC interactive TABLE]$
#singularity exec /hpc-home/cheemaj/BUILD/PYBEDTOOLS/pybed.simg bedtools  intersect \
#    -a 210s.M_persicae.onlySNPs.vcf.gz \
#    -b MYZPE13164_O_EIv2.1.annotation_genes_scaff1-6.gff3  \
#    -header > mperc-210-genic-regions.vcf
### sweed
##installation
#[cheemaj@NBI-HPC SweeD_v3.2.1_Linux]$ make -f Makefile.PTHREADS.gcc
#[cheemaj@NBI-HPC SweeD_v3.2.1_Linux]$ pwd
#/hpc-home/cheemaj/BUILD/PAML/SweeD_v3.2.1_Linux
#[cheemaj@NBI-HPC SweeD_v3.2.1_Linux]$ ./SweeD -h
###
#source package /tgac/software/testing/bin/bcftools-1.12
#bcftools stats mperc-210-genic-regions.vcf  > mperc-210-genic-regions_stats.txt
#bgzip -c file.vcf > file.vcf.gz
#11,870,914
# 7,216,570
#/hpc-home/cheemaj/BUILD/PAML/SweeD_v3.2.1_Linux/SweeD   -name mperc-210-genic-regions  -input mperc-210-genic-regions.vcf  -grid 100
#/hpc-home/cheemaj/BUILD/PAML/SweeD_v3.2.1_Linux/SweeD   -name mperc-210-genic-regionsF  -input Myzus_persicae_O_v2.0.scaffolds.fa -grid 100
```
```bash
cat resum1.sh #/hpc-home/cheemaj/BUILD/PAML/SweeD_v3.2.1_Linux/SweeD -name mperc-210-genic-regions-SweeD -input mperc-210-genic-regions.vcf -grid 100
#SweeD	(Sweep	Detector,	Pavlidis	&	Alachiotis) - https://cme.h-its.org/exelixis/resource/download/software/sweed3.0_manual.pdf - SweeD implements a composite likelihood ratio test which detects complete selective sweeps using Site Frequency Spectrum (SFS) patterns of single-nucleotide polymorphisms (SNPs).
SweeD_Info.mperc-210-genic-regions #output from SweeD
SweeD_Info.mperc-210-genic-regions-SweeD #output from SweeD
SweeD_Report.mperc-210-genic-regions #output from SweeD
SweeD_Report.mperc-210-genic-regions-SweeD #output from SweeD
sw.err #ouput from resum1.sh
sw.out #ouput from resum1.sh
sw.txt #ouput from resum1.sh
```
```bash
cat readme-tab.py # 41,304,161 mil snp s
```
```bash
subset-genes-gff3-scaffolds1-6.py #extract gene annotation information for chromosomes 1 to 6 and writes this to a new gff3 file
MYZPE13164_O_EIv2.1.annotation.gff3 #input for subset-genes-gff3-scaffolds1-6.py
MYZPE13164_O_EIv2.1.annotation_genes_scaff1-6.gff3 #output of subset-genes-gff3-scaffolds1-6.py

```
```bash
part-joiner.py #reads medelseq-allele-table-*.tab files and writes them to csv, checking that the no. of columns in correct
```
```bash
submitter.py #executes a singularity container to run vc2alle-part.py
vc2alle-part.py #polymorphic information content calculated
```
where do these come from/used for?:
```bash
mperc-210-genic-regions_stats.txt
mperc-210-genic-regions.vcf
```
#### /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/position-wise-snp-diff/TABLE/BCF
```bash
ls -l
-rwxrwx--- 1 cheemaj CBU   238964012 Apr 12 11:15 gene-hits-snp-O-vs-rest.tab
-rwxrwx--- 1 cheemaj CBU   222162123 Apr 12 11:15 gene-hits-snp.tab
-rwxrwx--- 1 cheemaj CBU         272 Apr 12 11:15 get_samples_names.py
-rwxrwx--- 1 cheemaj CBU 32927833880 Apr 12 11:24 mperc-210-genic-regions.vcf
-rwxrwx--- 1 cheemaj CBU  3676030269 Apr 12 11:25 mperc-210-genic-regions.vcf.gz
-rwxrwx--- 1 cheemaj CBU      269997 Apr 12 11:15 mperc-210-genic-regions.vcf.gz.tbi
-rwxrwx--- 1 cheemaj CBU     6054027 Apr 12 11:15 MYZPE13164_O_EIv2.1.annotation_genes_scaff1-6.gff3
-rwxrwx--- 1 cheemaj CBU        6886 Apr 12 11:15 perc-test.py
-rwxrwx--- 1 cheemaj CBU       10783 Apr 12 11:15 perc-test-sep-O-vs-rest.py
-rwxrwx--- 1 cheemaj CBU        7008 Apr 12 11:24 perc-test-sep.py
-rwxrwx--- 1 cheemaj CBU         405 Apr 12 11:15 README-manual.md
-rwxrwx--- 1 cheemaj CBU        2437 Apr 12 11:15 rest.py
-rwxrwx--- 1 cheemaj CBU        1702 Apr 12 11:15 samples.py
-rwxrwx--- 1 cheemaj CBU         346 Apr 12 11:38 scaffold-lengths.txt
-rwxrwx--- 1 cheemaj CBU        5345 Apr 12 12:09 snp-discover-saturation-curve-all-jic-bass.py
-rwxrwx--- 1 cheemaj CBU        2662 Apr 12 11:15 submitter.py
-rwxrwx--- 1 cheemaj CBU        2842 Apr 12 11:15 vc2alle-part.py

```

