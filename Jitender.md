# Jitender 

Notes from Jitender at \\jic-hpc-data\HPC-Home\jitender.md
#### /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/powerpoints

contains powerpoint presentations
## Data
References and filtered vcf files are from previous work by Roberto Beillo and others

#### /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/reference

Reference files from Mathers et al., 2020 Chromosome-scale genome assemblies of aphids reveal extensively rearranged autosomes and long-term conservation of the X chromosome (https://zenodo.org/record/3712089#.ZEqUPnbMJD8).
```bash
#Myzus_persicae_O_v2.0.scaffolds.braker2.gff3
#Myzus_persicae_O_v2.0.scaffolds.dict
#Myzus_persicae_O_v2.0.scaffolds.fa.amb
#Myzus_persicae_O_v2.0.scaffolds.fa.ann
#Myzus_persicae_O_v2.0.scaffolds.fa.bwt
#Myzus_persicae_O_v2.0.scaffolds.fa.fai
#Myzus_persicae_O_v2.0.scaffolds.fa.gz
#Myzus_persicae_O_v2.0.scaffolds.fa.pac
#Myzus_persicae_O_v2.0.scaffolds.fa.sa

find-tss-tes-from-exons-gen-bed-gene1.py # Generate the bed file from the GFF3 - by gathering exons for an isoform, input = MYZPE13164_O_EIv2.1.annotation.gff3
tss-test-from-exons-gene.bed #output from find-tss-tes-from-exons-gen-bed-gene1.py
MYZPE13164_O_EIv2.1.annotation.gff3 #source? EI?
mperc-gff3.db #
tair10-gff3.db #
```
#### /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/vcffilesafterfiltering
```bash
210s.M_persicae.onlySNPs_stats.txt
210s.M_persicae.onlySNPs.vcf.gz
```

## P_distance
p-distance  is the proportion (p) of nucleotide sites at which two sequences being compared are different. It is obtained by dividing the number of nucleotide differences by the total number of nucleotides compared.
#### /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia
```bash
color-combi.py #unknown function
```
readme.md describes:
```bash
VCF2Dis #VCF2Dis: A new simple and efficient software to calculate p-distance matrix based Variant Call Format. 
210s.M_persicae.onlySNPs.vcf.gz #input to VCF2Dis - presumably
p_dis.mat #output from VCF2Dis

bcftools stats 210s.M_persicae.onlySNPs.vcf.gz > 210s.M_persicae.onlySNPs_stats.txt
```
distmat calculates the evolutionary distance between every pair of sequences in a multiple sequence alignment. A variety of methods to estimate distance may be selected, and differ in how they correct the observed substitution rates to more accurately reflect the true evolutionary distance. An output file containing a distance matrix for the set of sequences is written. The distances are expressed in terms of the number of substitutions per 100 bases or amino acids.
```bash
generate-mperc-pdis-csv.py #Generate an ID based p-distmat in for R (p-distance matrix)
p_dis.mat #input for generate-mperc-pdis-csv.py
p_dis_mperc.csv #output of generate-mperc-pdis-csv.py
```

## corehunter
p_dis_mperc.csv as an input
#### step-size-one and steps-faster
step-size-one is mor complete? Has output files including satutration curve files
```bash
nano batcher-core-mperc-stepsize-1.py #makes resum and .R scripts for each like submitter does in other folder
```
```bash
nano resum11.sh
		#!/bin/bash
		#SBATCH --job-name=12
		#SBATCH -o  11.out
		#SBATCH -e  11.err
		#SBATCH --mem 40gb
		#SBATCH --nodes=1
		#SBATCH --ntasks-per-node=8
		#SBATCH -p jic-short

		# on dmatrix
		singularity exec ~/BUILD/ABEL/corehunter.cif  Rscript 11.R

		echo DONE
```
```bash
nano 11.R:
				options(java.parameters = "-Xmx50G")
				library(corehunter)
				options(java.parameters = "-Xmx50G")

				set.seed(1234)
				dist <- distances(file = "p_dis_mperc.csv")
				core <- sampleCore(dist, size = 11)
				write.table(as.list(core$sel), file = "11_core_sel.csv", row.names=FALSE, col.names=FALSE, sep=",") # output

#Core Hunter is a tool to sample diverse, representative subsets from large germplasm collections, with minimum redundancy. Such so-called core collections have applications in plant breeding and genetic resource management in general. Core Hunter can construct cores based on genetic marker data, phenotypic traits or precomputed distance matrices, optimizing one of many provided evaluation measures depending on the precise purpose of the core (e.g. high diversity, representativeness, or allelic richness).
```
```bash
nano melt-mperc-core-table.py #input = 1l_core_sel.csv #output="lablab_core_sel_stepsize1.csv" or similar 
```

## NJ tree
NJ tree was generated from p_distance combined with corehunter? So doesn't contain all samples just core?
#### /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/NJ
Contains neighbourhood joining phylogeny tree work and PCA work
README.md - using https://bioinformatics.phylolab.net/form/tree-format-conversion to convert from newick to nexus formatting
```bash
mperc_210_newick.nexus.tre #nexus formatted tree
mperc_210_newick.tre #newick formatted tree
```

## MDS
#### /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/NJ
Figures and images moved to own folder, MDS = multidimensional scaling
Principal component analysis starts with a correlation matrix, while multidimensional scaling can start with an inter-subject distance matrix or a correlation matrix. 

```bash
mds-3d-geographical.py #input = PCA_file_host.csv #output p_dis_mperc-MDS-geo.html

mds-3d.py #input p_dis_mperc-non-lig.csv #output p_dis_mperc-MDS.html

mds-3d-recent-bass.py #input PCA_file_host.csv, p_dis_mperc-non-lig.csv #output p_dis_mperc-non-lig-A015-MDS.csv, p_dis_mperc-MDS-recent.html

mds.py #input p_dis_mperc-non-lig.csv #output test.pdf

melt-mperc-core-table.py #input mperc_core_sel.csv #output

PCA_file_host.csv #table containing host and geography information for each samples as well as EV1-4

nano Genome_samples_location_host.txt #list of 106 aphid samples, country and host of origin
#Genome_samples_location_host.xlsx is the same as the above but an excel

RW_labmeeting_19_08_21.pptx #slides on the project unknown author, no script

/CORE/style-gen-core-mperc.py #input 0_core_sel.csv
/CORE #also contains tree images and other files
```

## Nei's gene diversity
#### /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/diversity

README.md - performing nei's diversity
```bash
nei-from-vcf.py --infile 210s.M_persicae.onlySNPs.vcf --outfile mperc_210.tab
```
```bash
x.py #print ['scaffold_' + str(c) for c in range(1,7)]
top500-lines.txt # top 500 lines from bcftoolsCommand=mpileup -Ou --threads 16 --annotate AD,DP -b ./txt/bam_files_fin.txt -R ./txt/scaff01.txt -f /jic/scratch/groups/Saskia-Hogenhout/roland/Myzus_periscae_popgen/reference/Myzus_persicae_O_v2.0.scaffolds.fa
```

## P450
#### /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/P450

- Look for the statUs of cyctochrome P450 in genomic and mitochondrial
- Find the copy number variation
- Clonal group "O" marked in the splitree plot
- Comparison(A): Group O versus  Prunus, GRC, G006
- O group is a cluster in the PCA plt witj mostly JIC and with some Bass
- Comparison(B): decade old Bass( prefixes S1, S2,..)  and and a year old JIC collection
- Note US1L and UK_SB are on the top of each other in the PCA plot

```bash
makeblastdb -in Myzus_persicae_O_v2.0.scaffolds.fa   -input_type fasta -dbtype nucl  -title MP  -parse_seqids -out MP

blastn -db MP  -query p450-mperc.fa  -out p450-hits.txt -outfmt "6 qseqid sseqid sstart send qstart qend length mismatch gapopen gaps evalue" -max_target_seqs 5
```






































## Genic SNPs

#### /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/position-wise-snp-diff

This first level of the directory structure seems to just contain groups of samples and scripts to generate these

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
# 
# 7,216,966 #of a total 11,870,914 SNPs
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
jic-bass-samples-extract.py #find samples from JIC and from Bass group, input = PCA_file_host.csv
PCA_file_host.csv #input for jic-bass-samples-extract.py
```
What am I looking at with these?:
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

#### /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/position-wise-snp-diff/TABLE
This level of the directory structure looks to be where genic SNPs are identified and a program called SweeD is run

```bash
subset-genes-gff3-scaffolds1-6.py #extract gene annotation information for chromosomes 1 to 6 and writes this to a new gff3 file
MYZPE13164_O_EIv2.1.annotation.gff3 #input for subset-genes-gff3-scaffolds1-6.py
MYZPE13164_O_EIv2.1.annotation_genes_scaff1-6.gff3 #output of subset-genes-gff3-scaffolds1-6.py
```
Finding genic SNPs:
```bash
cat protocol.md

#Finding genic SNPs:

singularity exec /hpc-home/cheemaj/BUILD/PYBEDTOOLS/pybed.simg bedtools  intersect \
    -a 210s.M_persicae.onlySNPs.vcf.gz \
    -b MYZPE13164_O_EIv2.1.annotation_genes_scaff1-6.gff3  \
    -header > mperc-210-genic-regions.vcf
```
### SweeD
Sweed is to do with genic regions as the input is genic-regions.vcf

SweeD implements a composite likelihood ratio test which detects complete selective sweeps using Site Frequency Spectrum (SFS) patterns of single-nucleotide polymorphisms (SNPs). - The site frequency spectrum summarizes the distribution of allele frequencies throughout the genome, and it is widely used as a summary statistic to infer demographic parameters and to detect signals of natural selection.

SweeD	(Sweep	Detector,	Pavlidis	&	Alachiotis) - https://cme.h-its.org/exelixis/resource/download/software/sweed3.0_manual.pdf
```bash
cat protocol.md

### sweed
##installation
make -f Makefile.PTHREADS.gcc
pwd
/hpc-home/cheemaj/BUILD/PAML/SweeD_v3.2.1_Linux
./SweeD -h

source package /tgac/software/testing/bin/bcftools-1.12
bcftools stats mperc-210-genic-regions.vcf  > mperc-210-genic-regions_stats.txt
bgzip -c file.vcf > file.vcf.gz
# 11,870,914
# 7,216,570

cat resum1.sh #runs SweeD
	/hpc-home/cheemaj/BUILD/PAML/SweeD_v3.2.1_Linux/SweeD   -name mperc-210-genic-regions  -input mperc-210-genic-regions.vcf  -grid 100

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
part-joiner.py #reads medelseq-allele-table-*.tab files and writes them to csv, checking that the no. of columns in correct
```
I think these are files copied from another project and edited versions for this project are one level down in the directory structure as I do not recognise the .vcf file given as an input in vc2alle-part.py or the sample names.
```bash
submitter.py #executes a singularity container to run vc2alle-part.py, the file to run this will be called resum[].sh
	vc2alle-part.py #polymorphic information content calculated
```
#### /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/position-wise-snp-diff/TABLE/BCF

README-manual.md 
```bash
#Again I do not recoginse this .vcf file, is this copied from another project? why is a singularity container needed to run a python script to collect sample names alone?
#---------------------------
bcftools view mperc-210-genic-regions.vcf  -Oz -o mperc-210-genic-regions.vcf.gz
bcftools index -t mperc-210-genic-regions.vcf.gz


singularity exec ~/BUILD/ABEL/abel.simg  python3 get_samples_names.py


       vcf_file = "../SNP.Missing-unphasing.ID.ann.PASS.allele2.DP5_30.MAF001.Missing01.InbreedingCoeff_retain.vcf.gz"
       vcf =  VCF(vcf_file)
#-----------------------
```
```bash
nano get_samples_names.py
#presumably does what it says on the tin
```

```bash
mperc-210-genic-regions.vcf
mperc-210-genic-regions.vcf.gz
mperc-210-genic-regions.vcf.gz.tbi
MYZPE13164_O_EIv2.1.annotation_genes_scaff1-6.gff3
```
Where has the printout at the end of these documents come from?:
```bash
nano perc-test.py
	#input=mperc-210-genic-regions.vcf.gz, MYZPE13164_O_EIv2.1.annotation_genes_scaff1-6.gff3
	#output=gene-hits-snp.tab
	#Python script that reads in a VCF file and extracts information for a subset of samples defined by the samples list and reads in a GFF3 file containing gene annotations for specific genomic regions and checks for the presence of variants within those regions in the selected subset of samples. If a variant is present, it records the gene name, chromosome, position, reference allele, alternate allele, and the number of individuals with a homozygous alternate genotype.
	#Genotype codes are converted to allele codes: 0: 'HOM_REF', 1: 'HET', 3: 'UNKNOWN', 2: 'HOM_VAR'

nano perc-test-sep.py
	#input=mperc-210-genic-regions.vcf.gz, MYZPE13164_O_EIv2.1.annotation_genes_scaff1-6.gff3
	#output=gene-hits-snp-O-vs-rest.tab

nano perc-test-sep-O-vs-rest.py
	#input=mperc-210-genic-regions.vcf.gz, MYZPE13164_O_EIv2.1.annotation_genes_scaff1-6.gff3
	#output=gene-hits-snp-O-vs-rest.tab
```
```bash
ls -l
nano gene-hits-snp.tab
nano gene-hits-snp-O-vs-rest.tab
#gene: the name of the gene where the variant was found
#CHROM: the chromosome where the variant was found
#POS: the position of the variant on the chromosome
#REF: the reference allele at the variant position
#ALT: the alternate allele(s) at the variant position
#alts: the number of homozygous alternate alleles in the samples for this variant
#extra 2 columns in gene-hits-snp-O-vs-rest.tab; combined, group 0, rest
```
```bash
nano samples.py
#Defines a python list of all the sample names - nothing else

nano rest.py
#finds all samples which are not in group 0, (s = element, so script finds all elements in sample list which are not in group 0 ist)


nano scaffold-lengths.txt #must be the output of something, but what?

nano snp-discover-saturation-curve-all-jic-bass.py #must make the discussed curves
	#input = mperc-210-genic-regions.vcf.gz
	#output = snp-alt-richness.tab -> this file isnt in this folder
	#This script seems to be extracting SNP genotypes from a VCF file for a set of target samples, and calculating the number of SNPs in which each sample has a homozygous variant genotype.
```


#### vcf2allele
```bash
nano /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/position-wise-snp-diff/TABLE/BCF/submitter.py
	#output= resum(N).sh -> there are no files of this name in this folder
	#submits individual samples slurm for vc2alle-part.py with custom wrapper.
	#eg. resum11.sh

#import sys
#import os
#import re
#import itertools
#import glob
#import StringIO
#from collections import defaultdict

#def strMUT(text, dic):
#    pat = "(%s)" % "|".join(map(re.escape, dic.keys()))
#    return re.sub(pat, lambda m: dic[m.group()], text)

#template='''#!/bin/bash
##SBATCH --job-name=SER
##SBATCH -o  BASE.out
##SBATCH -e  BASE.err
##SBATCH --mem 40gb
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH -p jic-medium
##SBATCH -t 1-23:00:00


#singularity exec ~/BUILD/ABEL/abel.simg  python3 vc2alle-part.py   --num NUM  --sample  SAMPLE

#echo DONE

#'''

#samples =  ['S98', 'K43', 'K40', 'A166', 'I1', 'A105', 'K16', 'A102', 'NL1', 'Q1200', 'S2B01', 'S2196G', 'NL21', 'S152', 'S232', 'UKW3', 'NL5', 'UK2', 'S204', 'A138', 'Crespys', 'S472', 'S196', 'MG1107', 'S481', 'Lierid$


#kodaf=open('ALL-BATCH.sh','w')
#for t, base in enumerate( samples, 1):
#  #if i not in [4]:## oparse error files
#  #   continue
#  sample = base
#  with open('resum'+ str(t) +'.sh','w') as out:
#    newst=strMUT(template, {'SER': str(t),  'BASE': base, 'SAMPLE': sample, 'NUM' : str(t) })
#    out.write(newst)
#  kodaf.write('sbatch   ' + 'resum'+ str(t) +'.sh' + '\n')

#print 'Done'

```
```bash
nano vc2alle-part.py #appears to be generating an allele table for a single sample from a VCF file
	#input = mperc-210-genic-regions.vcf.gz
	#output = mperg-allele-table-{num}.tab -> there are no files of this name in this folder

	#The output file would be a tab-separated text file with two columns: the first column would contain the input number (num) and the second column would contain the sample name (sample). The remaining columns would contain the allele information extracted from the VCF file.

	#The number of remaining columns would be equal to the number of positions in the VCF file that are within the genic regions (as specified in vcf_file). The values in the remaining columns would be either 0, 2, or an empty string. The 0 indicates a homozygous reference genotype, the 2 indicates a heterozygous or homozygous alternate genotype, and the empty string indicates that the genotype is unknown.

	#The number of rows in the output file would be equal to the number of samples specified in the samples list (which is ['S98', 'K43', 'K40'] in the current implementation). Each row would correspond to one sample, and the columns would contain the allele information for that sample at each position.
```






