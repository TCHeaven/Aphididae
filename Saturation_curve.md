# Exploring the diversity of the aphid samples as part of the M. persicae population genomics project

Contains SNP filtering, network analysis and rarefaction curve plotting.

## Contents
1. [Data quality control](#1)
        1. [Coverage](#2)
        2. [Missingness](#3)
                1. [VCFstats](#4)
2. [Distribution of SNPs](#5)
        1. [Sliding window](#6)
3. [Extract genic SNPs](#7)
        1. [Genic regions](#8)
        2. [CDS/exon/UTR regions](#9)
        3. [CDS regions](#10)
        4. [SNPEff non-synonymous SNPs](#11)
4. [P_distance](#12)
        1. [Whole genome SNPs](#13)
                1. [Splitstree](#14)
        2. [Genic and CDS SNPs](#15)
                1. [Splitstree](#16)
        3. [Non-synonymous SNPs](#17)
                1. [Splitstree](#18)
        4. [Synonymous SNPs](#19)
                1. [Splitstree](#20)
5. [SNP rarefaction curve](#21)
        1. [Corehunter (redundant)](#22)
        2. [True Rarefaction curves](#23)
                1. [Whole genome SNPs](#24)
                2. [Genic SNPs](#25)
                3. [For Bass and JIC samples seperately - Whole genome and genic](#26)
                4. [Fitting a curve](#27)
6. [vcftools plot of saturation curve](#28)
7. [Nei gene diversity](#29)
8. [Check for SNPs of clone O versus itself to control for diploidy](#30)

Unless stated otherwise work was performed from the directory: /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae

SNP calling has previously been performed by R.Wouters and R.biello.

Collecting data:
```bash
mkdir -p snp_calling/Myzus/persicae/biello/gatk/filtered
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/vcffilesafterfiltering/210s.M_persicae.onlySNPs.vcf.gz snp_calling/Myzus/persicae/biello/gatk/filtered/.

ln -s /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3

ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa.gz
```
## Data quality control <a name="1"></a>
The SNP data has already been filtered by R.Wouters and R.biello with vcftools etc.

### Coverage <a name="2"></a>
Coverage of each asample was explored:
```bash
mkdir snp_calling/Myzus/persicae/biello/gatk/reports
#JIC samples
for sample in $(ls /jic/scratch/groups/Saskia-Hogenhout/roland/Myzus_periscae_popgen/Output/realigned_BAM_files/*.sorted.mark_dups.realigned_stats/genome_results.txt); do
    SampleName=$(echo $sample | cut -d '/' -f10| cut -d '.' -f1)
    Coverage=$(grep 'mean coverageData' $sample)
    echo $SampleName $Coverage >> snp_calling/Myzus/persicae/biello/gatk/reports/JIC_coverage.txt
done
```
Some samples have <15X coverage according the Rolands results, but are included in the "final" .vcf file?

### Missingness <a name="3"></a>

The number of SNPs per sample was explored.
```bash
source package /nbi/software/testing/bin/vcftools-0.1.15
vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/209s.M_persicae.onlySNPs.vcf --missing-indv
#After filtering, kept 209 out of 209 Individuals
#Outputting Individual Missingness
#After filtering, kept 11870914 out of a possible 11870914 Sites
mv out.imiss /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/209s.M_persicae.out.imiss

source package /tgac/software/testing/bin/gnuplot-4.6.7
awk '!/IN/' snp_calling/Myzus/persicae/biello/gatk/filtered/209s.M_persicae.out.imiss | cut -f5 > totalmissing

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
rm totalmissing

                                        Histogram of % missing data per individual

  120 ++---------+----------+----------+----------+----------+----------+----------+----------+----------+---------++
      +          +          +          +          +          + 'totalmissing' using (bin($1,binwidth)):(1.0) ****** +
      **                                                                                                            |
      **                                                                                                            |
  100 **                                                                                                           ++
      **                                                                                                            |
      **                                                                                                            |
      **                                                                                                            |
   80 **                                                                                                           ++
      **                                                                                                            |
      **                                                                                                            |
   60 **                                                                                                           ++
      **                                                                                                            |
      **                                                                                                            |
      **                                                                                                            |
   40 ***                                                                                                          ++
      ***                                                                                                           |
      ***                                                                                                           |
      ***                                                                                                           |
   20 ***                                                                                                          ++
      *****                                                                                                         |
      *****                                                                                                         |
      ********** +          +          +          +          +          +          +          +          +          +
    0 ***************************************************************************************************************
      0         0.1        0.2        0.3        0.4        0.5        0.6        0.7        0.8        0.9         1
                                                     % of missing data
```
Plot the missingness level of samples:
```python
import matplotlib.pyplot as plt

# Read the F_MISS data from a file, skipping the header line
data = []
with open('snp_calling/Myzus/persicae/biello/gatk/filtered/209s.M_persicae.out.imiss', 'r') as file:
    next(file)  # Skip the header line
    for line in file:
        values = line.strip().split('\t')
        f_miss = float(values[-1])
        data.append(f_miss)

# Create the scatter plot
plt.scatter(range(len(data)), data, marker='o', color='blue')
plt.xlabel('Individual')
plt.ylabel('Missing Level (F_MISS)')
plt.title('Missing Level per Individual')
plt.savefig('snp_calling/Myzus/persicae/biello/gatk/filtered/209s_output_plot.png')
```
Remove samples with more than 10% of SNPs missing from the vcf:
```bash
awk -F'\t' '$NF > 0.1' snp_calling/Myzus/persicae/biello/gatk/filtered/209s.M_persicae.out.imiss
#INDV    N_DATA  N_GENOTYPES_FILTERED    N_MISS  F_MISS
#S98     11870914        0       3140218     0.26453
#A105    11870914        0       11848667    0.998126
#MISC47  11870914        0       8365879     0.704738
#UK20    11870914        0       1417373     0.119399
#A151    11870914        0       2908455     0.245007
#A156    11870914        0       4558493     0.384005
#S55     11870914        0       3173971     0.267374
#S95     11870914        0       1738082     0.146415
#S99     11870914        0       1516941     0.127786
#S123    11870914        0       2092827     0.176299
#S119    11870914        0       5572471     0.469422
#S104    11870914        0       1679059     0.141443
#S102    11870914        0       4697437     0.39571
#S97     11870914        0       4556016     0.383797
#S103    11870914        0       5892768     0.496404
#S116    11870914        0       2300217     0.193769

awk -F'\t' '$NF > 0.1' snp_calling/Myzus/persicae/biello/gatk/filtered/209s.M_persicae.out.imiss | wc -l
#17

patterns=("ligustri" "S98" "A105" "MISC47" "UK20" "A151" "A156" "S55" "S95" "S99" "S123" "S119" "S104" "S102" "S97" "S103" "S116")
for pattern in "${patterns[@]}"; do
    cat snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.vcf | awk -v pattern="$pattern" -F'\t' '{for (i=1; i<=NF; i++) if ($i ~ pattern) {print i; exit}}'
done


cut -f219,10,15,65,70,88,89,174,186,189,198,199,204,207,211,213,217 --complement 210s.M_persicae.onlySNPs.vcf > 193s.M_persicae.onlySNPs.vcf
head -n 396 192s.M_persicae.onlySNPs.vcf
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae

patterns=("ligustri" "S98" "A105" "MISC47" "UK20" "A151" "A156" "S55" "S95" "S99" "S123" "S119" "S104" "S102" "S97" "S103" "S116")
columns_to_remove=()

for pattern in "${patterns[@]}"; do
    column=$(cut -f1- --complement snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.vcf | awk -v pattern="$pattern" -F'\t' '{for (i=1; i<=NF; i++) if ($i ~ pattern) {print i; exit}}')
    columns_to_remove+=("$column")
done

cut -f"${columns_to_remove[*]}" --complement snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.vcf > snp_calling/Myzus/persicae/biello/gatk/filtered/193.M_persicae.onlySNPs.vcf


gzip -cd /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs.vcf.gz > snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.vcf
cut -f10,15,65,70,88,89,174,186,189,198,199,204,207,211,213,217 --complement snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.vcf > snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs.vcf

zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz | wc -l #11,871,310
```
Remove SNPs which have no minor allele now that the high missingness samples have been removed:
```bash
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz --mac 1 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1 #After filtering, kept 1,851,865 out of a possible 11,870,914 Sites
```
#### VCFstats <a name="4"></a>
```bash
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz); do
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity
    OutFile=193_genic
    Exclusion_list=
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_vcfstats.sh $InFile $OutDir $OutFile $Exclusion_list
done #55244504
```
Plot VCFStats
```R
setwd("C:/Users/did23faz/OneDrive - Norwich Bioscience Institutes/Desktop/R")

install.packages("tidyverse")
install.packages("ggplot2")

library(tidyverse)
library(ggplot2)

var_qual <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_site_qual.lqual", delim = "\t",
                       col_names = c("CHROM", "POS", "QUAL"), skip = 1)

var_depth <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_site-mean-depth.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

var_miss <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_missing_sites.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

var_freq <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_allele_freq.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

var_count <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_allele_count.frq.count", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "count", "a2"), skip = 1)

ind_depth <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_depth.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

ind_miss  <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_missing_ind.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

ind_het <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_het.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)


p1 <- ggplot(var_qual, aes(QUAL)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Per Site Quality")
p1 <- p1 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_var_qual.pdf", width=11, height=7)
plot(p1, pdf=T)
dev.off()

p2 <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + labs(title = "Genic SNPs: Mean Depth Per Site")
p2 <- p2 + theme(plot.title = element_text(hjust = 0.5))
summary(var_depth$mean_depth)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  8.249  16.212  19.529  18.710  21.285  30.829 
pdf("193_genic_var_depth.pdf", width=11, height=7)
plot(p2, pdf=T)
dev.off()

p3 <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + xlim(0, 100) + labs(title = "Genic SNPs: Mean Depth Per Site")
p3 <- p3 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_var_depth_x100.pdf", width=11, height=7)
plot(p3, pdf=T)
dev.off()

p4 <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Missingness Per Site")
p4 <- p4 + theme(plot.title = element_text(hjust = 0.5))
summary(var_miss$fmiss)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.005181 0.010363 0.016646 0.025907 0.108808 
pdf("193_genic_var_miss.pdf", width=11, height=7)
plot(p4, pdf=T)
dev.off()

var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
summary(var_freq$maf)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.00000 0.00000 0.01121 0.00000 0.50000 
p5 <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Allele Frequency for Each Site")
p5 <- p5 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_var_freq.pdf", width=11, height=7)
plot(p5, pdf=T)
dev.off()
p51 <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Allele Frequency for Each Site") + xlim(0, 0.1)
p51 <- p51 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_var_freq_0.1.pdf", width=11, height=7)
plot(p51, pdf=T)
dev.off()

p6 <- ggplot(var_count, aes(x = count)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Raw Allele Counts for Each Site")
p6 <- p6 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_allele_count_density.pdf", width = 11, height = 7)
plot(p6, pdf=T)
dev.off()

p7 <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Mean Depth Per Individual")
p7 <- p7 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_ind_depth.pdf", width=11, height=7)
plot(p7, pdf=T)
dev.off()

p8 <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Missingness Per Individual")
p8 <- p8 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_ind_miss.pdf", width=11, height=7)
plot(p8, pdf=T)
dev.off()

p9 <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Heterozygosity Per Individual")
p9 <- p9 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_ind_het.pdf", width=11, height=7)
plot(p9, pdf=T)
dev.off()
```
## Distribution of SNPs <a name="5"></a>

Find the number of SNPs per chromosome:
```bash
vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf --positions scaff1-6.bed --window-pi 10000 --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.snp_counts

bcftools query -f '%CHROM\t%POS\n' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf | \
    awk -F '\t' '{if ($1 ~ /^scaffold_[1-6]$/) count[$1]++;} END {for (i=1; i<=6; i++) print "scaffold_"i ": " count["scaffold_"i] " SNPs";}'

nano scaff1-6.bed
##CHROM  START   END
#scaffold1   1   105178091
#scaffold2   1   86073209
#scaffold3   1   69480500
#scaffold4   1   62328371
#scaffold5   1   30612500
#scaffold6   1   29865500
```
#### Sliding window <a name="6"></a>

Investigate the distribution of SNPs across the M. persicae genome, by counting the number of SNPs in sliding windows and plotting:
```python
import pysam
import matplotlib.pyplot as plt

# Define the window size
window_size = 10000
window_size = 50000
window_size = 100000

# Specify the output file path
output_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_10000.txt'
output_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_50000.txt'
output_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_100000.txt'

# Output image file name
output_image_template = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_10000_scaffold_{}_plot.png'
output_image_template = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_50000_scaffold_{}_plot.png'
output_image_template = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_100000_scaffold_{}_plot.png'

# Initialize a dictionary to store SNP counts for each scaffold
scaffold_snp_counts = {}

# Open the VCF file
with pysam.VariantFile('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf') as vcf_file:
    for record in vcf_file:
        scaffold = record.chrom  # Get the scaffold name
        pos = record.pos  # Get the position of the SNP
        # Check if the scaffold is in the range 1 to 6
        if scaffold.startswith('scaffold') and int(scaffold.split('_')[1]) in range(1, 7):
            if scaffold not in scaffold_snp_counts:
                scaffold_snp_counts[scaffold] = []
            # Calculate the window index for the current SNP
            window_index = pos // window_size
            if window_index >= len(scaffold_snp_counts[scaffold]):
                scaffold_snp_counts[scaffold].extend([0] * (window_index - len(scaffold_snp_counts[scaffold]) + 1))
            # Increment the SNP count for the corresponding window
            scaffold_snp_counts[scaffold][window_index] += 1

# Write SNP counts in sliding windows to the output file
with open(output_file, 'w') as outfile:
    for scaffold in [f'scaffold_{i}' for i in range(1, 7)]:
        if scaffold in scaffold_snp_counts:
            snp_counts = scaffold_snp_counts[scaffold]
            for i, count in enumerate(snp_counts):
                start = i * window_size
                end = start + window_size
                outfile.write(f'{scaffold}, Window {i+1} ({start}-{end} bp): {count} SNPs\n')

# Plot SNP counts in sliding windows for each scaffold and save the plots
for scaffold in [f'scaffold_{i}' for i in range(1, 7)]:
    if scaffold in scaffold_snp_counts:
        snp_counts = scaffold_snp_counts[scaffold]
        window_positions = [i * window_size for i in range(len(snp_counts))] 
        plt.figure(figsize=(20, 5))
        plt.plot(window_positions, snp_counts, marker='.', linestyle='', markersize=1)
        plt.xlabel('Window Position (bp)')
        plt.ylabel('Number of SNPs')
        plt.title(f'Scaffold {scaffold}')
        plt.grid(True)
        # Save the plot as an image file
        output_image = output_image_template.format(scaffold.split('_')[1])
        plt.savefig(output_image)
        plt.close()
```
Number of SNPs in sliding windows omitting windows with poor mappability
```bash
#Mappability masking is looking at individual nucleotides
WD=/jic/scratch/groups/Saskia-Hogenhout/roberto/m_persicae/popgen/VCF_filtering/mappability

cd $WD

bedtools intersect -a 212ind.M_persicae.bed -b Myzus_persicae_O.v2.gem.150_004.filt.tab.bed > Myzus_persicae_O.v2.mappable.bed

bedtools subtract -a Myzus_persicae_O.v2.mappable.bed -b Myzus_persicae_O.v2.masked.filt.tab.out.bed > Myzus_persicae_O.v2.callable.bed

#212ind.M_persicae.bed - bed file with all sites in vcf
#Myzus_persicae_O.v2.gem.150_004.filt.tab.bed - sites that have been mapped
#Myzus_persicae_O.v2.masked.filt.tab.out.bed intervals to subtract - poor mappability
#Myzus_persicae_O.v2.callable.bed - sites that have good mappability and can be kept

cp /jic/scratch/groups/Saskia-Hogenhout/roberto/m_persicae/popgen/VCF_filtering/mappability/Myzus_persicae_O.v2.masked.filt.tab.out.bed /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/.
```
Total length of poor mapping regions in Myzus_persicae_O.v2.masked.filt.tab.out.bed is 93,097,267 bp -error? this is the no. of bases masked from repeatmodeler
```python
import pysam
import matplotlib.pyplot as plt

# Define the window size
window_size = 100000

# Specify the output image file template
output_image_template = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_100000_<5_scaffold_{}_plot.png'

# Read the BED file containing regions to be excluded
exclude_regions = []
with open('/jic/scratch/groups/Saskia-Hogenhout/roberto/m_persicae/popgen/VCF_filtering/mappability/Myzus_persicae_O.v2.masked.filt.tab.out.bed', 'r') as bed_file:
    for line in bed_file:
        fields = line.strip().split('\t')
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        exclude_regions.append((chrom, start, end))

# Initialize a dictionary to store SNP counts for each scaffold
scaffold_snp_counts = {}

# Open the VCF file
with pysam.VariantFile('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf') as vcf_file:
    for record in vcf_file:
        scaffold = record.chrom  # Get the scaffold name
        pos = record.pos  # Get the position of the SNP
        # Check if the scaffold is in the range 1 to 6
        if scaffold.startswith('scaffold') and int(scaffold.split('_')[1]) in range(1, 7):
            if scaffold not in scaffold_snp_counts:
                scaffold_snp_counts[scaffold] = []
            # Calculate the window index for the current SNP
            window_index = pos // window_size
            if window_index >= len(scaffold_snp_counts[scaffold]):
                scaffold_snp_counts[scaffold].extend([0] * (window_index - len(scaffold_snp_counts[scaffold]) + 1))
            # Increment the SNP count for the corresponding window
            scaffold_snp_counts[scaffold][window_index] += 1

# Function to check if a window overlaps with excluded regions with more than 10% overlap
def is_overlapping_with_exclude_regions(start, end, regions, overlap_threshold=0.05):
    window_size = end - start
    for _, region_start, region_end in regions:
        overlap = max(0, min(end, region_end) - max(start, region_start))
        overlap_percentage = overlap / window_size
        if overlap_percentage > overlap_threshold:
            return True
    return False

# Plot SNP counts in sliding windows for each scaffold and save the plots
for scaffold in [f'scaffold_{i}' for i in range(1, 7)]:
    if scaffold in scaffold_snp_counts:
        snp_counts = scaffold_snp_counts[scaffold]
        window_positions = [i * window_size for i in range(len(snp_counts))]
        # Filter out windows overlapping with excluded regions
        filtered_window_positions = []
        filtered_snp_counts = []
        for i, start in enumerate(window_positions):
            end = start + window_size
            if not is_overlapping_with_exclude_regions(start, end, exclude_regions):
                filtered_window_positions.append(start)
                filtered_snp_counts.append(snp_counts[i])     
        plt.figure(figsize=(20, 5))
        plt.plot(filtered_window_positions, filtered_snp_counts, marker='.', linestyle='', markersize=4)
        plt.xlabel('Window Position (bp)')
        plt.ylabel('Number of SNPs')
        plt.title(f'Scaffold {scaffold} <5% overlap with poor mapping regions')
        plt.grid(True)
        # Save the plot as an image file
        output_image = output_image_template.format(scaffold.split('_')[1])
        plt.savefig(output_image)
        plt.close()
```

## Extract genic SNPs <a name="7"></a>

#### Genic regions <a name="8"></a>

Extract genic SNPs to new .vcf files:
```bash
wc -l snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3
#902,603

#the reduced set of 193 SNPs after removing high missingness samples:
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 \
-header > snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf

#the reduced set of 193 SNPs after removing high missingness samples, + M. ligustri:
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 \
-header > snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-genic-regions.vcf

#The original 210 samples from wouters+biello:
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 \
-header > /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf

source package 01ef5a53-c149-4c9e-b07d-0b9a46176cc0
bgzip -c snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf > snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz
zcat snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz | wc -l #31,637,658 there are many duplicate SNPs
bgzip -c /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf > /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf.gz
zcat /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf.gz | wc -l #31,637,658

#deduplicate files with:
Files=("snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz" "/jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf.gz")
for file in "${Files[@]}"; do
    singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 /hpc-home/did23faz/git_repos/Scripts/NBI/remove_duplicate_snps.py $file 
done

#check intersect of original and deduplicated files: 
bcftools isec  -p temp_dir /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf.gz /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf.gz_2 
#CONFIRMS THERE WERE DUPLICATED ENTRIES IN THE ORIGINAL AND THAT THESE HAVE BEEN SUCESSFULLY REMOVED, WHILST RETAINING ONE COPY

for file in "${Files[@]}"; do
    mv ${file}_2 $file 
done

zcat snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz | wc -l #6,783,065 -> 6,782,669 
zcat /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf.gz | wc -l #6,783,065

vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz --mac 1 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-mac1-regions #After filtering, kept 944,154 out of a possible 6,782,669 Sites

awk '{ if ($3 == "gene") sum += $5 - $4 + 1 } END { print sum }' snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 #247,120,340 = the number of base pairs covered by gene features - possibly an overestimate as some may overlap?
```
Gene regions cover over half of the genome and over half of SNPs are retained when extracting genic SNPs - as would be expected.

#### CDS/exon/UTR regions <a name="9"></a>

Extract CDS/exon/UTR regions region SNPs:
```bash
grep '#\|CDS\|exon\|five_prime_UTR\|three_prime_UTR' snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 > snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_exon_UTR_annotation.gff3
wc -l snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_exon_UTR_annotation.gff3
#806,274

#the reduced set of 193 SNPs after removing high missingness samples:
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_exon_UTR_annotation.gff3 \
-header > snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf

#the reduced set of 193 SNPs after removing high missingness samples, + M. ligustri:
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_exon_UTR_annotation.gff3 \
-header > snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf

#The original 210 samples from wouters+biello:
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_exon_UTR_annotation.gff3 \
-header > /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf

source package 01ef5a53-c149-4c9e-b07d-0b9a46176cc0
bgzip -c snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf > snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz
bgzip -c snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf > snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz
bgzip -c /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf > /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz
zcat snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz | wc -l #8,581,242 many duplicates
zcat snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz | wc -l #8,581,242 many duplicates
zcat /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz | wc -l #8,581,242 many duplicates

rm snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf 

#deduplicate files with:
Files=("snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz" "snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz" "/jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz")
for file in "${Files[@]}"; do
    singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 /hpc-home/did23faz/git_repos/Scripts/NBI/remove_duplicate_snps.py $file 
done

#check intersect of original and deduplicated files: 
bcftools isec  -p temp_dir /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz_2 
#CONFIRMS THERE WERE DUPLICATED ENTRIES IN THE ORIGINAL AND THAT THESE HAVE BEEN SUCESSFULLY REMOVED, WHILST RETAINING ONE COPY

for file in "${Files[@]}"; do
    mv ${file}_2 $file 
done

zcat snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz | wc -l #1,925,747
zcat snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz | wc -l #1,925,747
zcat /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf.gz | wc -l #1,925,747
```
#### CDS regions <a name="10"></a>

Extract CDS region SNPs:
```bash
grep '#\|CDS' snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 > snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.gff3
wc -l snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.gff3
#336,709

#the reduced set of 193 SNPs after removing high missingness samples:
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.gff3 \
-header > snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic-regions.vcf

#the reduced set of 193 SNPs after removing high missingness samples, + M. ligustri:
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.gff3 \
-header > snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_genic-regions.vcf

#The original 210 samples from wouters+biello:
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.gff3 \
-header > /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_genic-regions.vcf

source package 01ef5a53-c149-4c9e-b07d-0b9a46176cc0
bgzip -c snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic-regions.vcf > snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz
bgzip -c snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_genic-regions.vcf > snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz
bgzip -c /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_genic-regions.vcf > /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz

zcat snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz | wc -l #3,007,451 there do not appear to be duplicates, but I dnt know if I trust this now...
zcat snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz | wc -l #3,007,451
zcat /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz | wc -l #3,007,451

rm snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic-regions.vcf snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_genic-regions.vcf /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_genic-regions.vcf

#deduplicate files with:
Files=("snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz" "snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz" "/jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz")
for file in "${Files[@]}"; do
    singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 /hpc-home/did23faz/git_repos/Scripts/NBI/remove_duplicate_snps.py $file 
done

#check intersect of original and deduplicated files: 
bcftools isec  -p temp_dir /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz_2 
#CONFIRMS THERE WERE DUPLICATED ENTRIES IN THE ORIGINAL AND THAT THESE HAVE BEEN SUCESSFULLY REMOVED, WHILST RETAINING ONE COPY

for file in "${Files[@]}"; do
    mv ${file}_2 $file 
done

zcat snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz | wc -l #1,240,584 -> 1,240,188
zcat snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz | wc -l #1,240,584
zcat /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz | wc -l #1,240,584

#REMOVE SNPS THAT HAVE NO MINOR ALLELE COUNT NOW THAT HIGH MISSINGNESS SAMPLES HAVE BEEN REMOVED:
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz --mac 1 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions #After filtering, kept 126,162 out of a possible 1,240,188 Sites
bgzip snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.vcf
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz --mac 1 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_genic_mac1-regions #After filtering, kept 406,201 out of a possible 1,240,188 Sites
bgzip snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.vcf

awk '{ if ($3 == "CDS") sum += $5 - $4 + 1 } END { print sum }' snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.gff3 #69,433,195 = the number of base pairs covered by CDS features - possibly an overestimate as some may overlap?

#Check how many CDS genes have SNPs:
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.gff3 \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.vcf.gz \
-header > /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/MYZPE13164_O_EIv2.1.CDS_annotation_withsnpsonly.gff3
awk -F "\t" '{print $9}' /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/MYZPE13164_O_EIv2.1.CDS_annotation_withsnpsonly.gff3 | cut -d ';' -f2 | sort | uniq | wc -l #33,428
```
Many CDS contain SNPs

#### SNPEff non-synonymous SNPs <a name="11"></a>

Extract Non-synonymous CDS SNPs using Variant effect predictor (VEP):
```bash
source package /tgac/software/testing/bin/bedops-2.2.0
source package 4028d6e4-21a8-45ec-8545-90e4ed7e1a64
source package /tgac/software/production/bin/tabix-0.2.6
gff2bed < snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.gff3 > snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.bed

gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3
vcf=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.vcf.gz
bed=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.bed
genome=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa

grep -v "#" $gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > myData.gff.gz
tabix -p gff myData.gff.gz

grep -v "#" $bed | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > myData.bed.gz
tabix -p bed myData.bed.gz

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/vep.sif vep -i $vcf --gff myData.gff.gz --fasta $genome --vcf --output_file variant_effect_output.txt2

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/vep.sif vep -i $vcf --cache --output_file variant_effect_output.txt
/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/vep.sif

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/vep.sif vep --custom file=Filename,short_name=Short_name,format=File_type,type=Annotation_type,coords=Force_report_coordinates,fields=VCF_fields
```
Extract Non-synonymous CDS SNPs using and SNPEff:
```bash
source package 3e7beb4d-f08b-4d6b-9b6a-f99cc91a38f9
nano ~/snpEff/snpEff.config
## Myzus persicae genome, version O2_0
#m_persicae_O_2_0.genome : Peach_potato_aphid

mkdir -p ~/snpEff/data/m_persicae_O_2_0
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 ~/snpEff/data/m_persicae_O_2_0/genes.gff
cp /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.gtf ~/snpEff/data/m_persicae_O_2_0/genes.gtf
cp /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.aa.fa ~/snpEff/data/m_persicae_O_2_0/protein.fa
cp /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.CDS.fa ~/snpEff/data/m_persicae_O_2_0/cds.fa
mkdir -p ~/snpEff/data/genomes
cp /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa ~/snpEff/data/genomes/m_persicae_O_2_0.fa
cd ~/snpEff/data/m_persicae_O_2_0/
gzip genes.gff
gzip genes.gtf
cd ../..
java17 -jar snpEff.jar build -gff3 -v m_persicae_O_2_0
java17 -jar snpEff.jar build -gtf22 -v m_persicae_O_2_0
#56982461
#56982568
#56989016

java17 -Xmx8g -jar snpEff.jar m_persicae_O_2_0 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.vcf.gz > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.ann.vcf
#56989169

cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.ann.vcf | grep -v 'synonymous_variant' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions-nonsyn.recode.ann.vcf 
wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions-nonsyn.recode.ann.vcf #57,940 -> 57,544 non-synonymous SNPs
```
There are 57,544 non-synonymous SNPs

## P_distance <a name="12"></a>
A distance matrix was calculated for the samples.

p-distance  is the proportion (p) of nucleotide sites at which two sequences being compared are different. It is obtained by dividing the number of nucleotide differences by the total number of nucleotides compared.

distmat calculates the evolutionary distance between every pair of sequences in a multiple sequence alignment. A variety of methods to estimate distance may be selected, and differ in how they correct the observed substitution rates to more accurately reflect the true evolutionary distance. An output file containing a distance matrix for the set of sequences is written. The distances are expressed in terms of the number of substitutions per 100 bases or amino acids.

Generated for with and without ligustri as some analysis require an outgroup and others do not.
### Whole genome SNPs <a name="13"></a>

Generate a distance matrix from .vcf file:
```bash
source package /nbi/software/testing/bin/bcftools-1.8
interactive
bcftools stats snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.vcf.gz > snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs_stats.txt

#The original 210 samples from wouters+biello:
mkdir snp_calling/Myzus/persicae/biello/gatk/p_distance
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.vcf.gz); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=p_dis.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#Submitted batch job 54236748

#The original 210 samples from wouters+biello, without ligustri:
cd snp_calling/Myzus/persicae/biello/gatk/filtered
gunzip -c 210s.M_persicae.onlySNPs.vcf.gz > 210s.M_persicae.onlySNPs.vcf
cat 210s.M_persicae.onlySNPs.vcf | awk -F'\t' '{for (i=1; i<=NF; i++) if ($i ~ /ligustri/) {print i; exit}}'
#219
cut -f219 --complement 210s.M_persicae.onlySNPs.vcf > 209s.M_persicae.onlySNPs.vcf
head -n 396 209s.M_persicae.onlySNPs.vcf
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae

for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/209s.M_persicae.onlySNPs.vcf); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=p_dis_209.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#Submitted batch job 54322980, 54349296,54351157
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered
gunzip 209s.M_persicae.onlySNPs.vcf.gz 

#the reduced set of 193 SNPs after removing high missingness samples:
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=p_dis_193.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#54355967

#the reduced set of 193 SNPs after removing high missingness samples + ligustri:
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs.vcf.gz); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=p_dis_193+ligustri.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#54435632
```
Generate an ID based p-distmat input for R
```bash
source package /nbi/software/production/bin/python-2.7.11
python
```
```python
from collections import defaultdict
# load distance matrix
acc = defaultdict(list)
acc_order = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis.mat') as inp:
    next(inp)
    for line in inp:
        A = line.strip().split()
        acc[A[0]] = A[1:]
        acc_order.append(A[0])

assert len(acc) == 210

# write the distance matrix in a desired format
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_mperc.csv", "w") as outf:
    outf.write(",".join(["ID"] + acc_order) + "\n")
    for a in acc_order:
        outf.write(",".join([a] + acc[a]) + "\n")

print("Done")

############################################################################################

from collections import defaultdict
# load distance matrix
acc = defaultdict(list)
acc_order = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193.mat') as inp:
    next(inp)
    for line in inp:
        A = line.strip().split()
        acc[A[0]] = A[1:]
        acc_order.append(A[0])

assert len(acc) == 193

# write the distance matrix in a desired format
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_mperc.csv", "w") as outf:
    outf.write(",".join(["ID"] + acc_order) + "\n")
    for a in acc_order:
        outf.write(",".join([a] + acc[a]) + "\n")

print("Done")

################################################################################################

from collections import defaultdict
# load distance matrix
acc = defaultdict(list)
acc_order = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193+ligustri.mat') as inp:
    next(inp)
    for line in inp:
        A = line.strip().split()
        acc[A[0]] = A[1:]
        acc_order.append(A[0])

assert len(acc) == 194

# write the distance matrix in a desired format
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193+ligustri_mperc.csv", "w") as outf:
    outf.write(",".join(["ID"] + acc_order) + "\n")
    for a in acc_order:
        outf.write(",".join([a] + acc[a]) + "\n")

print("Done")


exit()
```
Reformat:
```python
import csv

def read_distance_matrix_from_csv(csv_file):
    distance_matrix = {}
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        taxa_names = reader.fieldnames[1:]
        for row in reader:
            taxon = row['ID']
            distances = {taxa_names[i]: float(row[taxa_names[i]]) for i in range(len(taxa_names))}
            distance_matrix[taxon] = distances
    return distance_matrix

def convert_to_phylip(distance_matrix, output_file):
    num_taxa = len(distance_matrix)
    taxa_names = list(distance_matrix.keys())
    with open(output_file, 'w') as f:
        f.write(str(num_taxa) + '\n')
        for i in range(num_taxa):
            taxon_distances = []
            for j in range(num_taxa):
                taxon_distances.append(str(distance_matrix[taxa_names[i]][taxa_names[j]]))
            f.write(taxa_names[i] + '\t' + '\t'.join(taxon_distances) + '\n')

csv_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193+ligustri_mperc.csv'
output_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193+ligustri_mperc.dist'

distance_matrix = read_distance_matrix_from_csv(csv_file)
convert_to_phylip(distance_matrix, output_file)
```
#### Splitstree <a name="14"></a>

Construct a phylogenetic network for the samples, based upon the SNP data.
```bash
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree
```

### Genic and CDS SNPs <a name="15"></a>

Generate a distance matrix from .vcf file, for genic and CDS SNPs:
```bash
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=p_dis_193_genic.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#54584280

for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-genic-regions.vcf.gz); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=p_dis_193+ligustri_genic.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#54584287


#True CDS only SNPS:
bcftools stats /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.vcf.gz > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions_stats.txt
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.vcf.gz); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=p_dis_193CDS.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#56697394

for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.vcf.gz); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=p_dis_193+ligustri_CDS.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#56700915

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
```
Generate an ID based p-distmat input for R
```bash
source package /nbi/software/production/bin/python-2.7.11
python
```
```python
from collections import defaultdict
# load distance matrix
acc = defaultdict(list)
acc_order = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_genic.mat') as inp:
    next(inp)
    for line in inp:
        A = line.strip().split()
        acc[A[0]] = A[1:]
        acc_order.append(A[0])

assert len(acc) == 193

# write the distance matrix in a desired format
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_genic.csv", "w") as outf:
    outf.write(",".join(["ID"] + acc_order) + "\n")
    for a in acc_order:
        outf.write(",".join([a] + acc[a]) + "\n")

# load distance matrix
acc = defaultdict(list)
acc_order = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193+ligustri_genic.mat') as inp:
    next(inp)
    for line in inp:
        A = line.strip().split()
        acc[A[0]] = A[1:]
        acc_order.append(A[0])

assert len(acc) == 194
# write the distance matrix in a desired format
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193+ligustri_genic.csv", "w") as outf:
    outf.write(",".join(["ID"] + acc_order) + "\n")
    for a in acc_order:
        outf.write(",".join([a] + acc[a]) + "\n")

print("Done")

#################################################################################################

from collections import defaultdict
# load distance matrix
acc = defaultdict(list)
acc_order = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193CDS.mat') as inp:
    next(inp)
    for line in inp:
        A = line.strip().split()
        acc[A[0]] = A[1:]
        acc_order.append(A[0])

assert len(acc) == 194

# write the distance matrix in a desired format
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193CDS.csv", "w") as outf:
    outf.write(",".join(["ID"] + acc_order) + "\n")
    for a in acc_order:
        outf.write(",".join([a] + acc[a]) + "\n")

print("Done")

# load distance matrix
acc = defaultdict(list)
acc_order = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193+ligustri_CDS.mat') as inp:
    next(inp)
    for line in inp:
        A = line.strip().split()
        acc[A[0]] = A[1:]
        acc_order.append(A[0])

assert len(acc) == 194

# write the distance matrix in a desired format
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193+ligustri_CDS.csv", "w") as outf:
    outf.write(",".join(["ID"] + acc_order) + "\n")
    for a in acc_order:
        outf.write(",".join([a] + acc[a]) + "\n")

print("Done")

exit()
```
Reformat:
```python
import csv

def read_distance_matrix_from_csv(csv_file):
    distance_matrix = {}
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        taxa_names = reader.fieldnames[1:]
        for row in reader:
            taxon = row['ID']
            distances = {taxa_names[i]: float(row[taxa_names[i]]) for i in range(len(taxa_names))}
            distance_matrix[taxon] = distances
    return distance_matrix

def convert_to_phylip(distance_matrix, output_file):
    num_taxa = len(distance_matrix)
    taxa_names = list(distance_matrix.keys())
    with open(output_file, 'w') as f:
        f.write(str(num_taxa) + '\n')
        for i in range(num_taxa):
            taxon_distances = []
            for j in range(num_taxa):
                taxon_distances.append(str(distance_matrix[taxa_names[i]][taxa_names[j]]))
            f.write(taxa_names[i] + '\t' + '\t'.join(taxon_distances) + '\n')

csv_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193CDS.csv'
output_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193CDS.dist'

distance_matrix = read_distance_matrix_from_csv(csv_file)
convert_to_phylip(distance_matrix, output_file)
```
#### Splitstree <a name="16"></a>

Construct a phylogenetic network for the samples, based upon the genic and CDS SNP data.
```bash
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree
```
### Non-synonymous SNPs <a name="17"></a>

Generate a distance matrix from .vcf file, for non-synonymous SNPs:
```bash
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions-nonsyn.recode.ann.vcf.gz); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=p_dis_193_nonsyn.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#57060954

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
```
Generate an ID based p-distmat input for R
```bash
source package /nbi/software/production/bin/python-2.7.11
python
```
```python
from collections import defaultdict
# load distance matrix
acc = defaultdict(list)
acc_order = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_nonsyn.mat') as inp:
    next(inp)
    for line in inp:
        A = line.strip().split()
        acc[A[0]] = A[1:]
        acc_order.append(A[0])

assert len(acc) == 193

# write the distance matrix in a desired format
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_nonsyn.csv", "w") as outf:
    outf.write(",".join(["ID"] + acc_order) + "\n")
    for a in acc_order:
        outf.write(",".join([a] + acc[a]) + "\n")
```
Reformat:
```python
import csv

def read_distance_matrix_from_csv(csv_file):
    distance_matrix = {}
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        taxa_names = reader.fieldnames[1:]
        for row in reader:
            taxon = row['ID']
            distances = {taxa_names[i]: float(row[taxa_names[i]]) for i in range(len(taxa_names))}
            distance_matrix[taxon] = distances
    return distance_matrix

def convert_to_phylip(distance_matrix, output_file):
    num_taxa = len(distance_matrix)
    taxa_names = list(distance_matrix.keys())
    with open(output_file, 'w') as f:
        f.write(str(num_taxa) + '\n')
        for i in range(num_taxa):
            taxon_distances = []
            for j in range(num_taxa):
                taxon_distances.append(str(distance_matrix[taxa_names[i]][taxa_names[j]]))
            f.write(taxa_names[i] + '\t' + '\t'.join(taxon_distances) + '\n')

csv_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_nonsyn.csv'
output_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_nonsyn.dist'

distance_matrix = read_distance_matrix_from_csv(csv_file)
convert_to_phylip(distance_matrix, output_file)

exit()
```
#### Splitstree <a name="18"></a>

Construct a phylogenetic network for the samples, based upon the non-synonymous SNP data.
```bash
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree
```
### Synonymous SNPs <a name="19"></a>

Extract synonymous SNPs to a seperate .vcf file:
```bash
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions-nonsyn.recode.ann.vcf.gz
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.vcf.gz
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf

# Assign unique identifiers and create compressed VCF files with identifiers
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -o file1_with_ids.vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -o file2_with_ids.vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions-nonsyn.recode.ann.vcf.gz
bgzip file1_with_ids.vcf
bgzip file2_with_ids.vcf
bcftools index file1_with_ids.vcf.gz
bcftools index file2_with_ids.vcf.gz

# Use bcftools isec to find variants unique to file1.vcf.gz
bcftools isec -n=1 -c none -o unique_variants -p output_directory file1_with_ids.vcf.gz file2_with_ids.vcf.gz
mv output_directory/0000.vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1-nononsyn.recode.vcf
rm -r output_directory
bgzip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1-nononsyn.recode.vcf

###########################################################################################################################################################
cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.ann.vcf | grep 'synonymous_variant\|#' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions-syn.recode.ann.vcf 
bgzip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions-syn.recode.ann.vcf 
```
Generate a distance matrix from .vcf file, for synonymous SNPs:
```bash
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1-nononsyn.recode.vcf.gz); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=p_dis_193_nononsyn.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done #57111795

for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions-syn.recode.ann.vcf.gz); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=p_dis_193_syn.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#57109676
```
Generate an ID based p-distmat input for R
```bash
source package /nbi/software/production/bin/python-2.7.11
python
```
```python
from collections import defaultdict
# load distance matrix
acc = defaultdict(list)
acc_order = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_nononsyn.mat') as inp:
    next(inp)
    for line in inp:
        A = line.strip().split()
        acc[A[0]] = A[1:]
        acc_order.append(A[0])

assert len(acc) == 193

# write the distance matrix in a desired format
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_nononsyn.csv", "w") as outf:
    outf.write(",".join(["ID"] + acc_order) + "\n")
    for a in acc_order:
        outf.write(",".join([a] + acc[a]) + "\n")


# load distance matrix
acc = defaultdict(list)
acc_order = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_syn.mat') as inp:
    next(inp)
    for line in inp:
        A = line.strip().split()
        acc[A[0]] = A[1:]
        acc_order.append(A[0])

assert len(acc) == 193

# write the distance matrix in a desired format
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_syn.csv", "w") as outf:
    outf.write(",".join(["ID"] + acc_order) + "\n")
    for a in acc_order:
        outf.write(",".join([a] + acc[a]) + "\n")
```
Reformat:
```python
import csv

def read_distance_matrix_from_csv(csv_file):
    distance_matrix = {}
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        taxa_names = reader.fieldnames[1:]
        for row in reader:
            taxon = row['ID']
            distances = {taxa_names[i]: float(row[taxa_names[i]]) for i in range(len(taxa_names))}
            distance_matrix[taxon] = distances
    return distance_matrix

def convert_to_phylip(distance_matrix, output_file):
    num_taxa = len(distance_matrix)
    taxa_names = list(distance_matrix.keys())
    with open(output_file, 'w') as f:
        f.write(str(num_taxa) + '\n')
        for i in range(num_taxa):
            taxon_distances = []
            for j in range(num_taxa):
                taxon_distances.append(str(distance_matrix[taxa_names[i]][taxa_names[j]]))
            f.write(taxa_names[i] + '\t' + '\t'.join(taxon_distances) + '\n')

csv_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_nononsyn.csv'
output_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_nononsyn.dist'

distance_matrix = read_distance_matrix_from_csv(csv_file)
convert_to_phylip(distance_matrix, output_file)

csv_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_syn.csv'
output_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_syn.dist'

distance_matrix = read_distance_matrix_from_csv(csv_file)
convert_to_phylip(distance_matrix, output_file)

exit()
```
#### Splitstree <a name="20"></a>

Construct a phylogenetic network for the samples, based upon the nsynonymous SNP data.
```bash
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree
```
## SNP rarefaction curve <a name="21"></a>
### Corehunter (redundant) <a name="22"></a>
Run corehunter with 2:209 samples as a proxy for a saturation curve:
```bash
wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_mperc.csv

for csv in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_mperc.csv); do
    echo $csv >> logs/corehunterlog.txt
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/corehunter
    mkdir $OutDir
    for i in {0..209}; do
        if [[ "$i" -ne "1" ]]; then
            echo $i
            echo $i >> logs/corehunterlog.txt
            Coreno=$(echo $i | grep -vw '1')
            echo $Coreno
            Outfile=$(echo $Coreno)_core_sel.csv
            sbatch $ProgDir/run_corehunter.sh $csv $Coreno $OutDir $Outfile 2>&1 >> logs/corehunterlog.txt
        fi
    done
done

sacct -j 54241155 --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode

#There are 210 samples in the dataset, corehunter cannot run with only 1 or all of them, therefore with 0 for default settings there should be 209 output files
ls snp_calling/Myzus/persicae/biello/gatk/corehunter/*_core_sel.csv | wc -l
#209

#########################################################################################

wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_mperc.csv

for csv in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_mperc.csv); do
    echo $csv >> logs/corehunterlog_193.txt
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/corehunter_193
    mkdir $OutDir
    for i in {0..192}; do
        if [[ "$i" -ne "1" ]]; then
            echo $i
            echo $i >> logs/corehunterlog.txt
            Coreno=$(echo $i | grep -vw '1')
            echo $Coreno
            Outfile=$(echo $Coreno)_core_sel.csv
            sbatch $ProgDir/run_corehunter.sh $csv $Coreno $OutDir $Outfile 2>&1 >> logs/corehunterlog.txt
        fi
    done
done

#There are 193 samples in the dataset, corehunter cannot run with only 1 or all of them, therefore with 0 for default settings there should be 192 output files
ls snp_calling/Myzus/persicae/biello/gatk/corehunter/*_core_sel.csv | wc -l
#192
```
Combine corehunter outputs into one file:
```bash
interactive
source package /nbi/software/production/bin/python-2.7.11
python
```
```python
import sys
import os

samples = ['0'] + [ str(t) for t in range(1,212, 1)]

output_file = 'snp_calling/Myzus/persicae/biello/gatk/corehunter/core_sel_stepsize1.csv'
input_folder = 'snp_calling/Myzus/persicae/biello/gatk/corehunter/'

with open(output_file, 'w') as outf:
    for base in samples:
        inp_file = os.path.join(input_folder, base + "_core_sel.csv")
        if os.path.isfile(inp_file):
            with open(inp_file) as inp:
                A = inp.read().strip().split(',')
                iso = [a.replace('"', '') for a in A]
                outf.write(','.join([base] + iso) + '\n')

print('Done')
```
```python
#As Excel:
import sys
import pandas as pd
import os

samples = ['0'] + [ str(t) for t in range(1,212, 1)]

output_file = 'snp_calling/Myzus/persicae/biello/gatk/corehunter/core_sel_stepsize1.xlsx'
input_folder = 'snp_calling/Myzus/persicae/biello/gatk/corehunter/'

data = []
for base in samples:
    inp_file = os.path.join(input_folder, base + "_core_sel.csv")
    if os.path.isfile(inp_file):
        with open(inp_file) as inp:
            A = inp.read().strip().split(',')
            iso = [a.replace('"', '') for a in A]
            data.append([base] + iso)

df = pd.DataFrame(data, columns=['Base'] + ['Iso' + str(i+1) for i in range(len(iso))])

df.to_excel(output_file, index=False)

print('Done')

#################################################################
import sys
import pandas as pd
import os

samples = ['0'] + [ str(t) for t in range(1,212, 1)]

output_file = 'snp_calling/Myzus/persicae/biello/gatk/corehunter_193/core_sel_stepsize1_193.xlsx'
input_folder = 'snp_calling/Myzus/persicae/biello/gatk/corehunter_193/'

data = []
for base in samples:
    inp_file = os.path.join(input_folder, base + "_core_sel.csv")
    if os.path.isfile(inp_file):
        with open(inp_file) as inp:
            A = inp.read().strip().split(',')
            iso = [a.replace('"', '') for a in A]
            data.append([base] + iso)

df = pd.DataFrame(data, columns=['Base'] + ['Iso' + str(i+1) for i in range(len(iso))])

df.to_excel(output_file, index=False)

print('Done')
```
Remove individual corehunter output files:
```bash
rm snp_calling/Myzus/persicae/biello/gatk/corehunter/*_core_sel.csv
rm snp_calling/Myzus/persicae/biello/gatk/corehunter_193/*_core_sel.csv
```
Collate data from corehunter results

```python
import pandas as pd
import openpyxl

# Open the input Excel file and read in the data from the first sheet
input_file = 'snp_calling/Myzus/persicae/biello/gatk/corehunter/core_sel_stepsize1.xlsx'
input_df = pd.read_excel(input_file, sheet_name='Sheet1')

# Create an empty DataFrame to store the output data
output_df = pd.DataFrame(columns=['No.', 'Expected', 'Actual'])

# Initialize an empty set to store the core accessions, and a counter for the total load of core accessions
core = set()
load = 0

# Iterate over each row in the input DataFrame
for i, row in input_df.iterrows():
    # Extract the base value (an integer) and the accession IDs (strings)
    base = row['Base']
    accs = row.iloc[2:]
    
    # If the base value is greater than 0, update the set of core accessions, increment the load counter, and add a new row to the output DataFrame
    if base > 0:
        core.update(accs)
        load += base
        output_df.loc[i] = [base, load, len(core)]

# Open the output Excel file and write the output DataFrame to a new sheet
output_file = 'snp_calling/Myzus/persicae/biello/gatk/corehunter/core_sel_stepsize1.xlsx'
with pd.ExcelWriter(output_file, engine='openpyxl', mode='a') as writer:
    writer.book = openpyxl.load_workbook(output_file)
    output_df.to_excel(writer, sheet_name='Sheet3', index=False)

print('Done')
#FutureWarning: Setting the `book` attribute is not part of the public API, usage can give unexpected or corrupted results and will be removed in a future version

########################################################################################################
import pandas as pd
import openpyxl

# Open the input Excel file and read in the data from the first sheet
input_file = 'snp_calling/Myzus/persicae/biello/gatk/corehunter_193/core_sel_stepsize1_193.xlsx'
input_df = pd.read_excel(input_file, sheet_name='Sheet1')

# Create an empty DataFrame to store the output data
output_df = pd.DataFrame(columns=['No.', 'Expected', 'Actual'])

# Initialize an empty set to store the core accessions, and a counter for the total load of core accessions
core = set()
load = 0

# Iterate over each row in the input DataFrame
for i, row in input_df.iterrows():
    # Extract the base value (an integer) and the accession IDs (strings)
    base = row['Base']
    accs = row.iloc[2:]
    
    # If the base value is greater than 0, update the set of core accessions, increment the load counter, and add a new row to the output DataFrame
    if base > 0:
        core.update(accs)
        load += base
        output_df.loc[i] = [base, load, len(core)]

# Open the output Excel file and write the output DataFrame to a new sheet
output_file = 'snp_calling/Myzus/persicae/biello/gatk/corehunter_193/core_sel_stepsize1_193.xlsx'
with pd.ExcelWriter(output_file, engine='openpyxl', mode='a') as writer:
    writer.book = openpyxl.load_workbook(output_file)
    output_df.to_excel(writer, sheet_name='Sheet2', index=False)

print('Done')
#FutureWarning: Setting the `book` attribute is not part of the public API, usage can give unexpected or corrupted results and will be removed in a future version
```
```python
import random 
random.seed(1234)
import pandas as pd

from random import sample 

universal = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_mperc.csv') as inp:
    for line in inp:
        A = line.strip().split(',')
        if 'ID' in A:
            universal.extend(A[1:])
 
accessions = list(set(universal))
total_accessions = len(accessions)
print ('accessions = ', accessions)
print ('total_accessions = ', total_accessions)

sampled = []
df = pd.read_excel('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/corehunter/core_sel_stepsize1.xlsx')
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/corehunter/core_sel_stepsize1_curvey.xlsx', 'w') as outf:
    core = set()
    r_core = set()
    load = 0  # total sampled thus far
    pipe = map(str, ['sample_size', 'sampled', 'unique_accession', 'unique_accession_randomized'])
    outf.write(','.join(pipe) + '\n')
    for index, row in df.iterrows():
        base = int(row[0])
        accs = row[1:]
        pick_size = len(accs)
        if base > 1:
            assert pick_size == base, base 
            # random sample 
            r_accs =  sample(accessions, pick_size)
            r_core.update(r_accs)
            core.update(accs)
            load += base 
            pipe = map(str, [pick_size, load, len(core), len(r_core)])         
            outf.write(','.join(pipe) + '\n') 
            sampled.append(base)
            print (sum(sampled), len(core), base, len(r_core)) 

```
### True Rarefaction curves <a name="23"></a>
#### Whole genome SNPs <a name="24"></a>
Perform rarefaction analysis:
```bash
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf);do
    Reference=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
    Replicates=100
    Steps=193
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction
    OutFile=rarefied-HomozygousALT-SNPS-genomic-mac1-193-100-193
    mkdir $OutDir
    sbatch $ProgDir/run_snprarefaction.sh $InFile $OutDir $OutFile $Reference $Replicates $Steps
done #57151199

#########################################################################################################################

for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/reference/Myzus_persicae_O_v2.0.scaffolds.fa.gz
    Replicates=100
    Steps=193
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction
    OutFile=rarefied-HomozygousALT-SNPS-genomic-193-100-193
    mkdir $OutDir
    sbatch $ProgDir/run_snprarefaction.sh $InFile $OutDir $OutFile $Reference $Replicates $Steps
done
#54359840, 54371114, 54373352, 54385547

#Column titles from csv =
head -n 1 $OutDir/$OutFile

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
```
Plot the rarefaction curve:
```python
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc

rc('mathtext', default='regular')

df = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genomic-193-100-193-curve.csv")

#Column titles from csv = sample_size,rep_polymorphisms_median,rep_polymorphisms_median_fraction

X = df['sample_size']
Y1 = df['rep_polymorphisms_median_fraction']*100.
Y2 = df['rep_polymorphisms_median']

fig, ax1 = plt.subplots(figsize=(10, 5)) 

color = 'tab:orange'
ax1.set_xlabel('Sampled Size') 
ax1.set_ylabel('Percentage', color = color) 
ax1.plot(X, Y1, color = color) 
ax1.tick_params(axis ='y', labelcolor = color) 

ax2 = ax1.twinx() 

color = 'tab:cyan'
ax2.set_ylabel('Unique SNPs', color = color) 
ax2.plot(X, Y2, color = 'tab:orange', linestyle='--', marker='o',) 
ax2.tick_params(axis ='y', labelcolor = color) 

# Add a global title to the plot
suptitle = plt.suptitle('Nuclear genomic SNP rarefaction curve with 193 samples using 100 replicates, step size 1', y=1.02)

plt.tight_layout()

# When you save the fig, add the suptitle text object as an extra artist
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genome-193-100-193-curve.pdf", bbox_extra_artists=(suptitle,), bbox_inches="tight") 
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genome-193-100-193-curve.png", bbox_extra_artists=(suptitle,), bbox_inches="tight") 
```
#### Genic SNPs <a name="25"></a>
Perform rarefaction analysis:
```bash
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz);do
    Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa.gz
    Replicates=100
    Steps=193
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction
    OutFile=rarefied-HomozygousALT-SNPS-genic-193-100-193
    mkdir $OutDir
    sbatch $ProgDir/run_snprarefaction.sh $InFile $OutDir $OutFile $Reference $Replicates $Steps
done
#54362194, 54371112, 54373403, 54379262

#Column titles from csv =
head -n 1 $OutDir/${OutFile}-curve.csv 

########################################################################################################################

bcftools view -v snps /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.vcf.gz | grep -v "^#" | wc -l #126,162
bcftools view -v snps /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.vcf.gz | grep -v "^#" | grep '1/1:'| wc -l  #28,091
#True CDS SNPS:
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.vcf.gz);do
    Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa.gz
    Replicates=100
    Steps=193
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction
    OutFile=rarefied-HomozygousALT-SNPS-CDS-193-100-193
    mkdir $OutDir
    sbatch $ProgDir/run_snprarefaction.sh $InFile $OutDir $OutFile $Reference $Replicates $Steps
done #56698263
```
Plot the rarefaction curve:
```python
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc

rc('mathtext', default='regular')

df = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-CDS-193-100-193-curve.csv")

#Column titles from csv = sample_size,rep_polymorphisms_median,rep_polymorphisms_median_fraction

X = df['sample_size']
Y1 = df['rep_polymorphisms_median_fraction']*100.
Y2 = df['rep_polymorphisms_median']

fig, ax1 = plt.subplots(figsize=(10, 5)) 

color = 'tab:orange'
ax1.set_xlabel('Sampled Size') 
ax1.set_ylabel('Percentage', color = color) 
ax1.plot(X, Y1, color = color) 
ax1.tick_params(axis ='y', labelcolor = color) 

ax2 = ax1.twinx() 

color = 'tab:cyan'
ax2.set_ylabel('Unique SNPs', color = color) 
ax2.plot(X, Y2, color = 'tab:orange', linestyle='--', marker='o',) 
ax2.tick_params(axis ='y', labelcolor = color) 

# Add a global title to the plot
suptitle = plt.suptitle('Nuclear CDS SNP rarefaction curve with 193 samples using 100 replicates, step size 1', y=1.02)

plt.tight_layout()

# When you save the fig, add the suptitle text object as an extra artist
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/CDS-193-100-193-curve2.pdf", bbox_extra_artists=(suptitle,), bbox_inches="tight") 
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/CDS-193-100-193-curve2.png", bbox_extra_artists=(suptitle,), bbox_inches="tight") 

df2 = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genomic-193-100-193-curve.csv")
X2 = df2['sample_size']
Y1_2 = df2['rep_polymorphisms_median_fraction']*100.
Y2_2 = df2['rep_polymorphisms_median']
ax1.plot(X2, Y1_2, color='tab:blue')  # Plot Y1_2 data on the same ax1
ax2.plot(X2, Y2_2, color='tab:blue', linestyle='--', marker='o')  # Plot Y2_2 data on the same ax2
suptitle = plt.suptitle('Comparison of SNP rarefaction curves with 100 reps, step size 1', y=1.02)
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genic+genomic-193-100-193-curve.pdf", bbox_extra_artists=(suptitle,), bbox_inches="tight")
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genic+genomic-193-100-193-curve.png", bbox_extra_artists=(suptitle,), bbox_inches="tight")
```
#### For Bass and JIC samples seperately - Whole genome and genic <a name="26"></a>

```bash
#High missingness samples:
awk -F'\t' '$NF > 0.1 {print $1}' snp_calling/Myzus/persicae/biello/gatk/filtered/209s.M_persicae.out.imiss | grep -v 'INDV' > temp.txt
#Ligustri
echo ligustri >> temp.txt
#JIC samples:
cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'JIC' | awk -F',' '{print $2}' >> temp.txt
wc -l temp.txt #96
for vcf in $(ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs.vcf.gz);do
    Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa.gz
    Replicates=100
    Steps=193
    Exclusion_list=temp.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction
    OutFile=rarefied-HomozygousALT-SNPS-genomic-Bass
    mkdir $OutDir
    sbatch $ProgDir/run_snprarefaction.sh $InFile $OutDir $OutFile $Reference $Replicates $Steps $Exclusion_list
done #54395430
for vcf in $(ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf.gz);do
    Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa.gz
    Replicates=100
    Steps=193
    Exclusion_list=temp.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction
    OutFile=rarefied-HomozygousALT-SNPS-genic-Bass
    mkdir $OutDir
    sbatch $ProgDir/run_snprarefaction.sh $InFile $OutDir $OutFile $Reference $Replicates $Steps $Exclusion_list
done #54395431

#High missingness samples:
awk -F'\t' '$NF > 0.1 {print $1}' snp_calling/Myzus/persicae/biello/gatk/filtered/209s.M_persicae.out.imiss | grep -v 'INDV' > temp.txt
#Ligustri
echo ligustri >> temp2.txt
#Bass samples:
cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Bass' | awk -F',' '{print $2}' >> temp.txt
wc -l temp.txt #147
for vcf in $(ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs.vcf.gz);do
    Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa.gz
    Replicates=100
    Steps=193
    Exclusion_list=temp.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction
    OutFile=rarefied-HomozygousALT-SNPS-genomic-JIC
    mkdir $OutDir
    sbatch $ProgDir/run_snprarefaction.sh $InFile $OutDir $OutFile $Reference $Replicates $Steps $Exclusion_list
done #54394326
for vcf in $(ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf.gz);do
    Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa.gz
    Replicates=100
    Steps=193
    Exclusion_list=temp2.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction
    OutFile=rarefied-HomozygousALT-SNPS-genic-JIC
    mkdir $OutDir
    sbatch $ProgDir/run_snprarefaction.sh $InFile $OutDir $OutFile $Reference $Replicates $Steps $Exclusion_list
done #54394327

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
```
Plot the rarefaction curves for Bass and JIC samples seperately - Whole genome and genic:
```python
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc

rc('mathtext', default='regular')

df = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-JIC-curve.csv")

#Column titles from csv = sample_size,rep_polymorphisms_median,rep_polymorphisms_median_fraction

X = df['sample_size']
Y1 = df['rep_polymorphisms_median_fraction']*100.
Y2 = df['rep_polymorphisms_median']

fig, ax1 = plt.subplots(figsize=(10, 5)) 

color = 'tab:orange'
ax1.set_xlabel('Sampled Size') 
ax1.set_ylabel('Percentage', color = color) 
ax1.plot(X, Y1, color = color) 
ax1.tick_params(axis ='y', labelcolor = color) 

ax2 = ax1.twinx() 

color = 'tab:cyan'
ax2.set_ylabel('Unique SNPs', color = color) 
ax2.plot(X, Y2, color = 'tab:orange', linestyle='--', marker='o',) 
ax2.tick_params(axis ='y', labelcolor = color) 

# Add a global title to the plot
suptitle = plt.suptitle('Nuclear genic SNP rarefaction curve with JIC samples using 100 replicates, step size 1', y=1.02)

plt.tight_layout()

# When you save the fig, add the suptitle text object as an extra artist
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-JIC-curve.pdf", bbox_extra_artists=(suptitle,), bbox_inches="tight") 
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-JIC-curve.png", bbox_extra_artists=(suptitle,), bbox_inches="tight") 

df2 = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genomic-JIC-curve.csv")
X2 = df2['sample_size']
Y1_2 = df2['rep_polymorphisms_median_fraction']*100.
Y2_2 = df2['rep_polymorphisms_median']
ax1.plot(X2, Y1_2, color='tab:blue')  # Plot Y1_2 data on the same ax1
ax2.plot(X2, Y2_2, color='tab:blue', linestyle='--', marker='o')  # Plot Y2_2 data on the same ax2
suptitle = plt.suptitle('Comparison of JIC SNP rarefaction curves with 100 reps, step size 1', y=1.02)
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genic+genomic-JIC-curve.pdf", bbox_extra_artists=(suptitle,), bbox_inches="tight")
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genic+genomic-JIC-curve.png", bbox_extra_artists=(suptitle,), bbox_inches="tight")

###################################################################################################

df = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-Bass-curve.csv")

#Column titles from csv = sample_size,rep_polymorphisms_median,rep_polymorphisms_median_fraction

X = df['sample_size']
Y1 = df['rep_polymorphisms_median_fraction']*100.
Y2 = df['rep_polymorphisms_median']

fig, ax1 = plt.subplots(figsize=(10, 5)) 

color = 'tab:orange'
ax1.set_xlabel('Sampled Size') 
ax1.set_ylabel('Percentage', color = color) 
ax1.plot(X, Y1, color = color) 
ax1.tick_params(axis ='y', labelcolor = color) 

ax2 = ax1.twinx() 

color = 'tab:cyan'
ax2.set_ylabel('Unique SNPs', color = color) 
ax2.plot(X, Y2, color = 'tab:orange', linestyle='--', marker='o',) 
ax2.tick_params(axis ='y', labelcolor = color) 

# Add a global title to the plot
suptitle = plt.suptitle('Nuclear genic SNP rarefaction curve with Bass samples using 100 replicates, step size 1', y=1.02)

plt.tight_layout()

# When you save the fig, add the suptitle text object as an extra artist
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-Bass-curve.pdf", bbox_extra_artists=(suptitle,), bbox_inches="tight") 
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-Bass-curve.png", bbox_extra_artists=(suptitle,), bbox_inches="tight") 

df2 = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genomic-Bass-curve.csv")
X2 = df2['sample_size']
Y1_2 = df2['rep_polymorphisms_median_fraction']*100.
Y2_2 = df2['rep_polymorphisms_median']
ax1.plot(X2, Y1_2, color='tab:blue')  # Plot Y1_2 data on the same ax1
ax2.plot(X2, Y2_2, color='tab:blue', linestyle='--', marker='o')  # Plot Y2_2 data on the same ax2
suptitle = plt.suptitle('Comparison of Bass SNP rarefaction curves with 100 reps, step size 1', y=1.02)
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genic+genomic-Bass-curve.pdf", bbox_extra_artists=(suptitle,), bbox_inches="tight")
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genic+genomic-Bass-curve.png", bbox_extra_artists=(suptitle,), bbox_inches="tight")
```
Plot the rarefaction curves for Bass and JIC samples whole genome and genic SNPs together:
With % and absolute values plotted:

```python
df = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-Bass-curve.csv")

#Column titles from csv = sample_size,rep_polymorphisms_median,rep_polymorphisms_median_fraction

X = df['sample_size']
Y1 = df['rep_polymorphisms_median_fraction']*100.
Y2 = df['rep_polymorphisms_median']

fig, ax1 = plt.subplots(figsize=(10, 5)) 

color = 'tab:orange'
ax1.set_xlabel('Sampled Size') 
ax1.set_ylabel('Percentage', color = color) 
ax1.plot(X, Y1, color = color) 
ax1.tick_params(axis ='y', labelcolor = color) 

ax2 = ax1.twinx() 

color = 'tab:cyan'
ax2.set_ylabel('Unique SNPs', color = color) 
ax2.plot(X, Y2, color = 'tab:orange', linestyle='--', marker='o',) 
ax2.tick_params(axis ='y', labelcolor = color) 

plt.tight_layout()

df2 = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genomic-Bass-curve.csv")
X2 = df2['sample_size']
Y1_2 = df2['rep_polymorphisms_median_fraction']*100.
Y2_2 = df2['rep_polymorphisms_median']
ax1.plot(X2, Y1_2, color='tab:blue')  # Plot Y1_2 data on the same ax1
ax2.plot(X2, Y2_2, color='tab:blue', linestyle='--', marker='o')  # Plot Y2_2 data on the same ax2

df3 = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-JIC-curve.csv")
X3 = df3['sample_size']
Y1_3 = df3['rep_polymorphisms_median_fraction']*100.
Y2_3 = df3['rep_polymorphisms_median']
ax1.plot(X3, Y1_3, color='tab:green')  
ax2.plot(X3, Y2_3, color='tab:green', linestyle='--', marker='o')  

df4 = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genomic-JIC-curve.csv")
X4 = df4['sample_size']
Y1_4 = df4['rep_polymorphisms_median_fraction']*100.
Y2_4 = df4['rep_polymorphisms_median']
ax1.plot(X4, Y1_4, color='tab:red')  # Plot Y1_2 data on the same ax1
ax2.plot(X4, Y2_4, color='tab:red', linestyle='--', marker='o')  # Plot Y2_2 data on the same ax2

suptitle = plt.suptitle('Comparison of Bass and JIC percentage SNP rarefaction curves with 100 reps, step size 1', y=1.02)
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genic+genomic-Bass+JIC-percent-curve.pdf", bbox_extra_artists=(suptitle,), bbox_inches="tight")
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genic+genomic-Bass+JIC-percent-curve.png", bbox_extra_artists=(suptitle,), bbox_inches="tight")

```
With only absolute points plotted:
```python

df = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-Bass-curve.csv")

#Column titles from csv = sample_size,rep_polymorphisms_median,rep_polymorphisms_median_fraction

X = df['sample_size']
Y1 = df['rep_polymorphisms_median_fraction']*100.
Y2 = df['rep_polymorphisms_median']

fig, ax1 = plt.subplots(figsize=(10, 5)) 

color = 'tab:orange'
ax1.set_xlabel('Sampled Size') 
ax1.set_ylabel('Percentage', color = color) 
#ax1.plot(X, Y1, color = color) 
ax1.tick_params(axis ='y', labelcolor = color) 

ax2 = ax1.twinx() 

color = 'tab:cyan'
ax2.set_ylabel('Unique SNPs', color = color) 
ax2.plot(X, Y2, color = 'tab:orange', linestyle='--', marker='o',) 
ax2.tick_params(axis ='y', labelcolor = color) 

plt.tight_layout()

df2 = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genomic-Bass-curve.csv")
X2 = df2['sample_size']
Y1_2 = df2['rep_polymorphisms_median_fraction']*100.
Y2_2 = df2['rep_polymorphisms_median']
#ax1.plot(X2, Y1_2, color='tab:blue')  # Plot Y1_2 data on the same ax1
ax2.plot(X2, Y2_2, color='tab:blue', linestyle='--', marker='o')  # Plot Y2_2 data on the same ax2

df3 = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-JIC-curve.csv")
X3 = df3['sample_size']
Y1_3 = df3['rep_polymorphisms_median_fraction']*100.
Y2_3 = df3['rep_polymorphisms_median']
#ax1.plot(X3, Y1_3, color='tab:green')  
ax2.plot(X3, Y2_3, color='tab:green', linestyle='--', marker='o')  

df4 = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genomic-JIC-curve.csv")
X4 = df4['sample_size']
Y1_4 = df4['rep_polymorphisms_median_fraction']*100.
Y2_4 = df4['rep_polymorphisms_median']
#ax1.plot(X4, Y1_4, color='tab:red')  # Plot Y1_2 data on the same ax1
ax2.plot(X4, Y2_4, color='tab:red', linestyle='--', marker='o')  # Plot Y2_2 data on the same ax2

suptitle = plt.suptitle('Comparison of Bass and JIC SNP rarefaction curves with 100 reps, step size 1', y=1.02)
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genic+genomic-Bass+JIC-curve.pdf", bbox_extra_artists=(suptitle,), bbox_inches="tight")
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genic+genomic-Bass+JIC-curve.png", bbox_extra_artists=(suptitle,), bbox_inches="tight")
```
Plot again with % points only and scaled x axis to allow direct comparison of JIC and Bass datasets:
```python
df = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-Bass-curve.csv")
X = df['sample_size']
Y1 = df['rep_polymorphisms_median_fraction'] * 100.

df3 = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-JIC-curve.csv")
X3 = df3['sample_size']
Y1_3 = df3['rep_polymorphisms_median_fraction'] * 100.

fig, ax1 = plt.subplots(figsize=(10, 5))

color = 'tab:orange'
ax1.set_xlabel('Sampled Size')
ax1.set_ylabel('Percentage', color=color)
ax1.plot(X, Y1, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twiny()  # Create a twin x-axis
ax2.plot(X3, Y1_3, color='tab:green')

# Set the x-axis limits based on the maximum and minimum values of both datasets
min_x = min(min(X), min(X3))
max_x = max(max(X), max(X3))
ax1.set_xlim(0, max_x)
ax2.set_xlim(0, max(X3))

# Set the y-axis limits to start at zero
ax1.set_ylim(0, max(Y1))
ax2.set_ylim(0, max(Y1_3))

ax1.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

ax1.tick_params(axis='x', direction='out')
ax1.tick_params(axis='y', colors=color)
ax2.tick_params(axis='x', colors='tab:green')

ax1.set_xlabel('Sampled Size (Bass)', color=color)
ax2.set_xlabel('Sampled Size (JIC)', color='tab:green')

suptitle = plt.suptitle('Comparison of Bass and JIC percentage SNP rarefaction curves with 100 reps, step size 1', y=1.02)
plt.tight_layout()
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genic+genomic-Bass+JIC-percent-curve-scaled.png", bbox_extra_artists=(suptitle,), bbox_inches="tight")
```
### Fitting a curve <a name="27"></a>
```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def logarithmic_func(x, a, b):
    return a * np.log(x) + b

df = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-Bass-curve.csv")
X = df['sample_size']
Y1 = df['rep_polymorphisms_median_fraction'] * 100.

# Curve fitting
try:
    popt, pcov = curve_fit(logarithmic_func, X, Y1)
except RuntimeError:
    popt = [1, 1]
    pcov = None

# Find the value of the horizontal asymptote
if popt[0] == 0:
    horizontal_asymptote = popt[1]
else:
    horizontal_asymptote = None

# Print the value of the horizontal asymptote
print("Horizontal Asymptote:", horizontal_asymptote)

# Plot the original data and fitted curve
plt.scatter(X, Y1, label='Data')
plt.plot(X, logarithmic_func(X, *popt), color='red', label='Fitted Curve')
plt.xlabel('Sample Size')
plt.ylabel('Polymorphisms Median Fraction (%)')
plt.legend()
plt.show()

# Save the figure as an image file
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/curve_with_asymptote_bass.png", bbox_inches="tight")

##################################################################################################################
def logarithmic_func(x, a, b):
    return a * np.log(x) + b

df = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-JIC-curve.csv")
X = df['sample_size']
Y1 = df['rep_polymorphisms_median_fraction'] * 100.

# Curve fitting
try:
    popt, pcov = curve_fit(logarithmic_func, X, Y1)
except RuntimeError:
    popt = [1, 1]
    pcov = None

# Find the value of the horizontal asymptote
if popt[0] == 0:
    horizontal_asymptote = popt[1]
else:
    horizontal_asymptote = None

# Print the value of the horizontal asymptote
print("Horizontal Asymptote:", horizontal_asymptote)

# Plot the original data and fitted curve
plt.scatter(X, Y1, label='Data')
plt.plot(X, logarithmic_func(X, *popt), color='red', label='Fitted Curve')
plt.xlabel('Sample Size')
plt.ylabel('Polymorphisms Median Fraction (%)')
plt.legend()
plt.show()

# Save the figure as an image file
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/curve_with_asymptote_jic.png", bbox_inches="tight")
##################################################################################################################
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def logarithmic_func(x, a, b):
    return a * np.log(x) + b

df = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-193-100-193-curve.csv")
X = df['sample_size']
Y1 = df['rep_polymorphisms_median']

# Curve fitting
try:
    popt, pcov = curve_fit(logarithmic_func, X, Y1)
except RuntimeError:
    popt = [1, 1]
    pcov = None

# Find the value of the horizontal asymptote
x_value = 1000000
estimated_y_value = logarithmic_func(x_value, *popt)
print("Estimated y value at x = 1000000:", estimated_y_value)
#243.8990609631632

# Plot the original data and fitted curve
plt.scatter(X, Y1, label='Data')
plt.plot(X, logarithmic_func(X, *popt), color='red', label='Fitted Curve')
plt.xlabel('Sample Size')
plt.ylabel('Genic Polymorphisms Median')
plt.legend()

# Save the figure as an image file
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/curve_with_asymptote-logarithmic.png", bbox_inches="tight")
```
```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def exponential_func(x, a, b):
    return a * np.exp(-b * x)

df = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-193-100-193-curve.csv")
X = df['sample_size']
Y1 = df['rep_polymorphisms_median']

# Curve fitting
X_scaled = X / 1000

try:
    popt, pcov = curve_fit(exponential_func, X_scaled, Y1)
except RuntimeError:
    popt = [1, 1, 1]
    pcov = None


# Find the value of the curve at a specific x value
x_value = 1000000
estimated_y_value = exponential_func(x_value, *popt)
print("Estimated y value at x = 1000000:", estimated_y_value)

# Plot the original data and fitted curve
plt.scatter(X, Y1, label='Data')
plt.plot(X, exponential_func(X, *popt), color='red', label='Fitted Curve')
plt.xlabel('Sample Size')
plt.ylabel('Genic Polymorphisms Median')
plt.legend()

# Save the figure as an image file
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/curve_with_asymptote-exponential.png", bbox_inches="tight")

```
## vcftools plot of saturation curve <a name="28"></a>

```bash
gzip -cd snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.vcf.gz | head -n 2000 > temp.vcf
rm temp.vcf

#Remove SNPs without variant:
source package /nbi/software/testing/bin/vcftools-0.1.15
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.vcf.gz --mac 1 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.vcf
#After filtering, kept 7,873,274 out of a possible 11,870,914 Sites = 1/3 removed - this file contains both ligustri and those high missingness samples that we remove later, I dont understand how there are ~4 million snps with no minor alleles present in this file

vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.vcf --mac 1 --remove-filtered-all --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf

#Using zlib version: 1.2.11
#After filtering, kept 210 out of 210 Individuals
#Outputting VCF file...
#After filtering, kept 7873274 out of a possible 7873274 Sites
#Run Time = 3819.00 seconds

#Using zlib version: 1.2.11
#After filtering, kept 210 out of 210 Individuals
#Outputting VCF file...
#After filtering, kept 7873274 out of a possible 7873274 Sites
#Run Time = 3910.00 seconds

vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --max-missing 1 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.nomissing.vcf
#After filtering, kept 96 out of a possible 7873274 Sites
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --max-missing 0.99 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.nomissing99.vcf
#After filtering, kept 206,656 out of a possible 7873274 Sites
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --max-missing 0.98 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.nomissing98.vcf
#After filtering, kept 1,304,612 out of a possible 7873274 Sites
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --max-missing 0.97 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.nomissing97.vcf
#After filtering, kept 2,943,804 out of a possible 7873274 Sites
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --max-missing 0.96 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.nomissing96.vcf
#After filtering, kept 4,452,503 out of a possible 7873274 Sites
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --max-missing 0.95 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.nomissing95.vcf
#After filtering, kept 5674625 out of a possible 7873274 Sites
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --max-missing 0.94 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.nomissing94.vcf
#After filtering, kept 6558942 out of a possible 7873274 Sites
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --max-missing 0.93 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.nomissing93.vcf
#After filtering, kept 7120848 out of a possible 7873274 Sites
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --max-missing 0.92 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.nomissing92.vcf
#After filtering, kept 7458687 out of a possible 7873274 Sites
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --max-missing 0.91 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.nomissing91.vcf
#After filtering, kept 7669588 out of a possible 7873274 Sites
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --max-missing 0.9 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.nomissing90.vcf
#After filtering, kept 7873274 out of a possible 7873274 Sites
echo finished


for vcf in $(ls snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.nomissing99.vcf.recode.vcf); do
    InFile=$vcf
    OutDir=snp_calling/Myzus/persicae/biello/gatk/filtered
    OutFile=snpsaturation_99.tsv
    SampleNo=209
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_SNPsaturation.sh $InFile $OutDir $OutFile $SampleNo
done
#54315804, 28mins elapsed

for vcf in $(ls snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.nomissing97.vcf.recode.vcf); do
    InFile=$vcf
    OutDir=snp_calling/Myzus/persicae/biello/gatk/filtered
    OutFile=snpsaturation_97.tsv
    SampleNo=209
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_SNPsaturation.sh $InFile $OutDir $OutFile $SampleNo
done
#54316707,  mins elapsed

vcftools --vcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --thin 1000 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs_thinned1000.vcf
#After filtering, kept 275,413 out of a possible 7873274 Sites
vcftools --vcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --thin 500 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs_thinned500.vcf
#After filtering, kept 493302 out of a possible 7873274 Sites
vcftools --vcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --thin 250 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs_thinned250.vcf
#After filtering, kept 863346 out of a possible 7873274 Sites
vcftools --vcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --thin 150 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs_thinned150.vcf
#After filtering, kept 1282244 out of a possible 7873274 Sites
vcftools --vcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --thin 100 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs_thinned100.vcf
#After filtering, kept 1726133 out of a possible 7873274 Sites
vcftools --vcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.removefiltered.vcf.recode.vcf --thin 50 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs_thinned50.vcf
#After filtering, kept 2730798 out of a possible 7873274 Sites

for vcf in $(ls snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs_thinned50.vcf.recode.vcf); do
    InFile=$vcf
    OutDir=snp_calling/Myzus/persicae/biello/gatk/filtered
    OutFile=snpsaturation_thinned50.tsv
    SampleNo=209
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_SNPsaturation.sh $InFile $OutDir $OutFile $SampleNo
done
#Submitted batch job 54316712, hours

for vcf in $(ls snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.vcf); do
    InFile=$vcf
    OutDir=snp_calling/Myzus/persicae/biello/gatk/filtered
    OutFile=snpsaturation.tsv
    SampleNo=209
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_SNPsaturation.sh $InFile $OutDir $OutFile $SampleNo
done
#54309843, 18.5 hours
```

## Nei gene diversity <a name="29"></a>

Calculate Nei's gene diversity
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/nei
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz); do
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tomF_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/nei
    OutFile=193s.M_persicae.onlySNPs.tab
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch  $ProgDir/run_nei-from-vcf.sh $InFile $OutDir $OutFile
done
#54436237
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz); do
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/nei
    OutFile=193s.M_persicae.onlySNPs-genic-regions.tab
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch  $ProgDir/run_nei-from-vcf.sh $InFile $OutDir $OutFile
done
#54436239
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions.recode.vcf.gz); do
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/nei
    OutFile=193s.M_persicae.onlySNPs-CDS-regions.tab
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch  $ProgDir/run_nei-from-vcf.sh $InFile $OutDir $OutFile
done
#56698277
```
## Check for SNPs of clone O versus itself to control for diploidy <a name="30"></a>

```bash
source package /nbi/software/testing/bin/vcftools-0.1.15
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz --indv O --counts
grep -v "^#" out.frq.count| cut -f3 | grep -cv "0"

source package /nbi/software/testing/bin/bcftools-1.8
bcftools view -c1 -s O -Ou /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz | bcftools query -e 'REF=ALT' -f '%CHROM\t%POS\n' | wc -l
#272,953
```
Remove SNPs that also occu between clone O and itself:
```bash
bcftools view -c1 -s O -Ou /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz | bcftools query -e 'REF=ALT' -f '%CHROM\t%POS\n' > temp.vcf
source package 01ef5a53-c149-4c9e-b07d-0b9a46176cc0
bgzip -cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz > 193s.M_persicae.onlySNPs.vcf
```
```python
# Read the first file and store the values in a set
first_file_lines = set()
with open('temp.vcf', 'r') as f1:
    for line in f1:
        columns = line.strip().split('\t')
        if len(columns) >= 2:
            first_file_lines.add((columns[0], columns[1]))

# Read the second file and write lines that don't match the criteria to a new file
with open('193s.M_persicae.onlySNPs.vcf', 'r') as f2, open('192s.M_persicae.onlySNPs.vcf', 'w') as output:
    for line in f2:
        columns = line.strip().split('\t')
        if len(columns) >= 2 and (columns[0], columns[1]) not in first_file_lines:
            output.write(line)

exit()
```
```bash
grep '#' 193s.M_persicae.onlySNPs.vcf >> 193s-no0.M_persicae.onlySNPs.vcf
grep -v '#' 192s.M_persicae.onlySNPs.vcf >> 193s-no0.M_persicae.onlySNPs.vcf
rm 192s.M_persicae.onlySNPs.vcf 193s.M_persicae.onlySNPs.vcf
bgzip 193s-no0.M_persicae.onlySNPs.vcf
mv 193s-no0.M_persicae.onlySNPs.vcf.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/.


bcftools view -s O -Ou /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz | bcftools query -e 'REF=ALT' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT\t%O[\t]\n' > temp.vcf
bcftools query -f '%CHROM\t%POS\n' -i 'ALT!="."' -s O /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz > snps_to_remove.txt
source package 01ef5a53-c149-4c9e-b07d-0b9a46176cc0
bgzip -cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz > 193s.M_persicae.onlySNPs.vcf
bcftools view -s O -v snps -f 'GT="0/0"' 193s.M_persicae.onlySNPs.vcf > filtered_file.vcf
rm 193s.M_persicae.onlySNPs.vcf
bgzip filtered_file.vcf
tabix -p vcf filtered_file.vcf.gz
```




