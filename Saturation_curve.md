# Exploring the diversity of the aphid samples
Unless stated otherwise work was performed from the directory: /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae

Collecting data:
```bash
mkdir -p snp_calling/Myzus/persicae/biello/gatk/filtered
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/vcffilesafterfiltering/210s.M_persicae.onlySNPs.vcf.gz snp_calling/Myzus/persicae/biello/gatk/filtered/.

ln -s /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3

ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa.gz
```
## Data quality control
The SNP data has already been filtered by R.Wouters and R.biello with vcftools etc.

### Coverage
Coverage:
```bash
mkdir snp_calling/Myzus/persicae/biello/gatk/reports
#JIC samples
for sample in $(ls /jic/scratch/groups/Saskia-Hogenhout/roland/Myzus_periscae_popgen/Output/realigned_BAM_files/*.sorted.mark_dups.realigned_stats/genome_results.txt); do
    SampleName=$(echo $sample | cut -d '/' -f10| cut -d '.' -f1)
    Coverage=$(grep 'mean coverageData' $sample)
    echo $SampleName $Coverage >> snp_calling/Myzus/persicae/biello/gatk/reports/JIC_coverage.txt
done

grep 'mean coverageData' /jic/scratch/groups/Saskia-Hogenhout/roland/Myzus_periscae_popgen/Output/realigned_BAM_files/UKW3.sorted.mark_dups.realigned_stats/genome_results.txt

\\jic-hpc-data\Group-Scratch\Saskia-Hogenhout\roland\Myzus_periscae_popgen\Output\realigned_BAM_files\A105.sorted.mark_dups.realigned_stats

switch-institute nbi; source jdk-7u45; /jic/research-groups/Saskia-Hogenhout/Roland/qualimap_v2.2.1/qualimap bamqc -nt 16 -bam ../A105.sorted.mark_dups.realigned.bam
```
### Missingness
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
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz --mac 1 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1 #After filtering, kept 1,851,865 out of a possible 11,870,914 Sites

nano scaff1-6.bed
##CHROM  START   END
#scaffold1   1   105178091
#scaffold2   1   86073209
#scaffold3   1   69480500
#scaffold4   1   62328371
#scaffold5   1   30612500
#scaffold6   1   29865500

vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf --positions scaff1-6.bed --window-pi 10000 --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.snp_counts

bcftools query -f '%CHROM\t%POS\n' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf | \
    awk -F '\t' '{if ($1 ~ /^scaffold_[1-6]$/) count[$1]++;} END {for (i=1; i<=6; i++) print "scaffold_"i ": " count["scaffold_"i] " SNPs";}'

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
```
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

## Extract genic SNPs
```bash
wc -l snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3
#902,603

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 \
-header > snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 \
-header > snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-genic-regions.vcf

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 \
-header > /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf

source package 01ef5a53-c149-4c9e-b07d-0b9a46176cc0
bgzip -c snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf > snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz
zcat snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz | wc -l #31,637,658 there are many duplicate SNPs
bgzip -c /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf > /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf.gz
zcat /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf.gz | wc -l #31,637,658

#deduplicated files with:
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


######################################################################################################################
# Mouse genome, version mm37.61
mm37.61.genome : Mouse

# Peach potato genome, version 2.1
M_persicae.genome : aphid



snpEff ann -i vcf -o gff -ud 0 myzus_persicae snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz > snp_effects.gff

```
Extract CDS/exon/UTR regions region SNPs:
```bash
grep '#\|CDS\|exon\|five_prime_UTR\|three_prime_UTR' snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 > snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_exon_UTR_annotation.gff3
wc -l snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_exon_UTR_annotation.gff3
#806,274

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_exon_UTR_annotation.gff3 \
-header > snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_exon_UTR_annotation.gff3 \
-header > snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_exon_UTR_genic-regions.vcf

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

#deduplicated files with:
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
Extract CDS region SNPs:
```bash
grep '#\|CDS' snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 > snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.gff3
wc -l snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.gff3
#336,709

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.gff3 \
-header > snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic-regions.vcf

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs.vcf.gz \
-b snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.CDS_annotation.gff3 \
-header > snp_calling/Myzus/persicae/biello/gatk/filtered/ligustri+193s.M_persicae.onlySNPs-CDS_genic-regions.vcf

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

#deduplicated files with:
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

#################################################################################################################################
#Convert .vcf to annovar format
convert2annovar.pl -format vcf4 snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz -outfile ex2.avinput -allsample -includeinfo #-withzyg -withfreq
annotate_variation.pl --regionanno --dbtype gff3 --gff3dbfile MYZPE13164_O_EIv2.1.annotation.gff3 ex2.avinput /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/ --out ex1


annotate_variation.pl --geneanno --dbtype gff3 --gff3dbfile MYZPE13164_O_EIv2.1.annotation.gff3 ex2.avinput /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/ --out ex1
```
Non-synonymous CDS SNPs
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

scaffold_1      MYZPE13164_O_EIv2.1     gene    8196    9175    545     -       .       ID=MYZPE13164_O_EIv2.1_0000010;Name=MYZPE13164_O_EIv2.1_0000010;biotype=transposable_element_gene;confidence=High
scaffold_1      MYZPE13164_O_EIv2.1     mRNA    8196    9175    545     -       .       ID=MYZPE13164_O_EIv2.1_0000010.1;Parent=MYZPE13164_O_EIv2.1_0000010;Name=MYZPE13164_O_EIv2.1_0000010.1;Note=MYZPE13164_O_EIv2.1_00000$
scaffold_1      MYZPE13164_O_EIv2.1     three_prime_UTR 8196    8329    .       -       .       ID=MYZPE13164_O_EIv2.1_0000010.1.three_prime_UTR1;Parent=MYZPE13164_O_EIv2.1_0000010.1
scaffold_1      MYZPE13164_O_EIv2.1     exon    8196    9175    .       -       .       ID=MYZPE13164_O_EIv2.1_0000010.1.exon1;Parent=MYZPE13164_O_EIv2.1_0000010.1
scaffold_1      MYZPE13164_O_EIv2.1     CDS     8330    9085    .       -       0       ID=MYZPE13164_O_EIv2.1_0000010.1.CDS1;Parent=MYZPE13164_O_EIv2.1_0000010.1
scaffold_1      MYZPE13164_O_EIv2.1     five_prime_UTR  9086    9175    .       -       .       ID=MYZPE13164_O_EIv2.1_0000010.1.five_prime_UTR1;Parent=MYZPE13164_O_EIv2.1_0000010.1


scaffold_1      MYZPE13164_O_EIv2.1     three_prime_UTR 8196    8329    .       -       .       ID=MYZPE13164_O_EIv2.1_0000010.1.three_prime_UTR1;Parent=MYZPE13164_O_EIv2.1_0000010.1
scaffold_1      MYZPE13164_O_EIv2.1     exon    8196    9175    .       -       .       ID=MYZPE13164_O_EIv2.1_0000010.1.exon1;Parent=MYZPE13164_O_EIv2.1_0000010.1
scaffold_1      MYZPE13164_O_EIv2.1     gene    8196    9175    545     -       .       ID=MYZPE13164_O_EIv2.1_0000010;Name=MYZPE13164_O_EIv2.1_0000010;biotype=transposable_element_gene;confidence=High
scaffold_1      MYZPE13164_O_EIv2.1     mRNA    8196    9175    545     -       .       ID=MYZPE13164_O_EIv2.1_0000010.1;Parent=MYZPE13164_O_EIv2.1_0000010;Name=MYZPE13164_O_EIv2.1_0000010.1;Note=MYZPE13164_O_EIv2.1_00000$
scaffold_1      MYZPE13164_O_EIv2.1     CDS     8330    9085    .       -       0       ID=MYZPE13164_O_EIv2.1_0000010.1.CDS1;Parent=MYZPE13164_O_EIv2.1_0000010.1
scaffold_1      MYZPE13164_O_EIv2.1     five_prime_UTR  9086    9175    .       -       .       ID=MYZPE13164_O_EIv2.1_0000010.1.five_prime_UTR1;Parent=MYZPE13164_O_EIv2.1_0000010.1


scaffold_1      MYZPE13164_O_EIv2.1     three_prime_UTR 8196    8329    .       -       .       ID=MYZPE13164_O_EIv2.1_0000010.1.three_prime_UTR1;Parent=MYZPE13164_O_EIv2.1_0000010.1
scaffold_1      MYZPE13164_O_EIv2.1     exon    8196    9175    .       -       .       ID=MYZPE13164_O_EIv2.1_0000010.1.exon1;Parent=MYZPE13164_O_EIv2.1_0000010.1
scaffold_1      MYZPE13164_O_EIv2.1     gene    8196    9175    545     -       .       ID=MYZPE13164_O_EIv2.1_0000010;Name=MYZPE13164_O_EIv2.1_0000010;biotype=transposable_element_gene;confidence=High
scaffold_1      MYZPE13164_O_EIv2.1     mRNA    8196    9175    545     -       .       ID=MYZPE13164_O_EIv2.1_0000010.1;Parent=MYZPE13164_O_EIv2.1_0000010;Name=MYZPE13164_O_EIv2.1_0000010.1;Note=MYZPE13164_O_EIv2.1_00000$
scaffold_1      MYZPE13164_O_EIv2.1     CDS     8330    9085    .       -       0       ID=MYZPE13164_O_EIv2.1_0000010.1.CDS1;Parent=MYZPE13164_O_EIv2.1_0000010.1
scaffold_1      MYZPE13164_O_EIv2.1     five_prime_UTR  9086    9175    .       -       .       ID=MYZPE13164_O_EIv2.1_0000010.1.five_prime_UTR1;Parent=MYZPE13164_O_EIv2.1_0000010.1

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



WARNING_FILE_NOT_FOUND: Rare Amino Acid analysis: Cannot read protein sequence file '/hpc-home/did23faz/snpEff/./data/m_persicae_O_2_0/protein.fa', nothing done.
ERROR: CDS check file '/hpc-home/did23faz/snpEff/./data/m_persicae_O_2_0/cds.fa' not found.
ERROR: Protein check file '/hpc-home/did23faz/snpEff/./data/m_persicae_O_2_0/protein.fa' not found.
ERROR: Database check failed.

source package /tgac/software/testing/bin/bedops-2.2.0
source package 4028d6e4-21a8-45ec-8545-90e4ed7e1a64
source package /tgac/software/testing/bin/gffread-0.11.4
for scaffold in $(ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/*.fa); do
name=$(echo $scaffold | cut -d '/' -f11 | sed 's@.fa@@g')
echo $name
#To extract whole gene sequences (NOTE: some are very long):
gff2bed < /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.gff > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.bed
grep "gene" /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.bed > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}_2.bed
bedtools getfasta -fi $scaffold -bed /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}_2.bed -fo /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.gff3.nt2 -name+
#To extract CDS only:
gffread -x /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.gff3.nt3 -g $scaffold /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.gff
done 

cat /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/scaffold_*.gff3.nt2 > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.gene.fa
cat /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/scaffold_*.gff3.nt3 > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.CDS.fa
rm -r /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds
```

## P_distance
A distance matrix was calculated for the samples.

p-distance  is the proportion (p) of nucleotide sites at which two sequences being compared are different. It is obtained by dividing the number of nucleotide differences by the total number of nucleotides compared.

distmat calculates the evolutionary distance between every pair of sequences in a multiple sequence alignment. A variety of methods to estimate distance may be selected, and differ in how they correct the observed substitution rates to more accurately reflect the true evolutionary distance. An output file containing a distance matrix for the set of sequences is written. The distances are expressed in terms of the number of substitutions per 100 bases or amino acids.

Generated for with and without ligustri as some analysis require an outgroup and others do not.
#### Whole genome SNPs
```bash
source package /nbi/software/testing/bin/bcftools-1.8
interactive
bcftools stats snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.vcf.gz > snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs_stats.txt

mkdir snp_calling/Myzus/persicae/biello/gatk/p_distance
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.vcf.gz); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=p_dis.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#Submitted batch job 54236748

#Without ligustri:
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

for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=p_dis_193.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#54355967

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
#### Genic SNPs
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
```bash
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree
```
#### Non-synonymous SNPs

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
```bash
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree
```
Collect genes that have non-synonymous SNPs in them:
```bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 \
-b /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions-nonsyn.recode.ann.vcf.gz \
-header > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation-nonsynonymous.gff3

for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions-nonsyn.recode.ann.vcf.gz); do
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/non_synonymous_snps_per_gene3
    OutFile=NA
    GffFile=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_extract_gene_snps.sh $InFile $OutDir $OutFile $GffFile
done #57083837, 57090045,57090519

cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/non_synonymous_snps_per_gene3
for file in *.gz; do bgzip -d "$file"; done
for x in $(ls *.vcf); do
snpno=$(cat $x | wc -l) 
snpno=$((snpno - 2))
echo $x $snpno >> gene_snp_report.txt
done
sed -i 's/ /,/g' gene_snp_report.txt
for file in *.vcf; do bgzip "$file"; done
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae

######################################################################################################################################################
#Second approach
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/non_synonymous_snps_per_gene2
grep 'gene' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation-nonsynonymous.gff3 > gene_lines.gff

while IFS= read -r line; do
    genename=$(echo $line | cut -d '=' -f2 | cut -d ';' -f1)
    grep '##gff-version\|##sequence-region' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation-nonsynonymous.gff3 > ${genename}.gff3
    sleep 3s
    grep "$genename" /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation-nonsynonymous.gff3 > temp_gene_line.gff3 >> ${genename}.gff3
    sleep 2s
    genename=$(echo $line | cut -d '=' -f2 | cut -d ';' -f1)
    echo $genename >> logs/bedtools_intersect.txt
    vcf=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic_mac1-regions-nonsyn.recode.ann.vcf.gz
    OutFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/non_synonymous_snps_per_gene2/${genename}_snps.vcf
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    Jobs=$(squeue -u did23faz | grep 'bedtools' | wc -l)
    while [ $Jobs -gt 99 ]; do
    sleep 15s
    printf "."
    Jobs=$(squeue -u did23faz | grep 'bedtools' | wc -l)
    done 
    sbatch $ProgDir/bedtools_intersect.sh $vcf ${genename}.gff3 $OutFile 2>&1 >> logs/bedtools_intersect.txt
done < gene_lines.gff
```
Plot cumulative frequency of SNPs in genes:
```R
df <- read.table(file = "//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/non_synonymous_snps_per_gene3/gene_snp_report.txt", sep = ',', header = FALSE)

# declaring data points
data_points <- df[,2]

# declaring the break points
break_points = seq(0, 1000, by=1)
# transforming the data
data_transform = cut(data_points, break_points, right=FALSE)
# creating the frequency table
freq_table = table(data_transform)
# printing the frequency table
print("Frequency Table")
print(freq_table)
# calculating cumulative frequency
cumulative_freq = c(0, cumsum(freq_table))
print("Cumulative Frequency")
print(cumulative_freq)
# plotting the data
plot(break_points, cumulative_freq,
     xlab="Data Points",
     ylab="Cumulative Frequency")
# creating line graph
lines(break_points, cumulative_freq)

##############################################################################################################

# Assuming your dataframe is called 'df' and the second column is named 'column2'
# Create breaks for the histogram
breaks <- seq(0, max(df$V2) + 10, by = 10)

# Plot the histogram
hist(df$V2, breaks = breaks, main = "Histogram of gene SNPs - all", xlab = "Values", ylab = "Frequency")
hist(df$V2, breaks = breaks, main = "Histogram of gene SNPs - all", xlab = "Values", ylab = "Frequency", ylim = c(0, 20))

hist(df$V2, breaks = breaks, main = "Histogram of gene SNPs - all", xlab = "Values", ylab = "Frequency", xlim = c(0, 1000), ylim = c(0, 100))
hist(df$V2, breaks = breaks, main = "Histogram of gene SNPs - all", xlab = "Values", ylab = "Frequency", xlim = c(0, 100), ylim = c(0, 1000))

# Filter values between 1 and 100
filtered_values <- df$V2[df$V2 >= 1 & df$V2 <= 100]

# Create breaks for the histogram
breaks <- seq(0, 100, by = 10)

# Plot the histogram
hist(filtered_values, breaks = breaks, main = "Histogram of gene SNP frequency (Values between 1 and 100)", xlab = "Values", ylab = "Frequency")
```
Plot gene length against SNP count
```bash
#find gene length:
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/non_synonymous_snps_per_gene3
for line in $(cat gene_snp_report.txt); do
gene=$(echo $line | cut -d '_' -f1,2,3,4)
length=$(sed -n '2p' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas/het_${gene}.1_CDS*.fa | wc -c)
echo ${line},${length} >> gene_snp_report2.txt
done
#NOTE: this will find the length of one spliec variant of each gene only
sed -i 's/_snps.vcf//g' gene_snp_report2.txt

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3

MYZPE13164_O_EIv2.1_0002220
MYZPE13164_O_EIv2.1_0037470
MYZPE13164_O_EIv2.1_0037650

MYZPE13164_O_EIv2.1_0000020,1,1933
MYZPE13164_O_EIv2.1_0000050,5,3364
MYZPE13164_O_EIv2.1_0000070,1,853
```
```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the main CSV file into a pandas DataFrame
data = pd.read_csv('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/non_synonymous_snps_per_gene3/gene_snp_report2.txt', header=None, names=['ID', 'SNPs', 'Length'])

# Read the subset file into a pandas DataFrame
subset_data = pd.read_csv('/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/effector_candidates.txt', header=None, names=['ID'])

# Filter the main data based on the subset of IDs
subset_filtered_data = data[data['ID'].isin(subset_data['ID'])]

# Define the bins for histogram
bins = [0, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86, 91, 96, 101, float('inf')]

# Group data based on bins and count rows in each group for both filtered and original data
grouped_data = data.groupby(pd.cut(data['SNPs'], bins=bins)).size()
subset_grouped_data = subset_filtered_data.groupby(pd.cut(subset_filtered_data['SNPs'], bins=bins)).size()

# Plotting the histogram for both original and filtered data
plt.figure(figsize=(10, 6))
plt.bar(range(len(grouped_data)), grouped_data, color='skyblue', edgecolor='black', label='All genes with non-synonymous SNPs')
plt.bar(range(len(subset_grouped_data)), subset_grouped_data, color='orange', edgecolor='black', label='Candidate effector genes with non-synonymous SNPs')
plt.title('Number of Non-synonymous SNPs per gene')
plt.xlabel('Number of SNPs')
plt.ylabel('Number of genes')
plt.yscale('log')
plt.xticks(range(len(grouped_data)), [f'{x.left}-{x.right - 1}' if x.right != float('inf') else f'>{x.left - 1}' for x in grouped_data.index], rotation=45)
plt.legend()
plt.savefig('Non-synonymousSNPspergenebarplot.png')


############################################################################################

# Calculate the new column values based on the formula: 1 - log10(Value2 / Value3)
data['Ratio'] = data['Length'] / data['SNPs'] 

# Save the updated DataFrame back to the CSV file
data.to_csv('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/non_synonymous_snps_per_gene3/gene_snp_report3.txt', index=False, header=False)
############################################################################################

# Read the main CSV file into a pandas DataFrame
data = pd.read_csv('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/non_synonymous_snps_per_gene3/gene_snp_report3.txt', header=None, names=['ID', 'SNPs', 'Length','Ratio'])

# Read the subset file into a pandas DataFrame
subset_data = pd.read_csv('/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/effector_candidates.txt', header=None, names=['ID'])

# Filter the main data based on the subset of IDs
subset_filtered_data = data[data['ID'].isin(subset_data['ID'])]

# Define the bins for histogram
bins = [0, 101, 201, 301, 401, 501, 601, 701, 801, 901, 1001, 1101, 1201, 1301, 1401, 1501, 1601, 1701, 1801, 1901, 2001, 2101, 2201, 2301, 2401, 2501, 2601, 2701, 2801, 2901, 3001, 4001, 5001, 6001, 7001, 8001, float('inf')]

# Group data based on bins and count rows in each group for both filtered and original data
grouped_data = data.groupby(pd.cut(data['Ratio'], bins=bins)).size()
subset_grouped_data = subset_filtered_data.groupby(pd.cut(subset_filtered_data['Ratio'], bins=bins)).size()

# Plotting the histogram for both original and filtered data
plt.figure(figsize=(10, 6))
plt.bar(range(len(grouped_data)), grouped_data, color='skyblue', edgecolor='black', label='All genes with non-synonymous SNPs')
plt.bar(range(len(subset_grouped_data)), subset_grouped_data, color='orange', edgecolor='black', label='Candidate effector genes with non-synonymous SNPs')
plt.title('Number of Non-synonymous SNPs per gene')
plt.xlabel('Gene Length / Number of SNPs')
plt.ylabel('Number of genes')
plt.yscale('log')
plt.xticks(range(len(grouped_data)), [f'{x.left}-{x.right - 1}' if x.right != float('inf') else f'>{x.left - 1}' for x in grouped_data.index], rotation=45)
plt.legend()
plt.savefig('Non-synonymousSNPspergenebarplot-ratio.png')

```
```python
import matplotlib.pyplot as plt
import csv
import numpy as np  # Import NumPy for log scaling

# Define the CSV file name
csv_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/non_synonymous_snps_per_gene3/gene_snp_report2.txt"

# Define the file containing values to highlight
highlight_file = "/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/effector_candidates.txt"

# Read the values to highlight from the file
with open(highlight_file, 'r') as highlight_file:
    highlight_values = [line.strip() for line in highlight_file]

# Initialize empty lists to store data from the CSV file
x_values = []
y_values = []
filenames = []

# Read data from the CSV file
with open(csv_file, 'r') as file:
    csv_reader = csv.reader(file)
    for row in csv_reader:
        # Assuming the CSV file format is consistent with the provided example
        filename, x, y = row
        filenames.append(filename)
        x_values.append(int(y))
        y_values.append(int(x))

# Create a list to hold the colors for each point
colors = ['red' if filename in highlight_values else 'black' for filename in filenames]

# Create a scatter plot for blue points first (with a size parameter)
plt.scatter(x_values, y_values, c='black', label='Non-Effector', s=5)

# Create a scatter plot for red points on top (with a size parameter)
red_indices = [i for i, color in enumerate(colors) if color == 'red']
plt.scatter([x_values[i] for i in red_indices], [y_values[i] for i in red_indices], c='red', label='Effector candidates', s=5)

# Label the axes
plt.ylabel("Number of non-synonymous SNPs")
plt.xlabel("Protein coding length")
plt.title("Non-synonymous vs CDS length")

# Add a legend to distinguish between highlighted and non-highlighted points
plt.legend()

# Save the plot as an image file (e.g., PNG)
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/non-synonymous-vs-CDS-length.png")
```
```python
import matplotlib.pyplot as plt
import csv
import numpy as np  # Import NumPy for log scaling

# Define the CSV file name
csv_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/non_synonymous_snps_per_gene3/gene_snp_report2.txt"

# Define the file containing values to highlight
highlight_file = "/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/effector_candidates.txt"

# Read the values to highlight from the file
with open(highlight_file, 'r') as highlight_file:
    highlight_values = [line.strip() for line in highlight_file]

# Initialize empty lists to store data from the CSV file
x_values = []
y_values = []
filenames = []

# Read data from the CSV file
with open(csv_file, 'r') as file:
    csv_reader = csv.reader(file)
    for row in csv_reader:
        # Assuming the CSV file format is consistent with the provided example
        filename, x, y = row
        filenames.append(filename)
        x_values.append(int(y))
        y_values.append(int(x))

# Create a list to hold the colors for each point
colors = ['red' if filename in highlight_values else 'black' for filename in filenames]

# Create a scatter plot for blue points first (with a size parameter)
plt.scatter(x_values, y_values, c='black', label='Non-Effector', s=5)

# Create a scatter plot for red points on top (with a size parameter)
red_indices = [i for i, color in enumerate(colors) if color == 'red']
plt.scatter([x_values[i] for i in red_indices], [y_values[i] for i in red_indices], c='red', label='Effector candidates', s=5)

# Label the axes
plt.ylabel("Number of non-synonymous SNPs")
plt.xlabel("Protein coding length (log scale)")
plt.title("Non-synonymous vs CDS length")

# Use a log scale for both the x-axis and y-axis
plt.xscale('log')
#plt.yscale('log')

# Add a legend to distinguish between highlighted and non-highlighted points
plt.legend()

# Save the plot as an image file (e.g., PNG)
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/non-synonymous-vs-CDS-length2.png")
```
```python
import matplotlib.pyplot as plt
import csv
import numpy as np  # Import NumPy for log scaling

# Define the CSV file name
csv_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/non_synonymous_snps_per_gene3/gene_snp_report2.txt"

# Define the file containing values to highlight
highlight_file = "/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/effector_candidates.txt"

# Read the values to highlight from the file
with open(highlight_file, 'r') as highlight_file:
    highlight_values = [line.strip() for line in highlight_file]

# Initialize empty lists to store data from the CSV file
x_values = []
y_values = []
filenames = []

# Read data from the CSV file
with open(csv_file, 'r') as file:
    csv_reader = csv.reader(file)
    for row in csv_reader:
        # Assuming the CSV file format is consistent with the provided example
        filename, x, y = row
        filenames.append(filename)
        x_values.append(int(y))
        y_values.append(int(x))

# Create a list to hold the colors for each point
colors = ['red' if filename in highlight_values else 'black' for filename in filenames]

# Create a scatter plot for blue points first (with a size parameter)
plt.scatter(x_values, y_values, c='black', label='Non-Effector', s=5)

# Create a scatter plot for red points on top (with a size parameter)
red_indices = [i for i, color in enumerate(colors) if color == 'red']
plt.scatter([x_values[i] for i in red_indices], [y_values[i] for i in red_indices], c='red', label='Effector candidates', s=5)

# Label the axes
plt.ylabel("Number of non-synonymous SNPs (log scale)")
plt.xlabel("Protein coding length (log scale)")
plt.title("Non-synonymous vs CDS length")

# Use a log scale for both the x-axis and y-axis
plt.xscale('log')
plt.yscale('log')

# Add a legend to distinguish between highlighted and non-highlighted points
plt.legend()

# Save the plot as an image file (e.g., PNG)
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/non-synonymous-vs-CDS-length3.png")
```
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

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
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
```bash
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree
```
## SNP rarefaction curve
### Corehunter (redundant)
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
```bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
```
```bash
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
```bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
```
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
### True Rarefaction curves
#### Whole genome SNPs
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
#### Genic SNPs
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
#### For Bass and JIC samples seperately - Whole genome and genic

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
### Fitting a curve
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
#### Nei gene diversity
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
#### Gemma
```bash
perl ~/git_repos/Scripts/NBI/bcf2bbgeno.pl -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.bbgeno -p H-W -s -r
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/
gzip 193s.M_persicae.onlySNPs-genic-regions.bbgeno
#57527288
```
```bash
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz | awk 'BEGIN{OFS="\t"} {if ($3 == ".") $3 = $1"_"$2"_"$5; print}' | gzip > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.vcf.gz

zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.vcf.gz | sed 's/scaffold_//g' | gzip > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.tmp.vcf.gz && mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.tmp.vcf.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.vcf.gz

##################################################################################################################
source package 01ef5a53-c149-4c9e-b07d-0b9a46176cc0
bgzip -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf.gz

zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf.gz | awk 'BEGIN{OFS="\t"} {if ($3 == ".") $3 = $1"_"$2"_"$5; print}' | bgzip -c > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf.gz

#57451298

#there are 360 scaffolds in the assembly, however plink cannot process more than 95 chromosome, therefore will use only the large scaffolds 1-6 for admixture analysis.
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf.gz | head -n 5 > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf.gz | grep -w '##ALT\|##INFO\|##FORMAT\|#CHROM\|##bcftools_callVersion\|##bcftools_callCommand\|##bcftools_concatVersion\|##bcftools_concatCommand\|##bcftools_filterVersion\|##bcftools_filterCommand\|##bcftools_viewVersion\|##bcftools_viewCommand\|scaffold_1\|scaffold_2\|scaffold_3\|scaffold_4\|scaffold_5\|scaffold_6' | sed 's@scaffold_@chr@g' >> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/
bgzip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf
#57451608

source package /nbi/software/testing/bin/plink-1.9 
plink --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf.gz --make-bed --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod
#57452909
plink --bfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod --genome --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod
#57501018

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
```
```python
import pandas as pd
import numpy as np

ibd_data = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod.genome", delim_whitespace=True)
df = ibd_data[['FID1', 'IID1', 'FID2', 'IID2', 'PI_HAT']]
kinship_matrix = df.pivot_table(values='PI_HAT', index=['FID1', 'IID1'], columns=['FID2', 'IID2']).values
kinship_matrix = np.nan_to_num(kinship_matrix)
kinship_matrix -= np.mean(kinship_matrix)
np.savetxt("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/kinship_matrix.txt", kinship_matrix, fmt='%.6f', delimiter='\t')
```
Graphs were plotted in excel
#### vcftools plot of saturation curve
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





```bash
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/step-size-one

cp batcher-core-mperc-stepsize-1.py batcher-core-mperc-stepsize-1-2.py
nano batcher-core-mperc-stepsize-1-2.py

```
```
ls step-size-one/
ALL-BATCH.sh #submits each to HPC
batcher-core-mperc-stepsize-1.py
melt-mperc-core-table.py #compling .csv outputs into one
p_dis_mperc.csv #input

lablab_core_sel_stepsize1.csv #output
mperc_core_sel_stepsize1.csv
mperc_core_sel_stepsize1.xlsx
saturation/
mperc_core_sel_stepsize1.csv           mperc_core_sel_stepsize1_curve.xlsx
mperc_core_sel_stepsize1_curve.csv     saturation-curve.py



ls steps-faster/
ALL-BATCH.sh       
batcher-core-mperc-stepsize-1.py                             
melt-mperc-core-table.py          
p_dis_mperc.csv  

CORE/
0_core_sel.csv  style-gen-core-lablab.py 


import sys
import os
import re
import itertools
import glob
import jitu

```

## Breadth and Depth/Coverage
Download Singh data:
```bash
C:\Users\did23faz\Downloads\sratoolkit.current-win64\sratoolkit.3.0.5-win64\bin>prefetch --max-size 1t --option-file C:\Users\did23faz\Documents\accessions.txt --output-directory \\jic-hpc-data\Group-Scratch\Saskia-Hogenhout\tom_heaven\Aphididae\raw_data\Myzus\persicae\singh\download
```
```bash
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/singh/download
nano key.txt #SRR and sample number mapping
source package /nbi/software/testing/bin/sratoolkit-2.9.0

for dir in $(ls -d *); do
ID=$(grep $dir key.txt | cut -d $'\t' -f2 | sed 's@S110@NIC_23@g' | sed 's@S111@NIC_410G@g' | sed 's@S106@NIC_5191A@g'| sed 's@S114@NIC_57@g'| sed 's@S115@NIC_8124@g'| sed 's@S112@NIC_926B@g'| sed 's@S105@SUS_4106a@g'| sed 's@S108@SUS_4225A@g'| sed 's@S107@SUS_NS@g'| sed 's@S109@SUS_US1L@g')
if [ -e ../$ID/qualimap ]; then
echo $ID
fastq-dump --split-files --gzip -O ../$ID $dir/*.sra
else 
echo $dir
fi
done

for file in $(ls ../*/*_1.fastq.gz); do
echo $file
SRR=$(echo $file | cut -d '/' -f3 | cut -d '_' -f1)
rm -r $SRR
done
```
Singh raw data has subsequently been deleted from the HPC save space as it is backed up with NCBI.

```bash
for sample in $(cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep -v '#' | cut -d ',' -f2); do
    ls -d raw_data/Myzus/persicae/*/${sample}
done
awk '{print $2}' snp_calling/Myzus/persicae/biello/PCA_file_host.csv


ls: cannot access raw_data/Myzus/persicae/*/NIC_23: No such file or directory   Italy 14,565,309 SRR13326441, SRR13326434, SRR10199535, SRR10199534
ls: cannot access raw_data/Myzus/persicae/*/NIC_410G: No such file or directory Greece 16,030,966 SRR10199540, SRR10199536
ls: cannot access raw_data/Myzus/persicae/*/NIC_57: No such file or directory   Italy  16,512,008 SRR13326441, SRR13326434, SRR10199535, SRR10199534
ls: cannot access raw_data/Myzus/persicae/*/NIC_8124: No such file or directory Italy 14,252,463 SRR13326441, SRR13326434, SRR10199535, SRR10199534
ls: cannot access raw_data/Myzus/persicae/*/NIC_926B: No such file or directory Greece 15,408,286 SRR10199540, SRR10199536
ls: cannot access raw_data/Myzus/persicae/*/NIC_410R: No such file or directory Zimbabwe 15,834,607 SRR10199537
ls: cannot access raw_data/Myzus/persicae/*/NIC_5191A: No such file or directory    Zimbabwe 39,307,583 SRR10199537
ls: cannot access raw_data/Myzus/persicae/*/SUS_4106a: No such file or directory    UK 38,838,364 SRR10199546, SRR10199544, SRR10199542
ls: cannot access raw_data/Myzus/persicae/*/SUS_4225A: No such file or directory  Italy  14,657,289 SRR10199541
ls: cannot access raw_data/Myzus/persicae/*/SUS_US1L: No such file or directory UK R2=14,445,344 SRR10199546,  SRR10199544, SRR10199542
ls: cannot access raw_data/Myzus/persicae/*/SUS_1X: No such file or directory   Italy  14,657,289 SRR10199541 > Algeria SRR13326486
ls: cannot access raw_data/Myzus/persicae/*/SUS_NS: No such file or directory   Germany 15,833,996 SRR10199543

SRR13326486, SRR13326389, SRR13326456, ERX1223986

ls: cannot access raw_data/Myzus/persicae/*/NL_IRS: No such file or directory
ls: cannot access raw_data/Myzus/persicae/*/NL_WUR: No such file or directory

ls: cannot access raw_data/Myzus/persicae/*/A156: No such file or directory

for sample in $(ls -d raw_data/Myzus/persicae/*/* | rev | cut -d '/' -f1 | rev); do
x=$(cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep -w "$sample" | cut -d ',' -f2)
if [ -z "$x" ]; then
  echo "$sample doesnt exist"
fi
done

S100 doesnt exist
S105 doesnt exist
S106 doesnt exist
S107 doesnt exist
S108 doesnt exist
S109 doesnt exist
S110 doesnt exist
S111 doesnt exist
S112 doesnt exist
S114 doesnt exist
S115 doesnt exist
S122 doesnt exist
S21 doesnt exist
S35 doesnt exist
S41 doesnt exist
S82 doesnt exist
BE2_batch1 doesnt exist
BE44 doesnt exist
BE84 doesnt exist
BE85 doesnt exist
BE86 doesnt exist
BE87 doesnt exist
C11 doesnt exist
C13 doesnt exist
C16 doesnt exist
C18 doesnt exist
C19 doesnt exist
C25 doesnt exist
C4 doesnt exist
C5 doesnt exist
ES14 doesnt exist
FR44 doesnt exist
FR6 doesnt exist
FRC09 doesnt exist
HUN12 doesnt exist
IT102 doesnt exist
IT68 doesnt exist
MISC15 doesnt exist
MISC18 doesnt exist
MISC19 doesnt exist
MISC20 doesnt exist
MISC21 doesnt exist
MISC23 doesnt exist
MISC26 doesnt exist
NL13 doesnt exist
NL53 doesnt exist
RUS14 doesnt exist
S135 doesnt exist
S145 doesnt exist
S258 doesnt exist
UK32 doesnt exist
UK35 doesnt exist
UK37 doesnt exist

du -Lh /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/*/*/*.gz | sort -rh | head -n 1
#21G     /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/wouters/ligustri/M_lig_095_cat_run1_run2_R2_val_2.fq.gz
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/wouters/*); do
    if [ ! -e ${ReadDir}/qualimap/*genome_results_gff.txt ]; then
    echo Running for:
    ls ${ReadDir}/qualimap/*genome_results_gff.txt
    Fread=$(ls $ReadDir/*_1.fq.gz)
    Rread=$(ls $ReadDir/*_2.fq.gz)
    Fread2=$(ls $ReadDir/*_3.fq.gz)
    Rread2=$(ls $ReadDir/*_4.fq.gz)
    OutDir=$(echo $ReadDir)
    Reference_genome=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
    Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3
    ProgDir=~/git_repos/Wrappers/NBI
    Jobs=$(squeue -u did23faz| grep 'qualimap'  | wc -l)
    echo x
    while [ $Jobs -gt 19 ]; do
        sleep 300s
        printf "."
        Jobs=$(squeue -u did23faz| grep 'qualimap'  | wc -l)
    done
    echo $ReadDir >> logs/raw_qualimap_report.txt
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread $Fread2 $Rread2 2>&1 >> logs/raw_qualimap_report.txt
    else
    echo Already Done:
    ls ${ReadDir}/qualimap/*genome_results.txt
    fi
done

source package /tgac/software/testing/bin/seqkit-0.10.0
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/*/*); do
sample=$(echo $ReadDir | rev | cut -d '/' -f1 | rev)
coverage=$(grep 'mean coverageData' $ReadDir/qualimap/*genome_results.txt | rev | cut -d ' ' -f1 | rev | sed 's@X@@g')
if (( $(echo "$coverage < 15" | bc -l) )); then
    echo $sample has low coverage: $coverage >> Reports/Low_coverage_report2.txt
fi
target_coverage=15
subsampling_ratio=$(bc <<< "scale=2; $target_coverage / $coverage")
for reads in $(ls $ReadDir/*.fq.gz); do
OutFile=$(basename $reads | cut -d '.' -f1)_subsampled.fq.gz
OutDir=${ReadDir}/subsampled
OutDir2=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/M_persicae/${sample}
mkdir $OutDir
mkdir $OutDir2
echo $sample raw reads have average coverage of ${coverage}, these were subsampled with ratio $subsampling_ratio >> Reports/Subsampling_report.txt
seqkit sample -s 1234 -p $subsampling_ratio -o ${OutDir2}/${OutFile} $reads
ln -s ${OutDir2}/${OutFile} ${OutDir}/${OutFile}
done
if [ -e ${OutDir2}/${OutFile} ]; then
    echo ${sample} Done
#rm -r /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/M_persicae_Bass_SRA/$sample
else 
rm ${OutDir}/*_subsampled.fq.gz
for reads in $(ls $ReadDir/*.fq.gz); do
name=$(echo $reads | rev | cut -d '/' -f1 | rev)
cp $reads ${OutDir2}/${name}
ln -s ${OutDir2}/${name} ${OutDir}/${name}
#rm -r /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/M_persicae_Bass_SRA/$sample
done
fi
done
#NOTE: some samples used in the final SNP calling analysis have low coverage

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/*/*/subsampled); do
    if [ ! -e ${ReadDir}/qualimap/*genome_results.txt ]; then
    echo Running for:
    ls ${ReadDir}/qualimap/*genome_results.txt
    Fread=$(ls $ReadDir/*_1*fq.gz)
    Rread=$(ls $ReadDir/*_2*fq.gz)
    Fread2=$(ls $ReadDir/*_3*fq.gz)
    Rread2=$(ls $ReadDir/*_4*fq.gz)
    OutDir=$(echo $ReadDir)
    Reference_genome=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
    ProgDir=~/git_repos/Wrappers/NBI
    Jobs=$(squeue -u did23faz| grep 'qualimap'  | wc -l)
    echo x
    while [ $Jobs -gt 19 ]; do
        sleep 300s
        printf "."
        Jobs=$(squeue -u did23faz| grep 'qualimap'  | wc -l)
    done
    echo $ReadDir >> logs/raw_qualimap_report.txt
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Fread $Rread $Fread2 $Rread2 2>&1 >> logs/raw_qualimap_report.txt
    else
    echo Already Done:
    ls ${ReadDir}/qualimap/*genome_results.txt
    fi
done

for sample in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/*/*/subsampled/qualimap/*genome_results.txt); do
cov=$(grep 'mean coverageData' $sample | rev | cut -d ' ' -f1 | rev )
name=$(echo $sample | cut -d '/' -f12)
echo $name $cov
done

#NOTE: to save space the script has been edited to delete input files upon production of outputs.
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/*/*/subsampled); do
    sample=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
    Fread=$(ls $ReadDir/*_1*fq.gz)
    Rread=$(ls $ReadDir/*_2*fq.gz)
    Fread2=$(ls $ReadDir/*_3*fq.gz)
    Rread2=$(ls $ReadDir/*_4*fq.gz)
    OutDir=/jic/research-groups/Saskia-Hogenhout/TCHeaven/dna_qc/M_persicae/${sample}
    OutFile=${sample}_subsampled_trimmed
    Quality=20
    Length=50
    ProgDir=~/git_repos/Wrappers/NBI
    Jobs=$(squeue -u did23faz| grep 'trim_g'  | wc -l)
    echo x
    while [ $Jobs -gt 19 ]; do
        sleep 300s
        printf "."
        Jobs=$(squeue -u did23faz| grep 'trim_g'  | wc -l)
    done
    mkdir -p $OutDir
    echo $sample >> logs/trim_galore_report.txt
    sbatch $ProgDir/run_trim_galore.sh $OutDir $OutFile $Quality $Length $Fread $Rread $Fread2 $Rread2  2>&1 >> logs/trim_galore_report.txt
done 

if [ -e ${OutDir}/${OutFile}_1.fq.gz ] && [ -e ${OutDir}/${OutFile}_2.fq.gz ]; then
echo Outputs detected
else
echo Outputs not detected
fi
#Subsampled Trimmed read files have been deleted following successful generation of mitochondrial VCF file, can be regenertated from raw data if needed.

#Symlink files to the main Aphididae directory
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/*/*/subsampled); do
    sample=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
    Dir=/jic/research-groups/Saskia-Hogenhout/TCHeaven/dna_qc/M_persicae/${sample}
    OutDir=$(echo $ReadDir | sed 's@raw_data@dna_qc@g')/trim_galore
    mkdir -p $OutDir
    ln -s /jic/research-groups/Saskia-Hogenhout/TCHeaven/dna_qc/M_persicae/${sample}/* $OutDir/.
done

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/*/*/subsampled/trim_galore); do
    if [ ! -e ${ReadDir}/qualimap/*genome_results.txt ]; then
    echo Running for:
    ls ${ReadDir}/qualimap/*genome_results.txt
    Fread=$(ls $ReadDir/*_1*fq.gz)
    Rread=$(ls $ReadDir/*_2*fq.gz)
    Fread2=$(ls $ReadDir/*_3*fq.gz)
    Rread2=$(ls $ReadDir/*_4*fq.gz)
    OutDir=$(echo $ReadDir)
    Reference_genome=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
    ProgDir=~/git_repos/Wrappers/NBI
    Jobs=$(squeue -u did23faz| grep 'qualimap'  | wc -l)
    echo x
    while [ $Jobs -gt 19 ]; do
        sleep 300s
        printf "."
        Jobs=$(squeue -u did23faz| grep 'qualimap'  | wc -l)
    done
    echo $ReadDir >> logs/raw_qualimap_report.txt
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Fread $Rread $Fread2 $Rread2 2>&1 >> logs/raw_qualimap_report.txt
    else
    echo Already Done:
    ls ${ReadDir}/qualimap/*genome_results.txt
    fi
done

for sample in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/*/*/subsampled/trim_galore/qualimap/*genome_results.txt); do
cov=$(grep 'mean coverageData' $sample | rev | cut -d ' ' -f1 | rev )
name=$(echo $sample | cut -d '/' -f12)
echo $name $cov
done
```
MOVE FILES FROM /jic/research-groups/Saskia-Hogenhout/ TO /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae AS THIS IS THE CHEAP SPACE:


The M. persicae mitochondrial reference genome was downloaded from NCBI (Accession: NC_029727).

Reads were aligned to the mitochondrial reference genome:
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/*/*/subsampled/trim_galore); do
    Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/mito/Wang2016_NC_029727.fasta
    Fread=$(ls $ReadDir/*_1*fq.gz)
    Rread=$(ls $ReadDir/*_2*fq.gz)
    OutDir=$(echo $ReadDir | sed 's@subsampled/trim_galore@bwa@g' | sed 's@dna_qc@alignment@g' | sed 's@persicae@persicae/SNPs/mito@g')
    OutFile=$(echo $ReadDir | cut -d '/' -f11,12 | sed 's@/@_@g')_mito
    echo ${OutDir}/${OutFile}
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    Jobs=$(squeue -u did23faz| grep 'bwa'  | wc -l)
    echo x
    while [ $Jobs -gt 9 ]; do
        sleep 300s
        printf "."
        Jobs=$(squeue -u did23faz| grep 'bwa'  | wc -l)
    done
    mkdir -p $OutDir    
    sbatch $ProgDir/bwa-mem.sh $OutDir $OutFile $Reference $Fread $Rread 
done 
```
FYI, Reads previously aligned to the nuclear genome by Roland are here, I do not have the permissions to symlink them into my file system:
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/*/*/subsampled/trim_galore); do
    isolate=$(echo $ReadDir | cut -d '/' -f12)
    isolates=$(echo $ReadDir | cut -d '/' -f11,12)
    echo $isolates
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/SNPs/genomic/${isolates}
    mkdir -p $OutDir
    ln -s /jic/scratch/groups/Saskia-Hogenhout/roland/Myzus_periscae_popgen/Output/realigned_BAM_files/${isolate}.sorted.mark_dups.realigned.bam ${Outdir}/${isolate}.sorted.mark_dups.realigned.bam
    ln -s /jic/scratch/groups/Saskia-Hogenhout/roland/Myzus_periscae_popgen/Output/realigned_BAM_files/${isolate}.sorted.mark_dups.realigned.bam.bai ${Outdir}/${isolate}.sorted.mark_dups.realigned.bam.bai
    ln -s /jic/scratch/groups/Saskia-Hogenhout/roland/Myzus_periscae_popgen/Output/realigned_BAM_files/${isolate}.sorted.mark_dups.realigned_stats ${Outdir}/${isolate}.sorted.mark_dups.realigned_stats
done
```
Files were sorted by scaffold and coordinate level, duplicates were marked and removed and the files were re-indexed:
NOTE: this script does not remove the input file and .bam files take up a lot of space - was therefore actually run 1 step at a time
```bash
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/SNPs/mito/*/*/bwa/*_mito.bam); do
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_sort_bam.sh $file
done 

#Input files were subsequently deleted to save space
```
GATK
```bash
#The mitochondrial reference genome was indexed:
interactive
source switch-institute ei
source package 638df626-d658-40aa-80e5-14a275b7464b
source pilon-1.22
source package /tgac/software/testing/bin/picardtools-2.1.1
cd /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/mito
java -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar CreateSequenceDictionary R=Wang2016_NC_029727.fasta O=Wang2016_NC_029727.dict
samtools faidx Wang2016_NC_029727.fasta

#GATK was used to perform indel realignment:
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/SNPs/mito/*/*/bwa/*MarkDups.bam); do
Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/mito/Wang2016_NC_029727.fasta
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_realign.sh $file $Reference 
done #56494561-56494688
```
Combined into a .vcf and called variants with bcftools:
```bash
source package 638df626-d658-40aa-80e5-14a275b7464b
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/SNPs/mito/*/*/bwa/*realigned.bam > bamlist.txt
bcftools mpileup -b /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/bamlist.txt --annotate AD,DP --fasta-ref /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/mito/Wang2016_NC_029727.fasta -O z -o /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/255s.M_persicae.mito.vcf.gz
sbatch $ProgDir/temp.sh #56496243, 56496558

cp bamlist.txt bamlist2.txt
nano bamlist2.txt #edit to keep 193 samples only #NOTE: SUS_USIL is missing, so there are 192 samples
bcftools mpileup -b /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/bamlist2.txt --annotate AD,DP --fasta-ref /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/mito/Wang2016_NC_029727.fasta -O z -o /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/192s.M_persicae.mito.vcf.gz
sbatch $ProgDir/temp.sh #56609348

source package 638df626-d658-40aa-80e5-14a275b7464b
bcftools call -Oz -v -m -o /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/255.M_persicae.mito.vcf.gz /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/255s.M_persicae.mito.vcf.gz
bcftools call -Oz -v -m -o /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/192.M_persicae.mito.vcf.gz /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/192s.M_persicae.mito.vcf.gz
bcftools call --ploidy 1 -Oz -v -m -o /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/255.M_persicae.mito_hap.vcf.gz /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/255s.M_persicae.mito.vcf.gz
bcftools call --ploidy 1 -Oz -v -m -o /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/192.M_persicae.mito_hap.vcf.gz /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/192s.M_persicae.mito.vcf.gz
sbatch $ProgDir/temp.sh #56624720, 56626500, 56627364 (haploid calls)

bgzip -cd /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/255.M_persicae.mito.vcf.gz | grep -v '#' | wc -l #4,165
bgzip -cd /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/192.M_persicae.mito.vcf.gz | grep -v '#' | wc -l #603

#bam files were subsequently removed to save space
```
Filtering:

I have not done mappability masking for the SNPs, the programs used to do this with the nuclear SNPs do not seem to be installed on the JIC HPC, given that these SNPs are called against only the ~17.5k bp mitochondrial genome which is not expected to have lengthy repetative regions it does not seem worth it to install new software to perform this analysis.
```bash
#Keep only samples that were in the filtered nuclear SNPs .vcf:
bgzip -cd /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/255.M_persicae.mito.vcf.gz > 255.M_persicae.mito.vcf
bgzip -cd /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/255.M_persicae.mito_hap.vcf.gz > 255.M_persicae.mito_hap.vcf
sed -i 's/singh_//g' 255.M_persicae.mito.vcf
sed -i 's/wouters_//g' 255.M_persicae.mito.vcf
sed -i 's/singh_//g' 255.M_persicae.mito_hap.vcf
sed -i 's/wouters_//g' 255.M_persicae.mito_hap.vcf
patterns=("CHROM" "POS" "ID" "REF" "ALT" "QUAL" "FILTER" "INFO" "FORMAT" "K43" "K40" "A166" "I1" "K16" "A102" "NL1" "Q1200" "S2B01" "S2196G" "NL21" "S152" "S232" "UKW3" "NL5" "UK2" "S204" "A138" "Crespys" "S472" "S196" "MG1107" "S481" "Lierida" "ES01" "Generac" "FR15" "MISC5" "MISC10" "MISC14" "ES149" "MISC81" "ES146" "MISC48" "BE1" "ES92" "MISC34" "MISC27" "MISC46" "MISC42" "UK21" "MISC51" "MISC56" "K66" "MISC33" "MISC30" "MISC28" "MISC31" "BE23" "MISC53" "MISC4" "MISC38" "BE2_batch2" "ES88" "BE6" "K88" "UK25" "K87" "FR433" "a456BE" "MISC79" "MISC64" "MISC76" "MISC67" "BE49" "MISC65" "NL93" "MISC60" "UK1" "BE33A" "MISC78" "A161" "MISC69" "UK19" "NL_IRS" "NL_WUR" "SUS_4106a" "23" "4106a" "FRC" "G006" "NIC_410G" "NIC_410R" "NIC_5191A" "NIC_57" "NIC_8124" "NIC_926B" "O" "SUS_1X" "SUS_4225A" "SUS_NS" "SUS_US1L" "UKSB" "US1L" "S2" "S14" "S3" "S17" "S18" "S5" "S4" "S19" "S6" "S7" "S20" "S11" "S23" "S1" "S9" "S16" "S8" "S22" "S15" "S12" "S13" "S28" "S10" "S25" "S27" "S29" "S24" "S33" "S36" "S34" "S26" "S37" "S38" "S30" "S39" "S51" "S32" "S42" "S43" "S31" "S48" "S45" "S50" "S40" "S44" "S46" "S49" "S59" "S47" "S53" "S52" "S56" "S60" "S54" "S64" "S62" "S57" "S63" "S66" "S70" "S58" "S68" "S72" "S61" "S71" "S74" "S73" "S65" "S69" "S76" "S78" "S67" "S77" "S80" "S79" "S86" "S75" "S88" "S101" "S89" "S87" "S85" "S96" "S84" "S118" "S81" "S93" "S91" "S94" "S90" "S125" "S83" "S92" "S124" "S120" "S121" "S113" "S117" "HFH5FBCX2")
#^samples to keep added to keep.txt
source package /nbi/software/testing/bin/vcftools-0.1.15
source package /tgac/software/testing/bin/bcftools-1.9
vcftools --keep keep.txt --vcf 255.M_persicae.mito.vcf --mac 1 --recode --recode-INFO-all --out temp193 #After filtering, kept 571 out of a possible 4165 Sites, After filtering, kept 193 out of 255 Individuals
nano temp193.recode.vcf #edit wouters samples with same names as singh samples to be unique
vcftools --vcf temp193.recode.vcf --remove-indels --mac 1 --recode --recode-INFO-all --out temp193_2 #After filtering, kept 549 out of a possible 571 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193.recode.vcf --remove-indels --max-alleles 2 --mac 1 --recode --recode-INFO-all --out temp193_2 #After filtering, kept 541 out of a possible 571 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --recode --recode-INFO-all --out temp193_2 #After filtering, kept 536 out of a possible 571 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --recode --recode-INFO-all --out temp193_2 #After filtering, kept 496 out of a possible 571 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --recode --recode-INFO-all --out temp193_2 #After filtering, kept 493 out of a possible 571 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --recode --recode-INFO-all --out temp193_2 #After filtering, kept 493 out of a possible 571 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --recode --recode-INFO-all --out temp193_2 #After filtering, kept 0 out of a possible 571 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --max-meanDP 300 --recode --recode-INFO-all --out temp193_2 #After filtering, kept 8 out of a possible 571 Sites, After filtering, kept 193 out of 193 Individuals -> most site have very high depth coverage
bcftools query -f '%CHROM\t%POS\t%INFO/DP[\t%DP]\n' temp193.recode.vcf > ../temp.tsv
vcftools --vcf temp193.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 10 --min-meanDP 40 --maxDP 500 --max-meanDP 300 --recode --recode-INFO-all --out /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/192s.M_persicae.mito #After filtering, kept 490 out of a possible 571 Sites, and kept 193 out of 193 Individuals

vcftools --keep keep.txt --vcf 255.M_persicae.mito_hap.vcf --mac 1 --recode --recode-INFO-all --out temp193_hap #After filtering, kept 165 out of a possible 3709 Sites, After filtering, kept 193 out of 255 Individuals
nano temp193_hap.recode.vcf #edit wouters samples with same names as singh samples to be unique
vcftools --vcf temp193_hap.recode.vcf --remove-indels --mac 1 --recode --recode-INFO-all --out temp193_hap_2 #After filtering, kept 146 out of a possible 165 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193_hap.recode.vcf --remove-indels --max-alleles 2 --mac 1 --recode --recode-INFO-all --out temp193_hap #After filtering, kept 145 out of a possible 165 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193_hap.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --recode --recode-INFO-all --out temp193_hap_2 #After filtering, kept 140 out of a possible 165 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193_hap.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --recode --recode-INFO-all --out temp193_hap_2 #After filtering, kept 133 out of a possible 165 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193_hap.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --recode --recode-INFO-all --out temp193_hap_2 #After filtering, kept 130 out of a possible 165 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193_hap.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --recode --recode-INFO-all --out temp193_hap_2 #After filtering, kept 130 out of a possible 165 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193_hap.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --recode --recode-INFO-all --out temp193_hap_2 #After filtering, kept 0 out of a possible 165 Sites, After filtering, kept 193 out of 193 Individuals
vcftools --vcf temp193_hap.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 10 --min-meanDP 40 --maxDP 500 --max-meanDP 300 --recode --recode-INFO-all --out /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/192s.M_persicae.mito_hap #After filtering, kept 129 out of a possible 165 Sites, and kept 193 out of 193 Individuals 

bgzip /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/192s.M_persicae.mito.recode.vcf
bgzip /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/192s.M_persicae.mito_hap.recode.vcf

#remove g006 from analysis:
vcftools --keep keep2.txt --vcf 255.M_persicae.mito_hap.vcf --mac 1 --recode --recode-INFO-all --out temp192_hap #After filtering, kept 165 out of a possible 3709 Sites, After filtering, kept 192 out of 255 Individuals
nano temp192_hap.recode.vcf #edit wouters samples with same names as singh samples to be unique
vcftools --vcf temp192_hap.recode.vcf --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 10 --min-meanDP 40 --maxDP 500 --max-meanDP 300 --recode --recode-INFO-all --out /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/191s.M_persicae.mito_hap #After filtering, kept 129 out of a possible 165 Sites, and kept 192 out of 192 Individuals 
bgzip /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/191s.M_persicae.mito_hap.recode.vcf
```
```bash
for vcf in $(ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/191s.M_persicae.mito_hap.recode.vcf.gz); do
echo $vcf
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance
Outfile=mito_p_dis_-g006.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done #56651077, 56651731

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
```
```python
from collections import defaultdict
# load distance matrix
acc = defaultdict(list)
acc_order = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/mito_p_dis_-g006.mat') as inp:
    next(inp)
    for line in inp:
        A = line.strip().split()
        acc[A[0]] = A[1:]
        acc_order.append(A[0])

assert len(acc) == 193

# write the distance matrix in a desired format
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/mito_p_dis_-g006.csv", "w") as outf:
    outf.write(",".join(["ID"] + acc_order) + "\n")
    for a in acc_order:
        outf.write(",".join([a] + acc[a]) + "\n")
```
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

# Example usage:
csv_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/mito_p_dis_-g006.csv'
output_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/mito_p_dis_-g006.dist'

distance_matrix = read_distance_matrix_from_csv(csv_file)
convert_to_phylip(distance_matrix, output_file)
```
```bash
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree
#Import failed: unknown format or error in file (eg. illegal characters or multiple occurrences of same taxon name)
#NOTE: G006 is triggering the above error, the mitochondrial sequence used as reference is from G006 -> there are therefore no SNPs for this sample causing the error? G006 was therefore reomved during filtering
```

```bash
rm /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/SNPs/mito/*/*/bwa/*.sam
#################################################################################################

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/SNPs/mito/*/*/bwa/*_mito_sorted.bam); do
OutFile=/jic/research-groups/Saskia-Hogenhout/TCHeaven/temp/$(basename $file | sed 's@_mito_sorted.bam@_mitochondrial_sorted_MarkDups.bam@g')
if [ ! -e $OutFile ]; then
echo $OutFile
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/temp.sh $file
else
echo $OutFile exists
fi
done

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/SNPs/mito/*/*/bwa/*_mito.bam); do
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_sort_bam.sh $file
done 

interactive
source package 638df626-d658-40aa-80e5-14a275b7464b
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/SNPs/mito/*/*/bwa/*_mito.bam); do
name=$(basename $file | cut -d '.' -f1)
OutDir=$(dirname $file)
samtools sort -o ${OutDir}/${name}_sorted.sam -T ${name}_1234 $file
samtools view -bS ${OutDir}/${name}_sorted.sam > ${OutDir}/${name}_sorted.bam
rm ${OutDir}/${name}_sorted.sam
rm $file
done

source switch-institute ei
source pilon-1.22
source package /tgac/software/testing/bin/picardtools-2.1.1

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/SNPs/mito/*/*/bwa/*_sorted.bam); do
OutDir=$(dirname $file)
OutFile=$(basename $file | sed 's@mito@mitochondrial@g')
java -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar SortSam I=$file O=${OutDir}/${OutFile} SORT_ORDER=coordinate
if [ -e ${OutDir}/${OutFile} ]; then
rm $file
else
echo Output missing
fi
done

#!/bin/bash
#SBATCH --job-name=bam_sort
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=1
#SBATCH -p jic-medium,nbi-medium
#SBATCH --time=02-00:00

#Collect inputs
CurPath=$PWD
WorkDir=$PWD${TMPDIR}_${SLURM_JOB_ID}
InFile=$1

echo CurPth:
echo $CurPath
echo WorkDir:
echo $WorkDir
echo InFile:
echo $InFile

source switch-institute ei
source package 638df626-d658-40aa-80e5-14a275b7464b
source pilon-1.22
source package /tgac/software/testing/bin/picardtools-2.1.1

name=$(basename $InFile | cut -d '.' -f1)
OutDir=$(dirname $InFile)
samtools sort -o ${OutDir}/${name}_sorted.sam -T ${name}_1234 $file
samtools view -bS ${OutDir}/${name}_sorted.sam > ${OutDir}/${name}_sorted.bam
if [ -e ${OutDir}/${name}_sorted.bam ]; then
rm ${OutDir}/${name}_sorted.sam
rm $InFile
else
echo Output missing
fi


OutFile=$(basename ${OutDir}/${name}_sorted.bam | sed 's@mito@mitochondrial@g')
java -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar SortSam I=${OutDir}/${name}_sorted.bam O=${OutDir}/${OutFile} SORT_ORDER=coordinate
if [ -e ${OutDir}/${OutFile} ]; then
rm ${OutDir}/${name}_sorted.bam
else
echo Output2 missing
fi

OutFile2=$(echo ${OutDir}/${OutFile} | sed 's@.bam@_MarkDups.bam@g')
java -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar MarkDuplicates I=${OutDir}/${OutFile} O=${OutFile2} -M ${OutDir}/${name}_marked_dup_metrics.txt
if [ -e ${OutFile2} ]; then
rm ${OutDir}/${OutFile}
else
echo Output3 missing
fi

cd ${OutDir}
samtools index ${OutFile2}

echo DONE

```
#### Mitochondrial blast
```bash
source package d6092385-3a81-49d9-b044-8ffb85d0c446
mkdir M_persicae_mitochondria_blast_db/
cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mito1.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mito2.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mito3.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mito4.fasta > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mito.fasta
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mito.fasta -dbtype nucl -parse_seqids -out M_persicae_mitochondria_blast_db/blast_db
blastn -query /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa -db M_persicae_mitochondria_blast_db/blast_db -out results.txt -evalue 1e-5 -outfmt 6 -num_threads 1

```
#### Extract CathB SNPs
```bash
head -n 396 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf > CathB_snps.vcf
awk -F '\t' '$1 == "scaffold_1" && $2 >= 17900000 && $2 <= 17960000' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf >> CathB_snps.vcf
awk -F '\t' '$1 == "scaffold_1" && $2 >= 66140000 && $2 <= 66170000' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf >> CathB_snps.vcf
awk -F '\t' '$1 == "scaffold_1" && $2 >= 90330000 && $2 <= 90360000' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf >> CathB_snps.vcf
awk -F '\t' '$1 == "scaffold_1" && $2 >= 98550000 && $2 <= 98700000' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf >> CathB_snps.vcf
awk -F '\t' '$1 == "scaffold_2" && $2 >= 36300000 && $2 <= 36350000' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf >> CathB_snps.vcf
awk -F '\t' '$1 == "scaffold_2" && $2 >= 57500000 && $2 <= 57540000' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf >> CathB_snps.vcf
awk -F '\t' '$1 == "scaffold_2" && $2 >= 66850000 && $2 <= 66900000' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf >> CathB_snps.vcf
awk -F '\t' '$1 == "scaffold_4" && $2 >= 8200000 && $2 <= 8250000' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf >> CathB_snps.vcf
awk -F '\t' '$1 == "scaffold_4" && $2 >= 12550000 && $2 <= 12630000' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf >> CathB_snps.vcf
awk -F '\t' '$1 == "scaffold_4" && $2 >= 47700000 && $2 <= 47850000' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf >> CathB_snps.vcf
awk -F '\t' '$1 == "scaffold_134"' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf >> CathB_snps.vcf
```
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/*/*); do
    if [ ! -e ${ReadDir}/qualimap/*genome_results.txt ]; then
    echo Running for:
    ls ${ReadDir}/qualimap/*genome_results.txt
    Fread=$(ls $ReadDir/*_1.fq.gz)
    Rread=$(ls $ReadDir/*_2.fq.gz)
    Fread2=$(ls $ReadDir/*_3.fq.gz)
    Rread2=$(ls $ReadDir/*_4.fq.gz)
    OutDir=$(echo $ReadDir)
    Reference_genome=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
    ProgDir=~/git_repos/Wrappers/NBI
    Jobs=$(squeue -u did23faz| grep 'bwa'  | wc -l)
    echo x
    while [ $Jobs -gt 19 ]; do
        sleep 300s
        printf "."
        Jobs=$(squeue -u did23faz| grep 'bwa'  | wc -l)
    done
    echo $ReadDir >> logs/raw_bwamem_report.txt
    sbatch $ProgDir/run_bwa-mem.sh $OutDir $Reference_genome $Fread $Rread $Fread2 $Rread2 2>&1 >> logs/raw_bwamem_report.txt
    else
    echo Already Done:
    ls ${ReadDir}/qualimap/*genome_results.txt
    fi
done
```