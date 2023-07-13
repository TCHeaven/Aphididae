# Contructing phylogenies as part of the M. persicae population genomics project

## Extract gene sequences
```bash
srun -p jic-short --ntasks-per-node=8 --mem 60G --nodes=1 singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif gff3_to_fasta -g /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 -f /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa -st gene -d simple -o /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.fa

source package /tgac/software/testing/bin/bedops-2.2.0
source package 4028d6e4-21a8-45ec-8545-90e4ed7e1a64
gff2bed < /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 > output.bed
grep "gene" output.bed > genes.bed
bedtools getfasta -fi /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa -bed genes.bed -fo /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt2 -name+
#These methods extract the entire gene sequence, exons and introns, some of which are VERY large.
grep -A 1 'MYZPE13164_O_EIv2.1_0361080' /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 | grep 'gene'
#scaffold_6      MYZPE13164_O_EIv2.1     gene    20285247        20383537        1.05e+03        +       .       ID=MYZPE13164_O_EIv2.1_0361080;Name=MYZPE13164_O_EIv2.1_0361080;biotype=protein_coding_gene;confidence=High
#20383537-20285247=98290
grep -A 1 'MYZPE13164_O_EIv2.1_0361080' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.fa | wc -c #98,321
grep -A 1 'MYZPE13164_O_EIv2.1_0361080' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt2 | wc -c #98,351

#To extract CDS only:
source package /tgac/software/testing/bin/gffread-0.11.4
gffread -x /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.CDS.fa -g /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
```
### Seperate SNP data for each gene in the M. persicae genome
Create a vcf file for each gene:
#### The entire gene region:
```bash
for vcf in $(ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/193s.M_persicae.onlySNPs-genic-regions.vcf.gz); do
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene
    OutFile=NA
    GffFile=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_extract_gene_snps.sh $InFile $OutDir $OutFile $GffFile
done #54396579

cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene
for file in $(ls *.vcf); do 
    singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 /hpc-home/did23faz/git_repos/Scripts/NBI/remove_duplicate_snps.py $file 
    rm $file
done

for file in $(ls *.gz);do
    nongz=$(echo $file | sed 's/dedup_//g')
    if [ ! -f "$nongz" ]; then
    echo "File does not exist: $nongz"
    fi
done

for file in *.gz; do bgzip -d "$file"; done
for x in $(ls *.vcf); do
snpno=$(cat $x | wc -l) 
echo $x $snpno >> gene_snp_report.txt
done
sed -i 's/ /,/g' gene_snp_report.txt
for file in *; do bgzip "$file"; done
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae
```
Plot cumulative frequency of SNPs in genes:
```R
df <- read.table(file = "//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_snp_report.txt", sep = ',', header = FALSE)

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
#### The CDS gene region
```bash
for vcf in $(ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/193s.M_persicae.onlySNPs-genic-regions.vcf.gz); do
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS
    OutFile=NA
    GffFile=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_extract_CDS_snps.sh $InFile $OutDir $OutFile $GffFile
done #55768658
```
#### Extract gene information
Extract gene predictions for the 6 largest chromosomes, this is 35,552 of 37,720:
```bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3

chroms = ['scaffold_1', 'scaffold_2', 'scaffold_3', 'scaffold_4', 'scaffold_5', 'scaffold_6']

with open('/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3') as inp, open('/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation_genes_scaff1-6.gff3', 'w') as outf:
 for line in inp:
     A = line.strip().split('\t')
     if A[0] in chroms and A[2] == 'gene':
        outf.write(line.strip() + '\n')

print ('Done')
exit()
```
Extract information for each gene:
```python
import re

gff_file_path = '/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3'
output_file_path = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_info.txt'

def extract_gene_info_from_gff(gff_file, output_file):
    gene_info = []
    with open(gff_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):  # Ignore comment lines
                fields = line.split('\t')
                feature_type = fields[2]
                if feature_type == 'gene':  # Process only gene features
                    chromosome = fields[0]
                    attributes = fields[8]
                    gene_name = re.search(r'Name=(.*?)(;|$)', attributes).group(1)
                    start_pos = int(fields[3])
                    stop_pos = int(fields[4])
                    gene_info.append((chromosome, gene_name, start_pos, stop_pos))
    # Write gene information to output file
    with open(output_file, 'w') as output:
        for chromosome, gene_name, start_pos, stop_pos in gene_info:
            output.write(f"Chromosome: {chromosome}\tGene Name: {gene_name}\tStart Position: {start_pos}\tStop Position: {stop_pos}\n")
    print("Extraction complete. Gene information written to", output_file)

extract_gene_info_from_gff(gff_file_path, output_file_path)

gff_file_path = '/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation_genes_scaff1-6.gff3'
output_file_path = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_info1-6.txt'

def extract_gene_info_from_gff(gff_file, output_file):
    gene_info = []
    with open(gff_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):  # Ignore comment lines
                fields = line.split('\t')
                feature_type = fields[2]
                if feature_type == 'gene':  # Process only gene features
                    chromosome = fields[0]
                    attributes = fields[8]
                    gene_name = re.search(r'Name=(.*?)(;|$)', attributes).group(1)
                    start_pos = int(fields[3])
                    stop_pos = int(fields[4])
                    gene_info.append((chromosome, gene_name, start_pos, stop_pos))
    # Write gene information to output file
    with open(output_file, 'w') as output:
        for chromosome, gene_name, start_pos, stop_pos in gene_info:
            output.write(f"Chromosome: {chromosome}\tGene Name: {gene_name}\tStart Position: {start_pos}\tStop Position: {stop_pos}\n")
    print("Extraction complete. Gene information written to", output_file)

extract_gene_info_from_gff(gff_file_path, output_file_path)
```
Extract information for genes with SNPs in them:
```bash
for vcf in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/ -name "*_snps.vcf" -exec readlink -f {} \;); do
GeneName=$(echo $vcf | cut -d '/' -f15 | cut -d '_' -f2,3,4,5)
#echo $GeneName
grep $GeneName /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/gene_info.txt >> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_gene_info.txt
done
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene | wc -l #23,459
wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_gene_info.txt #23,457
rm /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_gene_info.txt
```
#### Create mutant gene multifastas
For homozygous mutant mutations:
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/homo_gene_fastas
for vcf in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/ -name "*_snps.vcf" -exec readlink -f {} \;); do
    Jobs=$(squeue -u did23faz | wc -l)
    while [ $Jobs -gt 99 ]; do
    sleep 15s
    printf "."
    Jobs=$(squeue -u did23faz | wc -l)
    done 
    gene_info_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_info.txt
    reference_fasta=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.fa
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/homo_gene_fastas
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
#    GeneName=$(echo $vcf | cut -d '/' -f15 | cut -d '_' -f2,3,4,5 | grep -e MYZPE13164_O_EIv2.1_0035290)
    GeneName=$(echo $vcf | cut -d '/' -f15 | cut -d '_' -f2,3,4,5)
    if [ -n "$GeneName" ]; then
     echo $GeneName
     echo "$GeneName" >> logs/run_create_sample_sequence_log.txt
     vcf_file=$vcf
     OutFile=${GeneName}.fa
     echo $vcf_file >> logs/run_create_sample_sequence_log.txt
     echo $OutFile >> logs/run_create_sample_sequence_log.txt
     sbatch $ProgDir/run_create_sample_sequence_files.sh $vcf_file $OutDir $OutFile $reference_fasta $gene_info_file 2>&1 >> logs/run_create_sample_sequence_log.txt
    fi
done
echo done
#NOTE: The above script will replace missing SNP data with Ns, the actual nucleotide could be reference or alternative, probably best not to use these positions for plotting trees etc., or exclude those samples with missing data from these.

for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/homo_gene_fastas -name "*hom*.fa" -exec readlink -f {} \;); do
    echo $file
    grep '>' $file | wc -l
done
awk '{print length}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/homo_gene_fastas2/hom_MYZPE13164_O_EIv2.1_0010900.fa #All samples are same length as expected
#ERROR: multifastas are all different lengths and do not contain all of the clones 54749041 0010900.fa -> fixed 12/07/2023?
```
For heterozygous mutations:
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hetero_gene_fastas
for vcf in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/ -name "*_snps.vcf" -exec readlink -f {} \;); do
    Jobs=$(squeue -u did23faz | wc -l)
    while [ $Jobs -gt 99 ]; do
    sleep 15s
    printf "."
    Jobs=$(squeue -u did23faz | wc -l)
    done 
    gene_info_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_info.txt
    reference_fasta=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.fa
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hetero_gene_fastas
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
#    GeneName=$(echo $vcf | cut -d '/' -f15 | cut -d '_' -f2,3,4,5 | grep -e MYZPE13164_O_EIv2.1_0035290)
    GeneName=$(echo $vcf | cut -d '/' -f15 | cut -d '_' -f2,3,4,5)
    if [ -n "$GeneName" ]; then
     echo "$GeneName" >> logs/run_create_sample_sequence_log.txt
     vcf_file=$vcf
     OutFile=${GeneName}.fa
     echo $vcf_file >> logs/run_create_sample_sequence_log.txt
     echo $OutFile >> logs/run_create_sample_sequence_log.txt
     sbatch $ProgDir/run_create_sample_sequence_files_het.sh $vcf_file $OutDir $OutFile $reference_fasta $gene_info_file 2>&1 >> logs/run_create_sample_sequence_log.txt
    fi
done
echo done
#NOTE: The above script will replace missing SNP data with Ns, the actual nucleotide could be reference or alternative, probably best not to use these positions for plotting trees etc., or exclude those samples with missing data from these.
#Reverse engineering the sequences in this way means that they are effectively already aligned and trimmed, so there is no need to use MAFFT on the multi fastas follwed by Trim-Al.
```
#### Create mutant CDS multifastas
For homozygous mutant mutations:
```bash
for gene_multifasta in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/homo_gene_fastas -name "hom_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \;); do  
Count=1
for ((i=Count; i<2; i+=1)); do
genename=$(basename $gene_multifasta | cut -d '_' -f2,3,4,5 | cut -d '.' -f1,2).${i}
gff=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
if grep -q "$genename" "$gff"; then
OutDir=$(echo $gene_multifasta | cut -d '/' -f1,2,3,4,5,6,7,8,9,10,11,12,13,14)/homo_CDS_fastas
OutFile=hom_$(echo $genename)_CDS.fa
echo ${OutDir}/${OutFile} >> logs/hom_splice_CDS.txt
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
Jobs=$(squeue -u did23faz| grep 'create_m'  | wc -l)
echo x
while [ $Jobs -gt 19 ]; do
sleep 60s
printf "."
Jobs=$(squeue -u did23faz| grep 'create_m'  | wc -l)
done
mkdir $OutDir
sbatch $ProgDir/run_splice_CDS.sh $gene_multifasta $gff $genename $OutDir $OutFile 2>&1 >> logs/hom_splice_CDS.txt
fi
done
done 
```
For heterozygous mutations:
```bash
for gene_multifasta in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hetero_gene_fastas -name "het_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \;); do  
Count=1
for ((i=Count; i<2; i+=1)); do
genename=$(basename $gene_multifasta | cut -d '_' -f2,3,4,5 | cut -d '.' -f1,2).${i}
gff=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
if grep -q "$genename" "$gff"; then
OutDir=$(echo $gene_multifasta | cut -d '/' -f1,2,3,4,5,6,7,8,9,10,11,12,13,14)/hetero_CDS_fastas
OutFile=het_$(echo $genename)_CDS.fa
echo ${OutDir}/${OutFile} >> logs/het_splice_CDS.txt
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
Jobs=$(squeue -u did23faz| grep 'create_m'  | wc -l)
echo x
while [ $Jobs -gt 19 ]; do
sleep 60s
printf "."
Jobs=$(squeue -u did23faz| grep 'create_m'  | wc -l)
done
mkdir $OutDir
sbatch $ProgDir/run_splice_CDS.sh $gene_multifasta $gff $genename $OutDir $OutFile 2>&1 >> logs/het_splice_CDS.txt
fi
done
done
```
#### Create mutant genomes
For homozygous mutant mutations:
```bash
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz); do
reference_fasta=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/mutant_genomes
OutFile=NA
mkdir $OutDir
sbatch $ProgDir/run_create_samples.sh $file $OutDir $OutFile $reference_fasta
done #55516570

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz); do
reference_fasta=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/mutant_genomes
OutFile=NA
mkdir $OutDir
sbatch $ProgDir/run_create_samples.sh $file $OutDir $OutFile $reference_fasta
done #55510639

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz); do
reference_fasta=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/mutant_genomes
OutFile=NA
mkdir $OutDir
sbatch $ProgDir/run_create_samples_lowmem.sh $file $OutDir $OutFile $reference_fasta
done #55516661

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz); do
reference_fasta=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/mutant_genomes
OutFile=NA
mkdir $OutDir
sbatch $ProgDir/run_create_samples_lowmem2.sh $file $OutDir $OutFile $reference_fasta
done #55588276

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz); do
reference_fasta=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/mutant_genomes
OutFile=NA
mkdir $OutDir
sbatch $ProgDir/run_create_samples_lowmem3.sh $file $OutDir $OutFile $reference_fasta
done #55611138 - out of memory with 8 cpus, 55612463
```


















































```bash
source package fa33234e-dceb-4a58-9a78-7bcf9809edd7
bwa index /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.fa
bwa mem /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.fa temp/MYZPE13164_O_EIv2.1_0035290.fa > temp/alignments.sam
```
```bash
source package /tgac/software/production/bin/tabix-0.2.6
tabix -p vcf /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/193s.M_persicae.onlySNPs-genic-regions.vcf.gz
```
```bash
#RAxML - problems
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/RAxML
for Alignment in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_fastas/ -name "homo_MYZPE13164_O_EIv2.1_0377200.fa" -exec readlink -f {} \;); do
#Jobs=$(squeue -u did23faz | wc -l)
#echo x
#while [ $Jobs -gt 2000000 ]; do
#sleep 15s
#printf "."
#Jobs=$(squeue -u did23faz | wc -l)
#done   
printf "\n" >> logs/raxmllog.txt
Prefix=$(basename $Alignment | cut -f1,2 -d '.')
echo $Prefix >> logs/raxmllog.txt
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/RAxML/$Prefix
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_RAxML.sh $Alignment $OutDir $Prefix 2>&1 >> logs/raxmllog.txt
done

#IQtree
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/IQtree
for Alignment in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_fastas/ -name "homo_MYZPE13164_O_EIv2.1_0377200.fa" -exec readlink -f {} \;); do
#Jobs=$(squeue -u did23faz | wc -l)
#echo x
#while [ $Jobs -gt 2000000 ]; do
#sleep 15s
#printf "."
#Jobs=$(squeue -u did23faz | wc -l)
#done   
printf "\n" >> logs/IQtreelog.txt
Prefix=$(basename $Alignment | cut -f1,2 -d '.')
echo $Prefix >> logs/IQtreelog.txt
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/IQtree/$Prefix
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_RAxML.sh $Alignment $OutDir $Prefix 2>&1 >> logs/IQtreelog.txt
done

source package /nbi/software/production/bin/openmpi-1.6.3

source package /tgac/software/testing/bin/fasttree-2.1.11
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
csv_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_genic.csv'
output_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_genic.dist'

distance_matrix = read_distance_matrix_from_csv(csv_file)
convert_to_phylip(distance_matrix, output_file)
```
```bash
sed "s/$(cat JIC.txt)/J/g" /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_genic.dist > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_genic_sampler.dist

```
```
source package 3e7beb4d-f08b-4d6b-9b6a-f99cc91a38f9
source package 7654f72b-1692-46bb-9a56-443406d03fd9
source package d37013e7-5691-40b6-8885-f029fe5fad54
source package /tgac/software/testing/bin/pgdspider-2.1.1.5
SplitsTree
```