# Contructing phylogenies as part of the M. persicae population genomics project

## Extract gene sequences
```bash
#Convert fasta to single line format:
awk ' {if (NR==1) {print $0} else {if ($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa > /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fasta

#Separate different scaffolds:
for scaffold in $(grep '>' /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fasta); do
name=$(echo $scaffold | sed 's@>@@g')
grep -w -A 1 "$scaffold" /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fasta > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.fa
grep -w "##gff-version 3\|$name" /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.gff
done 

#Generate gene sequences from scaffolds
for scaffold in $(ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/*.fa); do
name=$(echo $scaffold | cut -d '/' -f11 | sed 's@.fa@@g')
echo $name
#To extract whole gene sequences (NOTE: some are very long):
source package /tgac/software/testing/bin/bedops-2.2.0
source package 4028d6e4-21a8-45ec-8545-90e4ed7e1a64
gff2bed < /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.gff > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.bed
grep "gene" /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.bed > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}_2.bed
bedtools getfasta -fi $scaffold -bed /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}_2.bed -fo /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.gff3.nt2 -name+
#To extract CDS only:
source package /tgac/software/testing/bin/gffread-0.11.4
gffread -x /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.gff3.nt3 -g $scaffold /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/${name}.gff
done 

cat /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/scaffold_*.gff3.nt2 > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.gene.fa
cat /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds/scaffold_*.gff3.nt3 > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.CDS.fa
rm -r /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/scaffolds
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
#I have investigated a few individual genes and confirmed that the genes missing from /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene have no SNPs in them

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

######################################################################################################################################################
#Second appraoch
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene2
grep 'gene' /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 > gene_lines.gff

while IFS= read -r line; do
    genename=$(echo $line | cut -d '=' -f2 | cut -d ';' -f1)
    grep '##gff-version\|##sequence-region' /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 > ${genename}.gff3
    sleep 3s
    grep "$genename" /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 > temp_gene_line.gff3 >> ${genename}.gff3
    sleep 2s
    genename=$(echo $line | cut -d '=' -f2 | cut -d ';' -f1)
    echo $genename >> logs/bedtools_intersect.txt
    vcf=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz
    OutFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene2/${genename}_snps.vcf
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    Jobs=$(squeue -u did23faz | grep 'bedtools' | wc -l)
    while [ $Jobs -gt 99 ]; do
    sleep 15s
    printf "."
    Jobs=$(squeue -u did23faz | grep 'bedtools' | wc -l)
    done 
    sbatch $ProgDir/bedtools_intersect.sh $vcf ${genename}.gff3 $OutFile 2>&1 >> logs/bedtools_intersect.txt
done < gene_lines.gff
0000030
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
bgzip -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/dedup_MYZPE13164_O_EIv2.1_0213490_snps.vcf > dedup_MYZPE13164_O_EIv2.1_0213490_snps.vcf.gz

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS
grep 'gene' /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 | wc -l #37,720 this is the total number of genes
find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/ -name "*_snps.vcf" -exec readlink -f {} \; >> temp_conter.txt
wc -l temp_conter.txt #23,456 this is the number of genes with SNPs in them.

for vcf in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/ -name "*_snps.vcf" -exec readlink -f {} \;); do
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS
    OutFile=NA
    GffFile=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    Jobs=$(squeue -u did23faz | grep 'extract' | wc -l)
    while [ $Jobs -gt 99 ]; do
    sleep 15s
    printf "."
    Jobs=$(squeue -u did23faz | grep 'extract' | wc -l)
    done 
    echo $InFile >> logs/cds_vcfs.txt
    sbatch $ProgDir/run_extract_CDS_snps.sh $InFile $OutDir $OutFile $GffFile >> logs/cds_vcfs.txt
done 

ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS | wc -l #13,991 there are many missing

for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/ -name "*_snps.vcf.gz" -exec readlink -f {} \;); do
rm $file
done

find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas -name "het*.fa" -exec readlink -f {} \; | wc -l

######################################################################################################################
#Use the gene script but a .vcf with only CDS entries:
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS2
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-CDS_genic-regions.vcf.gz); do
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS2
    OutFile=NA
    GffFile=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    sbatch $ProgDir/run_extract_gene_snps.sh $InFile $OutDir $OutFile $GffFile #56425989
done 

ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS2 | wc -l #17,794 
find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/ -name "*dedup_MYZPE13164_O_EIv2.1_*_snps.vcf" -exec readlink -f {} \; | wc -l #23,456 - 5,662 have no CDS SNPs
#I have investigated a few individual genes and confirmed that the genes missing from /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS2 have no SNPs in their CDS regions
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
Extract information for each gene, for all and the largest 6 scaffolds:
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
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas

#As an array job:
vcf_dir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/
gene_info_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_info.txt
reference_fasta=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.gene.fa
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
log=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/logs/run_create_sample_sequence_log.txt
sbatch $ProgDir/run_create_sample_sequence_files_array100.sh $vcf_dir $gene_info_file $reference_fasta $OutDir $log
#MAX ARRAY SIXE PER USER IS 6000, THEREFORE CAN'T USE ARRAYS FOR MANY JOBS MAKING THEM POINTLESS...


for vcf in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/ -name "*MYZPE13164_O_EIv2.1_*_snps.vcf" -exec readlink -f {} \;); do
#for vcf in $(ls *MYZPE13164_O_EIv2.1_0317870_snps.vcf); do
Jobs=$(squeue -u did23faz | wc -l)
while [ $Jobs -gt 99 ]; do
sleep 15s
printf "."
Jobs=$(squeue -u did23faz | wc -l)
done 
gene_info_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_info.txt
reference_fasta=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.gene.fa
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
GeneName=$(echo $vcf | cut -d '/' -f15 | cut -d '_' -f2,3,4,5)
if [ -n "$GeneName" ]; then
vcf_file=$vcf
OutFile=hom_${GeneName}.fa
if [ ! -e ${OutDir}/${OutFile} ]; then
echo $GeneName
echo "$GeneName" >> logs/run_create_sample_sequence_log_log_28072023.txt
echo $vcf_file >> logs/run_create_sample_sequence_log_log_28072023.txt
echo $OutFile >> logs/run_create_sample_sequence_log_log_28072023.txt
#sbatch $ProgDir/run_create_sample_sequence_files.sh $vcf_file $OutDir $OutFile $reference_fasta $gene_info_file 2>&1 >> logs/run_create_sample_sequence_log.txt 
sbatch $ProgDir/run_create_sample_sequence_files.sh $vcf_file $OutDir $OutFile $reference_fasta $gene_info_file 2>&1 >> logs/run_create_sample_sequence_log_28072023.txt 
else
echo $GeneName already run
fi
fi
done
echo done
#NOTE: The above script will replace missing SNP data with Ns, the actual nucleotide could be reference or alternative, probably best not to use these positions for plotting trees etc., or exclude those samples with missing data from these.
#NOTE: Slurm files are output to logs/create_sample_sequence_files/. the script prints and error message if the reference base is not what is expected.
#NOTE: The script prefixes output files with the slurm job number

#Remove job no. prefixes:
for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas/ -name "*_hom_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \;); do
newname=$(basename $file | rev | cut -d '_' -f1,2,3,4,5 | rev)
#echo $file
#echo /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas/${newname}
mv $file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas/${newname}
done

#Check that none of the jobs failed:
for job in $(grep 'Submitted' logs/run_create_sample_sequence_log_28072023.txt | cut -d ' ' -f4); do
sacct -j $job --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode | grep -v 'COMPLETED\|----------\|State'
done #No jobs failed

#Check that the multifasta have all samples, and all genes are the expected lengths: 
#NOTE: insertions and deletions will not be captured as these are filtered out of the vcf - other than the passing SNPs the variant genes will be shown as the same as the reference.
for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas -name "hom*.fa" -exec readlink -f {} \;); do
    echo $file
    grep '>' $file | wc -l
done #All files have the expected length 194

awk '{print length}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas/hom_MYZPE13164_O_EIv2.1_0010900.fa #All samples are same length as expected
awk '{print length}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas/hom_MYZPE13164_O_EIv2.1_0110900.fa #All samples are same length as expected
awk '{print length}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas/hom_MYZPE13164_O_EIv2.1_0210900.fa #All samples are same length as expected
awk '{print length}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas/hom_MYZPE13164_O_EIv2.1_0310900.fa #All samples are same length as expected
```
For heterozygous mutations:
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas
mkdir logs/create_sample_sequence_files_het
for vcf in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/ -name "*_snps.vcf" -exec readlink -f {} \;); do
Jobs=$(squeue -u did23faz | wc -l)
while [ $Jobs -gt 99 ]; do
sleep 15s
printf "."
Jobs=$(squeue -u did23faz | wc -l)
done 
gene_info_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_info.txt
reference_fasta=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.gene.fa
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
GeneName=$(echo $vcf | cut -d '/' -f15 | cut -d '_' -f2,3,4,5)
if [ -n "$GeneName" ]; then
vcf_file=$vcf
OutFile=het_${GeneName}.fa
if [ ! -e ${OutDir}/${OutFile} ]; then
echo $GeneName
echo "$GeneName" >> logs/run_create_sample_sequence_het_log_01082023.txt
echo $vcf_file >> logs/run_create_sample_sequence_het_log_01082023.txt
echo $OutFile >> logs/run_create_sample_sequence_het_log_01082023.txt
sbatch $ProgDir/run_create_sample_sequence_files_het.sh $vcf_file $OutDir $OutFile $reference_fasta $gene_info_file 2>&1 >> logs/run_create_sample_sequence_het_log_01082023.txt 
else
echo $GeneName already run
fi
fi
done
echo done
#NOTE: The above script will replace missing SNP data with Ns, the actual nucleotide could be reference or alternative, probably best not to use these positions for plotting trees etc., or exclude those samples with missing data from these.
#NOTE: Slurm files are output to logs/create_sample_sequence_files_het/. the script prints and error message if the reference base is not what is expected.
#NOTE: The script prefixes output files with the slurm job number

#Remove job no. prefixes:
for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas/ -name "*_het_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \;); do
newname=$(basename $file | rev | cut -d '_' -f1,2,3,4,5 | rev)
mv $file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas/${newname}
done

#Reverse engineering the sequences in this way means that they are effectively already aligned and trimmed, so there is no need to use MAFFT on the multi fastas follwed by Trim-Al.

#Check that none of the jobs failed:
for job in $(grep 'Submitted' logs/run_create_sample_sequence_het_log_01082023.txt | cut -d ' ' -f4); do
sacct -j $job --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode | grep -v 'COMPLETED\|----------\|State'
done

#Check that the multifasta have all samples, and all genes are the expected lengths: 
#NOTE: insertions and deletions will not be captured as these are filtered out of the vcf - other than the passing SNPs the variant genes will be shown as the same as the reference.
for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas -name "het*.fa" -exec readlink -f {} \;); do
    echo $file
    grep '>' $file | wc -l
done #All files have the expected length 194

awk '{print length}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas/het_MYZPE13164_O_EIv2.1_0010900.fa #All samples are same length as expected
awk '{print length}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas/het_MYZPE13164_O_EIv2.1_0110900.fa #All samples are same length as expected
awk '{print length}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas/het_MYZPE13164_O_EIv2.1_0210900.fa #All samples are same length as expected
awk '{print length}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas/het_MYZPE13164_O_EIv2.1_0310900.fa #All samples are same length as expected
```

#### Create mutant CDS multifastas
For homozygous mutant mutations:
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas
for gene_multifasta in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas -name "hom_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \;); do  
Jobs=$(squeue -u did23faz| grep 'create_m' | wc -l)
while [ $Jobs -gt 19 ]; do
sleep 60s
printf "."
Jobs=$(squeue -u did23faz| grep 'create_m' | wc -l)
done
gff=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
Count=1
geneID=$(echo $gene_multifasta | rev | cut -d '_' -f1 | rev | cut -d '.' -f1)  
limit=$(grep "$geneID" $gff | grep 'mRNA' | wc -l)
for ((i=Count; i<$((limit + 1)); i+=1)); do
variant=$(basename $gene_multifasta | cut -d '_' -f2,3,4,5 | cut -d '.' -f1,2).${i}
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas
OutFile=hom_$(echo $variant)_CDS.fa 
echo ${OutDir}/${OutFile} >> logs/hom_splice_CDS.txt
echo $OutFile  
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI 
sbatch $ProgDir/run_splice_CDS.sh $gene_multifasta $gff $variant $OutDir $OutFile 2>&1 >> logs/hom_splice_CDS.txt
done
done 

#CHECK HAS WORKED AS INTENDED
for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom*.fa" -exec readlink -f {} \;); do
    echo $file
    grep '>' $file | wc -l
done #All samples present

awk '{print length}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/hom_MYZPE13164_O_EIv2.1_0010900.1_CDS.fa #same lengths and multiple of 3
awk '{print length}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/hom_MYZPE13164_O_EIv2.1_0110900.1_CDS.fa #same lengths and multiple of 3
awk '{print length}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/hom_MYZPE13164_O_EIv2.1_0210900.1_CDS.fa #missing - this has no mRNA and is ncRNA therefore the loop did not work
awk '{print length}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/hom_MYZPE13164_O_EIv2.1_0310900.1_CDS.fa #same lengths and multiple of 3

grep 'mRNA' /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 | wc -l #47,508
for gene_multifasta in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas -name "hom_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \;); do  
gff=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
Count=1
geneID=$(echo $gene_multifasta | rev | cut -d '_' -f1 | rev | cut -d '.' -f1)  
limit=$(grep "$geneID" $gff | grep 'mRNA' | wc -l)
for ((i=Count; i<$((limit + 1)); i+=1)); do
variant=$(basename $gene_multifasta | cut -d '_' -f2,3,4,5 | cut -d '.' -f1,2).${i}
OutFile=hom_$(echo $variant)_CDS.fa 
echo $OutFile >> temp_count_hom.txt
done
done 
wc -l temp_count_hom.txt #36,535
find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "*" -exec readlink -f {} \; >> temp_hom.txt
wc -l temp_hom.txt #54,381 <- this makes no sense there should not be more CDS fasta files than there are mRNA in the .gff, there should be less as not all genes have SNPs
find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom*.fa" -exec readlink -f {} \; >> temp_hom2.txt
find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom*+.fa" -exec readlink -f {} \; >> temp_hom3.txt
find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom*-.fa" -exec readlink -f {} \; >> temp_hom3.txt
wc -l temp_hom3.txt #36,434
find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas -name "hom_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \; >> temp_hom4.txt
wc -l temp_hom4.txt #23,455
for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas -name "het_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \;); do
name=$(basename $file | sed 's@het_@@g')
if ! grep -q "$name" temp_hom3.txt ; then
echo $name >> temp_missing_hom.txt
fi
done



for gene_multifasta in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas -name "hom_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \;); do  
name=$(basename $gene_multifasta | cut -d '.' -f1,2)
if grep -q "$name" temp.txt ; then
sleep 0
else
echo $gene_multifasta >> temp_erro.txt
fi
done

wc -l temp_erro.txt #4,406

for gene in $(cat temp_erro.txt); do
name=$(basename $gene | cut -d '.' -f1,2 | cut -d '_' -f2,3,4,5)
if ! grep -q "$name" /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 || ! grep -q "ncRNA" /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 ; then
echo $gene is missng and not a ncRNA
fi
done
#This confirms that all missing genes are ncRNAs
```
Correct for strandedness:
```bash
for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom*S.fa" -exec readlink -f {} \;); do
gff=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
name=$(basename $file | cut -d '_' -f2,3,4,5 | cut -d '.' -f1,2)
echo $name
strand=$(grep "$name" $gff | grep 'gene' | cut -f7)
echo $strand
if [ "$strand" = "-" ]; then
OutFile=$(dirname $file)/$(basename $file| sed 's@CDS.fa@CDS-.fa@g')
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 /hpc-home/did23faz/git_repos/Scripts/NBI/reverse_compliment.py $file $OutFile
if [ ! -e $OutFile ]; then
rm $file
fi
else
OutFile=$(dirname $file)/$(basename $file| sed 's@CDS.fa@CDS+.fa@g')
mv $file $OutFile
fi
done

for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom*S.fa" -exec readlink -f {} \;); do
rm $file
done

for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom_MYZPE13164_O_EIv2.1_*_CDS-.fa" -exec readlink -f {} \;); do
awk ' {if (NR==1) {print $0} else {if ($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $file > temp.fa
cat temp.fa > $file
done

for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom_MYZPE13164_O_EIv2.1_*_CDS*.fa" -exec readlink -f {} \;); do
cp $file /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/hom_CDS_fastas/.
done

find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom*.fa" -exec readlink -f {} \; | wc -l #36,480
```
For heterozygous mutations:
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas
for gene_multifasta in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas -name "het_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \;); do  
Jobs=$(squeue -u did23faz| grep 'create_m' | wc -l)
while [ $Jobs -gt 19 ]; do
sleep 60s
printf "."
Jobs=$(squeue -u did23faz| grep 'create_m' | wc -l)
done
gff=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
Count=1
geneID=$(echo $gene_multifasta | rev | cut -d '_' -f1 | rev | cut -d '.' -f1)  
limit=$(grep "$geneID" $gff | grep 'mRNA' | wc -l)
for ((i=Count; i<$((limit + 1)); i+=1)); do
variant=$(basename $gene_multifasta | cut -d '_' -f2,3,4,5 | cut -d '.' -f1,2).${i}
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas
OutFile=het_$(echo $variant)_CDS.fa 
echo ${OutDir}/${OutFile} >> logs/het_splice_CDS.txt
echo $OutFile  
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI 
sbatch $ProgDir/run_splice_CDS.sh $gene_multifasta $gff $variant $OutDir $OutFile 2>&1 >> logs/het_splice_CDS.txt
done
done 

for gene_multifasta in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas -name "het_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \;); do  
gff=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
Count=1
geneID=$(echo $gene_multifasta | rev | cut -d '_' -f1 | rev | cut -d '.' -f1)  
limit=$(grep "$geneID" $gff | grep 'mRNA' | wc -l)
for ((i=Count; i<$((limit + 1)); i+=1)); do
variant=$(basename $gene_multifasta | cut -d '_' -f2,3,4,5 | cut -d '.' -f1,2).${i}
OutFile=het_$(echo $variant)_CDS.fa 
echo $OutFile >> temp_count_het.txt
done
done 
find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas -name "*" -exec readlink -f {} \; >> temp_het.txt
wc -l temp_het.txt #36,536

find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas -name "*" -exec readlink -f {} \; >> temp2.txt
for gene_multifasta in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas -name "het_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \;); do  
name=$(basename $gene_multifasta | cut -d '.' -f1)
if grep -q "$name" temp2.txt ; then
else
echo $gene_multifasta > temp_erro.txt
done

56399073
```
Correct for strandedness:
```bash
for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas -name "het_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \;); do
gff=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3
name=$(basename $file | cut -d '_' -f2,3,4,5 | cut -d '.' -f1,2)
echo $name
strand=$(grep "$name" $gff | grep 'gene' | cut -f7)
echo $strand
if [ "$strand" = "-" ]; then
OutFile=$(dirname $file)/$(basename $file| sed 's@CDS.fa@CDS-.fa@g')
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 /hpc-home/did23faz/git_repos/Scripts/NBI/reverse_compliment.py $file $OutFile
rm $file
if [ ! -e $OutFile ]; then
rm $file
fi
else
OutFile=$(dirname $file)/$(basename $file| sed 's@CDS.fa@CDS+.fa@g')
mv $file $OutFile
fi
done

for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas -name "het_MYZPE13164_O_EIv2.1_*_CDS-.fa" -exec readlink -f {} \;); do
awk ' {if (NR==1) {print $0} else {if ($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $file > temp.fa
cat temp.fa > $file
done

for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas -name "het_MYZPE13164_O_EIv2.1_*_CDS*.fa" -exec readlink -f {} \;); do
cp $file /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/het_CDS_fastas/.
done

find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas -name "het_MYZPE13164_O_EIv2.1_*_CDS*.fa" -exec readlink -f {} \; | wc -l #37,769
```
#### RAxML
```bash
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas/het_MYZPE13164_O_EIv2.1_0363660.1_CDS+.fa temp77.fa
sed -i 's@MYZPE13164_O_EIv2.1_@_@g' temp77.fa
sed -i 's@_CDS@@g' temp77.fa
sed -i 's@_@@g' temp77.fa

for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas -name "het_MYZPE13164_O_EIv2.1_*_CDS*.fa" -exec readlink -f {} \;); do
sed -i '/^$/d' $file
multifasta=$file
OutFile=$(basename $multifasta | sed 's@.fa@@g')
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas/RAxML/$OutFile
ProgDir=~/git_repos/Wrappers/NBI
Jobs=$(squeue -u did23faz| grep 'RAxML' | wc -l)
while [ $Jobs -gt 99 ]; do
sleep 60s
printf "."
Jobs=$(squeue -u did23faz| grep 'RAxML' | wc -l)
done
if [ ! -e ${OutDir}/${OutFile}.raxml.bestTree ]; then
mkdir $OutDir
echo $OutFile  
echo $OutFile >> logs/raxml_het.txt
sbatch $ProgDir/run_RAxML_msa.sh $multifasta $OutDir $OutFile 2>&1 >> logs/raxml_het.txt
else
echo Done for $OutFile
fi
done

#Check that none of the jobs failed:
for job in $(grep 'Submitted' logs/raxml_het.txt | cut -d ' ' -f4); do
sacct -j $job --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode | grep -v 'COMPLETED\|----------\|State'
done

multifasta=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas/het_MYZPE13164_O_EIv2.1_0363660.1_CDS+.fa

#######################################################################################################################

for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom_MYZPE13164_O_EIv2.1_*_CDS*.fa" -exec readlink -f {} \;); do
sed -i '/^$/d' $file
multifasta=$file
OutFile=$(basename $multifasta | sed 's@.fa@@g')
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/RAxML/$OutFile
ProgDir=~/git_repos/Wrappers/NBI
Jobs=$(squeue -u did23faz| grep 'RAxML' | wc -l)
while [ $Jobs -gt 99 ]; do
sleep 60s
printf "."
Jobs=$(squeue -u did23faz| grep 'RAxML' | wc -l)
done
if [ ! -e ${OutDir}/${OutFile}.raxml.bestTree ]; then
mkdir $OutDir
echo $OutFile  
echo $OutFile >> logs/raxml_hom.txt
sbatch $ProgDir/run_RAxML_msa.sh $multifasta $OutDir $OutFile 2>&1 >> logs/raxml_hom.txt
else
echo Done for $OutFile
fi
done

#Check that none of the jobs failed:
for job in $(grep 'Submitted' logs/raxml_hom.txt | cut -d ' ' -f4); do
sacct -j $job --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode | grep -v 'COMPLETED\|----------\|State'
done

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/RAxML/ | grep -v 'CDS+\|CDS-'); do
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/RAxML/$file
done
```
Missing:
```bash
for Seqfile in $(tac temp_files.txt); do
    TreeFile=$(dirname $Seqfile)/RAxML/$(basename $Seqfile | sed 's@.fa@@g')/$(basename $Seqfile | sed 's@.fa@@g').raxml.bestTree
    if [ ! -e "${TreeFile}" ] || [ ! -s "${TreeFile}" ]; then
        echo $Seqfile >> temp_missing.txt
    fi
done

for x in $(cat temp_missing.txt); do
    y=$(echo $x | rev | cut -d '/' -f1 | rev | sed 's@.fa@@g')
    z=$(grep -A 1 "$y" logs/raxml_hom.txt | awk 'NR==2' | cut -d ' ' -f4)
    sacct -j $z --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode
done
```
#### PAML - Omega DN/DS
```bash
find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom_MYZPE13164_O_EIv2.1_*_CDS*.fa" -exec readlink -f {} \; > temp_files.txt
wc -l temp_files.txt #36480


function is_valid_time {
    current_hour=$(date +"%H")
    day_of_week=$(date +"%u")  # 1 (Monday) through 7 (Sunday)
    if [ "$day_of_week" -ge 1 ] && [ "$day_of_week" -le 5 ]; then
        if [ "$current_hour" -ge 18 ] || [ "$current_hour" -lt 4 ]; then
            return 0  # Valid time on weekdays
        else
            return 1  # Invalid time on weekdays
        fi
    else
        return 0  # Valid time on weekends
    fi
}

#for Seqfile in $(cat temp_csep_files.txt); do

for Seqfile in $(cat temp_files.txt); do
Jobs=$(squeue -u did23faz| grep 'paml'  | wc -l)
echo $Jobs 1
TreeFile=$(dirname $Seqfile)/RAxML/$(basename $Seqfile | sed 's@.fa@@g')/$(basename $Seqfile | sed 's@.fa@@g').raxml.bestTree
OutDir=$(dirname $TreeFile)/paml
OutFile=$(basename $Seqfile | sed 's@_CDS-.fa@@' | sed 's@_CDS+.fa@@').out
ProgDir=~/git_repos/Wrappers/NBI
if is_valid_time; then
    if [ ! -e "${OutDir}/${OutFile}" ] || [ ! -s "${OutDir}/${OutFile}" ]; then
        while [ $Jobs -gt 189 ]; do
            sleep 300s
            printf "."
            Jobs=$(squeue -u did23faz| grep 'paml'| wc -l)
        done
        ls $TreeFile
        mkdir $OutDir
        ls $TreeFile 2>&1 >> logs/pamllog.txt
        sbatch $ProgDir/run_paml_omega.sh $Seqfile $TreeFile $OutDir $OutFile 2>&1 >> logs/pamllog.txt
    else
        echo Already run for ${OutFile}
    fi 
else
    if [ ! -e "${OutDir}/${OutFile}" ] || [ ! -s "${OutDir}/${OutFile}" ]; then
        while [ $Jobs -gt 189 ]; do
            sleep 300s
            printf "."
            Jobs=$(squeue -u did23faz| grep 'paml'| wc -l)
        done
        ls $TreeFile
        mkdir $OutDir
        ls $TreeFile 2>&1 >> logs/pamllog.txt
        sbatch $ProgDir/run_paml_omega.sh $Seqfile $TreeFile $OutDir $OutFile 2>&1 >> logs/pamllog.txt
    else
        echo Already run for ${OutFile}
    fi 
fi
done

/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/hom_MYZPE13164_O_EIv2.1_0002220*
MYZPE13164_O_EIv2.1_0037470
MYZPE13164_O_EIv2.1_0037650

find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom_MYZPE13164_O_EIv2.1_*_CDS*.fa" -exec readlink -f {} \; > temp_files.txt

for gene in $(cat /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/effector_candidates.txt); do
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/hom_${gene}* >> temp_csep_files.txt
done


find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom_MYZPE13164_O_EIv2.1_*_CDS*.fa" -exec readlink -f {} \; > temp_files_all.txt
for Seqfile in $(cat temp_files_all.txt); do
TreeFile=$(dirname $Seqfile)/RAxML/$(basename $Seqfile | sed 's@.fa@@g')/$(basename $Seqfile | sed 's@.fa@@g').raxml.bestTree
OutDir=$(dirname $TreeFile)/paml
OutFile=$(basename $Seqfile | sed 's@_CDS-.fa@@' | sed 's@_CDS+.fa@@').out 
if [ -e "${OutDir}/${OutFile}" ] && [ -s "${OutDir}/${OutFile}" ]; then
    ls ${OutDir}/${OutFile} >> temp_count.txt
fi
done

for file in $(cat temp_count.txt); do
grep 'omega (dN/dS)' $file
done















tail -n 4000 temp_files.txt | tac > temp_temp_files.txt
for Seqfile in $(cat temp_temp_files.txt); do
Jobs=$(squeue -u did23faz| grep 'paml'  | wc -l)
echo $Jobs 1
TreeFile=$(dirname $Seqfile)/RAxML/$(basename $Seqfile | sed 's@.fa@@g')/$(basename $Seqfile | sed 's@.fa@@g').raxml.bestTree
OutDir=$(dirname $TreeFile)/paml
OutFile=$(basename $Seqfile | sed 's@_CDS-.fa@@' | sed 's@_CDS+.fa@@').out
ProgDir=~/git_repos/Wrappers/NBI
if [ ! -e "${OutDir}/${OutFile}" ] || [ ! -s "${OutDir}/${OutFile}" ]; then
while [ $Jobs -gt 189 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -u did23faz| grep 'paml'| wc -l)
done
ls $TreeFile
mkdir $OutDir
ls $TreeFile 2>&1 >> logs/pamllog.txt
sbatch $ProgDir/run_paml_omega.sh $Seqfile $TreeFile $OutDir $OutFile 2>&1 >> logs/pamllog.txt
else
echo Already run for ${OutFile}
fi 
done

head -n 1000 temp_temp_files.txt > temp_array.txt
sbatch $ProgDir/run_paml_omega_array.sh /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_array.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/logs/pamllog.txt

for file in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas -name "hom_MYZPE13164_O_EIv2.1_*_CDS*.fa" -exec readlink -f {} \;); do
TreeFile=$(dirname $Seqfile)/RAxML/$(basename $Seqfile | sed 's@.fa@@g')/$(basename $Seqfile | sed 's@.fa@@g').raxml.bestTree
OutDir=$(dirname $TreeFile)/paml
OutFile=$(basename $file | sed 's@_CDS-.fa@@' | sed 's@_CDS+.fa@@').out

done





for split in $(split -d temp_files.txt); do
Jobs=$(squeue -u did23faz| grep 'paml'  | wc -l)
while [ $Jobs -gt 190 ]; do
    sleep 1800s
    printf "."
    Jobs=$(squeue -u did23faz| grep 'paml'| wc -l)
done
Log=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/logs/pamllog.txt
sbatch $ProgDir/run_paml_omega_array.sh $split $Log
done

for split in $(split -d temp_files.txt); do
Jobs=$(squeue -u did23faz| grep 'paml'  | wc -l)
while [ $Jobs -gt 190 ]; do
    sleep 1800s
    printf "."
    Jobs=$(squeue -u did23faz| grep 'paml'| wc -l)
done
cat $split >> temp_test.txt
echo _ >> temp_test.txt
echo _ >> temp_test.txt
#Log=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/logs/pamllog.txt
#sbatch $ProgDir/run_paml_omega_array.sh $split $Log
done


#!/bin/bash

TempFile="temp_files.txt"
BatchSize=1000
TotalLines=$(wc -l < "$TempFile")
for ((start = 1; start <= TotalLines; start += BatchSize)); do
    end=$((start + BatchSize - 1))
    sed -n "$start,${end}p" "$TempFile" > batch_files.txt
    while read -r split; do
        Jobs=$(squeue -u did23faz | grep 'paml' | wc -l)
        while [ "$Jobs" -gt 189 ]; do
            sleep 1800s
            printf "."
            Jobs=$(squeue -u did23faz | grep 'paml' | wc -l)
        done
        Log="/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/logs/pamllog.txt"
        sbatch "$ProgDir/run_paml_omega_array.sh" "$split" "$Log"
    done < batch_files.txt
    rm batch_files.txt
done

```
Orthofinder dn/ds
```bash
  ProjDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae
  cd $ProjDir
  IsolateAbrv=persicae_v_ligustri
  WorkDir=analysis/orthology/orthofinder/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  

source package /nbi/software/production/bin/porthomcl-40f497e

Taxon_code=MYZ
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=LIG
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_ligustri/v1.1/Myzus_ligustri_v1.1.scaffolds.braker.filtered.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

for WorkDir in $(ls -d analysis/orthology/orthofinder/*); do
Input_dir=$WorkDir/formatted
Min_length=10
Max_percent_stops=20
Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
done

sbatch ~/git_repos/Wrappers/NBI/run_orthofinder.sh $WorkDir 
#57518180

OrthogroupsTxt=$WorkDir/formatted/OrthoFinder/*/Orthogroups/Orthogroups.txt
GoodProts=$WorkDir/goodProteins/goodProteins.fasta
OutDir=$WorkDir/orthogroups_fasta
mkdir -p $OutDir
source package /tgac/software/production/bin/python-2.7.10
python ~/git_repos/Scripts/NBI/orthoMCLgroups2fasta.py --orthogroups $OrthogroupsTxt --fasta $GoodProts --out_dir $OutDir
#57527510
mv slurm.57527510.out $WorkDir/orthoreport.txt

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/orthology/orthofinder/persicae_v_ligustri/orthogroups_fasta/paired
for fasta in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/orthology/orthofinder/persicae_v_ligustri/orthogroups_fasta -name "orthogroupOG*.fa" -exec readlink -f {} \;); do
    if grep -q '^>MYZ' "$fasta" && grep -q '^>LIG' "$fasta"; then
    OutFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/orthology/orthofinder/persicae_v_ligustri/orthogroups_fasta/paired/$(basename $fasta | sed 's@.fa@_paired.fa@g')
    singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/find_longest_myzlig.py $fasta $OutFile
    fi
done
```
```python
fasta_file_path = 'analysis/orthology/orthofinder/persicae_v_ligustri/orthogroups_fasta/orthogroupOG0003376.fa'  # Replace 'your_file.txt' with the actual path to your file

import sys
import os

fasta_file_path = sys.argv[1]
output_path = sys.argv[2]

def process_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    max_length_myz = 0
    max_sequence_myz = ""
    max_length_lig = 0
    max_sequence_lig = ""
    current_header = ""
    for line in lines:
        line = line.strip()
        if line.startswith('>MYZ') or line.startswith('>LIG'):
            current_header = line
        else:
            sequence_length = len(line)
            if current_header.startswith('>MYZ') and sequence_length > max_length_myz:
                max_length_myz = sequence_length
                max_sequence_myz = line
                max_header_myz = current_header
            elif current_header.startswith('>LIG') and sequence_length > max_length_lig:
                max_length_lig = sequence_length
                max_sequence_lig = line
                max_header_lig = current_header
    with open(output_path, 'w') as output_file:
        output_file.write(f"{max_header_myz}")
        output_file.write(f"{max_sequence_myz}")
        output_file.write(f"{max_header_lig}")
        output_file.write(f"{max_sequence_lig}")

process_fasta(fasta_file_path)
```
Orthofinder dn/ds
```bash
  ProjDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae
  cd $ProjDir
  IsolateAbrv=all_aphid
  WorkDir=analysis/orthology/orthofinder/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  

source package /nbi/software/production/bin/porthomcl-40f497e

Taxon_code=ACYpis
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Acyrthosiphon_pisum/JIC1_v1/Acyrthosiphon_pisum_JIC1_v1.0.scaffolds.braker2.gff.filtered.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=APHfab
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Aphis_fabae/JIC1_v2/Aphis_fabae_JIC1_v2.scaffolds.braker.filtered.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=APHgly
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Aphis_glycines/biotype_4_v3/Aphis_glycines_4.v3.scaffolds.braker.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=APHgos
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Aphis_gossypii/JIC1_v1/Aphis_gossypii_JIC1_v1.scaffolds.braker.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=APHrum
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Aphis_rumicis/v1/Aphis_rumicis_v1.scaffolds.braker_gthr.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=APHtha
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Aphis_thalictri/v1/Aphis_thalictri_v1.scaffolds.braker.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=BRAcar
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Brachycaudus_cardui/v1.1/Brachycaudus_cardui_v1.1.scaffolds.braker.filtered.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=BRAhel
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Brachycaudus_helichrysi/v1.1/Brachycaudus_helichrysi_v1.1.scaffolds.braker.filtered.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=BRAklu
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Brachycaudus_klugkisti/v1.1/Brachycaudus_klugkisti_v1.1.scaffolds.braker.filtered.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=BREbra
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Brevicoryne_brassicae/v2/Brevicoryne_brassicae_v2.scaffolds.braker.filtered.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=CINced
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Cinara_cedri/v1/cinced3A.pep.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=DAKvit
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Daktulosphaira_vitifoliae/INRAPcf7_v5/Daktulosphaira_vitifoliae_INRAPcf7_v5.scaffolds.braker.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=DIUnox
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Diuraphis_noxia/SAM_v1.1/Diuraphis_noxia_SAM.v1.1.scaffolds.braker.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

#Taxon_code=DREpla
#Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Drepanosiphum_platanoidis/)
#Id_field=1
#orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
#mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=ERIlan
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Eriosoma_lanigerum/v1/Eriosoma_lanigerum_v1.0.scaffolds.braker.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=HORcor
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Hormaphis_cornu/v1/Augustus.updated_w_annots.21Aug20.gff3.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=MACalb
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Macrosiphum_albifrons/v1/Macrosiphum_albifrons_v1.scaffolds.braker.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=METdir
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Metopolophium_dirhodum/JIC1_v1.1/Metopolophium_dirhodum_UK035_v1.1.scaffolds.braker.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=MYZcer
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_cerasi/Thorpe_v1.2/Myzus_cerasi_v1.2.scaffolds.braker.filtered.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=MYZlig
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_ligustri/v1.1/Myzus_ligustri_v1.1.scaffolds.braker.filtered.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=MYZlyt
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_lythri/v1.1/Myzus_lythri_v1.1.scaffolds.braker.filtered.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=MYZper
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=MYZvar
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_varians/v1.1/Myzus_varians_v1.1.scaffolds.braker.filtered.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

#Taxon_code=PEMphi
#Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Pemphigus_spyrothecae/)
#Id_field=1
#orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
#mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=PENnig
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Pentalonia_nigronervosa/v1/Pentalonia_nigronervosa.v1.scaffolds.braker.filtered.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=PHOcan
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Phorodon_cannabis/v1/Phorodon_cannabis_v1.scaffolds.braker.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=PHOhum
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Phorodon_humuli/v2/Phorodon_humuli_v2_scaffolds.braker.filtered.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

#Taxon_code=RHOmai #####this crashes orthofinder for some reason
#Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Rhopalosiphum_maidis/v1/rmaidis_v2.gff3.prot.fa)
#Id_field=1
#orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
#mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=RHOpad
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Rhopalosiphum_padi/JIC1_v1/Rhopalosiphum_padi_JIC1_v1.scaffolds.braker.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=SCHchi
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Schlechtendalia_chinensis/v1/proteins.fasta)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=SITave
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Sitobion_avenae/JIC1_v2.1/Sitobion_avenae_JIC1_v2.1.scaffolds.braker.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=SITfra
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Sitobion_fragariae/v1/Sitobion_fragariae_v1.scaffolds.braker.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=SITmis
Fasta_file=$(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Sitobion_miscanthi/v2/Sitobion_miscanthi_v2.scaffolds.braker.aa.fa)
Id_field=1
orthomclAdjustFasta  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta


for Dir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/orthology/orthofinder/$IsolateAbrv); do
Input_dir=$Dir/formatted
Min_length=10
Max_percent_stops=20
Good_proteins_file=$Dir/goodProteins/goodProteins.fasta
Poor_proteins_file=$Dir/badProteins/poorProteins.fasta
orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
done

sbatch ~/git_repos/Wrappers/NBI/run_orthofinder.sh $WorkDir 
#57528443,57528511

OrthogroupsTxt=$WorkDir/formatted/OrthoFinder/*/Orthogroups/Orthogroups.txt
GoodProts=$WorkDir/goodProteins/goodProteins.fasta
OutDir=$WorkDir/orthogroups_fasta
mkdir -p $OutDir
source package /tgac/software/production/bin/python-2.7.10
python ~/git_repos/Scripts/NBI/orthoMCLgroups2fasta.py --orthogroups $OrthogroupsTxt --fasta $GoodProts --out_dir $OutDir 2>&1 | tee -a $WorkDir/orthoreport.txt

mkdir ${WorkDir}/orthogroups_fasta/paired
for fasta in $(find ${WorkDir}/orthogroups_fasta -name "orthogroupOG*.fa" -exec readlink -f {} \;); do
    codes=">SITmis,>SITfra,>SITave,>SCHchi,>RHOpad,>RHOmai,>PHOhum,>PHOcan,>PENnig,>MYZvar,>MYZlyt,>MYZlig,>MYZcer,>METdir,>MACalb,>HORcor,>ERIlan,>DIUnox,>DAKvit,>CINced,>BREbra,>BRAklu,>BRAhel,>BRAcar,>APHtha,>APHrum,>APHgos,>APHgly,>APHfab,>ACYpis"
    count=$(grep -E "^($(echo $codes | tr ',' '|'))" "$fasta" | wc -l)
    if grep -q '^>MYZ' "$fasta" && [ "$count" -ge 10 ]; then
    OutFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/orthology/orthofinder/persicae_v_ligustri/orthogroups_fasta/paired/$(basename $fasta | sed 's@.fa@_paired.fa@g')
    singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/find_longest_myzlig.py $fasta $OutFile
    fi
done
#57528412, 57528547

MLSRLNSKYGLDVILVGNEAIKNARYMGKIKIEMVVTASEFYVYGKYDLDFNNKSDSKYKRQIISTEALT
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
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas/effector_candidates
for gene_multifasta in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas -name "het_MYZPE13164_O_EIv2.1_*.fa" -exec readlink -f {} \; | grep 'MYZPE13164_O_EIv2.1_0002220\|MYZPE13164_O_EIv2.1_0037470\|MYZPE13164_O_EIv2.1_0037650\|MYZPE13164_O_EIv2.1_0037670\|MYZPE13164_O_EIv2.1_0043160\|MYZPE13164_O_EIv2.1_0043750\|MYZPE13164_O_EIv2.1_0045990\|MYZPE13164_O_EIv2.1_0049040\|MYZPE13164_O_EIv2.1_0054380\|MYZPE13164_O_EIv2.1_0055470\|MYZPE13164_O_EIv2.1_0056440\|MYZPE13164_O_EIv2.1_0059180\|MYZPE13164_O_EIv2.1_0061940\|MYZPE13164_O_EIv2.1_0063800\|MYZPE13164_O_EIv2.1_0063860\|MYZPE13164_O_EIv2.1_0063890\|MYZPE13164_O_EIv2.1_0064380\|MYZPE13164_O_EIv2.1_0067290\|MYZPE13164_O_EIv2.1_0075010\|MYZPE13164_O_EIv2.1_0080430\|MYZPE13164_O_EIv2.1_0080440\|MYZPE13164_O_EIv2.1_0080950\|MYZPE13164_O_EIv2.1_0082260\|MYZPE13164_O_EIv2.1_0082270\|MYZPE13164_O_EIv2.1_0082380\|MYZPE13164_O_EIv2.1_0082580\|MYZPE13164_O_EIv2.1_0082750\|MYZPE13164_O_EIv2.1_0083010\|MYZPE13164_O_EIv2.1_0083400\|MYZPE13164_O_EIv2.1_0084880\|MYZPE13164_O_EIv2.1_0087020\|MYZPE13164_O_EIv2.1_0087460\|MYZPE13164_O_EIv2.1_0087920\|MYZPE13164_O_EIv2.1_0087930\|MYZPE13164_O_EIv2.1_0087950\|MYZPE13164_O_EIv2.1_0088410\|MYZPE13164_O_EIv2.1_0092760\|MYZPE13164_O_EIv2.1_0092990\|MYZPE13164_O_EIv2.1_0094080\|MYZPE13164_O_EIv2.1_0097540\|MYZPE13164_O_EIv2.1_0098320\|MYZPE13164_O_EIv2.1_0098420\|MYZPE13164_O_EIv2.1_0098460\|MYZPE13164_O_EIv2.1_0099380\|MYZPE13164_O_EIv2.1_0100510\|MYZPE13164_O_EIv2.1_0100520\|MYZPE13164_O_EIv2.1_0106960\|MYZPE13164_O_EIv2.1_0107140\|MYZPE13164_O_EIv2.1_0107850\|MYZPE13164_O_EIv2.1_0111080\|MYZPE13164_O_EIv2.1_0111090\|MYZPE13164_O_EIv2.1_0111120\|MYZPE13164_O_EIv2.1_0123610\|MYZPE13164_O_EIv2.1_0123730\|MYZPE13164_O_EIv2.1_0124690\|MYZPE13164_O_EIv2.1_0124700\|MYZPE13164_O_EIv2.1_0133690\|MYZPE13164_O_EIv2.1_0134410\|MYZPE13164_O_EIv2.1_0135120\|MYZPE13164_O_EIv2.1_0135180\|MYZPE13164_O_EIv2.1_0135320\|MYZPE13164_O_EIv2.1_0135620\|MYZPE13164_O_EIv2.1_0135720\|MYZPE13164_O_EIv2.1_0135840\|MYZPE13164_O_EIv2.1_0136290\|MYZPE13164_O_EIv2.1_0136330\|MYZPE13164_O_EIv2.1_0136390\|MYZPE13164_O_EIv2.1_0136470\|MYZPE13164_O_EIv2.1_0136520\|MYZPE13164_O_EIv2.1_0136530\|MYZPE13164_O_EIv2.1_0136540\|MYZPE13164_O_EIv2.1_0136550\|MYZPE13164_O_EIv2.1_0137070\|MYZPE13164_O_EIv2.1_0138160\|MYZPE13164_O_EIv2.1_0138430\|MYZPE13164_O_EIv2.1_0138570\|MYZPE13164_O_EIv2.1_0138680\|MYZPE13164_O_EIv2.1_0138800\|MYZPE13164_O_EIv2.1_0138820\|MYZPE13164_O_EIv2.1_0139180\|MYZPE13164_O_EIv2.1_0139910\|MYZPE13164_O_EIv2.1_0139920\|MYZPE13164_O_EIv2.1_0140600\|MYZPE13164_O_EIv2.1_0140720\|MYZPE13164_O_EIv2.1_0140880\|MYZPE13164_O_EIv2.1_0140910\|MYZPE13164_O_EIv2.1_0141310\|MYZPE13164_O_EIv2.1_0141320\|MYZPE13164_O_EIv2.1_0141340\|MYZPE13164_O_EIv2.1_0141350\|MYZPE13164_O_EIv2.1_0141360\|MYZPE13164_O_EIv2.1_0142140\|MYZPE13164_O_EIv2.1_0143240\|MYZPE13164_O_EIv2.1_0144130\|MYZPE13164_O_EIv2.1_0144560\|MYZPE13164_O_EIv2.1_0144610\|MYZPE13164_O_EIv2.1_0145380\|MYZPE13164_O_EIv2.1_0145510\|MYZPE13164_O_EIv2.1_0145520\|MYZPE13164_O_EIv2.1_0145590\|MYZPE13164_O_EIv2.1_0146220\|MYZPE13164_O_EIv2.1_0146520\|MYZPE13164_O_EIv2.1_0146530\|MYZPE13164_O_EIv2.1_0146650\|MYZPE13164_O_EIv2.1_0147220\|MYZPE13164_O_EIv2.1_0147240\|MYZPE13164_O_EIv2.1_0147400\|MYZPE13164_O_EIv2.1_0147410\|MYZPE13164_O_EIv2.1_0148500\|MYZPE13164_O_EIv2.1_0148540\|MYZPE13164_O_EIv2.1_0148810\|MYZPE13164_O_EIv2.1_0149720\|MYZPE13164_O_EIv2.1_0151640\|MYZPE13164_O_EIv2.1_0151650\|MYZPE13164_O_EIv2.1_0151660\|MYZPE13164_O_EIv2.1_0151700\|MYZPE13164_O_EIv2.1_0151720\|MYZPE13164_O_EIv2.1_0151730\|MYZPE13164_O_EIv2.1_0152360\|MYZPE13164_O_EIv2.1_0152530\|MYZPE13164_O_EIv2.1_0152540\|MYZPE13164_O_EIv2.1_0152630\|MYZPE13164_O_EIv2.1_0153210\|MYZPE13164_O_EIv2.1_0153400\|MYZPE13164_O_EIv2.1_0153420\|MYZPE13164_O_EIv2.1_0154480\|MYZPE13164_O_EIv2.1_0154750\|MYZPE13164_O_EIv2.1_0154880\|MYZPE13164_O_EIv2.1_0155220\|MYZPE13164_O_EIv2.1_0155260\|MYZPE13164_O_EIv2.1_0159250\|MYZPE13164_O_EIv2.1_0159280\|MYZPE13164_O_EIv2.1_0160140\|MYZPE13164_O_EIv2.1_0163620\|MYZPE13164_O_EIv2.1_0163630\|MYZPE13164_O_EIv2.1_0163650\|MYZPE13164_O_EIv2.1_0163800\|MYZPE13164_O_EIv2.1_0163810\|MYZPE13164_O_EIv2.1_0164300\|MYZPE13164_O_EIv2.1_0164440\|MYZPE13164_O_EIv2.1_0164630\|MYZPE13164_O_EIv2.1_0164660\|MYZPE13164_O_EIv2.1_0164710\|MYZPE13164_O_EIv2.1_0165190\|MYZPE13164_O_EIv2.1_0165320\|MYZPE13164_O_EIv2.1_0165390\|MYZPE13164_O_EIv2.1_0166420\|MYZPE13164_O_EIv2.1_0166450\|MYZPE13164_O_EIv2.1_0167360\|MYZPE13164_O_EIv2.1_0167600\|MYZPE13164_O_EIv2.1_0169070\|MYZPE13164_O_EIv2.1_0169320\|MYZPE13164_O_EIv2.1_0169370\|MYZPE13164_O_EIv2.1_0169450\|MYZPE13164_O_EIv2.1_0169480\|MYZPE13164_O_EIv2.1_0171150\|MYZPE13164_O_EIv2.1_0171450\|MYZPE13164_O_EIv2.1_0171820\|MYZPE13164_O_EIv2.1_0171940\|MYZPE13164_O_EIv2.1_0172040\|MYZPE13164_O_EIv2.1_0172050\|MYZPE13164_O_EIv2.1_0173370\|MYZPE13164_O_EIv2.1_0173410\|MYZPE13164_O_EIv2.1_0174480\|MYZPE13164_O_EIv2.1_0174560\|MYZPE13164_O_EIv2.1_0175570\|MYZPE13164_O_EIv2.1_0177210\|MYZPE13164_O_EIv2.1_0178420\|MYZPE13164_O_EIv2.1_0178600\|MYZPE13164_O_EIv2.1_0178740\|MYZPE13164_O_EIv2.1_0179580\|MYZPE13164_O_EIv2.1_0179840\|MYZPE13164_O_EIv2.1_0181430\|MYZPE13164_O_EIv2.1_0181660\|MYZPE13164_O_EIv2.1_0181980\|MYZPE13164_O_EIv2.1_0182210\|MYZPE13164_O_EIv2.1_0182260\|MYZPE13164_O_EIv2.1_0182520\|MYZPE13164_O_EIv2.1_0183050\|MYZPE13164_O_EIv2.1_0183710\|MYZPE13164_O_EIv2.1_0183960\|MYZPE13164_O_EIv2.1_0183970\|MYZPE13164_O_EIv2.1_0184720\|MYZPE13164_O_EIv2.1_0185380\|MYZPE13164_O_EIv2.1_0185430\|MYZPE13164_O_EIv2.1_0186140\|MYZPE13164_O_EIv2.1_0186640\|MYZPE13164_O_EIv2.1_0186650\|MYZPE13164_O_EIv2.1_0187090\|MYZPE13164_O_EIv2.1_0187420\|MYZPE13164_O_EIv2.1_0188030\|MYZPE13164_O_EIv2.1_0188050\|MYZPE13164_O_EIv2.1_0188720\|MYZPE13164_O_EIv2.1_0189060\|MYZPE13164_O_EIv2.1_0190050\|MYZPE13164_O_EIv2.1_0190480\|MYZPE13164_O_EIv2.1_0191000\|MYZPE13164_O_EIv2.1_0193570\|MYZPE13164_O_EIv2.1_0195210\|MYZPE13164_O_EIv2.1_0195220\|MYZPE13164_O_EIv2.1_0195600\|MYZPE13164_O_EIv2.1_0196170\|MYZPE13164_O_EIv2.1_0199840\|MYZPE13164_O_EIv2.1_0200020\|MYZPE13164_O_EIv2.1_0200040\|MYZPE13164_O_EIv2.1_0200140\|MYZPE13164_O_EIv2.1_0200400\|MYZPE13164_O_EIv2.1_0200950\|MYZPE13164_O_EIv2.1_0200980\|MYZPE13164_O_EIv2.1_0201680\|MYZPE13164_O_EIv2.1_0204830\|MYZPE13164_O_EIv2.1_0206480\|MYZPE13164_O_EIv2.1_0206630\|MYZPE13164_O_EIv2.1_0206880\|MYZPE13164_O_EIv2.1_0206970\|MYZPE13164_O_EIv2.1_0207040\|MYZPE13164_O_EIv2.1_0207270\|MYZPE13164_O_EIv2.1_0207370\|MYZPE13164_O_EIv2.1_0207560\|MYZPE13164_O_EIv2.1_0207890\|MYZPE13164_O_EIv2.1_0207900\|MYZPE13164_O_EIv2.1_0208060\|MYZPE13164_O_EIv2.1_0208180\|MYZPE13164_O_EIv2.1_0209200\|MYZPE13164_O_EIv2.1_0209240\|MYZPE13164_O_EIv2.1_0209300\|MYZPE13164_O_EIv2.1_0209730\|MYZPE13164_O_EIv2.1_0210180\|MYZPE13164_O_EIv2.1_0210270\|MYZPE13164_O_EIv2.1_0210920\|MYZPE13164_O_EIv2.1_0210940\|MYZPE13164_O_EIv2.1_0210990\|MYZPE13164_O_EIv2.1_0211040\|MYZPE13164_O_EIv2.1_0211490\|MYZPE13164_O_EIv2.1_0211990\|MYZPE13164_O_EIv2.1_0212120\|MYZPE13164_O_EIv2.1_0213140\|MYZPE13164_O_EIv2.1_0213480\|MYZPE13164_O_EIv2.1_0213680\|MYZPE13164_O_EIv2.1_0213980\|MYZPE13164_O_EIv2.1_0216020\|MYZPE13164_O_EIv2.1_0216250\|MYZPE13164_O_EIv2.1_0216310\|MYZPE13164_O_EIv2.1_0216320\|MYZPE13164_O_EIv2.1_0216330\|MYZPE13164_O_EIv2.1_0216740\|MYZPE13164_O_EIv2.1_0217290\|MYZPE13164_O_EIv2.1_0217510\|MYZPE13164_O_EIv2.1_0217890\|MYZPE13164_O_EIv2.1_0218140\|MYZPE13164_O_EIv2.1_0219900\|MYZPE13164_O_EIv2.1_0220560\|MYZPE13164_O_EIv2.1_0220790\|MYZPE13164_O_EIv2.1_0220980\|MYZPE13164_O_EIv2.1_0221360\|MYZPE13164_O_EIv2.1_0221370\|MYZPE13164_O_EIv2.1_0221890\|MYZPE13164_O_EIv2.1_0222540\|MYZPE13164_O_EIv2.1_0223090\|MYZPE13164_O_EIv2.1_0223940\|MYZPE13164_O_EIv2.1_0224920\|MYZPE13164_O_EIv2.1_0225090\|MYZPE13164_O_EIv2.1_0226260\|MYZPE13164_O_EIv2.1_0226530\|MYZPE13164_O_EIv2.1_0227980\|MYZPE13164_O_EIv2.1_0228650\|MYZPE13164_O_EIv2.1_0228810\|MYZPE13164_O_EIv2.1_0228840\|MYZPE13164_O_EIv2.1_0229070\|MYZPE13164_O_EIv2.1_0229570\|MYZPE13164_O_EIv2.1_0230870\|MYZPE13164_O_EIv2.1_0231260\|MYZPE13164_O_EIv2.1_0231810\|MYZPE13164_O_EIv2.1_0232600\|MYZPE13164_O_EIv2.1_0233810\|MYZPE13164_O_EIv2.1_0234500\|MYZPE13164_O_EIv2.1_0234520\|MYZPE13164_O_EIv2.1_0235230\|MYZPE13164_O_EIv2.1_0235660\|MYZPE13164_O_EIv2.1_0236470\|MYZPE13164_O_EIv2.1_0236800\|MYZPE13164_O_EIv2.1_0238080\|MYZPE13164_O_EIv2.1_0241310\|MYZPE13164_O_EIv2.1_0241990\|MYZPE13164_O_EIv2.1_0242890\|MYZPE13164_O_EIv2.1_0242910\|MYZPE13164_O_EIv2.1_0242970\|MYZPE13164_O_EIv2.1_0243080\|MYZPE13164_O_EIv2.1_0243090\|MYZPE13164_O_EIv2.1_0243160\|MYZPE13164_O_EIv2.1_0245060\|MYZPE13164_O_EIv2.1_0246590\|MYZPE13164_O_EIv2.1_0246720\|MYZPE13164_O_EIv2.1_0248170\|MYZPE13164_O_EIv2.1_0248740\|MYZPE13164_O_EIv2.1_0249490\|MYZPE13164_O_EIv2.1_0249680\|MYZPE13164_O_EIv2.1_0250270\|MYZPE13164_O_EIv2.1_0250940\|MYZPE13164_O_EIv2.1_0251130\|MYZPE13164_O_EIv2.1_0253480\|MYZPE13164_O_EIv2.1_0253730\|MYZPE13164_O_EIv2.1_0254200\|MYZPE13164_O_EIv2.1_0254240\|MYZPE13164_O_EIv2.1_0254250\|MYZPE13164_O_EIv2.1_0254760\|MYZPE13164_O_EIv2.1_0255000\|MYZPE13164_O_EIv2.1_0256470\|MYZPE13164_O_EIv2.1_0258410\|MYZPE13164_O_EIv2.1_0258780\|MYZPE13164_O_EIv2.1_0259530\|MYZPE13164_O_EIv2.1_0259540\|MYZPE13164_O_EIv2.1_0259560\|MYZPE13164_O_EIv2.1_0259700\|MYZPE13164_O_EIv2.1_0260190\|MYZPE13164_O_EIv2.1_0260300\|MYZPE13164_O_EIv2.1_0260310\|MYZPE13164_O_EIv2.1_0261410\|MYZPE13164_O_EIv2.1_0261420\|MYZPE13164_O_EIv2.1_0262120\|MYZPE13164_O_EIv2.1_0266490\|MYZPE13164_O_EIv2.1_0266500\|MYZPE13164_O_EIv2.1_0266760\|MYZPE13164_O_EIv2.1_0267430\|MYZPE13164_O_EIv2.1_0267480\|MYZPE13164_O_EIv2.1_0268430\|MYZPE13164_O_EIv2.1_0268970\|MYZPE13164_O_EIv2.1_0269400\|MYZPE13164_O_EIv2.1_0270740\|MYZPE13164_O_EIv2.1_0271340\|MYZPE13164_O_EIv2.1_0272280\|MYZPE13164_O_EIv2.1_0272350\|MYZPE13164_O_EIv2.1_0272650\|MYZPE13164_O_EIv2.1_0272700\|MYZPE13164_O_EIv2.1_0273350\|MYZPE13164_O_EIv2.1_0273540\|MYZPE13164_O_EIv2.1_0275060\|MYZPE13164_O_EIv2.1_0275340\|MYZPE13164_O_EIv2.1_0275790\|MYZPE13164_O_EIv2.1_0276040\|MYZPE13164_O_EIv2.1_0276660\|MYZPE13164_O_EIv2.1_0277210\|MYZPE13164_O_EIv2.1_0277330\|MYZPE13164_O_EIv2.1_0279010\|MYZPE13164_O_EIv2.1_0279090\|MYZPE13164_O_EIv2.1_0279280\|MYZPE13164_O_EIv2.1_0279710\|MYZPE13164_O_EIv2.1_0279720\|MYZPE13164_O_EIv2.1_0280140\|MYZPE13164_O_EIv2.1_0281260\|MYZPE13164_O_EIv2.1_0281280\|MYZPE13164_O_EIv2.1_0281380\|MYZPE13164_O_EIv2.1_0281480\|MYZPE13164_O_EIv2.1_0282080\|MYZPE13164_O_EIv2.1_0282100\|MYZPE13164_O_EIv2.1_0282500\|MYZPE13164_O_EIv2.1_0282560\|MYZPE13164_O_EIv2.1_0282590\|MYZPE13164_O_EIv2.1_0282600\|MYZPE13164_O_EIv2.1_0283700\|MYZPE13164_O_EIv2.1_0283790\|MYZPE13164_O_EIv2.1_0284230\|MYZPE13164_O_EIv2.1_0284250\|MYZPE13164_O_EIv2.1_0284700\|MYZPE13164_O_EIv2.1_0285760\|MYZPE13164_O_EIv2.1_0285960\|MYZPE13164_O_EIv2.1_0287290\|MYZPE13164_O_EIv2.1_0287910\|MYZPE13164_O_EIv2.1_0287930\|MYZPE13164_O_EIv2.1_0289020\|MYZPE13164_O_EIv2.1_0289480\|MYZPE13164_O_EIv2.1_0290630\|MYZPE13164_O_EIv2.1_0292910\|MYZPE13164_O_EIv2.1_0292940\|MYZPE13164_O_EIv2.1_0292960\|MYZPE13164_O_EIv2.1_0293740\|MYZPE13164_O_EIv2.1_0294670\|MYZPE13164_O_EIv2.1_0295400\|MYZPE13164_O_EIv2.1_0295790\|MYZPE13164_O_EIv2.1_0295860\|MYZPE13164_O_EIv2.1_0295880\|MYZPE13164_O_EIv2.1_0296460\|MYZPE13164_O_EIv2.1_0297030\|MYZPE13164_O_EIv2.1_0297270\|MYZPE13164_O_EIv2.1_0298300\|MYZPE13164_O_EIv2.1_0298620\|MYZPE13164_O_EIv2.1_0300870\|MYZPE13164_O_EIv2.1_0300970\|MYZPE13164_O_EIv2.1_0301280\|MYZPE13164_O_EIv2.1_0302000\|MYZPE13164_O_EIv2.1_0302010\|MYZPE13164_O_EIv2.1_0303250\|MYZPE13164_O_EIv2.1_0304300\|MYZPE13164_O_EIv2.1_0304970\|MYZPE13164_O_EIv2.1_0306100\|MYZPE13164_O_EIv2.1_0306680\|MYZPE13164_O_EIv2.1_0309020\|MYZPE13164_O_EIv2.1_0309030\|MYZPE13164_O_EIv2.1_0309330\|MYZPE13164_O_EIv2.1_0310090\|MYZPE13164_O_EIv2.1_0311040\|MYZPE13164_O_EIv2.1_0312120\|MYZPE13164_O_EIv2.1_0317580\|MYZPE13164_O_EIv2.1_0318530\|MYZPE13164_O_EIv2.1_0318960\|MYZPE13164_O_EIv2.1_0320460\|MYZPE13164_O_EIv2.1_0320470\|MYZPE13164_O_EIv2.1_0323070\|MYZPE13164_O_EIv2.1_0325490\|MYZPE13164_O_EIv2.1_0325510\|MYZPE13164_O_EIv2.1_0326160\|MYZPE13164_O_EIv2.1_0326560\|MYZPE13164_O_EIv2.1_0328570\|MYZPE13164_O_EIv2.1_0328670\|MYZPE13164_O_EIv2.1_0329710\|MYZPE13164_O_EIv2.1_0332210\|MYZPE13164_O_EIv2.1_0333430\|MYZPE13164_O_EIv2.1_0333880\|MYZPE13164_O_EIv2.1_0333970\|MYZPE13164_O_EIv2.1_0336300\|MYZPE13164_O_EIv2.1_0336990\|MYZPE13164_O_EIv2.1_0337240\|MYZPE13164_O_EIv2.1_0337250\|MYZPE13164_O_EIv2.1_0337260\|MYZPE13164_O_EIv2.1_0337270\|MYZPE13164_O_EIv2.1_0337580\|MYZPE13164_O_EIv2.1_0337590\|MYZPE13164_O_EIv2.1_0338150\|MYZPE13164_O_EIv2.1_0338940\|MYZPE13164_O_EIv2.1_0339430\|MYZPE13164_O_EIv2.1_0339930\|MYZPE13164_O_EIv2.1_0340140\|MYZPE13164_O_EIv2.1_0341080\|MYZPE13164_O_EIv2.1_0342070\|MYZPE13164_O_EIv2.1_0343870\|MYZPE13164_O_EIv2.1_0343990\|MYZPE13164_O_EIv2.1_0345790\|MYZPE13164_O_EIv2.1_0345820\|MYZPE13164_O_EIv2.1_0345840\|MYZPE13164_O_EIv2.1_0345860\|MYZPE13164_O_EIv2.1_0345880\|MYZPE13164_O_EIv2.1_0346220\|MYZPE13164_O_EIv2.1_0346420\|MYZPE13164_O_EIv2.1_0347830\|MYZPE13164_O_EIv2.1_0347890\|MYZPE13164_O_EIv2.1_0348170\|MYZPE13164_O_EIv2.1_0348220\|MYZPE13164_O_EIv2.1_0348550\|MYZPE13164_O_EIv2.1_0348900\|MYZPE13164_O_EIv2.1_0349310\|MYZPE13164_O_EIv2.1_0349540\|MYZPE13164_O_EIv2.1_0349610\|MYZPE13164_O_EIv2.1_0350040\|MYZPE13164_O_EIv2.1_0350050\|MYZPE13164_O_EIv2.1_0350350\|MYZPE13164_O_EIv2.1_0350610\|MYZPE13164_O_EIv2.1_0350640\|MYZPE13164_O_EIv2.1_0351370\|MYZPE13164_O_EIv2.1_0351790\|MYZPE13164_O_EIv2.1_0353780\|MYZPE13164_O_EIv2.1_0353810\|MYZPE13164_O_EIv2.1_0354030\|MYZPE13164_O_EIv2.1_0357250\|MYZPE13164_O_EIv2.1_0357590\|MYZPE13164_O_EIv2.1_0358240\|MYZPE13164_O_EIv2.1_0358520\|MYZPE13164_O_EIv2.1_0358840\|MYZPE13164_O_EIv2.1_0358890\|MYZPE13164_O_EIv2.1_0359910\|MYZPE13164_O_EIv2.1_0359970\|MYZPE13164_O_EIv2.1_0361760\|MYZPE13164_O_EIv2.1_0362100\|MYZPE13164_O_EIv2.1_0362450\|MYZPE13164_O_EIv2.1_0362830\|MYZPE13164_O_EIv2.1_0362930\|MYZPE13164_O_EIv2.1_0362940\|MYZPE13164_O_EIv2.1_0363870\|MYZPE13164_O_EIv2.1_0364310\|MYZPE13164_O_EIv2.1_0364360\|MYZPE13164_O_EIv2.1_0365150\|MYZPE13164_O_EIv2.1_0365510\|MYZPE13164_O_EIv2.1_0365860\|MYZPE13164_O_EIv2.1_0365920\|MYZPE13164_O_EIv2.1_0366100\|MYZPE13164_O_EIv2.1_0366700\|MYZPE13164_O_EIv2.1_0366730\|MYZPE13164_O_EIv2.1_0366970\|MYZPE13164_O_EIv2.1_0367250\|MYZPE13164_O_EIv2.1_0367640\|MYZPE13164_O_EIv2.1_0367680\|MYZPE13164_O_EIv2.1_0367740\|MYZPE13164_O_EIv2.1_0367750\|MYZPE13164_O_EIv2.1_0367760\|MYZPE13164_O_EIv2.1_0367810\|MYZPE13164_O_EIv2.1_0367980\|MYZPE13164_O_EIv2.1_0369920\|sequence-region\|gff-version'); do 
cp $gene_multifasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas/effector_candidates/.
done

 snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/effector_candidates.gff3
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