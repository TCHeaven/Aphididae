# Reproducing the work done by Jitender
Unless stated otherwise work was performed from the directory: /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae

Collecting data:
```bash
mkdir -p snp_calling/Myzus/persicae/biello/gatk/filtered
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mperc-analysis-jitender/saskia/vcffilesafterfiltering/210s.M_persicae.onlySNPs.vcf.gz snp_calling/Myzus/persicae/biello/gatk/filtered/.

ln -s /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3

ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa.gz
```
#### P_distance
p-distance  is the proportion (p) of nucleotide sites at which two sequences being compared are different. It is obtained by dividing the number of nucleotide differences by the total number of nucleotides compared.

distmat calculates the evolutionary distance between every pair of sequences in a multiple sequence alignment. A variety of methods to estimate distance may be selected, and differ in how they correct the observed substitution rates to more accurately reflect the true evolutionary distance. An output file containing a distance matrix for the set of sequences is written. The distances are expressed in terms of the number of substitutions per 100 bases or amino acids.
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
Missing data?:
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

219
10
15
65
70
88
89
174
186
189
198
199
204
207
211
213
217

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
```
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
Generate an ID based p-distmat in for R
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
#### Corehunter
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
Rarefaction
```bash
head -n 2000 snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf > snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs_1000.vcf

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
Genic regions:
```bash
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
bgzip -c /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf > /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/210s.M_persicae.onlySNPs-genic-regions.vcf.gz

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


singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
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
exit()

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

```
```python
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc

rc('mathtext', default='regular')

df = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/rarefied-HomozygousALT-SNPS-genic-193-100-193-curve.csv")

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
suptitle = plt.suptitle('Nuclear genic SNP rarefaction curve with 193 samples using 100 replicates, step size 1', y=1.02)

plt.tight_layout()

# When you save the fig, add the suptitle text object as an extra artist
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genic-193-100-193-curve.pdf", bbox_extra_artists=(suptitle,), bbox_inches="tight") 
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/rarefaction/genic-193-100-193-curve.png", bbox_extra_artists=(suptitle,), bbox_inches="tight") 

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
For Bass and JIC samples seperately
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
Bass vs JIC
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

################################################################################################################

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

######################################################################################
#Scaled x axis:

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
```
#### Gemma
```bash
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz | awk 'BEGIN{OFS="\t"} {if ($3 == ".") $3 = $1"_"$2"_"$5; print}' | gzip > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.vcf.gz

zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.vcf.gz | sed 's/scaffold_//g' | gzip > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.tmp.vcf.gz && mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.tmp.vcf.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.vcf.gz

source package /tgac/software/testing/bin/plink-1.90

plink --vcf temp.vcf  --recode --out temp
plink --allow-extra-chr --file temp --make-bed --out temp
```

Graphs were plotted in excel
#### vcftools plot of saturation curve
```bash
gzip -cd snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.vcf.gz | head -n 2000 > temp.vcf
rm temp.vcf

#Remove SNPs without variant:
source package /nbi/software/testing/bin/vcftools-0.1.15
vcftools --gzvcf snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.vcf.gz --mac 1 --recode --recode-INFO-all --out snp_calling/Myzus/persicae/biello/gatk/filtered/210s.M_persicae.onlySNPs.mac1.vcf
#After filtering, kept 7,873,274 out of a possible 11,870,914 Sites = 1/3 removed - presumably most from from lugustri

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
