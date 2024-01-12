# Mitochondrial SNP calling as part of the M. persicae population genomics project

Contains SNP analysis of the M. persicae mitochondria

## Contents
1. [Call SNPs against the M. persicae mitochondria](#1)
        1. [Collecting data](#2)
        2. [QC](#3)
                1. [FASTQC and qualimap](#4)
                2. [Seqtk Subsampling](#5)
        3. [Trimming](#6)
                1. [Trim_galore](#7)
        4. [Alignment](#8)
                1. [bwa-mem alignment](#9)
                2. [Picard and GATK](#10)
        5. [Variant calling](#11)
        6. [Filtering](#12)
                1. [Filter for SNP quality and missingness](#13)
2. [Network analysis](#14)
        1. [Distance matrix](#15)
        2. [Splitstree](#16)

Look for evidence of speciation between the aphid samples by looking at SNPs in the mitochondria:

## Call SNPs against the M. persicae mitochondria <a name="1"></a>
### Collecting data <a name="2"></a>
Raw data was collected by R. Wouters and Singh et al. and used to perform genomic SNP calling previously. However, many intermediate files have since been deleted to save space, eg. the Singh data has been removed entirely as it is backed up with NCBI. Additionally there is no clear record of breadth + depth of coverage at each stage of the data processing pipeline, which would be ideal to have. Trimmed subsampled reads will need to be regenerated.

Download Singh data from NCBI to local machine:
```bash
C:\Users\did23faz\Downloads\sratoolkit.current-win64\sratoolkit.3.0.5-win64\bin>prefetch --max-size 1t --option-file C:\Users\did23faz\Documents\accessions.txt --output-directory \\jic-hpc-data\Group-Scratch\Saskia-Hogenhout\tom_heaven\Aphididae\raw_data\Myzus\persicae\singh\download
```
Copy Singh data from local machine to the HPC and extract, note that some samples have a different naming scheme with Singh and on the HPC inherited from R. Wouters.
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

Wouters raw data is already on the HPC.

Check whether the mitochondrial genome has been removed from the O_v2 assemlby:
```bash
source package d6092385-3a81-49d9-b044-8ffb85d0c446
mkdir M_persicae_mitochondria_blast_db/
cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mito1.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mito2.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mito3.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mito4.fasta > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mito.fasta
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/mito.fasta -dbtype nucl -parse_seqids -out M_persicae_mitochondria_blast_db/blast_db
blastn -query /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa -db M_persicae_mitochondria_blast_db/blast_db -out results.txt -evalue 1e-5 -outfmt 6 -num_threads 1
```
The mitochondrial genome has been removed from the O_v2 assemlby.

The M. persicae mitochondrial reference genome was downloaded from NCBI (Accession: NC_029727).

### QC <a name="3"></a>
#### FASTQC and qualimap <a name="4"></a>

Quality and coverage of the raw reads was assessed
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/*/*); do
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
```
#### Seqtk Subsampling <a name="5"></a>

Raw reads were subsampled so that samples did not have greater than 15X coverage:
Note: this is 15X coverage vs. the nuclear genome, not the mitochondria
```bash
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
```
The quality and coverage was reassessed for the subsampled reads:
```bash
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
```
Subsampled reads have been deleted to save space - can be regenerated from raw reads if needed.

### Trimming <a name="6"></a>
#### Trim_galore <a name="7"></a>
Reads were trimmed to removed poor quality regions and adapters using trim_galore.

#NOTE: to save space the script has been edited to delete input files upon production of outputs.
```bash
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
```
The quality and coverage was reassessed for the trimmed reads:
```bash
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

### Alignment <a name="8"></a>
#### bwa-mem alignment <a name="9"></a>

Subsampled, trimmed reads were aligned to the mitochondrial reference genome:
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
Note, Reads previously aligned to the nuclear genome by Roland are here, I do not have the permissions to symlink them into my file system:
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
#### Picard and GATK <a name="10"></a>

Files were sorted by scaffold and coordinate level, duplicates were marked and removed and the files were re-indexed:
NOTE: this script does not remove the input file and .bam files take up a lot of space - was therefore actually run 1 step at a time
```bash
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/SNPs/mito/*/*/bwa/*_mito.bam); do
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_sort_bam.sh $file
done 
```
Input files were subsequently deleted to save space
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
### Variant calling <a name="11"></a>

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
```
bam files were subsequently removed to save space

### Filtering <a name="12"></a>

I have not done mappability masking for the mitochondrial SNPs, the programs used to do this with the nuclear SNPs (GEM-Mapper) do not seem to be installed on the JIC HPC, given that these SNPs are called against only the ~17.5k bp mitochondrial genome which is not expected to have lengthy repetative regions it does not seem worth it to install new software to perform this filtering.

#### Filter for SNP quality and missingness <a name="13"></a>

Samples which have been removed from the genomic SNP analysis due to high missingness etc. were removed from the mitochondrial data.

SNPs were filtered to keep only bi-allelic SNPs, with minimum depth of 40, minimum quality score of 30, and maximum depth of 500. The depth of sequencing for the mitochondria is far higher than that for the nuclear genome.
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
We are somehow missing a sample from the analysis :(, it doesnt seem worth repeating it all to fix this.

## Network analysis <a name="14"></a>

Investigate the relatedness of the samples, based upon the mitochondrial SNPs.

#### Distance matrix <a name="15"></a>

Convert from a VCF file to a distance matrix:
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
Reformat:
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

# Example usage:
csv_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/mito_p_dis_-g006.csv'
output_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/mito_p_dis_-g006.dist'

distance_matrix = read_distance_matrix_from_csv(csv_file)
convert_to_phylip(distance_matrix, output_file)
```
#### Splitstree <a name="16"></a>

Construct a phylogenetic network for the samples, based upon the mitochondrial SNP data.
```bash
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree
#Import failed: unknown format or error in file (eg. illegal characters or multiple occurrences of same taxon name)
#NOTE: G006 is triggering the above error, the mitochondrial sequence used as reference is from G006 -> there are therefore no SNPs for this sample causing the error? G006 was therefore reomved during filtering
<<<<<<< HEAD
```
=======
```
>>>>>>> 9bb487a11dfffd64de11be5c3e2043aee326f5d0
