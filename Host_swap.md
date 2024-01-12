# Three host swap 

Contains analysis of the 3 host swap experiment performed by M. Gravino and others

## Contents
1. [Collect data](#1)
2. [Check mutations in the aphid genomes over the course of the experiment with WGS](#2)
        1. [QC](#3)
                1. [fastqc and qualimap](#4)
        2. [Trimming](#5)
                1. [Trim galore](#6)
        3. [Alignment](#7)
                1. [Picard and GATK](#8)
        4. [Variant calling](#9)
        5. [Filter](#10)
                1. [Genmap filter for genome mappability](#11)
                2. [Filter samples for missingness:](#12)
                3. [Filter for SNP quality:](#13)
3. [Check for epigenetic changes over the course of the experiment with WGBS](#14)
                1. [fastqc and qualimap](#15)
        2. [Trimming](#16)
                1. [Trim galore](#17)
        3. [BSmooth](#18)
                1. [Bismark alignment](#19)
                2. [BSsmooth](#20)
        4. [Methylkit](#21)
                1. [BSmap alignment](#22)
                2. [Filter cytosine postions](#23)
                3. [Methylkit](#24)
                4. [Genomation](#25)
                5. [Sliding window](#26)
        5. [Custom](#27)
4. [Check for transcriptional changes over the course of the experiments 1 & 2](#28)
        1. [QC](#29)
                1. [fastqc](#30)
        2. [Trimming](#31)
        3. [Mapping](#32)
                1. [Alignment free approach](#33)
                        1. [Trinity](#34)
                2. [Alignment approach](#35)
                        1. [STAR](#36)
                        2. [Stringtie](#37)

## Collect data <a name="1"></a>
Data was copied from storage with the Swarbreck group to the Hogenhout scratch space (not backed up as already backed up with Swarbreck group).

All files that look like analysis not raw data have been removed/not copied over in order to save space.
```bash
#Host Swap Data:
#WGS Raw data:
cp -r /ei/projects/b/bda78a5f-7de9-4e0c-813f-1b1105bd6c24/WGS_Feb2021/* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/.
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021
unzip X201SC21011755-Z01-F001.zip
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/BR25_E2/BR25_E2_FDSW210041715-1a_HNTVNDSXY_L3_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/BR25_E2/BR25_E2_FDSW210041715-1a_HNTVNDSXY_L3_3.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/BR25_E2/BR25_E2_FDSW210041715-1a_HNTVNDSXY_L3_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/BR25_E2/BR25_E2_FDSW210041715-1a_HNTVNDSXY_L3_4.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/BR39_E1/BR39_E1_FDSW210041711-1a_HKJGKDSXY_L2_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/BR39_E1/BR39_E1_FDSW210041711-1a_HKJGKDSXY_L2_3.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/BR39_E1/BR39_E1_FDSW210041711-1a_HKJGKDSXY_L2_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/BR39_E1/BR39_E1_FDSW210041711-1a_HKJGKDSXY_L2_4.fq.gz

#WGBS Raw data:
cp /ei/projects/b/bda78a5f-7de9-4e0c-813f-1b1105bd6c24/WGBS_Mar2021/X201SC21011757-Z01-F001_*.zip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGBS/Archana_Mar2021/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGBS/Archana_Mar2021/At/* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGBS/Archana_Mar2021/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGBS/Archana_Mar2021/Nb/* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGBS/Archana_Mar2021/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGBS/Archana_Mar2021/Br/* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGBS/Archana_Mar2021/.
for dir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGBS/Archana_Mar2021/*/); do
    counter=1
    for file in $(ls ${dir}*_1.fq.gz); do
        sample_name=$(basename ${file} _1.fq.gz)
        echo ${dir}${sample_name}_1.fq.gz
        echo ${dir}${sample_name}_${counter}.fq.gz
        echo ${dir}${sample_name}_2.fq.gz 
        echo ${dir}${sample_name}_$((counter + 1)).fq.gz
        mv ${dir}${sample_name}_1.fq.gz ${dir}${sample_name}_${counter}.fq.gz
        mv ${dir}${sample_name}_2.fq.gz ${dir}${sample_name}_$((counter + 1)).fq.gz
        counter=$((counter + 2))
    done
    echo _
    echo _
done

#RNA data:
#Host swap experiment 1:
cp -r /ei/projects/b/bda78a5f-7de9-4e0c-813f-1b1105bd6c24/host_adjustment/Data_Dec2020/* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/.
#All files unzipped
mkdir At
mkdir Br
mkdir Nb
mv AT* At/.
mv BR* Br/.
mv NB* Nb/.

#Host swap experiment 2:
cp -r /ei/projects/b/bda78a5f-7de9-4e0c-813f-1b1105bd6c24/host_adjustment/HostAdaptation/raw_data/X201SC20052390-Z01-F001_0*_20200727.zip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/.
#All files unzipped
mkdir At
mkdir Br
mkdir Nb
mv AT* At/.
mv BR* Br/.
mv NB* Nb/.
#Experimet 2, week 9 results appear to have been stored in the wrong place
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/At/AT9_1 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/At/. 
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/At/AT9_2 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/At/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/At/AT9_3 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/At/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/At/AT9_4 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/At/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/At/AT9_5 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/At/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/Br/BR9_1 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/Br/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/Br/BR9_2 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/Br/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/Br/BR9_3 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/Br/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/Br/BR9_4 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/Br/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/Br/BR9_5 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/Br/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/Nb/NB9_1 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/Nb/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/Nb/NB9_2 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/Nb/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/Nb/NB9_3 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/Nb/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/Nb/NB9_4 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/Nb/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/Nb/NB9_5 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/Nb/.

#Host swap on Knock out plants:
cp -r /ei/projects/b/bda78a5f-7de9-4e0c-813f-1b1105bd6c24/RNA_Mar2021/* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Mar2021/.
#All files unzipped
mkdir BR
mkdir NB
mv BR_* BR/.
mv NB_* NB/.

########################################################################################################################
#Organ data raw data:
cp -r /ei/projects/a/a93e4b69-58e3-495b-ac4d-04e978fed5d1/data/organ_data/raw_data/for_archana/* /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/archana_organ_data/.

#This is analysis - removed to make room for all raw data
cp -r analysis /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/Archana_hostadaptation_analysis/.
cp -r Effector /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/Archana_hostadaptation_analysis/.
cp -r * /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/Archana_hostadaptation_analysis/alternative_splicing/.
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/Archana_hostadaptation_analysis
```
## Check mutations in the aphid genomes over the course of the experiment with WGS <a name="2"></a>
This is the purpose of the WGS data which was collected at the start and end of each experimental repeat from each host, ie.: generation 1 of aphids exposed to Br in experiments 1 and 2 (E1, E2), the generation 39 of aphids exposed to Br, At, Nb in E1, and the generation 25 of aphids exposed to Br, At, Nb in E2. 

### QC  <a name="3"></a>
#### fastqc and qualimap <a name="4"></a>
Raw reads were assessed for quality and coverage of the clone O_v2 reference assembly:
```bash
for ReadDir in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/ -mindepth 1 -type d); do
    if [ ! -e ${ReadDir}/qualimap/*genome_results_gff.txt ]; then
    echo Running for:
    ls ${ReadDir}/qualimap/*genome_results_gff.txt
    Fread=$(ls ${ReadDir}/*_1.fq.gz)
    Rread=$(ls ${ReadDir}/*_2.fq.gz)
    Fread2=$(ls ${ReadDir}/*_3.fq.gz)
    Rread2=$(ls ${ReadDir}/*_4.fq.gz)
    OutDir=$(echo ${ReadDir})
    Reference_genome=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
    Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread $Fread2 $Rread2 
    else
    echo Already Done:
    ls ${ReadDir}/qualimap/*genome_results.txt
    fi
done
#57483656-63

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/*/); do
sample=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
coverage=$(grep 'mean coverageData' ${ReadDir}qualimap/*_genome_results_gff.txt | rev | cut -d ' ' -f1 | rev | sed 's@X@@g')
echo $sample raw reads have average coverage of ${coverage}
done
```
BR1_E1 raw reads have average coverage of 23.8422
BR39_E1 raw reads have average coverage of 23.2599
AT39_E1 raw reads have average coverage of 20.6892
NB39_E1 raw reads have average coverage of 21.0765

BR1_E2 raw reads have average coverage of 20.0134
BR25_E2 raw reads have average coverage of 23.7613
AT25_E2 raw reads have average coverage of 19.7724
NB25_E2 raw reads have average coverage of 21.6007

### Trimming <a name="5"></a>
#### Trim galore <a name="6"></a>
Adapters and low quality regions were trimmed from raw reads via trim-galore:
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/*_E*); do
    sample=$(echo $ReadDir | rev | cut -d '/' -f1 | rev)
    Fread=$(ls $ReadDir/*_1*fq.gz)
    Rread=$(ls $ReadDir/*_2*fq.gz)
    Fread2=$(ls $ReadDir/*_3*fq.gz)
    Rread2=$(ls $ReadDir/*_4*fq.gz)
    OutDir=$(echo $ReadDir | sed 's@raw_data@dna_qc@g')/trim_galore
    OutFile=${sample}_trimmed
    Quality=20
    Length=50
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir -p $OutDir
    sbatch $ProgDir/run_trim_galore.sh $OutDir $OutFile $Quality $Length $Fread $Rread $Fread2 $Rread2 
done 
#57483703-10
```
Trimmed reads were re-assessed for quality and coverage of the clone O_v2 reference assembly:
```bash
for ReadDir in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/WGS/Archana_Feb2021/ -mindepth 1 -type d); do
    if [ ! -e ${ReadDir}/qualimap/*genome_results_gff.txt ]; then
    echo Running for:
    ls ${ReadDir}/qualimap/*genome_results_gff.txt
    Fread=$(ls ${ReadDir}/*_1.fq.gz)
    Rread=$(ls ${ReadDir}/*_2.fq.gz)
    Fread2=$(ls ${ReadDir}/*_3.fq.gz)
    Rread2=$(ls ${ReadDir}/*_4.fq.gz)
    OutDir=$(echo ${ReadDir})
    Reference_genome=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
    Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread $Fread2 $Rread2 
    else
    echo Already Done:
    ls ${ReadDir}/qualimap/*genome_results.txt
    fi
done
#57490888-57490895

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/WGS/Archana_Feb2021/*/trim_galore); do
sample=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
coverage=$(grep 'mean coverageData' ${ReadDir}qualimap/*_genome_results_gff.txt | rev | cut -d ' ' -f1 | rev | sed 's@X@@g')
echo $sample raw reads have average coverage of ${coverage}
done
```
BR1_E1 raw reads have average coverage of 23.7446
BR39_E1 raw reads have average coverage of 23.161
AT39_E1 raw reads have average coverage of 20.6012
NB39_E1 raw reads have average coverage of 20.9847

BR1_E2 raw reads have average coverage of 19.9298
BR25_E2 raw reads have average coverage of 23.6622
AT25_E2 raw reads have average coverage of 19.69
NB25_E2 raw reads have average coverage of 21.5096

### Alignment <a name="7"></a>

Alignment to reference genome is performed as part of run_raw_read_qc.sh with bwa-mem (previous step).

#### Picard and GATK <a name="8"></a>

Files were sorted by scaffold and coordinate level, duplicates were marked and removed and the files were re-indexed:
```bash
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/WGS/Archana_Feb2021/*/trim_galore/bwa-mem/*_trimmed_1_sorted.bam); do
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_sort_bam.sh $file
done 
#57500981-88

source package 3e7beb4d-f08b-4d6b-9b6a-f99cc91a38f9
source package 638df626-d658-40aa-80e5-14a275b7464b
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/WGS/Archana_Feb2021/*/trim_galore/bwa-mem/*_trimmed_1_sorted.bam); do
    OutDir=$(dirname $file)
    OutFile=$(basename $file | sed s'@_sorted.bam@_sorted_MarkDups.bam@g')
    metrics=$(basename $file | sed s'@_trimmed_1_sorted.bam@@g')_marked_dup_metrics.txt
    java17 -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar SortSam I=${file} O=${OutDir}/temp.bam SORT_ORDER=coordinate
    java17 -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${OutDir}/temp.bam O=${OutDir}/${OutFile} M=${OutDir}/${metrics}
    cd ${OutDir}
    rm temp.bam
    samtools index ${OutFile}
    cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae
done
#57510702

#Input files were subsequently deleted to save space
```

Reads near detected indels realigned to remove alignment artifacts:
```bash
#The reference genome was indexed and a dictionary created:
interactive
source switch-institute ei
source package 638df626-d658-40aa-80e5-14a275b7464b
source pilon-1.22
source package /tgac/software/testing/bin/picardtools-2.1.1
cd /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/
java -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar CreateSequenceDictionary R=Myzus_persicae_O_v2.0.scaffolds.fa O=Myzus_persicae_O_v2.0.scaffolds.dict
samtools faidx Myzus_persicae_O_v2.0.scaffolds.fa
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae

#GATK indel realignment:
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/WGS/Archana_Feb2021/*/trim_galore/bwa-mem/*MarkDups.bam); do
Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_realign.sh $file $Reference 
done #57519857-64
```

### Variant calling <a name="9"></a>

Combined into a .vcf and called variants with bcftools:
```bash
source package 638df626-d658-40aa-80e5-14a275b7464b
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/WGS/Archana_Feb2021/*/trim_galore/bwa-mem/gatk/*realigned.bam > bamlist.txt
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/
bcftools mpileup -b /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/bamlist.txt --annotate AD,DP --fasta-ref /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa -O z -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/BR_AT_NB_hostswaps.vcf.gz

bcftools call --ploidy 2 -Oz -v -m -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/BR_AT_NB_hostswap.vcf.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/BR_AT_NB_hostswaps.vcf.gz
#57520372
```
### Filter  <a name="10"></a>
#### Genmap filter for genome mappability <a name="11"></a>

We further refined our input files by calculating mappability of the genome with GenMap (v1.3.0) with the parameters -K 100 -E 2. This estimates k-mer uniqueness and identifies regions of the genome where Illumina reads are unable to map uniquely. We masked all regions larger than 100 bp with less than 1 i.e., max mappability. 10.1038/s41586-021-04269-6 and 10.1038/s41467-023-43383-z use k=100 for 150bp paired reads.
```bash
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/BR_AT_NB_hostswap.vcf.gz | wc -l #977,192

VCF=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/BR_AT_NB_hostswap.vcf.gz
OutDir=$(dirname $VCF)/genmap
OutFile=$(basename $VCF | sed 's@.vcf.gz@@g')
Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
Repeatmodeller=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/repeatmodeller/Myzus_persicae_O_v2.0.scaffolds.fa.out
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_genmap_mappability_masking.sh $OutDir $OutFile $Reference $Repeatmodeller $VCF
#57695712

zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable.vcf.gz | wc -l #644,700
```
Mappable bases (genmap): 328,601,147
Masked bases (repeatmodeler): 93,097,267
Callable bases (overlap of VCF and mappable genmap regions, minus repeatmodeler masked regions): 2,417,154

Plot genemap mappability to compare with qualimap mappability plots. These are not clear as the mappability regions are too dense to properly distinguish at whole genome scale, nonetheless it appears that genmap successfully identifies poor mappability at the start of Chromosome 1.
```python
import matplotlib.pyplot as plt

def plot_bedgraph(bedgraph_file, chromosome):
    starts, values = [], []
    with open(bedgraph_file, 'r') as f:
        for line in f:
            if line.startswith('track') or line.startswith('browser'):
                continue
            fields = line.strip().split('\t')
            if fields[0] == chromosome:
                starts.append(int(fields[1]))
                values.append(float(fields[3]))
    plt.plot(starts, values, label=f'{chromosome}')
    plt.xlabel('Position')
    plt.ylabel('Mappability')
    plt.title(f'Genmap Mappability Plot - {chromosome}')
    plt.legend()

bedgraph_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_genmap.bedgraph'
chromosomes_to_plot = ['scaffold_1', 'scaffold_2', 'scaffold_3', 'scaffold_4', 'scaffold_5', 'scaffold_6']
plt.figure(figsize=(96, 48))

for idx, chromosome in enumerate(chromosomes_to_plot, 1):
    plt.subplot(6, 1, idx)
    plot_bedgraph(bedgraph_file, chromosome)

plt.tight_layout()

plt.savefig('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/genemap_mappability_plot.png')
```
#### Filter samples for missingness: <a name="12"></a>
```bash
source package /nbi/software/testing/bin/vcftools-0.1.15
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable.vcf.gz --missing-indv
mv out.imiss /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable.out.imiss                                    
```
Plot the missingness level of samples:
```python
import matplotlib.pyplot as plt

# Read the F_MISS data from a file, skipping the header line
data = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable.out.imiss', 'r') as file:
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
plt.savefig('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable.out.imiss_plot.png')
```
Missingness for all eight samples is very low - they are all genuine M.persicae samples (as expected).

#### Filter for SNP quality: <a name="13"></a>

SNPs were filtered to keep only bi-allelic SNPs, with minimum depth of 5, minimum quality score of 30, and maximum missingness of SNPs of 10% (with 8 samples this means that SNP positions must be present in all samples). 
```bash
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable_filtered
bgzip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable_filtered.recode.vcf
#After filtering, kept 359,245 out of a possible 644,310 Sites
```
Split the two experiments:
```bash
bcftools view -s AT39_E1_trimmed_1_sorted_MarkDups,BR39_E1_trimmed_1_sorted_MarkDups,NB39_E1_trimmed_1_sorted_MarkDups,BR1_E1_trimmed_1_sorted_MarkDups -Oz -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable_filtered_expt1.vcf.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable_filtered.recode.vcf.gz
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable_filtered_expt1.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable_filtered_expt1
bgzip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable_filtered_expt1.recode.vcf
#After filtering, kept 358,369 out of a possible 359,245 Sites - 876 removed

bcftools view -s NB25_E2_trimmed_1_sorted_MarkDups,BR25_E2_trimmed_1_sorted_MarkDups,AT25_E2_trimmed_1_sorted_MarkDups,BR1_E2_trimmed_1_sorted_MarkDups -Oz -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable_filtered_expt2.vcf.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable_filtered.recode.vcf.gz
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable_filtered_expt2.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable_filtered_expt2
bgzip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable_filtered_expt2.recode.vcf
#After filtering, kept 358,377 out of a possible 359,245 Sites - 868 removed
```
The 1,744 SNPs that are specific to one experiment are unlikely to be responsible for adaptation to a new host.

**Question:** how many SNPs are there between BR1_E1 and BR39_E1, and BR1_E2 and BR25_E2? This is the number of SNPs without a host swap.
```bash
Directory=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap
zcat ${Directory}/BR_AT_NB_hostswap_callable_filtered_expt1.recode.vcf.gz > temp.vcf
perl ~/git_repos/Scripts/NBI/similarity-matrix.pl -c off -f temp.vcf -o matrix.txt
#sample_name     AT39_E1_trimmed_1_sorted_MarkDups       BR1_E1_trimmed_1_sorted_MarkDups        BR39_E1_trimmed_1_sorted_MarkDups       NB39_E1_trimmed_1_sorted_MarkDups
#AT39_E1_trimmed_1_sorted_MarkDups       358369/358369   344102/358369   344670/358369   343692/358369
#BR1_E1_trimmed_1_sorted_MarkDups        344102/358369   358369/358369   345270/358369   344237/358369
#BR39_E1_trimmed_1_sorted_MarkDups       344670/358369   345270/358369   358369/358369   345223/358369
#NB39_E1_trimmed_1_sorted_MarkDups       343692/358369   344237/358369   345223/358369   358369/358369

zcat ${Directory}/BR_AT_NB_hostswap_callable_filtered_expt2.recode.vcf.gz > temp.vcf
perl ~/git_repos/Scripts/NBI/similarity-matrix.pl -c off -f temp.vcf -o matrix.txt
#sample_name     AT25_E2_trimmed_1_sorted_MarkDups       BR1_E2_trimmed_1_sorted_MarkDups        BR25_E2_trimmed_1_sorted_MarkDups       NB25_E2_trimmed_1_sorted_MarkDups
#AT25_E2_trimmed_1_sorted_MarkDups       358377/358377   344954/358377   346106/358377   346064/358377
#BR1_E2_trimmed_1_sorted_MarkDups        344954/358377   358377/358377   346319/358377   346122/358377
#BR25_E2_trimmed_1_sorted_MarkDups       346106/358377   346319/358377   358377/358377   347062/358377
#NB25_E2_trimmed_1_sorted_MarkDups       346064/358377   346122/358377   347062/358377   358377/358377
```
Most SNPs are identical across all four samples from an experiment, and presumably reflect heterozygous positions in clone O called against the unphased reference genome or SNP mutations between the sequencing of the reference clone O genome and the start of the host swap experiment.

At the end of experiment 1 there were 13,099 SNPs between samples kept on Brassica rapa versus samples from the start of the experiment, compared to 14,267 and 14,132 in samples switched to Arabidopsis and N.benthamiana resectively.

At the end of experiment 2 there were 12,058 SNPs between samples kept on Brassica rapa versus samples from the start of the experiment, compared to 13,423 and 12,255 in samples switched to Arabidopsis and N.benthamiana resectively.

**Question:** did any SNPs occur in NB or AT host swaps in both experiments, but not in BR samples? If so then these would be candidates for causing host adaptation and could be investigated further - probably overkill to check this though.

```bash
Directory=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap
bcftools view -s BR1_E1_trimmed_1_sorted_MarkDups -Oz -o ${Directory}/BR1_E1.vcf.gz ${Directory}/BR_AT_NB_hostswap_callable_filtered_expt1.recode.vcf.gz
vcftools --gzvcf ${Directory}/BR1_E1.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out ${Directory}/BR1_E1
bgzip ${Directory}/BR1_E1.recode.vcf && mv ${Directory}/BR1_E1.recode.vcf.gz ${Directory}/BR1_E1.vcf.gz

bcftools view -s BR39_E1_trimmed_1_sorted_MarkDups -Oz -o ${Directory}/BR39_E1.vcf.gz ${Directory}/BR_AT_NB_hostswap_callable_filtered_expt1.recode.vcf.gz
vcftools --gzvcf ${Directory}/BR39_E1.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out ${Directory}/BR39_E1
bgzip ${Directory}/BR39_E1.recode.vcf && mv ${Directory}/BR39_E1.recode.vcf.gz ${Directory}/BR39_E1.vcf.gz

bcftools view -s AT39_E1_trimmed_1_sorted_MarkDups -Oz -o ${Directory}/AT39_E1.vcf.gz ${Directory}/BR_AT_NB_hostswap_callable_filtered_expt1.recode.vcf.gz
vcftools --gzvcf ${Directory}/AT39_E1.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out ${Directory}/AT39_E1
bgzip ${Directory}/AT39_E1.recode.vcf && mv ${Directory}/AT39_E1.recode.vcf.gz ${Directory}/AT39_E1.vcf.gz

bcftools view -s NB39_E1_trimmed_1_sorted_MarkDups -Oz -o ${Directory}/NB39_E1.vcf.gz ${Directory}/BR_AT_NB_hostswap_callable_filtered_expt1.recode.vcf.gz
vcftools --gzvcf ${Directory}/NB39_E1.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out ${Directory}/NB39_E1
bgzip ${Directory}/NB39_E1.recode.vcf && mv ${Directory}/NB39_E1.recode.vcf.gz ${Directory}/NB39_E1.vcf.gz

bcftools index ${Directory}/BR1_E1.vcf.gz -t -o ${Directory}/BR1_E1.vcf.gz.tbi
bcftools index ${Directory}/BR39_E1.vcf.gz -t -o ${Directory}/BR39_E1.vcf.gz.tbi
bcftools index ${Directory}/AT39_E1.vcf.gz -t -o ${Directory}/AT39_E1.vcf.gz.tbi
bcftools index ${Directory}/NB39_E1.vcf.gz -t -o ${Directory}/NB39_E1.vcf.gz.tbi

bcftools view -s BR1_E2_trimmed_1_sorted_MarkDups -Oz -o ${Directory}/BR1_E2.vcf.gz ${Directory}/BR_AT_NB_hostswap_callable_filtered_expt2.recode.vcf.gz
vcftools --gzvcf ${Directory}/BR1_E2.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out ${Directory}/BR1_E2
bgzip ${Directory}/BR1_E2.recode.vcf && mv ${Directory}/BR1_E2.recode.vcf.gz ${Directory}/BR1_E2.vcf.gz

bcftools view -s BR25_E2_trimmed_1_sorted_MarkDups -Oz -o ${Directory}/BR25_E2.vcf.gz ${Directory}/BR_AT_NB_hostswap_callable_filtered_expt2.recode.vcf.gz
vcftools --gzvcf ${Directory}/BR25_E2.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out ${Directory}/BR25_E2
bgzip ${Directory}/BR25_E2.recode.vcf && mv ${Directory}/BR25_E2.recode.vcf.gz ${Directory}/BR25_E2.vcf.gz

bcftools view -s AT25_E2_trimmed_1_sorted_MarkDups -Oz -o ${Directory}/AT25_E2.vcf.gz ${Directory}/BR_AT_NB_hostswap_callable_filtered_expt2.recode.vcf.gz
vcftools --gzvcf ${Directory}/AT25_E2.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out ${Directory}/AT25_E2
bgzip ${Directory}/AT25_E2.recode.vcf && mv ${Directory}/AT25_E2.recode.vcf.gz ${Directory}/AT25_E2.vcf.gz

bcftools view -s NB25_E2_trimmed_1_sorted_MarkDups -Oz -o ${Directory}/NB25_E2.vcf.gz ${Directory}/BR_AT_NB_hostswap_callable_filtered_expt2.recode.vcf.gz
vcftools --gzvcf ${Directory}/NB25_E2.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out ${Directory}/NB25_E2
bgzip ${Directory}/NB25_E2.recode.vcf && mv ${Directory}/NB25_E2.recode.vcf.gz ${Directory}/NB25_E2.vcf.gz

bcftools index ${Directory}/BR1_E2.vcf.gz -t -o ${Directory}/BR1_E2.vcf.gz.tbi
bcftools index ${Directory}/BR25_E2.vcf.gz -t -o ${Directory}/BR25_E2.vcf.gz.tbi
bcftools index ${Directory}/AT25_E2.vcf.gz -t -o ${Directory}/AT25_E2.vcf.gz.tbi
bcftools index ${Directory}/NB25_E2.vcf.gz -t -o ${Directory}/NB25_E2.vcf.gz.tbi

source package 4028d6e4-21a8-45ec-8545-90e4ed7e1a64
bedtools subtract -a ${Directory}/AT39_E1.vcf.gz -b ${Directory}/BR1_E1.vcf.gz ${Directory}/BR39_E1.vcf.gz > ${Directory}/AT39_E1_swaponly.vcf #1,637
bedtools subtract -a ${Directory}/NB39_E1.vcf.gz -b ${Directory}/BR1_E1.vcf.gz ${Directory}/BR39_E1.vcf.gz > ${Directory}/NB39_E1_swaponly.vcf #1,637
awk '$0 !~ /^#/ {print $1"\t"$2"\t"$2}' ${Directory}/AT39_E1_swaponly.vcf > ${Directory}/AT39_E1_swaponly.bed
awk '$0 !~ /^#/ {print $1"\t"$2"\t"$2}' ${Directory}/NB39_E1_swaponly.vcf > ${Directory}/NB39_E1_swaponly.bed

bedtools subtract -a ${Directory}/AT25_E2.vcf.gz -b ${Directory}/BR1_E2.vcf.gz ${Directory}/BR25_E2.vcf.gz > ${Directory}/AT25_E2_swaponly.vcf #1,935
bedtools subtract -a ${Directory}/NB25_E2.vcf.gz -b ${Directory}/BR1_E2.vcf.gz ${Directory}/BR25_E2.vcf.gz > ${Directory}/NB25_E2_swaponly.vcf #2,065
awk '$0 !~ /^#/ {print $1"\t"$2"\t"$2}' ${Directory}/AT25_E2_swaponly.vcf > ${Directory}/AT25_E2_swaponly.bed
awk '$0 !~ /^#/ {print $1"\t"$2"\t"$2}' ${Directory}/NB25_E2_swaponly.vcf > ${Directory}/NB25_E2_swaponly.bed

bedtools intersect -a ${Directory}/AT39_E1_swaponly.bed -b ${Directory}/AT25_E2_swaponly.bed > ${Directory}/AT25_swaponly.bed #186
bedtools intersect -a ${Directory}/NB39_E1_swaponly.bed -b ${Directory}/NB25_E2_swaponly.bed > ${Directory}/NB25_swaponly.bed #239
```
186 SNPs are found in both experiments from samples swapped onto Arabidopsis plants but not in samples kept on the same host, these could be investigated further to determine if they fall within coding regions and could be resposible for host adaptation, or are likely erroneous but missed by previous filters.

239 SNPs are found in both experiments from samples swapped onto N.benthamiana plants but not in samples kept on the same host, these could be investigated further to determine if they fall within coding regions and could be resposible for host adaptation, or are likely erroneous but missed by previous filters.

## Check for epigenetic changes over the course of the experiment with WGBS <a name="14"></a>
#### fastqc and qualimap <a name="15"></a>
Bisulfite treatment followed by PCR amplification specifically converts unmethylated cytosines to thymine, standard illumina sequencing is then performed. Fastqc can be used for QC, a specialist aligner is needed to take account of C2T conversions when aligning to the genome.

Raw reads were assessed for quality and coverage of the clone O_v2 reference assembly:
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGBS/Archana_Mar2021/*/); do
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    Fread2=$(ls ${ReadDir}*_3.fq.gz)
    Rread2=$(ls ${ReadDir}*_4.fq.gz)
    Fread3=$(ls ${ReadDir}*_5.fq.gz)
    Rread3=$(ls ${ReadDir}*_6.fq.gz)
    Fread4=$(ls ${ReadDir}*_7.fq.gz)
    Rread4=$(ls ${ReadDir}*_8.fq.gz)
    Fread5=$(ls ${ReadDir}*_9.fq.gz)
    Rread5=$(ls ${ReadDir}*_10.fq.gz)
    OutDir=$(echo ${ReadDir})
    Reference_genome=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
    Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread $Fread2 $Rread2 $Fread3 $Rread3 $Fread4 $Rread4 $Fread5 $Rread5
done
#57551367-57551414
```
Raw read folder has been compressed to save space.
### Trimming <a name="16"></a>
#### Trim galore <a name="17"></a>

Adapters and low quality regions were trimmed from raw reads via trim-galore:
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGBS/Archana_Mar2021/*/); do
    sample=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    Fread2=$(ls ${ReadDir}*_3.fq.gz)
    Rread2=$(ls ${ReadDir}*_4.fq.gz)
    Fread3=$(ls ${ReadDir}*_5.fq.gz)
    Rread3=$(ls ${ReadDir}*_6.fq.gz)
    Fread4=$(ls ${ReadDir}*_7.fq.gz)
    Rread4=$(ls ${ReadDir}*_8.fq.gz)
    Fread5=$(ls ${ReadDir}*_9.fq.gz)
    Rread5=$(ls ${ReadDir}*_10.fq.gz)
    OutDir=$(echo $ReadDir | sed 's@raw_data@dna_qc@g')trim_galore
    OutFile=${sample}_trimmed
    Quality=20
    Length=50
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir -p $OutDir
    sbatch $ProgDir/run_trim_galore.sh $OutDir $OutFile $Quality $Length $Fread $Rread $Fread2 $Rread2 $Fread3 $Rread3 $Fread4 $Rread4 $Fread5 $Rread5
done 
#57551302-57551349
```
Trimmed reads were re-assessed for quality, the C2T conversion aware aligner bsmap was used to allow assessment of coverage of the clone O_v2 reference assembly, although the script does not keep these .bam files as written:

```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/WGBS/Archana_Mar2021/*/trim_galore/); do
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    OutDir=$(echo ${ReadDir})
    Reference_genome=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
    Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch $ProgDir/run_raw_bs_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread 
done
#57795800-57795847

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/WGS/Archana_Feb2021/*/trim_galore); do
sample=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
coverage=$(grep 'mean coverageData' ${ReadDir}qualimap/*_genome_results_gff.txt | rev | cut -d ' ' -f1 | rev | sed 's@X@@g')
echo $sample raw reads have average coverage of ${coverage}
done
```
### BSmooth <a name="18"></a>
Common approaches for differential methylation analysis are Bsmooth, Methylkit or a custom approach to define differentially methylated regions (DMRs) DOI:10.1093/bib/bbx077. Here bsmooth:

#### Bismark alignment <a name="19"></a>

The best C2T aware aligner seems to be bsmap based upon the literature (eg.: DOI:10.1016/j.csbj.2022.08.051), however the bsmooth package does not have pre-built support for a bsmap input but does for bismark. Bismark alignment, assessment and 'genome wide cytosine report' were prepared, .bam files from this alignment were not saved as we have bsmap generated files already:
```bash
#Create bismark index
source package 33c48798-0827-4add-8153-909c1bd83e89
source package 29a74b59-88fc-4453-a30b-1310b34910b9
mkdir /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/bismark
ln -s /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/bismark/.
bismark_genome_preparation --bowtie2 --verbose /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/bismark

#Bismark Alignment
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/WGBS/Archana_Mar2021/*/trim_galore/); do
    Reference_dir=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/bismark
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    OutDir=$(echo ${ReadDir} | sed 's@trim_galore@bismark@g' | sed 's@dna_qc@alignment@g')
    OutFile=$(basename $Fread | sed 's@_trimmed_1.fq.gz@bismark@g')
    ProgDir=~/git_repos/Wrappers/NBI
    #if [ ! -d "$OutDir" ]; then
    sbatch $ProgDir/run_bismark0.24.1.sh $OutDir $OutFile $Reference_dir $Fread $Rread 
    #fi
done
#57768086-133
```
#### BSsmooth <a name="20"></a>

BSmooth is part of the BSseq R package, the package need to load the files into working memory, for our files this is more memory than a local machine has available, therefore R will have to be run on the HPC (which is a pain :().
```bash
#Make directory for output files:
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/bsseq1.38.0.sif R
#57758169
```
Differential methylation between week 1 samples:
```R
library(bsseq)
library(stats)
library(BiocParallel)

#Read in data:
col_names <- c("treatment", "replicate", "col")
row_names <- c("BR1_E2_1", "BR1_E2_2", "BR1_E2_3", "NB1_E2_1", "NB1_E2_2", "NB1_E2_3", "AT1_E2_1", "AT1_E2_2", "AT1_E2_3")
data <- matrix(c("control", "control", "control", "NB", "NB", "NB", "AT", "AT", "AT", 1, 2, 3, 1, 2, 3, 1, 2, 3, "blue", "blue", "blue", "red", "red", "red", "green", "green", "green"), nrow = length(row_names), ncol = length(col_names))
df <- data.frame(data, row.names = row_names)
colnames(df) <- col_names
print("Input dataframe:")
print(df)

bsseq <- read.bismark(files = c("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_1/bismark/BR1_E2_1_bismark.deduplicated.CpG_report.txt",
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_2/bismark/BR1_E2_2_bismark.deduplicated.CpG_report.txt",
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_3/bismark/BR1_E2_3_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_1/bismark/NB1_E2_1_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_2/bismark/NB1_E2_2_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_3/bismark/NB1_E2_3_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT1_E2_1/bismark/AT1_E2_1_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT1_E2_2/bismark/AT1_E2_2_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT1_E2_3/bismark/AT1_E2_3_bismark.deduplicated.CpG_report.txt"),
    colData = df,
    rmZeroCov = FALSE,
    strandCollapse = TRUE,
    verbose = TRUE)
save(bsseq, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/1_bsseq.RData")

sapply(assays(bsseq, withDimnames = FALSE), class)
print("Bsseq object:")
bsseq
pData(bsseq)

print("The average coverage of CpGs:")
round(colMeans(getCoverage(bsseq)), 1) 
#BR1_E2_1 BR1_E2_2 BR1_E2_3 NB1_E2_1 NB1_E2_2 NB1_E2_3 AT1_E2_1 AT1_E2_2 AT1_E2_3
#    19.3     26.0     22.2     22.0     23.5     24.7     23.0     33.6    20.5
print("The number of CpGs:")
length(bsseq) 
#10,962,492
print("Number of CpGs which are covered by at least 1 read in all samples:")
sum(rowSums(getCoverage(bsseq) >= 1) == 7)
#65,854
print("Number of CpGs with 0 coverage in all samples:")
sum(rowSums(getCoverage(bsseq)) == 0)
#265,410

#Perform smoothing:
#"ns is the minimum number of CpGs contained in each window, h is half the minimum window with (the actual window width is either 2 times h or wide enough to contain ns covered CpGs, whichever is greater). Note that the window width is different at each position in the genome and may also be different for different samples at the same position, since it depends on how many nearby CpGs with non-zero coverage. Per default, a smoothing cluster is a whole chromosome. By “cluster” we mean a set of CpGs which are processed together. This means that even if there is a large distance between two CpGs, we borrow strength between them. By setting maxGap this can be prevented since the argument describes the longest distance between two CpGs before a cluster is broken up into two clusters." - all default:

bssmooth <- BSmooth(
    BSseq = bsseq, 
    ns = 70,
    h = 1000,
    maxGap = 10^8,
    BPPARAM = MulticoreParam(workers = 1), 
    verbose = TRUE)
bssmooth
save(bssmooth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/1_bsmooth.RData")

#Remove CpGs with coverage below 4 in all samples:
BS.cov <- getCoverage(bssmooth)
keepLoci.ex <- which(rowSums(BS.cov[, bsseq$treatment == "control"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "NB"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "AT"] >= 4) >= 3)
print("The number of CpGs with coverage >=4 in all samples:")
length(keepLoci.ex)
#9,711,258
bssmooth <- bssmooth[keepLoci.ex,]


#Compute t-statistics based on smoothed whole-genome bisulfite sequencing data:
print("Compute t-stats vs AT:")
AT.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("AT1_E2_1", "AT1_E2_2", "AT1_E2_3"),  
                                    group2 = c("BR1_E2_1", "BR1_E2_2", "BR1_E2_3"),
                                    estimate.var = "same",
                                    local.correct = TRUE,
                                    verbose = TRUE)
AT.tstat
#9,711,258 methylation loci
temp <- AT.tstat@stats
temp2 <- data.frame(temp)
temp2 <- na.omit(temp2) #9,710,695 methylation loci remain
quantile(temp2$tstat, c(0.025, 0.975))
#     2.5%     97.5% 
#-1.575393  1.690201
print("Compute t-stats vs NB:")
NB.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("NB1_E2_1", "NB1_E2_2", "NB1_E2_3"),  
                                    group2 = c("BR1_E2_1", "BR1_E2_2", "BR1_E2_3"),
                                    estimate.var = "same",
                                    local.correct = TRUE,
                                    verbose = TRUE)
NB.tstat
#9,711,258 methylation loci
temp <- NB.tstat@stats
temp2 <- data.frame(temp)
temp2 <- na.omit(temp2) #9,710,588 methylation loci remain
quantile(temp2$tstat, c(0.025, 0.975))
#     2.5%     97.5% 
#-1.597031  1.661042


#Finding DMRs
print("Find DMRs vs AT:")
dmrs0 <- dmrFinder(AT.tstat, cutoff = c(-1.575393, 1.690201), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation of at least 0.1.
AT.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(AT.dmrs) #116
head(AT.dmrs, n = 3)
write.table(AT.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT1_dmrs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

print("Find DMRs vs NB:")
dmrs0 <- dmrFinder(NB.tstat, cutoff = c(-1.597031, 1.661042), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylationof at least 0.1.
NB.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(NB.dmrs) #109
head(NB.dmrs, n = 3)
write.table(NB.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB1_dmrs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Plot the top DMRs
#blue=BR,red=NB,green=AT
print("Plot top AT DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT1_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, AT.dmrs[1:116,], extend = 5000, 
                addRegions = AT.dmrs)
dev.off()
print("Plot top NB DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB1_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, NB.dmrs[1:109,], extend = 5000, 
                addRegions = NB.dmrs)
dev.off()
```
Differential methylation between week 3 samples:
```R
library(bsseq)
library(stats)
library(BiocParallel)

#Read in data:
col_names <- c("treatment", "replicate", "col")
row_names <- c("BR3_E2_1", "BR3_E2_2", "BR3_E2_3", "NB3_E2_1", "NB3_E2_2", "NB3_E2_3", "AT3_E2_1", "AT3_E2_2", "AT3_E2_3")
data <- matrix(c("control", "control", "control", "NB", "NB", "NB", "AT", "AT", "AT", 1, 2, 3, 1, 2, 3, 1, 2, 3, "blue", "blue", "blue", "red", "red", "red", "green", "green", "green"), nrow = length(row_names), ncol = length(col_names))
df <- data.frame(data, row.names = row_names)
colnames(df) <- col_names
print("Input dataframe:")
print(df)

bsseq <- read.bismark(files = c("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_1/bismark/BR3_E2_1_bismark.deduplicated.CpG_report.txt",
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_2/bismark/BR3_E2_2_bismark.deduplicated.CpG_report.txt",
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_3/bismark/BR3_E2_3_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB3_E2_1/bismark/NB3_E2_1_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB3_E2_2/bismark/NB3_E2_2_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB3_E2_3/bismark/NB3_E2_3_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT3_E2_1/bismark/AT3_E2_1_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT3_E2_2/bismark/AT3_E2_2_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT3_E2_3/bismark/AT3_E2_3_bismark.deduplicated.CpG_report.txt"),
    colData = df,
    rmZeroCov = FALSE,
    strandCollapse = TRUE,
    verbose = TRUE)
save(bsseq, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/3_bsseq.RData")

sapply(assays(bsseq, withDimnames = FALSE), class)
print("Bsseq object:")
bsseq
pData(bsseq)

print("The average coverage of CpGs:")
round(colMeans(getCoverage(bsseq)), 1)
#BR3_E2_1 BR3_E2_2 BR3_E2_3 NB3_E2_1 NB3_E2_2 NB3_E2_3 AT3_E2_1 AT3_E2_2 AT3_E2_3
#    30.3     26.1     24.2     28.7     23.0     23.1     24.9     24.9    22.3
print("The number of CpGs:")
length(bsseq)
#10,962,492
print("Number of CpGs which are covered by at least 1 read in all samples:")
sum(rowSums(getCoverage(bsseq) >= 1) == 7)
#63,718
print("Number of CpGs with 0 coverage in all samples:")
sum(rowSums(getCoverage(bsseq)) == 0)
#271,066


#Perform smoothing:
#"ns is the minimum number of CpGs contained in each window, h is half the minimum window with (the actual window width is either 2 times h or wide enough to contain ns covered CpGs, whichever is greater). Note that the window width is different at each position in the genome and may also be different for different samples at the same position, since it depends on how many nearby CpGs with non-zero coverage. Per default, a smoothing cluster is a whole chromosome. By “cluster” we mean a set of CpGs which are processed together. This means that even if there is a large distance between two CpGs, we borrow strength between them. By setting maxGap this can be prevented since the argument describes the longest distance between two CpGs before a cluster is broken up into two clusters." - all default:

bssmooth <- BSmooth(
    BSseq = bsseq, 
    ns = 70,
    h = 1000,
    maxGap = 10^8,
    BPPARAM = MulticoreParam(workers = 1), 
    verbose = TRUE)
bssmooth
save(bssmooth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/3_bsmooth.RData")

#Remove CpGs with coverage below 4 in all samples:
BS.cov <- getCoverage(bssmooth)
keepLoci.ex <- which(rowSums(BS.cov[, bsseq$treatment == "control"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "NB"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "AT"] >= 4) >= 3)
print("The number of CpGs with coverage >=4 in all samples:")
length(keepLoci.ex)
#9,764,358
bssmooth <- bssmooth[keepLoci.ex,]


#Compute t-statistics based on smoothed whole-genome bisulfite sequencing data:
print("Compute t-stats vs AT:")
AT.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("AT3_E2_1", "AT3_E2_2", "AT3_E2_3"),  
                                    group2 = c("BR3_E2_1", "BR3_E2_2", "BR3_E2_3"),
                                    estimate.var = "same",
                                    local.correct = TRUE,
                                    verbose = TRUE)
AT.tstat
#9,764,358 methylation loci
temp <- AT.tstat@stats
temp2 <- data.frame(temp)
temp2 <- na.omit(temp2) #9,763,691 methylation loci remain
quantile(temp2$tstat, c(0.025, 0.975))
#     2.5%     97.5% 
#-1.496376  1.784706
print("Compute t-stats vs NB:")
NB.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("NB3_E2_1", "NB3_E2_2", "NB3_E2_3"),  
                                    group2 = c("BR3_E2_1", "BR3_E2_2", "BR3_E2_3"),
                                    estimate.var = "same",
                                    local.correct = TRUE,
                                    verbose = TRUE)
NB.tstat
#9,764,358 methylation loci
temp <- NB.tstat@stats
temp2 <- data.frame(temp)
temp2 <- na.omit(temp2) #9,763,600 methylation loci remain
quantile(temp2$tstat, c(0.025, 0.975))
#     2.5%     97.5% 
#-1.529655  1.735956


#Finding DMRs
print("Find DMRs vs AT:")
dmrs0 <- dmrFinder(AT.tstat, cutoff = c(-1.496376, 1.784706), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation of at least 0.1.
AT.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(AT.dmrs) #130
head(AT.dmrs, n = 3)
write.table(AT.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT3_dmrs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

print("Find DMRs vs NB:")
dmrs0 <- dmrFinder(NB.tstat, cutoff = c(-1.529655, 1.735956), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylationof at least 0.1.
NB.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(NB.dmrs) #179
head(NB.dmrs, n = 3)
write.table(NB.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB3_dmrs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Plot the top DMRs
#blue=BR,red=NB,green=AT
print("Plot top AT DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT3_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, AT.dmrs[1:130,], extend = 5000, 
                addRegions = AT.dmrs)
dev.off()
print("Plot top NB DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB3_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, NB.dmrs[1:179,], extend = 5000, 
                addRegions = NB.dmrs)
dev.off()
```
Differential methylation between week 6 samples:
```R
library(bsseq)
library(stats)
library(BiocParallel)

#Read in data:
col_names <- c("treatment", "replicate", "col")
row_names <- c("BR6_E2_1", "BR6_E2_2", "BR6_E2_3", "NB6_E2_1", "NB6_E2_2", "NB6_E2_3", "AT6_E2_1", "AT6_E2_2", "AT6_E2_3")
data <- matrix(c("control", "control", "control", "NB", "NB", "NB", "AT", "AT", "AT", 1, 2, 3, 1, 2, 3, 1, 2, 3, "blue", "blue", "blue", "red", "red", "red", "green", "green", "green"), nrow = length(row_names), ncol = length(col_names))
df <- data.frame(data, row.names = row_names)
colnames(df) <- col_names
print("Input dataframe:")
print(df)

bsseq <- read.bismark(files = c("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_1/bismark/BR6_E2_1_bismark.deduplicated.CpG_report.txt",
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_2/bismark/BR6_E2_2_bismark.deduplicated.CpG_report.txt",
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_3/bismark/BR6_E2_3_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_1/bismark/NB6_E2_1_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_2/bismark/NB6_E2_2_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_3/bismark/NB6_E2_3_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_1/bismark/AT6_E2_1_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_2/bismark/AT6_E2_2_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_3/bismark/AT6_E2_3_bismark.deduplicated.CpG_report.txt"),
    colData = df,
    rmZeroCov = FALSE,
    strandCollapse = TRUE,
    verbose = TRUE)
save(bsseq, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/6_bsseq.RData")

sapply(assays(bsseq, withDimnames = FALSE), class)
print("Bsseq object:")
bsseq
pData(bsseq)

print("The average coverage of CpGs:")
round(colMeans(getCoverage(bsseq)), 1)
#BR6_E2_1 BR6_E2_2 BR6_E2_3 NB6_E2_1 NB6_E2_2 NB6_E2_3 AT6_E2_1 AT6_E2_2 AT6_E2_3
#    20.6     18.7     20.3     23.2     19.6     33.5     21.9     23.9    22.9
print("The number of CpGs:")
length(bsseq)
#10,962,492
print("Number of CpGs which are covered by at least 1 read in all samples:")
sum(rowSums(getCoverage(bsseq) >= 1) == 7)
#67,436
print("Number of CpGs with 0 coverage in all samples:")
sum(rowSums(getCoverage(bsseq)) == 0)
#261,473


#Perform smoothing:
#"ns is the minimum number of CpGs contained in each window, h is half the minimum window with (the actual window width is either 2 times h or wide enough to contain ns covered CpGs, whichever is greater). Note that the window width is different at each position in the genome and may also be different for different samples at the same position, since it depends on how many nearby CpGs with non-zero coverage. Per default, a smoothing cluster is a whole chromosome. By “cluster” we mean a set of CpGs which are processed together. This means that even if there is a large distance between two CpGs, we borrow strength between them. By setting maxGap this can be prevented since the argument describes the longest distance between two CpGs before a cluster is broken up into two clusters." - all default:

bssmooth <- BSmooth(
    BSseq = bsseq, 
    ns = 70,
    h = 1000,
    maxGap = 10^8,
    BPPARAM = MulticoreParam(workers = 1), 
    verbose = TRUE)
bssmooth
save(bssmooth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/6_bsmooth.RData")

#Remove CpGs with coverage below 4 in all samples:
BS.cov <- getCoverage(bssmooth)
keepLoci.ex <- which(rowSums(BS.cov[, bsseq$treatment == "control"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "NB"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "AT"] >= 4) >= 3)
print("The number of CpGs with coverage >=4 in all samples:")
length(keepLoci.ex)
#9,673,446
bssmooth <- bssmooth[keepLoci.ex,]


#Compute t-statistics based on smoothed whole-genome bisulfite sequencing data:
print("Compute t-stats vs AT:")
AT.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("AT6_E2_1", "AT6_E2_2", "AT6_E2_3"),  
                                    group2 = c("BR6_E2_1", "BR6_E2_2", "BR6_E2_3"),
                                    estimate.var = "same",
                                    local.correct = TRUE,
                                    verbose = TRUE)
AT.tstat
#9,673,446 methylation loci
temp <- AT.tstat@stats
temp2 <- data.frame(temp)
temp2 <- na.omit(temp2) #9,672,758 methylation loci remain
quantile(temp2$tstat, c(0.025, 0.975))
#     2.5%     97.5% 
#-1.821612  1.508577
print("Compute t-stats vs NB:")
NB.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("NB6_E2_1", "NB6_E2_2", "NB6_E2_3"),  
                                    group2 = c("BR6_E2_1", "BR6_E2_2", "BR6_E2_3"),
                                    estimate.var = "same",
                                    local.correct = TRUE,
                                    verbose = TRUE)
NB.tstat
#9,673,446 methylation loci
temp <- NB.tstat@stats
temp2 <- data.frame(temp)
temp2 <- na.omit(temp2) #9,672,792 methylation loci remain
quantile(temp2$tstat, c(0.025, 0.975))
#     2.5%     97.5% 
#-1.926534  1.424541


#Finding DMRs
print("Find DMRs vs AT:")
dmrs0 <- dmrFinder(AT.tstat, cutoff = c(-1.821612, 1.508577), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation of at least 0.1.
AT.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(AT.dmrs) #240
head(AT.dmrs, n = 3)
write.table(AT.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT6_dmrs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

print("Find DMRs vs NB:")
dmrs0 <- dmrFinder(NB.tstat, cutoff = c(-1.926534, 1.424541), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylationof at least 0.1.
NB.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(NB.dmrs) #353
head(NB.dmrs, n = 3)
write.table(NB.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB6_dmrs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Plot the top DMRs
#blue=BR,red=NB,green=AT
print("Plot top AT DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT6_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, AT.dmrs[1:240,], extend = 5000, 
                addRegions = AT.dmrs)
dev.off()
print("Plot top NB DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB6_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, NB.dmrs[1:353,], extend = 5000, 
                addRegions = NB.dmrs)
dev.off()
```
Differential methylation between week 9 samples:
```R
library(bsseq)
library(stats)
library(BiocParallel)

#Read in data:
col_names <- c("treatment", "replicate", "col")
row_names <- c("BR9_E2_1", "BR9_E2_2", "BR9_E2_3", "NB9_E2_1", "NB9_E2_2", "NB9_E2_3", "AT9_E2_1", "AT9_E2_2", "AT9_E2_3")
data <- matrix(c("control", "control", "control", "NB", "NB", "NB", "AT", "AT", "AT", 1, 2, 3, 1, 2, 3, 1, 2, 3, "blue", "blue", "blue", "red", "red", "red", "green", "green", "green"), nrow = length(row_names), ncol = length(col_names))
df <- data.frame(data, row.names = row_names)
colnames(df) <- col_names
print("Input dataframe:")
print(df)

bsseq <- read.bismark(files = c("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_1/bismark/BR9_E2_1_bismark.deduplicated.CpG_report.txt",
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_2/bismark/BR9_E2_2_bismark.deduplicated.CpG_report.txt",
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_3/bismark/BR9_E2_3_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB9_E2_1/bismark/NB9_E2_1_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB9_E2_2/bismark/NB9_E2_2_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB9_E2_3/bismark/NB9_E2_3_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT9_E2_1/bismark/AT9_E2_1_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT9_E2_2/bismark/AT9_E2_2_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT9_E2_3/bismark/AT9_E2_3_bismark.deduplicated.CpG_report.txt"),
    colData = df,
    rmZeroCov = FALSE,
    strandCollapse = TRUE,
    verbose = TRUE)
save(bsseq, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/9_bsseq.RData")

sapply(assays(bsseq, withDimnames = FALSE), class)
print("Bsseq object:")
bsseq
pData(bsseq)

print("The average coverage of CpGs:")
round(colMeans(getCoverage(bsseq)), 1)
#BR9_E2_1 BR9_E2_2 BR9_E2_3 NB9_E2_1 NB9_E2_2 NB9_E2_3 AT9_E2_1 AT9_E2_2 AT9_E2_3
#    23.1     22.8     23.3     23.5     28.3     19.0     23.6     22.5    23.3
print("The number of CpGs:")
length(bsseq)
#10,962,492
print("Number of CpGs which are covered by at least 1 read in all samples:")
sum(rowSums(getCoverage(bsseq) >= 1) == 7)
#65,102
print("Number of CpGs with 0 coverage in all samples:")
sum(rowSums(getCoverage(bsseq)) == 0)
#266,298


#Perform smoothing:
#"ns is the minimum number of CpGs contained in each window, h is half the minimum window with (the actual window width is either 2 times h or wide enough to contain ns covered CpGs, whichever is greater). Note that the window width is different at each position in the genome and may also be different for different samples at the same position, since it depends on how many nearby CpGs with non-zero coverage. Per default, a smoothing cluster is a whole chromosome. By “cluster” we mean a set of CpGs which are processed together. This means that even if there is a large distance between two CpGs, we borrow strength between them. By setting maxGap this can be prevented since the argument describes the longest distance between two CpGs before a cluster is broken up into two clusters." - all default:

bssmooth <- BSmooth(
    BSseq = bsseq, 
    ns = 70,
    h = 1000,
    maxGap = 10^8,
    BPPARAM = MulticoreParam(workers = 1), 
    verbose = TRUE)
bssmooth
save(bssmooth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/9_bsmooth.RData")

#Remove CpGs with coverage below 4 in all samples:
BS.cov <- getCoverage(bssmooth)
keepLoci.ex <- which(rowSums(BS.cov[, bsseq$treatment == "control"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "NB"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "AT"] >= 4) >= 3)
print("The number of CpGs with coverage >=4 in all samples:")
length(keepLoci.ex)
#9,717,924
bssmooth <- bssmooth[keepLoci.ex,]


#Compute t-statistics based on smoothed whole-genome bisulfite sequencing data:
print("Compute t-stats vs AT:")
AT.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("AT9_E2_1", "AT9_E2_2", "AT9_E2_3"),  
                                    group2 = c("BR9_E2_1", "BR9_E2_2", "BR9_E2_3"),
                                    estimate.var = "same",
                                    local.correct = TRUE,
                                    verbose = TRUE)
AT.tstat
#9,717,924 methylation loci
temp <- AT.tstat@stats
temp2 <- data.frame(temp)
temp2 <- na.omit(temp2) #9,717,196 methylation loci remain
quantile(temp2$tstat, c(0.025, 0.975))
#     2.5%     97.5% 
#-1.610949  1.665140
print("Compute t-stats vs NB:")
NB.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("NB9_E2_1", "NB9_E2_2", "NB9_E2_3"),  
                                    group2 = c("BR9_E2_1", "BR9_E2_2", "BR9_E2_3"),
                                    estimate.var = "same",
                                    local.correct = TRUE,
                                    verbose = TRUE)
NB.tstat
#9,717,924 methylation loci
temp <- NB.tstat@stats
temp2 <- data.frame(temp)
temp2 <- na.omit(temp2) #9,717,157 methylation loci remain
quantile(temp2$tstat, c(0.025, 0.975))
#     2.5%     97.5% 
#-1.575890  1.750382


#Finding DMRs
print("Find DMRs vs AT:")
dmrs0 <- dmrFinder(AT.tstat, cutoff = c(-1.610949, 1.665140), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation of at least 0.1.
AT.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(AT.dmrs) #208
head(AT.dmrs, n = 3)
write.table(AT.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT9_dmrs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

print("Find DMRs vs NB:")
dmrs0 <- dmrFinder(NB.tstat, cutoff = c(-1.575890, 1.750382), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylationof at least 0.1.
NB.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(NB.dmrs) #405
head(NB.dmrs, n = 3)
write.table(NB.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB9_dmrs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Plot the top DMRs
#blue=BR,red=NB,green=AT
print("Plot top AT DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT9_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, AT.dmrs[1:208,], extend = 5000, 
                addRegions = AT.dmrs)
dev.off()
print("Plot top NB DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB9_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, NB.dmrs[1:405,], extend = 5000, 
                addRegions = NB.dmrs)
dev.off()
```
Differential methylation between week 25 samples:
```R
library(bsseq)
library(stats)
library(BiocParallel)

#Read in data:
col_names <- c("treatment", "replicate", "col")
row_names <- c("BR25_E2_1", "BR25_E2_2", "BR25_E2_3", "NB25_E2_1", "NB25_E2_2", "NB25_E2_3", "AT25_E2_1", "AT25_E2_2", "AT25_E2_3")
data <- matrix(c("control", "control", "control", "NB", "NB", "NB", "AT", "AT", "AT", 1, 2, 3, 1, 2, 3, 1, 2, 3, "blue", "blue", "blue", "red", "red", "red", "green", "green", "green"), nrow = length(row_names), ncol = length(col_names))
df <- data.frame(data, row.names = row_names)
colnames(df) <- col_names
print("Input dataframe:")
print(df)

bsseq <- read.bismark(files = c("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_1/bismark/BR25_E2_1_bismark.deduplicated.CpG_report.txt",
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_2/bismark/BR25_E2_2_bismark.deduplicated.CpG_report.txt",
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_3/bismark/BR25_E2_3_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB25_E2_1/bismark/NB25_E2_1_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB25_E2_2/bismark/NB25_E2_2_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB25_E2_3/bismark/NB25_E2_3_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT25_E2_1/bismark/AT25_E2_1_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT25_E2_2/bismark/AT25_E2_2_bismark.deduplicated.CpG_report.txt", 
     "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT25_E2_3/bismark/AT25_E2_3_bismark.deduplicated.CpG_report.txt"),
    colData = df,
    rmZeroCov = FALSE,
    strandCollapse = TRUE,
    verbose = TRUE)
save(bsseq, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/25_bsseq.RData")

sapply(assays(bsseq, withDimnames = FALSE), class)
print("Bsseq object:")
bsseq
pData(bsseq)

print("The average coverage of CpGs:")
round(colMeans(getCoverage(bsseq)), 1)
#BR25_E2_1 BR25_E2_2 BR25_E2_3 NB25_E2_1 NB25_E2_2 NB25_E2_3 AT25_E2_1 AT25_E2_2 AT25_E2_3
#     20.3      19.4      22.2      22.5      23.5      21.7      25.0      20.8     30.0
print("The number of CpGs:")
length(bsseq)
#10,962,492
print("Number of CpGs which are covered by at least 1 read in all samples:")
sum(rowSums(getCoverage(bsseq) >= 1) == 7)
#65,610
print("Number of CpGs with 0 coverage in all samples:")
sum(rowSums(getCoverage(bsseq)) == 0)
#262,245


#Perform smoothing:
#"ns is the minimum number of CpGs contained in each window, h is half the minimum window with (the actual window width is either 2 times h or wide enough to contain ns covered CpGs, whichever is greater). Note that the window width is different at each position in the genome and may also be different for different samples at the same position, since it depends on how many nearby CpGs with non-zero coverage. Per default, a smoothing cluster is a whole chromosome. By “cluster” we mean a set of CpGs which are processed together. This means that even if there is a large distance between two CpGs, we borrow strength between them. By setting maxGap this can be prevented since the argument describes the longest distance between two CpGs before a cluster is broken up into two clusters." - all default:

bssmooth <- BSmooth(
    BSseq = bsseq, 
    ns = 70,
    h = 1000,
    maxGap = 10^8,
    BPPARAM = MulticoreParam(workers = 1), 
    verbose = TRUE)
bssmooth
save(bssmooth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/25_bsmooth.RData")

#Remove CpGs with coverage below 4 in all samples:
BS.cov <- getCoverage(bssmooth)
keepLoci.ex <- which(rowSums(BS.cov[, bsseq$treatment == "control"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "NB"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "AT"] >= 4) >= 3)
print("The number of CpGs with coverage >=4 in all samples:")
length(keepLoci.ex)
#9,706,324
bssmooth <- bssmooth[keepLoci.ex,]


#Compute t-statistics based on smoothed whole-genome bisulfite sequencing data:
print("Compute t-stats vs AT:")
AT.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("AT25_E2_1", "AT25_E2_2", "AT25_E2_3"),  
                                    group2 = c("BR25_E2_1", "BR25_E2_2", "BR25_E2_3"),
                                    estimate.var = "same",
                                    local.correct = TRUE,
                                    verbose = TRUE)
AT.tstat
#9706324 methylation loci
temp <- AT.tstat@stats
temp2 <- data.frame(temp)
temp2 <- na.omit(temp2) #9705661 methylation loci remain
quantile(temp2$tstat, c(0.025, 0.975))
#     2.5%     97.5% 
#-1.552562  1.728513 

print("Compute t-stats vs NB:")
NB.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("NB25_E2_1", "NB25_E2_2", "NB25_E2_3"),  
                                    group2 = c("BR25_E2_1", "BR25_E2_2", "BR25_E2_3"),
                                    estimate.var = "same",
                                    local.correct = TRUE,
                                    verbose = TRUE)
NB.tstat
#9706324 methylation loci
temp <- NB.tstat@stats
temp2 <- data.frame(temp)
temp2 <- na.omit(temp2) #9705661 methylation loci remain
quantile(temp2$tstat, c(0.025, 0.975))
#     2.5%     97.5% 
#-1.496220  1.789964

#Finding DMRs
print("Find DMRs vs AT:")
dmrs0 <- dmrFinder(AT.tstat, cutoff = c(-1.552562, 1.728513), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation of at least 0.1.
AT.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(AT.dmrs) #347
head(AT.dmrs, n = 3)
write.table(AT.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT25_dmrs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

print("Find DMRs vs NB:")
dmrs0 <- dmrFinder(NB.tstat, cutoff = c(-1.496220, 1.789964), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation of at least 0.1.
NB.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(NB.dmrs) #521
head(NB.dmrs, n = 3)
write.table(NB.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB25_dmrs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Plot the DMRs
#blue=BR,red=NB,green=AT
print("Plot top AT DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT25_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, AT.dmrs[1:347,], extend = 5000, 
                addRegions = AT.dmrs)
dev.off()
print("Plot top NB DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB25_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, NB.dmrs[1:521,], extend = 5000, 
                addRegions = NB.dmrs)
dev.off()
```
Convert from bsmap format to bsseq input format. - ERROR
```R
library(bsseq)

read.bsmap <- function(file) {
    dat <- read.table(
        file,
        skip = 1,
        row.names = NULL,
        col.names = c("chr", "pos", "strand", "context", "ratio", "eff_CT_count", "C_count", "CT_count", "rev_G_count", "rev_GA_count", "CI_lower", "CI_upper"),
        colClasses = c("character", "integer", "character", "character", "numeric", "numeric", "integer", "integer", "integer", "integer", "numeric", "numeric"))
    #remove all non-CpG calls.  This includes SNPs
    dat <- dat[dat$context == "CG", ]
    dat$context <- NULL
    dat$chr <- paste("chr", dat$chr, sep = "")
    #join separate lines for each strand
    tmp <- dat[dat$strand == "+", ]
    BS.forward <- BSseq(
        pos = tmp$pos,
        chr = tmp$chr,
        M = as.matrix(tmp$C_count, ncol = 1),
        Cov = as.matrix(tmp$CT_count, ncol = 1),
        sampleNames = "forward")
    tmp <- dat[dat$strand == "-", ]
    BS.reverse <- BSseq(
        pos = tmp$pos - 1L,
        chr = tmp$chr,
        M = as.matrix(tmp$C_count, ncol = 1),
        Cov = as.matrix(tmp$CT_count, ncol = 1),
        sampleNames = "reverse")
    BS <- combine(BS.forward, BS.reverse)
    BS <- collapseBSseq(BS, group = c("a", "a"), type = "integer")
    BS
}

# List all the files
file_paths <- list.files("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/", pattern = "bsmap_ratios_filtered.txt", recursive = TRUE, full.names = TRUE)

# Loop through each file and combine
bs_objects <- list()

for (file_path in file_paths) {
  sample_name <- sub("_bsmap_ratios_filtered", "", tools::file_path_sans_ext(basename(file_path)))
  assign(paste0("BS.", sample_name), read.bsmap(file_path))
  bs_objects[[sample_name]] <- get(paste0("BS.", sample_name))
  sampleNames(bs_objects[[sample_name]]) <- sample_name
}

BS.hostswap <- do.call(combine, bs_objects)

#Add replicate information
pData(BS.hostswap)$Rep <- c("replicate1", "replicate2", "replicate1", "replicate2", "replicate3", "replicate1", "replicate2", "replicate3", "replicate1", "replicate2", "replicate3", "replicate1", "replicate2", "replicate3", "replicate1", "replicate1", "replicate2", "replicate3", "replicate1", "replicate2", "replicate3", "replicate1", "replicate2", "replicate3", "replicate1", "replicate2", "replicate3", "replicate1", "replicate2", "replicate3", "replicate1", "replicate2", "replicate3", "replicate1", "replicate2", "replicate3", "replicate1", "replicate2", "replicate3", "replicate1", "replicate2", "replicate3", "replicate1", "replicate2")
validObject(BS.hostswap)
pData(BS.hostswap)

#Save to file
save(BS.hostswap, file = "BS.hostswap.rda")
tools::resaveRdaFiles("BS.hostswap.rda")
```

### Methylkit <a name="21"></a>
Common approaches for differential methylation analysis are Bsmooth, Methylkit or a custom approach to define differentially methylated regions (DMRs) DOI:10.1093/bib/bbx077. Here methylkit:

#### BSmap alignment <a name="22"></a>

The best C2T aware aligner seems to be bsmap based upon the literature (eg.: DOI:10.1016/j.csbj.2022.08.051). Bsmap alignment and cytosine report was performed:
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/WGBS/Archana_Mar2021/*/trim_galore/); do
    Jobs=$(squeue -u did23faz| grep 'bsmap'| wc -l)
    Reference_genome=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    OutDir=$(echo ${ReadDir} | sed 's@trim_galore@bsmap@g' | sed 's@dna_qc@alignment@g')
    OutFile=$(basename $Fread | sed 's@_trimmed_1.fq.gz@bsmap@g')
    ProgDir=~/git_repos/Wrappers/NBI
    if [ ! -e "${OutDir}/${OutFile}_ratios.txt" ] || [ ! -s "${OutDir}/${OutFile}_ratios.txt" ] || [ ! -e "${OutDir}/${OutFile}.bam" ] || [ ! -s "${OutDir}/${OutFile}.bam" ] || [ ! -e "${OutDir}/${OutFile}.bam.bai" ] || [ ! -s "${OutDir}/${OutFile}.bam.bai" ]; then
    while [ $Jobs -gt 5 ]; do
        sleep 300s
        printf "."
        Jobs=$(squeue -u did23faz| grep 'bsmap'| wc -l)
    done
    sbatch $ProgDir/run_bsmap.sh $OutDir $OutFile $Reference_genome $Fread $Rread 
    else
    echo Already run for $ReadDir  
    fi
done
```
**Cytosine report Output format**: tab delimited txt file with the following columns:
chr: Chromosome or scaffold name where the cytosine is located.
pos: Position of the cytosine in the chromosome or scaffold.
strand: Strand of DNA (either "+" or "-") where the cytosine is located.
context: The sequence context of the cytosine. In this case, it is "CHH," which means the cytosine is followed by two non-cytosine bases in the 3' to 5' direction.
ratio: Methylation ratio at the given cytosine position. It represents the proportion of methylated cytosines out of the total observed cytosines.
eff_CT_count: Effective count of cytosines considered for calculating the methylation ratio. This count excludes certain cytosines based on specific criteria (e.g., filtering low-quality reads).
C_count: Count of methylated cytosines.
CT_count: Total count of cytosines (both methylated and unmethylated).
rev_G_count: Count of guanines on the reverse strand corresponding to the cytosine.
rev_GA_count: Count of guanine-adenine pairs on the reverse strand corresponding to the cytosine.
CI_lower: Lower bound of the confidence interval for the methylation ratio.
CI_upper: Upper bound of the confidence interval for the methylation ratio.

#### Filter cytosine postions <a name="23"></a>
```bash
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/*/bsmap/*bsmap_ratios.txt); do
echo $file >> logs/bsmap_report.txt
cat $file | wc -l >> logs/bsmap_report.txt
done
```
Samples BR0_E2_2 and BR0_E2_3 have significantly fewer methylation sites than other samples - ommitted from further analysis.
```bash
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/*/bsmap/*bsmap_ratios.txt | grep -v 'BR0_E2_2\|BR0_E2_3' > temp_file_list.txt
```
Only cytosines with a depth of at least four in all libraries were considered:
```bash
#remove entries with depth below 4
for file in $(cat temp_file_list.txt); do
    awk '!/^chr/ && $6 >= 4' $file > $(dirname $file)/gooddepth_$(basename $file)
done

#get common line IDs across all files
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/*/bsmap/gooddepth_*bsmap_ratios.txt > temp_file_list2.txt
for file in $(cat temp_file_list2.txt); do
    echo $file
    LC_COLLATE=C sort -k1,1 -k2,2n $file > $(dirname $file)/$(basename $file | sed 's@.txt@_sorted.txt@g') && rm $file
done

files=($(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/*/bsmap/gooddepth_*bsmap_ratios_sorted.txt))
cut -f1,2 "${files[0]}" > temp_common_ids.txt
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/*/bsmap/gooddepth_*bsmap_ratios_sorted.txt > temp_file_list3.txt
cat temp_file_list3.txt | wc -l
for file in $(cat temp_file_list3.txt); do
    echo $file 
    comm --nocheck-order -12 <(cut -f1,2 $file) temp_common_ids.txt > temp_common_ids2.txt && mv temp_common_ids2.txt temp_common_ids.txt
    cat temp_common_ids.txt | wc -l 
done
cp temp_common_ids.txt common_ids.txt
#This reduces the number of postions: 104,274,111 -> 44,355,294, this seems like a lot given that all samples have >101,000,000 positions with 4x depth individually? 
cp temp_common_ids_3.txt common_ids3.txt
#With depth 3; this reduces the number of postions: 106220180 -> 60736387
cp temp_common_ids_2.txt common_ids2.txt
#With depth 2; this reduces the number of postions: 108190661 -> 84108985
cp temp_common_ids_1.txt common_ids1.txt
#With depth 1, ie. sites that overlap between all samples; this reduces the number of postions: 110393989 -> 92906702

#Investigate each sample for low overlap in depth:
files=($(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/*/bsmap/gooddepth_*bsmap_ratios_sorted.txt))
cut -f1,2 "${files[0]}" > temp_common_ids_00.txt
for file in $(tac temp_file_list3.txt); do
    echo $file >> logs/bsmap_4depth_report.txt
    comm --nocheck-order -12 <(cut -f1,2 $file) temp_common_ids_00.txt > temp.txt 
    cat temp.txt  | wc -l >> logs/bsmap_4depth_report.txt
done
#Also run for depths of 2 & 3: logs/bsmap_2depth_report.txt,logs/bsmap_3depth_report.txt,logs/bsmap_4depth_report.txt
#AT1_E2_3 has substantially fewer sites with depth >2, NB9_E2_3 has substantially fewer sites with depth 4. BR0_E2_1 also starts to drop off vs other samples at depth >3 although not to the same degree as AT1_E2_3, BR0_E2_1 is the only BR0 sample after discarding BR0_E2_2 and BR0_E2_3.

#Without samples with low depth:
files=($(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/*/bsmap/gooddepth_*bsmap_ratios_sorted.txt))
cut -f1,2 "${files[0]}" > temp_common_ids_000.txt
for file in $(tac temp_file_list3.txt | grep -v 'AT1_E2_3\|NB9_E2_3'); do
    echo $file 
    comm --nocheck-order -12 <(cut -f1,2 $file) temp_common_ids_000.txt > temp0.txt && mv temp0.txt temp_common_ids_000.txt
    cat temp_common_ids_000.txt | wc -l 
done
cp temp_common_ids_000.txt common_ids_000.txt
#This reduces the number of postions: 104,274,111 -> 72,842,711. 
cp temp_common_ids_003.txt common_ids_003.txt
#With depth 3; this reduces the number of postions: 100,041,850 -> 78,904,427

#Extract lines with common IDs from each file
for file in $(cat temp_file_list.txt); do
  head -n 1 $file > $(dirname $file)/$(basename $file | sed 's@.txt@_filtered.txt@g')
  grep -F -w -f common_ids_000.txt $file >> $(dirname $file)/$(basename $file | sed 's@.txt@_filtered.txt@g')
done
```
Therefore we have 72,842,711 cytosine positions with at least 4x coverage across all samples, but have dropped samples BR0_E2_2, BR0_E2_3, AT1_E2_3 & NB9_E2_3. This leaves only one 0 week sample remaining. No samples were dropped in the bsmooth analysis.

Renomve universally un-methylated sites (to save processing time later):
```bash
#Find sites where ratio is always zero
for file in $(cat temp_file_list.txt); do
    temp=$(dirname $file)/0_$(basename $file)
    awk '!/^chr/ && $5 == 0.000' $file > $temp
    LC_COLLATE=C sort -k1,1 -k2,2n $temp > $(dirname $temp)/$(basename $temp | sed 's@.txt@_sorted.txt@g') && rm $temp
done

files=($(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/*/bsmap/0_*bsmap_ratios_sorted.txt))
cut -f1,2 "${files[0]}" > temp_common_ids_0.txt
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/*/bsmap/0_*bsmap_ratios_sorted.txt > temp_file_list0.txt
cat temp_file_list0.txt | wc -l
for file in $(cat temp_file_list0.txt | grep -v 'AT1_E2_3\|NB9_E2_3'); do
    echo $file >> logs/unmethylated_report.txt
    comm --nocheck-order -12 <(cut -f1,2 $file) temp_common_ids_0.txt > temp_common_ids_20.txt && mv temp_common_ids_20.txt temp_common_ids_0.txt
    cat temp_common_ids_0.txt | wc -l >> logs/unmethylated_report.txt
done
cp temp_common_ids_0.txt common_ids_0.txt
```
The number of sites which are totally unmethylated across all samples is only 690. Given this low number I will not bother to remove these sites.

#### Methylkit <a name="24"></a>
```bash
#Seperate CpG methylation context Cs.
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/*/bsmap/*filtered.txt); do
    awk '{OFS="\t";if ($4 ~ /CG/) {print $1, $2, $3, $4, $5, int($6);}}' $file > $(dirname $file)/CpG_$(basename $file)
    cat $(dirname $file)/CpG_$(basename $file) | wc -l
done
#There are 13,682,374 CpG sites

#Create output directory for output files:
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/methylkit1.28.0.sif R
```
Differential methylation between CpG sites, week 1 samples:
```R
library(methylKit)
#Read in the methylation ratio files
NB.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_1/bsmap/CpG_BR1_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_2/bsmap/CpG_BR1_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_3/bsmap/CpG_BR1_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_1/bsmap/CpG_NB1_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_2/bsmap/CpG_NB1_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_3/bsmap/CpG_NB1_E2_3_bsmap_ratios_filtered.txt")
NB.week1=methRead(NB.file.list, 
    sample.id=list("BR1_E2_1","BR1_E2_2","BR1_E2_3","NB1_E2_1","NB1_E2_2","NB1_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
NB.week1
head(NB.week1[[1]])
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR1_E2_1_CpGhistogram.png")
getMethylationStats(NB.week1[[1]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR1_E2_2_CpGhistogram.png")
getMethylationStats(NB.week1[[2]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR1_E2_3_CpGhistogram.png")
getMethylationStats(NB.week1[[3]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_E2_1_CpGhistogram.png")
getMethylationStats(NB.week1[[4]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_E2_2_CpGhistogram.png")
getMethylationStats(NB.week1[[5]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_E2_3_CpGhistogram.png")
getMethylationStats(NB.week1[[6]], plot=TRUE, both.strands=FALSE)
dev.off()

#Normalisation and filtering
NB.week1.filt <- filterByCoverage(NB.week1,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

NB.week1.filt.norm <- normalizeCoverage(NB.week1.filt, method = "median")
NB.week1.meth <- unite(NB.week1.filt.norm, destrand=FALSE)
NB.week1.meth
# get percent methylation matrix
NB.pm=percMethylation(NB.week1.meth)
# calculate standard deviation of CpGs
NB.sds=matrixStats::rowSds(NB.pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_sdshistogram.png")
hist(NB.sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
NB.week1.meth <- NB.week1.meth[NB.sds > 2]
# This leaves us with this number of CpG sites
nrow(NB.week1.meth)

#Plot data structure
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_correlation.png")
getCorrelation(NB.week1.meth,plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_dendogram.png")
clusterSamples(NB.week1.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_PCA.png")
PCASamples(NB.week1.meth)
dev.off()

save(NB.week1.meth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1meth.RData")

#################################################################################################################################
#Read in the methylation ratio files
AT.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_1/bsmap/CpG_BR1_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_2/bsmap/CpG_BR1_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_3/bsmap/CpG_BR1_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT1_E2_1/bsmap/CpG_AT1_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT1_E2_2/bsmap/CpG_AT1_E2_2_bsmap_ratios_filtered.txt")
AT.week1=methRead(AT.file.list, 
    sample.id=list("BR1_E2_1","BR1_E2_2","BR1_E2_3","AT1_E2_1","AT1_E2_2"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
AT.week1
head(AT.week1[[1]])
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_E2_1_CpGhistogram.png")
getMethylationStats(AT.week1[[4]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_E2_2_CpGhistogram.png")
getMethylationStats(AT.week1[[5]], plot=TRUE, both.strands=FALSE)
dev.off()

#Normalisation and filtering
AT.week1.filt <- filterByCoverage(AT.week1,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

AT.week1.filt.norm <- normalizeCoverage(AT.week1.filt, method = "median")
AT.week1.meth <- unite(AT.week1.filt.norm, destrand=FALSE)
AT.week1.meth
# get percent methylation matrix
AT.pm=percMethylation(AT.week1.meth)
# calculate standard deviation of CpGs
AT.sds=matrixStats::rowSds(AT.pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_sdshistogram.png")
hist(AT.sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
AT.week1.meth <- AT.week1.meth[AT.sds > 2]
# This leaves us with this number of CpG sites
nrow(AT.week1.meth)

#Plot data structure
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_correlation.png")
getCorrelation(AT.week1.meth,plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_dendogram.png")
clusterSamples(AT.week1.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_PCA.png")
PCASamples(AT.week1.meth)
dev.off()

save(AT.week1.meth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1meth.RData")

#################################################################################################################################

#Identify differential methylation
NB.Diff <- calculateDiffMeth(NB.week1.meth,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
NB.Diff
save(NB.Diff, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_diffmeth.RData")
<<<<<<< HEAD

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_Volcano.png")
#plot(NB.Diff$meth.diff, -log10(NB.Diff$qvalue))
#abline(v=0)
#dev.off()

=======

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_Volcano.png")
#plot(NB.Diff$meth.diff, -log10(NB.Diff$qvalue))
#abline(v=0)
#dev.off()

>>>>>>> 9bb487a11dfffd64de11be5c3e2043aee326f5d0
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_diffmethperchr.png")
diffMethPerChr(NB.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
NB.Diff.25p.hyper <- getMethylDiff(NB.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
NB.Diff.25p.hyper <- NB.Diff.25p.hyper[order(NB.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
NB.Diff.25p.hypo <- getMethylDiff(NB.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
NB.Diff.25p.hypo <- NB.Diff.25p.hypo[order(NB.Diff.25p.hypo$qvalue), ]

# get all differentially methylated bases and order by qvalue
NB.Diff.25p <- getMethylDiff(NB.Diff,
                        difference=25,
                        qvalue=0.01)
NB.Diff.25p <- NB.Diff.25p[order(NB.Diff.25p$qvalue), ]

write.table(NB.Diff.25p.hyper, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_Diff_25p_hyper.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(NB.Diff.25p.hypo, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_Diff_25p_hypo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(NB.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_Diff_25p.txt", sep = "\t", quote = FALSE, row.names = FALSE)
<<<<<<< HEAD

save(NB.Diff.25p.hyper, NB.Diff.25p.hyper, NB.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1diffmeth.RData")

#################################################################################################################################

#Identify differential methylation
AT.Diff <- calculateDiffMeth(AT.week1.meth,
                            treatment=c(0,0,0,1,1),
                            overdispersion = "MN",
                            adjust="BH")
AT.Diff
save(AT.Diff, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_diffmeth.RData")

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_Volcano.png")
#plot(AT.Diff$meth.diff, -log10(AT.Diff$qvalue))
#abline(v=0)
#dev.off()

png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_diffmethperchr.png")
diffMethPerChr(AT.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
AT.Diff.25p.hyper <- getMethylDiff(AT.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
AT.Diff.25p.hyper <- AT.Diff.25p.hyper[order(AT.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
AT.Diff.25p.hypo <- getMethylDiff(AT.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
AT.Diff.25p.hypo <- AT.Diff.25p.hypo[order(AT.Diff.25p.hypo$qvalue), ]

# get all differentially methylated bases and order by qvalue
AT.Diff.25p <- getMethylDiff(AT.Diff,
                        difference=25,
                        qvalue=0.01)
AT.Diff.25p <- AT.Diff.25p[order(AT.Diff.25p$qvalue), ]

write.table(AT.Diff.25p.hyper, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_Diff_25p_hyper.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.Diff.25p.hypo, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_Diff_25p_hypo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_Diff_25p.txt", sep = "\t", quote = FALSE, row.names = FALSE)

=======

save(NB.Diff.25p.hyper, NB.Diff.25p.hyper, NB.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1diffmeth.RData")

#################################################################################################################################

#Identify differential methylation
AT.Diff <- calculateDiffMeth(AT.week1.meth,
                            treatment=c(0,0,0,1,1),
                            overdispersion = "MN",
                            adjust="BH")
AT.Diff
save(AT.Diff, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_diffmeth.RData")

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_Volcano.png")
#plot(AT.Diff$meth.diff, -log10(AT.Diff$qvalue))
#abline(v=0)
#dev.off()

png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_diffmethperchr.png")
diffMethPerChr(AT.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
AT.Diff.25p.hyper <- getMethylDiff(AT.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
AT.Diff.25p.hyper <- AT.Diff.25p.hyper[order(AT.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
AT.Diff.25p.hypo <- getMethylDiff(AT.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
AT.Diff.25p.hypo <- AT.Diff.25p.hypo[order(AT.Diff.25p.hypo$qvalue), ]

# get all differentially methylated bases and order by qvalue
AT.Diff.25p <- getMethylDiff(AT.Diff,
                        difference=25,
                        qvalue=0.01)
AT.Diff.25p <- AT.Diff.25p[order(AT.Diff.25p$qvalue), ]

write.table(AT.Diff.25p.hyper, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_Diff_25p_hyper.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.Diff.25p.hypo, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_Diff_25p_hypo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_Diff_25p.txt", sep = "\t", quote = FALSE, row.names = FALSE)

>>>>>>> 9bb487a11dfffd64de11be5c3e2043aee326f5d0
save(AT.Diff.25p.hyper, AT.Diff.25p.hyper, AT.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1diffmeth.RData")
```
Differential methylation between week 3 samples:
```R
library(methylKit)
#Read in the methylation ratio files
NB.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_1/bsmap/CpG_BR3_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_2/bsmap/CpG_BR3_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_3/bsmap/CpG_BR3_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB3_E2_1/bsmap/CpG_NB3_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB3_E2_2/bsmap/CpG_NB3_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB3_E2_3/bsmap/CpG_NB3_E2_3_bsmap_ratios_filtered.txt")
NB.week1=methRead(NB.file.list, 
    sample.id=list("BR3_E2_1","BR3_E2_2","BR3_E2_3","NB3_E2_1","NB3_E2_2","NB3_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
NB.week1
head(NB.week1[[1]])
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR3_E2_1_CpGhistogram.png")
getMethylationStats(NB.week1[[1]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR3_E2_2_CpGhistogram.png")
getMethylationStats(NB.week1[[2]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR3_E2_3_CpGhistogram.png")
getMethylationStats(NB.week1[[3]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_E2_1_CpGhistogram.png")
getMethylationStats(NB.week1[[4]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_E2_2_CpGhistogram.png")
getMethylationStats(NB.week1[[5]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_E2_3_CpGhistogram.png")
getMethylationStats(NB.week1[[6]], plot=TRUE, both.strands=FALSE)
dev.off()

#Normalisation and filtering
NB.week1.filt <- filterByCoverage(NB.week1,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

NB.week1.filt.norm <- normalizeCoverage(NB.week1.filt, method = "median")
NB.week1.meth <- unite(NB.week1.filt.norm, destrand=FALSE)
NB.week1.meth
# get percent methylation matrix
NB.pm=percMethylation(NB.week1.meth)
# calculate standard deviation of CpGs
NB.sds=matrixStats::rowSds(NB.pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_sdshistogram.png")
hist(NB.sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
NB.week1.meth <- NB.week1.meth[NB.sds > 2]
# This leaves us with this number of CpG sites
nrow(NB.week1.meth)

#Plot data structure
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_correlation.png")
getCorrelation(NB.week1.meth,plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_dendogram.png")
clusterSamples(NB.week1.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_PCA.png")
PCASamples(NB.week1.meth)
dev.off()

save(NB.week1.meth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3meth.RData")

#################################################################################################################################
#Read in the methylation ratio files
AT.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_1/bsmap/CpG_BR3_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_2/bsmap/CpG_BR3_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_3/bsmap/CpG_BR3_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT3_E2_1/bsmap/CpG_AT3_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT3_E2_2/bsmap/CpG_AT3_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT3_E2_3/bsmap/CpG_AT3_E2_3_bsmap_ratios_filtered.txt")
AT.week1=methRead(AT.file.list, 
    sample.id=list("BR3_E2_1","BR3_E2_2","BR3_E2_3","AT3_E2_1","AT3_E2_2","AT3_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
AT.week1
head(AT.week1[[1]])
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_E2_1_CpGhistogram.png")
getMethylationStats(AT.week1[[4]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_E2_2_CpGhistogram.png")
getMethylationStats(AT.week1[[5]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_E2_3_CpGhistogram.png")
getMethylationStats(AT.week1[[6]], plot=TRUE, both.strands=FALSE)
dev.off()

#Normalisation and filtering
AT.week1.filt <- filterByCoverage(AT.week1,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

AT.week1.filt.norm <- normalizeCoverage(AT.week1.filt, method = "median")
AT.week1.meth <- unite(AT.week1.filt.norm, destrand=FALSE)
AT.week1.meth
# get percent methylation matrix
AT.pm=percMethylation(AT.week1.meth)
# calculate standard deviation of CpGs
AT.sds=matrixStats::rowSds(AT.pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_sdshistogram.png")
hist(AT.sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
AT.week1.meth <- AT.week1.meth[AT.sds > 2]
# This leaves us with this number of CpG sites
nrow(AT.week1.meth)

#Plot data structure
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_correlation.png")
getCorrelation(AT.week1.meth,plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_dendogram.png")
clusterSamples(AT.week1.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_PCA.png")
PCASamples(AT.week1.meth)
dev.off()
<<<<<<< HEAD

save(AT.week1.meth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3meth.RData")

#################################################################################################################################

#Identify differential methylation
NB.Diff <- calculateDiffMeth(NB.week1.meth,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
NB.Diff
save(NB.Diff, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_diffmeth.RData")

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_Volcano.png")
#plot(NB.Diff$meth.diff, -log10(NB.Diff$qvalue))
#abline(v=0)
#dev.off()

png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_diffmethperchr.png")
diffMethPerChr(NB.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
NB.Diff.25p.hyper <- getMethylDiff(NB.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
NB.Diff.25p.hyper <- NB.Diff.25p.hyper[order(NB.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
NB.Diff.25p.hypo <- getMethylDiff(NB.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
NB.Diff.25p.hypo <- NB.Diff.25p.hypo[order(NB.Diff.25p.hypo$qvalue), ]

# get all differentially methylated bases and order by qvalue
NB.Diff.25p <- getMethylDiff(NB.Diff,
                        difference=25,
                        qvalue=0.01)
NB.Diff.25p <- NB.Diff.25p[order(NB.Diff.25p$qvalue), ]

write.table(NB.Diff.25p.hyper, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_Diff_25p_hyper.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(NB.Diff.25p.hypo, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_Diff_25p_hypo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(NB.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_Diff_25p.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save(NB.Diff.25p.hyper, NB.Diff.25p.hyper, NB.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3diffmeth.RData")

#################################################################################################################################

#Identify differential methylation
AT.Diff <- calculateDiffMeth(AT.week1.meth,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
AT.Diff
save(AT.Diff, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_diffmeth.RData")

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_Volcano.png")
#plot(AT.Diff$meth.diff, -log10(AT.Diff$qvalue))
#abline(v=0)
#dev.off()

png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_diffmethperchr.png")
diffMethPerChr(AT.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
AT.Diff.25p.hyper <- getMethylDiff(AT.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
AT.Diff.25p.hyper <- AT.Diff.25p.hyper[order(AT.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
AT.Diff.25p.hypo <- getMethylDiff(AT.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
AT.Diff.25p.hypo <- AT.Diff.25p.hypo[order(AT.Diff.25p.hypo$qvalue), ]

# get all differentially methylated bases and order by qvalue
AT.Diff.25p <- getMethylDiff(AT.Diff,
                        difference=25,
                        qvalue=0.01)
AT.Diff.25p <- AT.Diff.25p[order(AT.Diff.25p$qvalue), ]

=======

save(AT.week1.meth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3meth.RData")

#################################################################################################################################

#Identify differential methylation
NB.Diff <- calculateDiffMeth(NB.week1.meth,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
NB.Diff
save(NB.Diff, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_diffmeth.RData")

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_Volcano.png")
#plot(NB.Diff$meth.diff, -log10(NB.Diff$qvalue))
#abline(v=0)
#dev.off()

png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_diffmethperchr.png")
diffMethPerChr(NB.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
NB.Diff.25p.hyper <- getMethylDiff(NB.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
NB.Diff.25p.hyper <- NB.Diff.25p.hyper[order(NB.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
NB.Diff.25p.hypo <- getMethylDiff(NB.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
NB.Diff.25p.hypo <- NB.Diff.25p.hypo[order(NB.Diff.25p.hypo$qvalue), ]

# get all differentially methylated bases and order by qvalue
NB.Diff.25p <- getMethylDiff(NB.Diff,
                        difference=25,
                        qvalue=0.01)
NB.Diff.25p <- NB.Diff.25p[order(NB.Diff.25p$qvalue), ]

write.table(NB.Diff.25p.hyper, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_Diff_25p_hyper.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(NB.Diff.25p.hypo, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_Diff_25p_hypo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(NB.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_Diff_25p.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save(NB.Diff.25p.hyper, NB.Diff.25p.hyper, NB.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3diffmeth.RData")

#################################################################################################################################

#Identify differential methylation
AT.Diff <- calculateDiffMeth(AT.week1.meth,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
AT.Diff
save(AT.Diff, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_diffmeth.RData")

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_Volcano.png")
#plot(AT.Diff$meth.diff, -log10(AT.Diff$qvalue))
#abline(v=0)
#dev.off()

png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_diffmethperchr.png")
diffMethPerChr(AT.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
AT.Diff.25p.hyper <- getMethylDiff(AT.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
AT.Diff.25p.hyper <- AT.Diff.25p.hyper[order(AT.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
AT.Diff.25p.hypo <- getMethylDiff(AT.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
AT.Diff.25p.hypo <- AT.Diff.25p.hypo[order(AT.Diff.25p.hypo$qvalue), ]

# get all differentially methylated bases and order by qvalue
AT.Diff.25p <- getMethylDiff(AT.Diff,
                        difference=25,
                        qvalue=0.01)
AT.Diff.25p <- AT.Diff.25p[order(AT.Diff.25p$qvalue), ]

>>>>>>> 9bb487a11dfffd64de11be5c3e2043aee326f5d0
write.table(AT.Diff.25p.hyper, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_Diff_25p_hyper.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.Diff.25p.hypo, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_Diff_25p_hypo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_Diff_25p.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save(AT.Diff.25p.hyper, AT.Diff.25p.hyper, AT.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3diffmeth.RData")
```
Differential methylation between week 6 samples:
```R
library(methylKit)
#Read in the methylation ratio files
NB.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_1/bsmap/CpG_BR6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_2/bsmap/CpG_BR6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_3/bsmap/CpG_BR6_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_1/bsmap/CpG_NB6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_2/bsmap/CpG_NB6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_3/bsmap/CpG_NB6_E2_3_bsmap_ratios_filtered.txt")
NB.week1=methRead(NB.file.list, 
    sample.id=list("BR6_E2_1","BR6_E2_2","BR6_E2_3","NB6_E2_1","NB6_E2_2","NB6_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
NB.week1
head(NB.week1[[1]])
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR6_E2_1_CpGhistogram.png")
getMethylationStats(NB.week1[[1]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR6_E2_2_CpGhistogram.png")
getMethylationStats(NB.week1[[2]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR6_E2_3_CpGhistogram.png")
getMethylationStats(NB.week1[[3]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_E2_1_CpGhistogram.png")
getMethylationStats(NB.week1[[4]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_E2_2_CpGhistogram.png")
getMethylationStats(NB.week1[[5]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_E2_3_CpGhistogram.png")
getMethylationStats(NB.week1[[6]], plot=TRUE, both.strands=FALSE)
dev.off()

#Normalisation and filtering
NB.week1.filt <- filterByCoverage(NB.week1,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

NB.week1.filt.norm <- normalizeCoverage(NB.week1.filt, method = "median")
NB.week1.meth <- unite(NB.week1.filt.norm, destrand=FALSE)
NB.week1.meth
# get percent methylation matrix
NB.pm=percMethylation(NB.week1.meth)
# calculate standard deviation of CpGs
NB.sds=matrixStats::rowSds(NB.pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_sdshistogram.png")
hist(NB.sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
NB.week1.meth <- NB.week1.meth[NB.sds > 2]
# This leaves us with this number of CpG sites
nrow(NB.week1.meth)

#Plot data structure
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_correlation.png")
getCorrelation(NB.week1.meth,plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_dendogram.png")
clusterSamples(NB.week1.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_PCA.png")
PCASamples(NB.week1.meth)
dev.off()

save(NB.week1.meth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6meth.RData")

#################################################################################################################################
#Read in the methylation ratio files
AT.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_1/bsmap/CpG_BR6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_2/bsmap/CpG_BR6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_3/bsmap/CpG_BR6_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_1/bsmap/CpG_AT6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_2/bsmap/CpG_AT6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_3/bsmap/CpG_AT6_E2_3_bsmap_ratios_filtered.txt")
AT.week1=methRead(AT.file.list, 
    sample.id=list("BR6_E2_1","BR6_E2_2","BR6_E2_3","AT6_E2_1","AT6_E2_2","AT6_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
AT.week1
head(AT.week1[[1]])
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_E2_1_CpGhistogram.png")
getMethylationStats(AT.week1[[4]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_E2_2_CpGhistogram.png")
getMethylationStats(AT.week1[[5]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_E2_3_CpGhistogram.png")
getMethylationStats(AT.week1[[6]], plot=TRUE, both.strands=FALSE)
dev.off()

#Normalisation and filtering
AT.week1.filt <- filterByCoverage(AT.week1,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

AT.week1.filt.norm <- normalizeCoverage(AT.week1.filt, method = "median")
AT.week1.meth <- unite(AT.week1.filt.norm, destrand=FALSE)
AT.week1.meth
# get percent methylation matrix
AT.pm=percMethylation(AT.week1.meth)
# calculate standard deviation of CpGs
AT.sds=matrixStats::rowSds(AT.pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_sdshistogram.png")
hist(AT.sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
AT.week1.meth <- AT.week1.meth[AT.sds > 2]
# This leaves us with this number of CpG sites
nrow(AT.week1.meth)

#Plot data structure
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_correlation.png")
getCorrelation(AT.week1.meth,plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_dendogram.png")
clusterSamples(AT.week1.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_PCA.png")
PCASamples(AT.week1.meth)
dev.off()

save(AT.week1.meth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6meth.RData")

#################################################################################################################################

#Identify differential methylation
NB.Diff <- calculateDiffMeth(NB.week1.meth,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
NB.Diff
save(NB.Diff, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_diffmeth.RData")
<<<<<<< HEAD

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_Volcano.png")
#plot(NB.Diff$meth.diff, -log10(NB.Diff$qvalue))
#abline(v=0)
#dev.off()

=======

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_Volcano.png")
#plot(NB.Diff$meth.diff, -log10(NB.Diff$qvalue))
#abline(v=0)
#dev.off()

>>>>>>> 9bb487a11dfffd64de11be5c3e2043aee326f5d0
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_diffmethperchr.png")
diffMethPerChr(NB.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
NB.Diff.25p.hyper <- getMethylDiff(NB.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
NB.Diff.25p.hyper <- NB.Diff.25p.hyper[order(NB.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
NB.Diff.25p.hypo <- getMethylDiff(NB.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
NB.Diff.25p.hypo <- NB.Diff.25p.hypo[order(NB.Diff.25p.hypo$qvalue), ]

# get all differentially methylated bases and order by qvalue
NB.Diff.25p <- getMethylDiff(NB.Diff,
                        difference=25,
                        qvalue=0.01)
NB.Diff.25p <- NB.Diff.25p[order(NB.Diff.25p$qvalue), ]

write.table(NB.Diff.25p.hyper, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_Diff_25p_hyper.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(NB.Diff.25p.hypo, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_Diff_25p_hypo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(NB.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_Diff_25p.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save(NB.Diff.25p.hyper, NB.Diff.25p.hyper, NB.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6diffmeth.RData")

#################################################################################################################################

#Identify differential methylation
AT.Diff <- calculateDiffMeth(AT.week1.meth,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
AT.Diff
save(AT.Diff, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_diffmeth.RData")

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_Volcano.png")
#plot(AT.Diff$meth.diff, -log10(AT.Diff$qvalue))
#abline(v=0)
#dev.off()
<<<<<<< HEAD

png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_diffmethperchr.png")
diffMethPerChr(AT.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
AT.Diff.25p.hyper <- getMethylDiff(AT.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
AT.Diff.25p.hyper <- AT.Diff.25p.hyper[order(AT.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
AT.Diff.25p.hypo <- getMethylDiff(AT.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
AT.Diff.25p.hypo <- AT.Diff.25p.hypo[order(AT.Diff.25p.hypo$qvalue), ]

=======

png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_diffmethperchr.png")
diffMethPerChr(AT.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
AT.Diff.25p.hyper <- getMethylDiff(AT.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
AT.Diff.25p.hyper <- AT.Diff.25p.hyper[order(AT.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
AT.Diff.25p.hypo <- getMethylDiff(AT.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
AT.Diff.25p.hypo <- AT.Diff.25p.hypo[order(AT.Diff.25p.hypo$qvalue), ]

>>>>>>> 9bb487a11dfffd64de11be5c3e2043aee326f5d0
# get all differentially methylated bases and order by qvalue
AT.Diff.25p <- getMethylDiff(AT.Diff,
                        difference=25,
                        qvalue=0.01)
AT.Diff.25p <- AT.Diff.25p[order(AT.Diff.25p$qvalue), ]

write.table(AT.Diff.25p.hyper, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_Diff_25p_hyper.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.Diff.25p.hypo, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_Diff_25p_hypo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_Diff_25p.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save(AT.Diff.25p.hyper, AT.Diff.25p.hyper, AT.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6diffmeth.RData")
```
Differential methylation between week 9 samples:
```R
library(methylKit)
#Read in the methylation ratio files
NB.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_1/bsmap/CpG_BR9_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_2/bsmap/CpG_BR9_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_3/bsmap/CpG_BR9_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB9_E2_1/bsmap/CpG_NB9_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB9_E2_2/bsmap/CpG_NB9_E2_2_bsmap_ratios_filtered.txt")
NB.week1=methRead(NB.file.list, 
    sample.id=list("BR9_E2_1","BR9_E2_2","BR9_E2_3","NB9_E2_1","NB9_E2_2"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
NB.week1
head(NB.week1[[1]])
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR9_E2_1_CpGhistogram.png")
getMethylationStats(NB.week1[[1]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR9_E2_2_CpGhistogram.png")
getMethylationStats(NB.week1[[2]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR9_E2_3_CpGhistogram.png")
getMethylationStats(NB.week1[[3]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_E2_1_CpGhistogram.png")
getMethylationStats(NB.week1[[4]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_E2_2_CpGhistogram.png")
getMethylationStats(NB.week1[[5]], plot=TRUE, both.strands=FALSE)
dev.off()

#Normalisation and filtering
NB.week1.filt <- filterByCoverage(NB.week1,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

NB.week1.filt.norm <- normalizeCoverage(NB.week1.filt, method = "median")
NB.week1.meth <- unite(NB.week1.filt.norm, destrand=FALSE)
NB.week1.meth
# get percent methylation matrix
NB.pm=percMethylation(NB.week1.meth)
# calculate standard deviation of CpGs
NB.sds=matrixStats::rowSds(NB.pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_sdshistogram.png")
hist(NB.sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
NB.week1.meth <- NB.week1.meth[NB.sds > 2]
# This leaves us with this number of CpG sites
nrow(NB.week1.meth)

#Plot data structure
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_correlation.png")
getCorrelation(NB.week1.meth,plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_dendogram.png")
clusterSamples(NB.week1.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_PCA.png")
PCASamples(NB.week1.meth)
dev.off()

save(NB.week1.meth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9meth.RData")

#################################################################################################################################
#Read in the methylation ratio files
AT.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_1/bsmap/CpG_BR9_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_2/bsmap/CpG_BR9_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_3/bsmap/CpG_BR9_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT9_E2_1/bsmap/CpG_AT9_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT9_E2_2/bsmap/CpG_AT9_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT9_E2_3/bsmap/CpG_AT9_E2_3_bsmap_ratios_filtered.txt")
AT.week1=methRead(AT.file.list, 
    sample.id=list("BR9_E2_1","BR9_E2_2","BR9_E2_3","AT9_E2_1","AT9_E2_2","AT9_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
AT.week1
head(AT.week1[[1]])
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_E2_1_CpGhistogram.png")
getMethylationStats(AT.week1[[4]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_E2_2_CpGhistogram.png")
getMethylationStats(AT.week1[[5]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_E2_3_CpGhistogram.png")
getMethylationStats(AT.week1[[6]], plot=TRUE, both.strands=FALSE)
dev.off()

#Normalisation and filtering
AT.week1.filt <- filterByCoverage(AT.week1,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

AT.week1.filt.norm <- normalizeCoverage(AT.week1.filt, method = "median")
AT.week1.meth <- unite(AT.week1.filt.norm, destrand=FALSE)
AT.week1.meth
# get percent methylation matrix
AT.pm=percMethylation(AT.week1.meth)
# calculate standard deviation of CpGs
AT.sds=matrixStats::rowSds(AT.pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_sdshistogram.png")
hist(AT.sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
AT.week1.meth <- AT.week1.meth[AT.sds > 2]
# This leaves us with this number of CpG sites
nrow(AT.week1.meth)

#Plot data structure
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_correlation.png")
getCorrelation(AT.week1.meth,plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_dendogram.png")
clusterSamples(AT.week1.meth, dist="correlation", method="ward", plot=TRUE)
<<<<<<< HEAD
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_PCA.png")
PCASamples(AT.week1.meth)
dev.off()
=======
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_PCA.png")
PCASamples(AT.week1.meth)
dev.off()
>>>>>>> 9bb487a11dfffd64de11be5c3e2043aee326f5d0

save(AT.week1.meth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9meth.RData")

#################################################################################################################################

#Identify differential methylation
NB.Diff <- calculateDiffMeth(NB.week1.meth,
                            treatment=c(0,0,0,1,1),
                            overdispersion = "MN",
                            adjust="BH")
NB.Diff
save(NB.Diff, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_diffmeth.RData")

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_Volcano.png")
#plot(NB.Diff$meth.diff, -log10(NB.Diff$qvalue))
#abline(v=0)
#dev.off()

png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_diffmethperchr.png")
diffMethPerChr(NB.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
NB.Diff.25p.hyper <- getMethylDiff(NB.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
NB.Diff.25p.hyper <- NB.Diff.25p.hyper[order(NB.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
NB.Diff.25p.hypo <- getMethylDiff(NB.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
NB.Diff.25p.hypo <- NB.Diff.25p.hypo[order(NB.Diff.25p.hypo$qvalue), ]

# get all differentially methylated bases and order by qvalue
NB.Diff.25p <- getMethylDiff(NB.Diff,
                        difference=25,
                        qvalue=0.01)
NB.Diff.25p <- NB.Diff.25p[order(NB.Diff.25p$qvalue), ]

write.table(NB.Diff.25p.hyper, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_Diff_25p_hyper.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(NB.Diff.25p.hypo, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_Diff_25p_hypo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(NB.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_Diff_25p.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save(NB.Diff.25p.hyper, NB.Diff.25p.hyper, NB.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9diffmeth.RData")

#################################################################################################################################

#Identify differential methylation
AT.Diff <- calculateDiffMeth(AT.week1.meth,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
AT.Diff
save(AT.Diff, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_diffmeth.RData")

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_Volcano.png")
#plot(AT.Diff$meth.diff, -log10(AT.Diff$qvalue))
#abline(v=0)
#dev.off()

png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_diffmethperchr.png")
diffMethPerChr(AT.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
AT.Diff.25p.hyper <- getMethylDiff(AT.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
AT.Diff.25p.hyper <- AT.Diff.25p.hyper[order(AT.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
AT.Diff.25p.hypo <- getMethylDiff(AT.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
AT.Diff.25p.hypo <- AT.Diff.25p.hypo[order(AT.Diff.25p.hypo$qvalue), ]

# get all differentially methylated bases and order by qvalue
AT.Diff.25p <- getMethylDiff(AT.Diff,
                        difference=25,
                        qvalue=0.01)
AT.Diff.25p <- AT.Diff.25p[order(AT.Diff.25p$qvalue), ]

write.table(AT.Diff.25p.hyper, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_Diff_25p_hyper.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.Diff.25p.hypo, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_Diff_25p_hypo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_Diff_25p.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save(AT.Diff.25p.hyper, AT.Diff.25p.hyper, AT.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9diffmeth.RData")
```
Differential methylation between week 25 samples:
```R
library(methylKit)
#Read in the methylation ratio files
NB.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_1/bsmap/CpG_BR25_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_2/bsmap/CpG_BR25_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_3/bsmap/CpG_BR25_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB25_E2_1/bsmap/CpG_NB25_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB25_E2_2/bsmap/CpG_NB25_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB25_E2_3/bsmap/CpG_NB25_E2_3_bsmap_ratios_filtered.txt")
NB.week1=methRead(NB.file.list, 
    sample.id=list("BR25_E2_1","BR25_E2_2","BR25_E2_3","NB25_E2_1","NB25_E2_2","NB25_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
NB.week1
head(NB.week1[[1]])
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR25_E2_1_CpGhistogram.png")
getMethylationStats(NB.week1[[1]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR25_E2_2_CpGhistogram.png")
getMethylationStats(NB.week1[[2]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/BR25_E2_3_CpGhistogram.png")
getMethylationStats(NB.week1[[3]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_E2_1_CpGhistogram.png")
getMethylationStats(NB.week1[[4]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_E2_2_CpGhistogram.png")
getMethylationStats(NB.week1[[5]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_E2_3_CpGhistogram.png")
getMethylationStats(NB.week1[[6]], plot=TRUE, both.strands=FALSE)
dev.off()

#Normalisation and filtering
NB.week1.filt <- filterByCoverage(NB.week1,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

NB.week1.filt.norm <- normalizeCoverage(NB.week1.filt, method = "median")
NB.week1.meth <- unite(NB.week1.filt.norm, destrand=FALSE)
NB.week1.meth
# get percent methylation matrix
NB.pm=percMethylation(NB.week1.meth)
# calculate standard deviation of CpGs
NB.sds=matrixStats::rowSds(NB.pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_sdshistogram.png")
hist(NB.sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
NB.week1.meth <- NB.week1.meth[NB.sds > 2]
# This leaves us with this number of CpG sites
nrow(NB.week1.meth)

#Plot data structure
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_correlation.png")
getCorrelation(NB.week1.meth,plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_dendogram.png")
clusterSamples(NB.week1.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_PCA.png")
PCASamples(NB.week1.meth)
dev.off()

save(NB.week1.meth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25meth.RData")

#################################################################################################################################
#Read in the methylation ratio files
AT.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_1/bsmap/CpG_BR25_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_2/bsmap/CpG_BR25_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_3/bsmap/CpG_BR25_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT25_E2_1/bsmap/CpG_AT25_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT25_E2_2/bsmap/CpG_AT25_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT25_E2_3/bsmap/CpG_AT25_E2_3_bsmap_ratios_filtered.txt")
AT.week1=methRead(AT.file.list, 
    sample.id=list("BR25_E2_1","BR25_E2_2","BR25_E2_3","AT25_E2_1","AT25_E2_2","AT25_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
AT.week1
head(AT.week1[[1]])
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_E2_1_CpGhistogram.png")
getMethylationStats(AT.week1[[4]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_E2_2_CpGhistogram.png")
getMethylationStats(AT.week1[[5]], plot=TRUE, both.strands=FALSE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_E2_3_CpGhistogram.png")
getMethylationStats(AT.week1[[6]], plot=TRUE, both.strands=FALSE)
dev.off()

#Normalisation and filtering
AT.week1.filt <- filterByCoverage(AT.week1,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

AT.week1.filt.norm <- normalizeCoverage(AT.week1.filt, method = "median")
AT.week1.meth <- unite(AT.week1.filt.norm, destrand=FALSE)
AT.week1.meth
# get percent methylation matrix
AT.pm=percMethylation(AT.week1.meth)
# calculate standard deviation of CpGs
AT.sds=matrixStats::rowSds(AT.pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_sdshistogram.png")
hist(AT.sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
AT.week1.meth <- AT.week1.meth[AT.sds > 2]
# This leaves us with this number of CpG sites
nrow(AT.week1.meth)

#Plot data structure
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_correlation.png")
getCorrelation(AT.week1.meth,plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_dendogram.png")
clusterSamples(AT.week1.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_PCA.png")
PCASamples(AT.week1.meth)
dev.off()

save(AT.week1.meth, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25meth.RData")

#################################################################################################################################

#Identify differential methylation
NB.Diff <- calculateDiffMeth(NB.week1.meth,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
NB.Diff
save(NB.Diff, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_diffmeth.RData")
<<<<<<< HEAD

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_Volcano.png")
#plot(NB.Diff$meth.diff, -log10(NB.Diff$qvalue))
#abline(v=0)
#dev.off()

=======

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_Volcano.png")
#plot(NB.Diff$meth.diff, -log10(NB.Diff$qvalue))
#abline(v=0)
#dev.off()

>>>>>>> 9bb487a11dfffd64de11be5c3e2043aee326f5d0
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_diffmethperchr.png")
diffMethPerChr(NB.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
NB.Diff.25p.hyper <- getMethylDiff(NB.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
NB.Diff.25p.hyper <- NB.Diff.25p.hyper[order(NB.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
NB.Diff.25p.hypo <- getMethylDiff(NB.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
NB.Diff.25p.hypo <- NB.Diff.25p.hypo[order(NB.Diff.25p.hypo$qvalue), ]

# get all differentially methylated bases and order by qvalue
NB.Diff.25p <- getMethylDiff(NB.Diff,
                        difference=25,
                        qvalue=0.01)
NB.Diff.25p <- NB.Diff.25p[order(NB.Diff.25p$qvalue), ]

write.table(NB.Diff.25p.hyper, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_Diff_25p_hyper.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(NB.Diff.25p.hypo, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_Diff_25p_hypo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(NB.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_Diff_25p.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save(NB.Diff.25p.hyper, NB.Diff.25p.hyper, NB.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25diffmeth.RData")

#################################################################################################################################

#Identify differential methylation
AT.Diff <- calculateDiffMeth(AT.week1.meth,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
AT.Diff
save(AT.Diff, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_diffmeth.RData")

#png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_Volcano.png")
#plot(AT.Diff$meth.diff, -log10(AT.Diff$qvalue))
#abline(v=0)
#dev.off()

png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_diffmethperchr.png")
diffMethPerChr(AT.Diff) 
dev.off()

# get hyper methylated bases and order by qvalue
AT.Diff.25p.hyper <- getMethylDiff(AT.Diff,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
AT.Diff.25p.hyper <- AT.Diff.25p.hyper[order(AT.Diff.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
AT.Diff.25p.hypo <- getMethylDiff(AT.Diff,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
AT.Diff.25p.hypo <- AT.Diff.25p.hypo[order(AT.Diff.25p.hypo$qvalue), ]

# get all differentially methylated bases and order by qvalue
AT.Diff.25p <- getMethylDiff(AT.Diff,
                        difference=25,
                        qvalue=0.01)
AT.Diff.25p <- AT.Diff.25p[order(AT.Diff.25p$qvalue), ]

write.table(AT.Diff.25p.hyper, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_Diff_25p_hyper.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.Diff.25p.hypo, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_Diff_25p_hypo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_Diff_25p.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save(AT.Diff.25p.hyper, AT.Diff.25p.hyper, AT.Diff.25p, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25diffmeth.RData")
```

#### Genomation <a name="25"></a>
<<<<<<< HEAD

Annotate differentially methylated sites to determine if they fall within promoter, intronic or exonic regions, or CpG islands.

=======

Annotate differentially methylated sites to determine if they fall within promoter, intronic or exonic regions, or CpG islands.

>>>>>>> 9bb487a11dfffd64de11be5c3e2043aee326f5d0
Collect annotation info:
```bash
#Create bed12 file for gene annotations
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/gtf2bed.py /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.gtf > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed12 

#Identify CpG islands 
source package a684a2ed-d23f-4025-aa81-b21e27e458df
cpgplot -sequence /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa -window 100 -minlen 200 -minoe 0.6 -minpc 50. -graph cps -cg -pc -obsexp -outfile /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpgplot -outfeat /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.gff
awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $4-1, $5, $9, ".", $7}' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.gff > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.bed 

#Create output directory for output files:
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/genomation/

singularity exec --overlay /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/overlays/genomation1.34.0-overlay.sif:ro /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/genomation1.34.0.sif R
```
Annotate Week 1 differentially methylated sites:
```R
library("genomation")
library("methylKit")

#Read in data:
load(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1diffmeth.RData")
load(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1diffmeth.RData")
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed12",remove.unusual=FALSE)

#Promoter/exon/introns:
#Annotate differentially methylated sites with promoter/exon/intron information.
AT.Diff.25p.anot <- annotateWithGeneParts(target = as(AT.Diff.25p,"GRanges"),
                                       feature = refseq_anot)
NB.Diff.25p.anot <- annotateWithGeneParts(target = as(NB.Diff.25p,"GRanges"),
                                       feature = refseq_anot)

#Report distance to the nearest Transcription Start Site, the target.row column in the output indicates the row number in the initial target set
AT.dist_tss <- getAssociationWithTSS(AT.Diff.25p.anot)
head(AT.dist_tss)

NB.dist_tss <- getAssociationWithTSS(NB.Diff.25p.anot)
head(NB.dist_tss)

#Report whether the differentially methylated CpGs are within promoters,introns or exons; the order is the same as the target set
ATpie <- getMembers(AT.Diff.25p.anot)
NBpie <- getMembers(NB.Diff.25p.anot)
AT.dist_tss <- cbind(AT.dist_tss, ATpie)
NB.dist_tss <- cbind(NB.dist_tss, NBpie)

#Summarize for all differentially methylated CpGs
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/CpGs_in_gene_region_AT1.png")
plotTargetAnnotation(AT.Diff.25p.anot, main = "Differential Methylation Annotation BR vs AT Week 1")
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/CpGs_in_gene_region_NB1.png")
plotTargetAnnotation(NB.Diff.25p.anot, main = "Differential Methylation Annotation BR vs NB Week 1")
dev.off()

#CpG islands:
#Load the CpG info
cpg_anot <- readFeatureFlank("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)
NB.diffCpGann <- annotateWithFeatureFlank(as(NB.Diff.25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")
AT.diffCpGann <- annotateWithFeatureFlank(as(AT.Diff.25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")

#Do CpG belong to a CpG Island or Shore
NBissh <- getMembers(NB.diffCpGann)
head(NBissh)

ATissh <- getMembers(AT.diffCpGann)
head(ATissh)

AT.dist_tss <- cbind(AT.dist_tss, ATissh)
NB.dist_tss <- cbind(NB.dist_tss, NBissh)

write.table(NB.dist_tss, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_NB1.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.dist_tss, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_AT1.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Summarize for all differentially methylated CpGs
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/Island_shore_NB1.png")
plotTargetAnnotation(NB.diffCpGann, main = "Differential Methylation Annotation NB1")
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/Island_shore_AT1.png")
plotTargetAnnotation(AT.diffCpGann, main = "Differential Methylation Annotation AT1")
dev.off()
```
Annotate Week 3 differentially methylated sites:
```R
library("genomation")
library("methylKit")

#Read in data:
load(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3diffmeth.RData")
load(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3diffmeth.RData")
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed12",remove.unusual=FALSE)

#Promoter/exon/introns:
#Annotate differentially methylated sites with promoter/exon/intron information.
AT.Diff.25p.anot <- annotateWithGeneParts(target = as(AT.Diff.25p,"GRanges"),
                                       feature = refseq_anot)
NB.Diff.25p.anot <- annotateWithGeneParts(target = as(NB.Diff.25p,"GRanges"),
                                       feature = refseq_anot)

#Report distance to the nearest Transcription Start Site, the target.row column in the output indicates the row number in the initial target set
AT.dist_tss <- getAssociationWithTSS(AT.Diff.25p.anot)
head(AT.dist_tss)

NB.dist_tss <- getAssociationWithTSS(NB.Diff.25p.anot)
head(NB.dist_tss)

#Report whether the differentially methylated CpGs are within promoters,introns or exons; the order is the same as the target set
ATpie <- getMembers(AT.Diff.25p.anot)
NBpie <- getMembers(NB.Diff.25p.anot)
AT.dist_tss <- cbind(AT.dist_tss, ATpie)
NB.dist_tss <- cbind(NB.dist_tss, NBpie)

#Summarize for all differentially methylated CpGs
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/CpGs_in_gene_region_AT3.png")
plotTargetAnnotation(AT.Diff.25p.anot, main = "Differential Methylation Annotation BR vs AT Week 3")
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/CpGs_in_gene_region_NB3.png")
plotTargetAnnotation(NB.Diff.25p.anot, main = "Differential Methylation Annotation BR vs NB Week 3")
dev.off()

#CpG islands:
#Load the CpG info
cpg_anot <- readFeatureFlank("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)
NB.diffCpGann <- annotateWithFeatureFlank(as(NB.Diff.25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")
AT.diffCpGann <- annotateWithFeatureFlank(as(AT.Diff.25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")

#Do CpG belong to a CpG Island or Shore
NBissh <- getMembers(NB.diffCpGann)
head(NBissh)

ATissh <- getMembers(AT.diffCpGann)
head(ATissh)

AT.dist_tss <- cbind(AT.dist_tss, ATissh)
NB.dist_tss <- cbind(NB.dist_tss, NBissh)

write.table(NB.dist_tss, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_NB3.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.dist_tss, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_AT3.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Summarize for all differentially methylated CpGs
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/Island_shore_NB3.png")
plotTargetAnnotation(NB.diffCpGann, main = "Differential Methylation Annotation NB3")
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/Island_shore_AT3.png")
plotTargetAnnotation(AT.diffCpGann, main = "Differential Methylation Annotation AT3")
dev.off()
```
Annotate Week 6 differentially methylated sites:
```R
library("genomation")
library("methylKit")

#Read in data:
load(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6diffmeth.RData")
load(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6diffmeth.RData")
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed12",remove.unusual=FALSE)

#Promoter/exon/introns:
#Annotate differentially methylated sites with promoter/exon/intron information.
AT.Diff.25p.anot <- annotateWithGeneParts(target = as(AT.Diff.25p,"GRanges"),
                                       feature = refseq_anot)
NB.Diff.25p.anot <- annotateWithGeneParts(target = as(NB.Diff.25p,"GRanges"),
                                       feature = refseq_anot)

#Report distance to the nearest Transcription Start Site, the target.row column in the output indicates the row number in the initial target set
AT.dist_tss <- getAssociationWithTSS(AT.Diff.25p.anot)
head(AT.dist_tss)

NB.dist_tss <- getAssociationWithTSS(NB.Diff.25p.anot)
head(NB.dist_tss)

#Report whether the differentially methylated CpGs are within promoters,introns or exons; the order is the same as the target set
ATpie <- getMembers(AT.Diff.25p.anot)
NBpie <- getMembers(NB.Diff.25p.anot)
AT.dist_tss <- cbind(AT.dist_tss, ATpie)
NB.dist_tss <- cbind(NB.dist_tss, NBpie)

#Summarize for all differentially methylated CpGs
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/CpGs_in_gene_region_AT6.png")
plotTargetAnnotation(AT.Diff.25p.anot, main = "Differential Methylation Annotation BR vs AT Week 6")
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/CpGs_in_gene_region_NB6.png")
plotTargetAnnotation(NB.Diff.25p.anot, main = "Differential Methylation Annotation BR vs NB Week 6")
dev.off()

#CpG islands:
#Load the CpG info
cpg_anot <- readFeatureFlank("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)
NB.diffCpGann <- annotateWithFeatureFlank(as(NB.Diff.25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")
AT.diffCpGann <- annotateWithFeatureFlank(as(AT.Diff.25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")

#Do CpG belong to a CpG Island or Shore
NBissh <- getMembers(NB.diffCpGann)
head(NBissh)

ATissh <- getMembers(AT.diffCpGann)
head(ATissh)

AT.dist_tss <- cbind(AT.dist_tss, ATissh)
NB.dist_tss <- cbind(NB.dist_tss, NBissh)

write.table(NB.dist_tss, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_NB6.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.dist_tss, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_AT6.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Summarize for all differentially methylated CpGs
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/Island_shore_NB6.png")
plotTargetAnnotation(NB.diffCpGann, main = "Differential Methylation Annotation NB6")
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/Island_shore_AT6.png")
plotTargetAnnotation(AT.diffCpGann, main = "Differential Methylation Annotation AT6")
dev.off()
```
Annotate Week 9 differentially methylated sites:
```R
library("genomation")
library("methylKit")

#Read in data:
load(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9diffmeth.RData")
load(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9diffmeth.RData")
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed12",remove.unusual=FALSE)

#Promoter/exon/introns:
#Annotate differentially methylated sites with promoter/exon/intron information.
AT.Diff.25p.anot <- annotateWithGeneParts(target = as(AT.Diff.25p,"GRanges"),
                                       feature = refseq_anot)
NB.Diff.25p.anot <- annotateWithGeneParts(target = as(NB.Diff.25p,"GRanges"),
                                       feature = refseq_anot)

#Report distance to the nearest Transcription Start Site, the target.row column in the output indicates the row number in the initial target set
AT.dist_tss <- getAssociationWithTSS(AT.Diff.25p.anot)
head(AT.dist_tss)

NB.dist_tss <- getAssociationWithTSS(NB.Diff.25p.anot)
head(NB.dist_tss)

#Report whether the differentially methylated CpGs are within promoters,introns or exons; the order is the same as the target set
ATpie <- getMembers(AT.Diff.25p.anot)
NBpie <- getMembers(NB.Diff.25p.anot)
AT.dist_tss <- cbind(AT.dist_tss, ATpie)
NB.dist_tss <- cbind(NB.dist_tss, NBpie)

#Summarize for all differentially methylated CpGs
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/CpGs_in_gene_region_AT9.png")
plotTargetAnnotation(AT.Diff.25p.anot, main = "Differential Methylation Annotation BR vs AT Week 9")
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/CpGs_in_gene_region_NB9.png")
plotTargetAnnotation(NB.Diff.25p.anot, main = "Differential Methylation Annotation BR vs NB Week 9")
dev.off()

#CpG islands:
#Load the CpG info
cpg_anot <- readFeatureFlank("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)
NB.diffCpGann <- annotateWithFeatureFlank(as(NB.Diff.25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")
AT.diffCpGann <- annotateWithFeatureFlank(as(AT.Diff.25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")

#Do CpG belong to a CpG Island or Shore
NBissh <- getMembers(NB.diffCpGann)
head(NBissh)

ATissh <- getMembers(AT.diffCpGann)
head(ATissh)

AT.dist_tss <- cbind(AT.dist_tss, ATissh)
NB.dist_tss <- cbind(NB.dist_tss, NBissh)

write.table(NB.dist_tss, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_NB9.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.dist_tss, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_AT9.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Summarize for all differentially methylated CpGs
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/Island_shore_NB9.png")
plotTargetAnnotation(NB.diffCpGann, main = "Differential Methylation Annotation NB9")
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/Island_shore_AT9.png")
plotTargetAnnotation(AT.diffCpGann, main = "Differential Methylation Annotation AT9")
dev.off()
```
Annotate Week 25 differentially methylated sites:
```R
library("genomation")
library("methylKit")

#Read in data:
load(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25diffmeth.RData")
load(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25diffmeth.RData")
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed12",remove.unusual=FALSE)

#Promoter/exon/introns:
#Annotate differentially methylated sites with promoter/exon/intron information.
AT.Diff.25p.anot <- annotateWithGeneParts(target = as(AT.Diff.25p,"GRanges"),
                                       feature = refseq_anot)
NB.Diff.25p.anot <- annotateWithGeneParts(target = as(NB.Diff.25p,"GRanges"),
                                       feature = refseq_anot)

#Report distance to the nearest Transcription Start Site, the target.row column in the output indicates the row number in the initial target set
AT.dist_tss <- getAssociationWithTSS(AT.Diff.25p.anot)
head(AT.dist_tss)

NB.dist_tss <- getAssociationWithTSS(NB.Diff.25p.anot)
head(NB.dist_tss)

#Report whether the differentially methylated CpGs are within promoters,introns or exons; the order is the same as the target set
ATpie <- getMembers(AT.Diff.25p.anot)
NBpie <- getMembers(NB.Diff.25p.anot)
AT.dist_tss <- cbind(AT.dist_tss, ATpie)
NB.dist_tss <- cbind(NB.dist_tss, NBpie)

#Summarize for all differentially methylated CpGs
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/CpGs_in_gene_region_AT25.png")
plotTargetAnnotation(AT.Diff.25p.anot, main = "Differential Methylation Annotation BR vs AT Week 25")
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/CpGs_in_gene_region_NB25.png")
plotTargetAnnotation(NB.Diff.25p.anot, main = "Differential Methylation Annotation BR vs NB Week 25")
dev.off()

#CpG islands:
#Load the CpG info
cpg_anot <- readFeatureFlank("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)
NB.diffCpGann <- annotateWithFeatureFlank(as(NB.Diff.25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")
AT.diffCpGann <- annotateWithFeatureFlank(as(AT.Diff.25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")

#Do CpG belong to a CpG Island or Shore
NBissh <- getMembers(NB.diffCpGann)
head(NBissh)

ATissh <- getMembers(AT.diffCpGann)
head(ATissh)

AT.dist_tss <- cbind(AT.dist_tss, ATissh)
NB.dist_tss <- cbind(NB.dist_tss, NBissh)

write.table(NB.dist_tss, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_NB25.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AT.dist_tss, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_AT25.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Summarize for all differentially methylated CpGs
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/Island_shore_NB25.png")
plotTargetAnnotation(NB.diffCpGann, main = "Differential Methylation Annotation NB25")
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/Island_shore_AT25.png")
plotTargetAnnotation(AT.diffCpGann, main = "Differential Methylation Annotation AT25")
dev.off()
```
Combine annotations with DMCs:
```bash
df1 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_AT1.txt", header = TRUE, sep = "\t")
df2 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_Diff_25p.txt", header = TRUE, sep = "\t")
df3 <- cbind(df2, df1)
write.table(df3, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_Diff_25p_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df1 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_NB1.txt", header = TRUE, sep = "\t")
df2 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_Diff_25p.txt", header = TRUE, sep = "\t")
df3 <- cbind(df2, df1)
write.table(df3, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_Diff_25p_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)


df1 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_AT3.txt", header = TRUE, sep = "\t")
df2 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_Diff_25p.txt", header = TRUE, sep = "\t")
df3 <- cbind(df2, df1)
write.table(df3, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_Diff_25p_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df1 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_NB3.txt", header = TRUE, sep = "\t")
df2 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_Diff_25p.txt", header = TRUE, sep = "\t")
df3 <- cbind(df2, df1)
write.table(df3, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_Diff_25p_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)


df1 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_AT6.txt", header = TRUE, sep = "\t")
df2 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_Diff_25p.txt", header = TRUE, sep = "\t")
df3 <- cbind(df2, df1)
write.table(df3, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_Diff_25p_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df1 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_NB6.txt", header = TRUE, sep = "\t")
df2 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_Diff_25p.txt", header = TRUE, sep = "\t")
df3 <- cbind(df2, df1)
write.table(df3, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_Diff_25p_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)


df1 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_AT9.txt", header = TRUE, sep = "\t")
df2 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_Diff_25p.txt", header = TRUE, sep = "\t")
df3 <- cbind(df2, df1)
write.table(df3, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_Diff_25p_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df1 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_NB9.txt", header = TRUE, sep = "\t")
df2 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_Diff_25p.txt", header = TRUE, sep = "\t")
df3 <- cbind(df2, df1)
write.table(df3, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_Diff_25p_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)


df1 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_AT25.txt", header = TRUE, sep = "\t")
df2 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_Diff_25p.txt", header = TRUE, sep = "\t")
df3 <- cbind(df2, df1)
write.table(df3, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_Diff_25p_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df1 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_NB25.txt", header = TRUE, sep = "\t")
df2 <- read.table("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_Diff_25p.txt", header = TRUE, sep = "\t")
df3 <- cbind(df2, df1)
write.table(df3, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_Diff_25p_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```
#### Sliding window <a name="26"></a>

DMRs - sliding window - week 1
```R
library("genomation")
library("methylKit")
#Read in the methylation ratio files
NB.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_1/bsmap/CpG_BR1_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_2/bsmap/CpG_BR1_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_3/bsmap/CpG_BR1_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_1/bsmap/CpG_NB1_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_2/bsmap/CpG_NB1_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_3/bsmap/CpG_NB1_E2_3_bsmap_ratios_filtered.txt")
NB.week1=methRead(NB.file.list, 
    sample.id=list("BR1_E2_1","BR1_E2_2","BR1_E2_3","NB1_E2_1","NB1_E2_2","NB1_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
NB.week1
AT.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_1/bsmap/CpG_BR1_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_2/bsmap/CpG_BR1_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_3/bsmap/CpG_BR1_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT1_E2_1/bsmap/CpG_AT1_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT1_E2_2/bsmap/CpG_AT1_E2_2_bsmap_ratios_filtered.txt")
AT.week1=methRead(AT.file.list, 
    sample.id=list("BR1_E2_1","BR1_E2_2","BR1_E2_3","AT1_E2_1","AT1_E2_2"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
AT.week1

#Read in annotation info
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed12",remove.unusual=FALSE)
cpg_anot <- readFeatureFlank("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)

##################################################################################################################

#Group methylation count by sliding window region:
NB.tiles <- tileMethylCounts(NB.week1,win.size=1000,step.size=1000,cov.bases = 10)

#Filter and normalise
NB.tiles.filt <- filterByCoverage(NB.tiles,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
NB.tiles.filt.norm <- normalizeCoverage(NB.tiles.filt, method = "median")
NB.meth.tiles <- unite(NB.tiles.filt.norm, destrand=FALSE)
NB.meth.tiles
NB.diff.tiles <- calculateDiffMeth(NB.meth.tiles,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
NB.diff.tiles
save(NB.diff.tiles, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB1_diffmeth_windowed.RData")

#################################################################################################################

#Group methylation count by sliding window region:
AT.tiles <- tileMethylCounts(AT.week1,win.size=1000,step.size=1000,cov.bases = 10)

#Filter and normalise
AT.tiles.filt <- filterByCoverage(AT.tiles,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
AT.tiles.filt.norm <- normalizeCoverage(AT.tiles.filt, method = "median")
AT.meth.tiles <- unite(AT.tiles.filt.norm, destrand=FALSE)
AT.meth.tiles
AT.diff.tiles <- calculateDiffMeth(AT.meth.tiles,
                            treatment=c(0,0,0,1,1),
                            overdispersion = "MN",
                            adjust="BH")
AT.diff.tiles
save(AT.diff.tiles, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT1_diffmeth_windowed.RData")

#################################################################################################################

# Rank by significance
NB.diff.tiles <- NB.diff.tiles[order(NB.diff.tiles$qvalue),]
# get all differentially methylated regions
NB.diff.tiles.25p <- getMethylDiff(NB.diff.tiles,
                        difference=25,
                        qvalue=0.01)

NB.diff.tiles.25p
#0 rows - no differentially methylated regions

##################################################################################################################

# Rank by significance
AT.diff.tiles <- AT.diff.tiles[order(AT.diff.tiles$qvalue),]
# get all differentially methylated regions
AT.diff.tiles.25p <- getMethylDiff(AT.diff.tiles,
                        difference=25,
                        qvalue=0.01)
AT.diff.tiles.25p
#0 rows - no differentially methylated regions
```
DMRs - sliding window - week 3
```R
library("genomation")
library("methylKit")
#Read in the methylation ratio files
NB.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_1/bsmap/CpG_BR3_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_2/bsmap/CpG_BR3_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_3/bsmap/CpG_BR3_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB3_E2_1/bsmap/CpG_NB3_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB3_E2_2/bsmap/CpG_NB3_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB3_E2_3/bsmap/CpG_NB3_E2_3_bsmap_ratios_filtered.txt")
NB.week1=methRead(NB.file.list, 
    sample.id=list("BR3_E2_1","BR3_E2_2","BR3_E2_3","NB3_E2_1","NB3_E2_2","NB3_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
NB.week1
AT.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_1/bsmap/CpG_BR3_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_2/bsmap/CpG_BR3_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_3/bsmap/CpG_BR3_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT3_E2_1/bsmap/CpG_AT3_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT3_E2_2/bsmap/CpG_AT3_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT3_E2_3/bsmap/CpG_AT3_E2_3_bsmap_ratios_filtered.txt")
AT.week1=methRead(AT.file.list, 
    sample.id=list("BR3_E2_1","BR3_E2_2","BR3_E2_3","AT3_E2_1","AT3_E2_2","AT3_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
AT.week1

#Read in annotation info
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed12",remove.unusual=FALSE)
cpg_anot <- readFeatureFlank("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)

##################################################################################################################

#Group methylation count by sliding window region:
NB.tiles <- tileMethylCounts(NB.week1,win.size=1000,step.size=1000,cov.bases = 10)

#Filter and normalise
NB.tiles.filt <- filterByCoverage(NB.tiles,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
NB.tiles.filt.norm <- normalizeCoverage(NB.tiles.filt, method = "median")
NB.meth.tiles <- unite(NB.tiles.filt.norm, destrand=FALSE)
NB.meth.tiles
NB.diff.tiles <- calculateDiffMeth(NB.meth.tiles,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
NB.diff.tiles
save(NB.diff.tiles, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB3_diffmeth_windowed.RData")

#################################################################################################################

#Group methylation count by sliding window region:
AT.tiles <- tileMethylCounts(AT.week1,win.size=1000,step.size=1000,cov.bases = 10)

#Filter and normalise
AT.tiles.filt <- filterByCoverage(AT.tiles,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
AT.tiles.filt.norm <- normalizeCoverage(AT.tiles.filt, method = "median")
AT.meth.tiles <- unite(AT.tiles.filt.norm, destrand=FALSE)
AT.meth.tiles
AT.diff.tiles <- calculateDiffMeth(AT.meth.tiles,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
AT.diff.tiles
save(AT.diff.tiles, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT3_diffmeth_windowed.RData")

#################################################################################################################

# Rank by significance
NB.diff.tiles <- NB.diff.tiles[order(NB.diff.tiles$qvalue),]
# get all differentially methylated regions
NB.diff.tiles.25p <- getMethylDiff(NB.diff.tiles,
                        difference=25,
                        qvalue=0.01)

NB.diff.tiles.25p
#0 rows - no differentially methylated regions

##################################################################################################################
<<<<<<< HEAD

# Rank by significance
AT.diff.tiles <- AT.diff.tiles[order(AT.diff.tiles$qvalue),]
# get all differentially methylated regions
AT.diff.tiles.25p <- getMethylDiff(AT.diff.tiles,
                        difference=25,
                        qvalue=0.01)
AT.diff.tiles.25p
#0 rows - no differentially methylated regions
```
DMRs - sliding window - week 6
```R
library("genomation")
library("methylKit")
#Read in the methylation ratio files
NB.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_1/bsmap/CpG_BR6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_2/bsmap/CpG_BR6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_3/bsmap/CpG_BR6_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_1/bsmap/CpG_NB6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_2/bsmap/CpG_NB6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_3/bsmap/CpG_NB6_E2_3_bsmap_ratios_filtered.txt")
NB.week1=methRead(NB.file.list, 
    sample.id=list("BR6_E2_1","BR6_E2_2","BR6_E2_3","NB6_E2_1","NB6_E2_2","NB6_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
NB.week1
AT.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_1/bsmap/CpG_BR6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_2/bsmap/CpG_BR6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_3/bsmap/CpG_BR6_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_1/bsmap/CpG_AT6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_2/bsmap/CpG_AT6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_3/bsmap/CpG_AT6_E2_3_bsmap_ratios_filtered.txt")
AT.week1=methRead(AT.file.list, 
    sample.id=list("BR6_E2_1","BR6_E2_2","BR6_E2_3","AT6_E2_1","AT6_E2_2","AT6_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
AT.week1

#Read in annotation info
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed12",remove.unusual=FALSE)
cpg_anot <- readFeatureFlank("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)

##################################################################################################################

#Group methylation count by sliding window region:
NB.tiles <- tileMethylCounts(NB.week1,win.size=1000,step.size=1000,cov.bases = 10)

#Filter and normalise
NB.tiles.filt <- filterByCoverage(NB.tiles,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
NB.tiles.filt.norm <- normalizeCoverage(NB.tiles.filt, method = "median")
NB.meth.tiles <- unite(NB.tiles.filt.norm, destrand=FALSE)
NB.meth.tiles
NB.diff.tiles <- calculateDiffMeth(NB.meth.tiles,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
NB.diff.tiles
save(NB.diff.tiles, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_diffmeth_windowed.RData")

#################################################################################################################

#Group methylation count by sliding window region:
AT.tiles <- tileMethylCounts(AT.week1,win.size=1000,step.size=1000,cov.bases = 10)

#Filter and normalise
AT.tiles.filt <- filterByCoverage(AT.tiles,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
AT.tiles.filt.norm <- normalizeCoverage(AT.tiles.filt, method = "median")
AT.meth.tiles <- unite(AT.tiles.filt.norm, destrand=FALSE)
AT.meth.tiles
AT.diff.tiles <- calculateDiffMeth(AT.meth.tiles,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
AT.diff.tiles
save(AT.diff.tiles, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_diffmeth_windowed.RData")

#################################################################################################################

# Rank by significance
NB.diff.tiles <- NB.diff.tiles[order(NB.diff.tiles$qvalue),]
# get all differentially methylated regions
NB.diff.tiles.25p <- getMethylDiff(NB.diff.tiles,
                        difference=25,
                        qvalue=0.01)

=======

# Rank by significance
AT.diff.tiles <- AT.diff.tiles[order(AT.diff.tiles$qvalue),]
# get all differentially methylated regions
AT.diff.tiles.25p <- getMethylDiff(AT.diff.tiles,
                        difference=25,
                        qvalue=0.01)
AT.diff.tiles.25p
#0 rows - no differentially methylated regions
```
DMRs - sliding window - week 6
```R
library("genomation")
library("methylKit")
#Read in the methylation ratio files
NB.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_1/bsmap/CpG_BR6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_2/bsmap/CpG_BR6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_3/bsmap/CpG_BR6_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_1/bsmap/CpG_NB6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_2/bsmap/CpG_NB6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_3/bsmap/CpG_NB6_E2_3_bsmap_ratios_filtered.txt")
NB.week1=methRead(NB.file.list, 
    sample.id=list("BR6_E2_1","BR6_E2_2","BR6_E2_3","NB6_E2_1","NB6_E2_2","NB6_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
NB.week1
AT.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_1/bsmap/CpG_BR6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_2/bsmap/CpG_BR6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_3/bsmap/CpG_BR6_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_1/bsmap/CpG_AT6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_2/bsmap/CpG_AT6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_3/bsmap/CpG_AT6_E2_3_bsmap_ratios_filtered.txt")
AT.week1=methRead(AT.file.list, 
    sample.id=list("BR6_E2_1","BR6_E2_2","BR6_E2_3","AT6_E2_1","AT6_E2_2","AT6_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
AT.week1

#Read in annotation info
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed12",remove.unusual=FALSE)
cpg_anot <- readFeatureFlank("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)

##################################################################################################################

#Group methylation count by sliding window region:
NB.tiles <- tileMethylCounts(NB.week1,win.size=1000,step.size=1000,cov.bases = 10)

#Filter and normalise
NB.tiles.filt <- filterByCoverage(NB.tiles,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
NB.tiles.filt.norm <- normalizeCoverage(NB.tiles.filt, method = "median")
NB.meth.tiles <- unite(NB.tiles.filt.norm, destrand=FALSE)
NB.meth.tiles
NB.diff.tiles <- calculateDiffMeth(NB.meth.tiles,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
NB.diff.tiles
save(NB.diff.tiles, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB6_diffmeth_windowed.RData")

#################################################################################################################

#Group methylation count by sliding window region:
AT.tiles <- tileMethylCounts(AT.week1,win.size=1000,step.size=1000,cov.bases = 10)

#Filter and normalise
AT.tiles.filt <- filterByCoverage(AT.tiles,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
AT.tiles.filt.norm <- normalizeCoverage(AT.tiles.filt, method = "median")
AT.meth.tiles <- unite(AT.tiles.filt.norm, destrand=FALSE)
AT.meth.tiles
AT.diff.tiles <- calculateDiffMeth(AT.meth.tiles,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
AT.diff.tiles
save(AT.diff.tiles, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT6_diffmeth_windowed.RData")

#################################################################################################################

# Rank by significance
NB.diff.tiles <- NB.diff.tiles[order(NB.diff.tiles$qvalue),]
# get all differentially methylated regions
NB.diff.tiles.25p <- getMethylDiff(NB.diff.tiles,
                        difference=25,
                        qvalue=0.01)

>>>>>>> 9bb487a11dfffd64de11be5c3e2043aee326f5d0
NB.diff.tiles.25p
#0 rows - no differentially methylated regions

##################################################################################################################

# Rank by significance
AT.diff.tiles <- AT.diff.tiles[order(AT.diff.tiles$qvalue),]
# get all differentially methylated regions
AT.diff.tiles.25p <- getMethylDiff(AT.diff.tiles,
                        difference=25,
                        qvalue=0.01)
AT.diff.tiles.25p
#0 rows - no differentially methylated regions
```
DMRs - sliding window - week 9
```R
library("genomation")
library("methylKit")
#Read in the methylation ratio files
NB.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_1/bsmap/CpG_BR9_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_2/bsmap/CpG_BR9_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_3/bsmap/CpG_BR9_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB9_E2_1/bsmap/CpG_NB9_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB9_E2_2/bsmap/CpG_NB9_E2_2_bsmap_ratios_filtered.txt")
NB.week1=methRead(NB.file.list, 
    sample.id=list("BR9_E2_1","BR9_E2_2","BR9_E2_3","NB9_E2_1","NB9_E2_2"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
NB.week1
AT.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_1/bsmap/CpG_BR9_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_2/bsmap/CpG_BR9_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_3/bsmap/CpG_BR9_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT9_E2_1/bsmap/CpG_AT9_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT9_E2_2/bsmap/CpG_AT9_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT9_E2_3/bsmap/CpG_AT9_E2_3_bsmap_ratios_filtered.txt")
AT.week1=methRead(AT.file.list, 
    sample.id=list("BR9_E2_1","BR9_E2_2","BR9_E2_3","AT9_E2_1","AT9_E2_2","AT9_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
AT.week1

#Read in annotation info
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed12",remove.unusual=FALSE)
cpg_anot <- readFeatureFlank("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)

##################################################################################################################

#Group methylation count by sliding window region:
NB.tiles <- tileMethylCounts(NB.week1,win.size=1000,step.size=1000,cov.bases = 10)

#Filter and normalise
NB.tiles.filt <- filterByCoverage(NB.tiles,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
NB.tiles.filt.norm <- normalizeCoverage(NB.tiles.filt, method = "median")
NB.meth.tiles <- unite(NB.tiles.filt.norm, destrand=FALSE)
NB.meth.tiles
NB.diff.tiles <- calculateDiffMeth(NB.meth.tiles,
                            treatment=c(0,0,0,1,1),
                            overdispersion = "MN",
                            adjust="BH")
NB.diff.tiles
save(NB.diff.tiles, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB9_diffmeth_windowed.RData")

#################################################################################################################

#Group methylation count by sliding window region:
AT.tiles <- tileMethylCounts(AT.week1,win.size=1000,step.size=1000,cov.bases = 10)

#Filter and normalise
AT.tiles.filt <- filterByCoverage(AT.tiles,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
AT.tiles.filt.norm <- normalizeCoverage(AT.tiles.filt, method = "median")
AT.meth.tiles <- unite(AT.tiles.filt.norm, destrand=FALSE)
AT.meth.tiles
AT.diff.tiles <- calculateDiffMeth(AT.meth.tiles,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
AT.diff.tiles
save(AT.diff.tiles, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT9_diffmeth_windowed.RData")

#################################################################################################################

# Rank by significance
NB.diff.tiles <- NB.diff.tiles[order(NB.diff.tiles$qvalue),]
# get all differentially methylated regions
NB.diff.tiles.25p <- getMethylDiff(NB.diff.tiles,
                        difference=25,
                        qvalue=0.01)

NB.diff.tiles.25p
#2 rows 

#Annotate
NB.diff.tiles.25p.ann <- annotateWithGeneParts(target = as(NB.diff.tiles.25p, "GRanges"), feature = refseq_anot)
NB.dist_tss <- getAssociationWithTSS(NB.Diff.25p.anot)
NBpie <- getMembers(NB.Diff.25p.anot)
NB.dist_tss <- cbind(NB.dist_tss, NBpie)
NB.diffCpGann <- annotateWithFeatureFlank(as(NB.Diff.25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")
NBissh <- getMembers(NB.diffCpGann)
NB.dist_tss <- cbind(NB.dist_tss, NBissh)

write.table(NB.dist_tss, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_NB9_windowed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##################################################################################################################

# Rank by significance
AT.diff.tiles <- AT.diff.tiles[order(AT.diff.tiles$qvalue),]
# get all differentially methylated regions
AT.diff.tiles.25p <- getMethylDiff(AT.diff.tiles,
                        difference=25,
                        qvalue=0.01)
AT.diff.tiles.25p
#0 rows - no differentially methylated regions
```
DMRs - sliding window - week 25
```R
library("genomation")
library("methylKit")
#Read in the methylation ratio files
NB.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_1/bsmap/CpG_BR25_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_2/bsmap/CpG_BR25_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_3/bsmap/CpG_BR25_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB25_E2_1/bsmap/CpG_NB25_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB25_E2_2/bsmap/CpG_NB25_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB25_E2_3/bsmap/CpG_NB25_E2_3_bsmap_ratios_filtered.txt")
NB.week1=methRead(NB.file.list, 
    sample.id=list("BR25_E2_1","BR25_E2_2","BR25_E2_3","NB25_E2_1","NB25_E2_2","NB25_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
NB.week1
AT.file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_1/bsmap/CpG_BR25_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_2/bsmap/CpG_BR25_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_3/bsmap/CpG_BR25_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT25_E2_1/bsmap/CpG_AT25_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT25_E2_2/bsmap/CpG_AT25_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT25_E2_3/bsmap/CpG_AT25_E2_3_bsmap_ratios_filtered.txt")
AT.week1=methRead(AT.file.list, 
    sample.id=list("BR25_E2_1","BR25_E2_2","BR25_E2_3","AT25_E2_1","AT25_E2_2","AT25_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
AT.week1

#Read in annotation info
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed12",remove.unusual=FALSE)
cpg_anot <- readFeatureFlank("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.cpg.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)

##################################################################################################################

#Group methylation count by sliding window region:
NB.tiles <- tileMethylCounts(NB.week1,win.size=1000,step.size=1000,cov.bases = 10)

#Filter and normalise
NB.tiles.filt <- filterByCoverage(NB.tiles,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
NB.tiles.filt.norm <- normalizeCoverage(NB.tiles.filt, method = "median")
NB.meth.tiles <- unite(NB.tiles.filt.norm, destrand=FALSE)
NB.meth.tiles
NB.diff.tiles <- calculateDiffMeth(NB.meth.tiles,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
NB.diff.tiles
save(NB.diff.tiles, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/NB25_diffmeth_windowed.RData")

#################################################################################################################

#Group methylation count by sliding window region:
AT.tiles <- tileMethylCounts(AT.week1,win.size=1000,step.size=1000,cov.bases = 10)

#Filter and normalise
AT.tiles.filt <- filterByCoverage(AT.tiles,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
AT.tiles.filt.norm <- normalizeCoverage(AT.tiles.filt, method = "median")
AT.meth.tiles <- unite(AT.tiles.filt.norm, destrand=FALSE)
AT.meth.tiles
AT.diff.tiles <- calculateDiffMeth(AT.meth.tiles,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
AT.diff.tiles
save(AT.diff.tiles, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/AT25_diffmeth_windowed.RData")

#################################################################################################################

# Rank by significance
NB.diff.tiles <- NB.diff.tiles[order(NB.diff.tiles$qvalue),]
# get all differentially methylated regions
NB.diff.tiles.25p <- getMethylDiff(NB.diff.tiles,
                        difference=25,
                        qvalue=0.01)

NB.diff.tiles.25p
#0 rows - no differentially methylated regions

##################################################################################################################

# Rank by significance
AT.diff.tiles <- AT.diff.tiles[order(AT.diff.tiles$qvalue),]
# get all differentially methylated regions
AT.diff.tiles.25p <- getMethylDiff(AT.diff.tiles,
                        difference=25,
                        qvalue=0.01)
AT.diff.tiles.25p
#1 row

#Annotate
AT.diff.tiles.25p.ann <- annotateWithGeneParts(target = as(AT.diff.tiles.25p, "GRanges"), feature = refseq_anot)
AT.dist_tss <- getAssociationWithTSS(AT.Diff.25p.anot)
ATpie <- getMembers(AT.Diff.25p.anot)
AT.dist_tss <- cbind(AT.dist_tss, NBpie)
AT.diffCpGann <- annotateWithFeatureFlank(as(AT.Diff.25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")
ATissh <- getMembers(AT.diffCpGann)
AT.dist_tss <- cbind(AT.dist_tss, NBissh)

write.table(AT.dist_tss, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/methylkit/tss_AT25_windowed.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```

### Custom <a name="27"></a>
Common approaches for differential methylation analysis are Bsmooth, Methylkit or a custom approach to define differentially methylated regions (DMRs) DOI:10.1093/bib/bbx077

Custom
```bash

```




## Check for transcriptional changes over the course of the experiments 1 & 2 <a name="28"></a>
### QC  <a name="29"></a>
#### fastqc  <a name="30"></a>

Experiment 1:
```bash
for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/*/*/*.fq.gz); do
    OutDir=$(dirname $Reads)/fastqc
    OutFile=$(basename $Reads | sed 's@.fq.gz@@g')
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_fastqc.sh $Reads $OutDir $OutFile
done
```
Experiment 2:
```bash
for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/*/*/*.fq.gz); do
    OutDir=$(dirname $Reads)/fastqc
    OutFile=$(basename $Reads | sed 's@.fq.gz@@g')
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_fastqc.sh $Reads $OutDir $OutFile
done
```
### Trimming <a name="31"></a>

Trimming ommitted for alignment to genome approach as per doi: 10.1093/nargab/lqaa068 and doi: 10.1186/s12859-016-0956-2 - 'soft clipping' when using alignment based mapping approach should make trimming unnessessary - use <--quantTranscriptomeBan Singleend> with STAR.

Trim for trinity:
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/*/*/); do
    sample=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    OutDir=$(echo $ReadDir | sed 's@raw_data@dna_qc@g')trim_galore
    OutFile=${sample}_trimmed
    Quality=20
    Length=50
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir -p $OutDir
    sbatch $ProgDir/run_trim_galore.sh $OutDir $OutFile $Quality $Length $Fread $Rread $Fread2 $Rread2 $Fread3 $Rread3 $Fread4 $Rread4 $Fread5 $Rread5
done 
#57715374-57715448

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/*/*/); do
    sample=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    OutDir=$(echo $ReadDir | sed 's@raw_data@dna_qc@g')trim_galore
    OutFile=${sample}_trimmed
    Quality=20
    Length=50
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir -p $OutDir
    sbatch $ProgDir/run_trim_galore.sh $OutDir $OutFile $Quality $Length $Fread $Rread $Fread2 $Rread2 $Fread3 $Rread3 $Fread4 $Rread4 $Fread5 $Rread5
done 
#57715449-57715528

#George has already trimmed his reads:
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Mpersicae_organ_RNAseq/trim_galore
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/George_Seddon/tissue_specific_RNAseq_all_data/trimmed_reads/*1.fq.gz); do
ID=$(basename $file | cut -d '_' -f1,2)
echo $ID
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Mpersicae_organ_RNAseq/trim_galore/$ID
ln -s $file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Mpersicae_organ_RNAseq/trim_galore/$ID/.
ln -s $(echo $file | sed 's@_1_val_1.fq.gz@_2_val_2.fq.gz@g') /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Mpersicae_organ_RNAseq/trim_galore/$ID/.
done

#The 9 host swap data has already been trimmed:
for dir in $(ls -d /jic/research-groups/Saskia-Hogenhout/reads/RNASeq/Myzus_persicae_O_9_species_host_swap/reads/trimmed/*/); do
    ID=$(echo $dir | rev | cut -d '/' -f2 | rev)
    host=$(echo $dir | rev | cut -d '/' -f2 | rev| sed 's@[0-9]@@g')
    mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Myzus_persicae_O_9_species_host_swap/${host}/${ID}/
    ln -s ${dir}/* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Myzus_persicae_O_9_species_host_swap/${host}/${ID}/.
done
```

### Mapping <a name="32"></a>
We don't have a reference transcriptome for M.persicae for alignment free approach.

#### Alignment free approach <a name="33"></a>

Could generate de novo transcriptome with trinity. Trinity could find novel transcritps from George's data however discovering novel transcripts is not the focus of this study.

##### Trinity <a name="34"></a>

It is unclear whether the 'soft clipping' logic applies to trinity, will therefore use trinity trimming function.

On trinity: The Inchworm and Chrysalis steps can be memory intensive. A basic recommendation is to have ~1G of RAM per ~1M pairs of Illumina reads. Simpler transcriptomes (lower eukaryotes) require less memory than more complex transcriptomes such as from vertebrates. Trinity can also require hundreds of GB of disk space available and can generate many thousands of intermediate files during the run. However, the final output files generated are few and are often relatively small (MB rather than many GB). It's good to have a temporary workspace available with sufficient disk space to use during the job execution. The entire process can require ~1/2 hour to one hour per million pairs of reads.

Each fastq file from the host swap experiments contains ~18M reads requiring ~1500GB or 32 days for one run of the experiment if the above is correct, unclear if this estimate is cpu time or real time. --normalize_reads parameter should improve speed - as of sept 2016 this is on by default: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Insilico-Normalization. The largest NBI node has 4030GB memory, this could maybe handle experiment 1, experiment 2 and George's data if runtime is not an issue. Tinity manual suggests there are certain steps in the pipeline that cannot be memory limited... I beleive the inchworm step is the problem.

To greatly lessen memory requirements, include the option --min_kmer_cov 2, in which case no uniquely occurring kmer will be assayed, and --normalize_by_read_set to perform the initial in silico normalization step on each pair of fastq files separately rather than combining them all into one large read set.

```bash
#Experiment 1
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Archana_Dec2020/*/*/trim_galore/*.fq.gz | wc -l #150
#Experiment 2
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Archana/*/*/trim_galore/*.fq.gz | wc -l #160 - this contains BR0, experiment 1 does not
#9 host swap reads:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Myzus_persicae_O_9_species_host_swap/*/*/*.fq.gz | wc -l #90
#George's organ reads:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Mpersicae_organ_RNAseq/trim_galore/*/*.fq.gz | wc -l #78

forward=$(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Archana_Dec2020/*/*/trim_galore/**1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Archana/*/*/trim_galore/*1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Myzus_persicae_O_9_species_host_swap/*/*/*1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Mpersicae_organ_RNAseq/trim_galore/*/*1.fq.gz | tr '\n' ',' | sed 's/,$//')
reverse=$(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Archana_Dec2020/*/*/trim_galore/*2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Archana/*/*/trim_galore/*2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Myzus_persicae_O_9_species_host_swap/*/*/*2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/RNA_Seq/Mpersicae_organ_RNAseq/trim_galore/*/*2.fq.gz | tr '\n' ',' | sed 's/,$//') 

#Trinity is using a huge amount of hard drive space + appears to be amking 3 trimmed files for every input?.
source package 09da5776-7777-44d1-9fac-4c372b38fd37
Trinity --full_cleanup --seqType fq --CPU 64 --max_memory 4030G --min_kmer_cov 2 --normalize_by_read_set --verbose \
--left $forward \
--right $reverse \
--output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/M_persicae_trinity_transcriptome 2>&1 >> logs/trinity_run_log.txt
#57695683,57722987

#Trinity reaches 99.994% completed, but then errors referencing Java heap space occur. It seem that this is likely due to contaminating or endosymbiont bacterial genome or a plasmids: https://github.com/trinityrnaseq/trinityrnaseq/issues/1220, https://github.com/trinityrnaseq/trinityrnaseq/issues/1325

#Find those transcripts which are not completing successfully (the 0.006%, 23/401049):
for file in $(cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/M_persicae_trinity_transcriptome/partitioned_reads.files.list); do
x=$(basename $file)
if ! grep -q "$x" /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/M_persicae_trinity_transcriptome/recursive_trinity.cmds.completed ; then
    echo $file >> temp_trinity_23.txt
fi
done

sleep 10800s
#Run phase 2 step for these 23 with increased heap space:
for file in $(cat temp_trinity_23.txt); do
#Trinity --single $file --output ${file}.out --CPU 1 --max_memory 1G --run_as_paired --seqType fa --trinity_complete --full_cleanup --min_kmer_cov 2 --verbose --bflyHeapSpaceMax 100G
sbatch ~/git_repos/Wrappers/NBI/temp4.sh $file
done
#57795755-77,57806588-604

#Trinity's checkpoints only work if the original command hasn't been modified, --FORCE must therefore be used to force inclusion of both the 401026 and 23 in the final transcriptome:
Trinity --full_cleanup --seqType fq --CPU 64 --max_memory 4030G --min_kmer_cov 2 --normalize_by_read_set --verbose --FORCE\
--left $forward \
--right $reverse \
--output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/M_persicae_trinity_transcriptome 2>&1 >> logs/trinity_run_log2.txt
#57765347 - check for errors after the 99.994% sticking point that will now be forced but have not been run in isolation?

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/M_persicae_trinity_transcriptome.Trinity.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/M_persicae_trinity_transcriptome.Trinity-force1.fasta
```

```bash





grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

cat gencode.vM23.transcripts.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome.fa.gz

salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index 

salmon index -t transcripts.fa -i transcripts_index --decoys decoys.txt -k 31

/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.CDS.fa
```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz

#### Alignment approach <a name="35"></a>

##### STAR <a name="36"></a>
```bash
#Make STAR index for M.persicae
mkdir /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/GenomeDir
source package 266730e5-6b24-4438-aecb-ab95f1940339
STAR --runMode genomeGenerate --runThreadN 32 \
--genomeDir /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/GenomeDir \
--genomeFastaFiles /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa \
--sjdbGTFfile /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.gtf --sjdbOverhang 150
#57690797

#Experiment 1:
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/*/*/); do
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    OutDir=$(echo $ReadDir | sed 's@raw_data@alignment@g')star
    OutFile=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
    IndexDir=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/GenomeDir
    GtfFile=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.gtf
    Stringency=Singleend
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch $ProgDir/run_rnastar.sh $OutDir $OutFile $Fread $Rread $IndexDir $GtfFile $Stringency 
done
#

#Experiment 2:
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/*/*/); do
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    OutDir=$(echo $ReadDir | sed 's@raw_data@alignment@g')star
    OutFile=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
    IndexDir=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/GenomeDir
    GtfFile=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.gtf
    Stringency=Singleend
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch $ProgDir/run_rnastar.sh $OutDir $OutFile $Fread $Rread $IndexDir $GtfFile $Stringency 
done
#
```



##### Stringtie <a name="37"></a>

Stringtie generation of transcriptome for alignment free appraoch requires alignments...

```bash

```
