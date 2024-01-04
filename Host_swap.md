## Collect data
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
## Check mutations in the aphid genomes over the course of the experiment with WGS
This is the purpose of the WGS data which was collected at the start and end of each experimental repeat from each host, ie.: generation 1 of aphids exposed to Br in experiments 1 and 2 (E1, E2), the generation 39 of aphids exposed to Br, At, Nb in E1, and the generation 25 of aphids exposed to Br, At, Nb in E2. 

#### QC - fastqc and qualimap
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

#BR1_E1 raw reads have average coverage of 23.8422
#BR39_E1 raw reads have average coverage of 23.2599
#AT39_E1 raw reads have average coverage of 20.6892
#NB39_E1 raw reads have average coverage of 21.0765

#BR1_E2 raw reads have average coverage of 20.0134
#BR25_E2 raw reads have average coverage of 23.7613
#AT25_E2 raw reads have average coverage of 19.7724
#NB25_E2 raw reads have average coverage of 21.6007
```
Trim
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
QC
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

#BR1_E1 raw reads have average coverage of 23.7446
#BR39_E1 raw reads have average coverage of 23.161
#AT39_E1 raw reads have average coverage of 20.6012
#NB39_E1 raw reads have average coverage of 20.9847

#BR1_E2 raw reads have average coverage of 19.9298
#BR25_E2 raw reads have average coverage of 23.6622
#AT25_E2 raw reads have average coverage of 19.69
#NB25_E2 raw reads have average coverage of 21.5096
```
Alignment to reference genome is performed as part of run_raw_read_qc.sh with bwa-mem

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
GATK indel realignment:
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

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/WGS/Archana_Feb2021/*/trim_galore/bwa-mem/*MarkDups.bam); do
Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_realign.sh $file $Reference 
done #57519857-64
```
Combined into a .vcf and called variants with bcftools:
```bash
source package 638df626-d658-40aa-80e5-14a275b7464b
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/WGS/Archana_Feb2021/*/trim_galore/bwa-mem/gatk/*realigned.bam > bamlist.txt
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/
bcftools mpileup -b /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/bamlist.txt --annotate AD,DP --fasta-ref /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa -O z -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/BR_AT_NB_hostswaps.vcf.gz

bcftools call --ploidy 2 -Oz -v -m -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/BR_AT_NB_hostswap.vcf.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/BR_AT_NB_hostswaps.vcf.gz
#57520372
```
Filter for genome mappability:

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

#Mappable bases (genmap): 328,601,147
#Masked bases (repeatmodeler): 93,097,267
#Callable bases (overlap of VCF and mappable genmap regions, minus repeatmodeler masked regions): 2,417,154

zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/Archana_Feb2021/gatk/genmap/BR_AT_NB_hostswap_callable.vcf.gz | wc -l #644,700
```
Plot genemap mappability to compare with qualimap mappability plots. These are not clear as the mappability regions are generally too dense to distinguish, nonetheless it appears that genmap successfully identifies poor mappability at the start of Chromosome 1.
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
Filter samples for missingness:
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

Filter for SNP quality:

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

Question: how many SNPs are there between BR1_E1 and BR39_E1, and BR1_E2 and BR25_E2? This is the number of SNPs without a host swap.
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

Question: did any SNPs occur in NB or AT host swaps in both experiments, but not in BR samples? If so then these would be candidates for causing host adaptation and could be investigated further - probably overkill to check this though.

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

## Check for epigenetic changes over the course of the experiment with WGBS
QC
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
Trim
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
QC
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/dna_qc/Myzus/persicae/WGBS/Archana_Mar2021/*/trim_galore/); do
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    OutDir=$(echo ${ReadDir})
    Reference_genome=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
    Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread 
done
#57554425-57554472
```
Alignment
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

#Output format: tab delimited txt file with the following columns:
#chr: Chromosome or scaffold name where the cytosine is located.
#pos: Position of the cytosine in the chromosome or scaffold.
#strand: Strand of DNA (either "+" or "-") where the cytosine is located.
#context: The sequence context of the cytosine. In this case, it is "CHH," which means the cytosine is followed by two non-cytosine bases in the 3' to 5' direction.
#ratio: Methylation ratio at the given cytosine position. It represents the proportion of methylated cytosines out of the total observed cytosines.
#eff_CT_count: Effective count of cytosines considered for calculating the methylation ratio. This count excludes certain cytosines based on specific criteria (e.g., filtering low-quality reads).
#C_count: Count of methylated cytosines.
#CT_count: Total count of cytosines (both methylated and unmethylated).
#rev_G_count: Count of guanines on the reverse strand corresponding to the cytosine.
#rev_GA_count: Count of guanine-adenine pairs on the reverse strand corresponding to the cytosine.
#CI_lower: Lower bound of the confidence interval for the methylation ratio.
#CI_upper: Upper bound of the confidence interval for the methylation ratio.

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/*/bsmap/*bsmap_ratios.txt); do
echo $file >> logs/bsmap_report.txt
cat $file | wc -l >> logs/bsmap_report.txt
done
```
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
#57758418-65
#57765410-13
#57768086-133
```
Samples BR0_E2_2 and BR0_E2_3 have significantly fewer methylation sites than other samples - ommitted from further analysis.
```bash
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/*/bsmap/*bsmap_ratios.txt | grep -v 'BR0_E2_2\|BR0_E2_3' > temp_file_list.txt
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
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
Universally un-methylated sites were removed:
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
#The number of sites which are totally unmethylated across all samples is only 690. Given this low number I will not bother to remove these sites.
```
Custom
```bash

```
Methylkit
```bash
#Seperate CpG methylation context Cs.
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/*/bsmap/*filtered.txt); do
    awk '{OFS="\t";if ($4 ~ /CG/) {print $1, $2, $3, $4, $5, int($6);}}' $file > $(dirname $file)/CpG_$(basename $file)
    cat $(dirname $file)/CpG_$(basename $file) | wc -l
done
#There are 13,682,374 CpG sites

#Create bed file for gene annotations
awk -F'\t' '$3 == "gene" {print $1 "\t" $4-1 "\t" $5 "\t" $9}' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3 > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/methylkit1.28.0.sif R
sbatch ~/git_repos/Wrappers/NBI/temp3.sh
#57758156-61
```
test:
```bash
for file in \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_1/bsmap/CpG_BR1_E2_1_bsmap_ratios_filtered.txt" \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_2/bsmap/CpG_BR1_E2_2_bsmap_ratios_filtered.txt" \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_3/bsmap/CpG_BR1_E2_3_bsmap_ratios_filtered.txt" \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_1/bsmap/CpG_NB1_E2_1_bsmap_ratios_filtered.txt" \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_2/bsmap/CpG_NB1_E2_2_bsmap_ratios_filtered.txt" \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_3/bsmap/CpG_NB1_E2_3_bsmap_ratios_filtered.txt" \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT1_E2_1/bsmap/CpG_AT1_E2_1_bsmap_ratios_filtered.txt" \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT1_E2_2/bsmap/CpG_AT1_E2_2_bsmap_ratios_filtered.txt"; do
    tail -n 100000 $file > ./temp_$(basename $file)
done
```
```R
library(methylKit)
#Read in the methylation ratio files
week1_file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_CpG_BR1_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_CpG_BR1_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_CpG_BR1_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_CpG_NB1_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_CpG_NB1_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_CpG_NB1_E2_3_bsmap_ratios_filtered.txt")
week1=methRead(week1_file.list, 
    sample.id=list("BR1_E2_1","BR1_E2_2","BR1_E2_3","NB1_E2_1","NB1_E2_2","NB1_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
week1
head(week1[[1]])
png("CpGhistogram1.png")
getMethylationStats(week1[[1]], plot=TRUE, both.strands=FALSE)
dev.off()

#Normalisation and filtering
week1.filt <- filterByCoverage(week1,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

week1.filt.norm <- normalizeCoverage(week1.filt, method = "median")
week1.meth <- unite(week1.filt.norm, destrand=FALSE)
week1.meth
# get percent methylation matrix
pm=percMethylation(week1.meth)
# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("sdshistogram1.png")
hist(sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
week1.meth <- week1.meth[sds > 2]
# This leaves us with this number of CpG sites
nrow(week1.meth)

#Plot data structure
png("correlation1.png")
getCorrelation(week1.meth,plot=TRUE)
dev.off()
png("dendogram1.png")
clusterSamples(week1.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/PCA1.png")
PCASamples(week1.meth)
dev.off()

#Identify differential methylation
#Subset for NB
Diff1 <- calculateDiffMeth(week1.meth,
                            treatment=c(0,0,0,1,1,1),
                            overdispersion = "MN",
                            adjust="BH")
Diff1
png("Volcano1.png")
plot(Diff1$meth.diff, -log10(Diff1$qvalue))
abline(v=0)
dev.off()

diffMethPerChr(Diff1) #Cannot plot figure, excluded all available chromosomes.

# get hyper methylated bases and order by qvalue
Diff1.25p.hyper <- getMethylDiff(Diff1,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
Diff1.25p.hyper <- Diff1.25p.hyper[order(Diff1.25p.hyper$qvalue), ]

# get hypo methylated bases and order by qvalue
Diff1.25p.hypo <- getMethylDiff(Diff1,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
Diff1.25p.hypo <- Diff1.25p.hypo[order(Diff1.25p.hypo$qvalue), ]

# get all differentially methylated bases and order by qvalue
Diff1.25p <- getMethylDiff(Diff1,
                        difference=25,
                        qvalue=0.01)
Diff1.25p <- Diff1.25p[order(Diff1.25p$qvalue), ]

#Annotate CpGs in genic regions
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed")
#Annotate hypermethylated CpGs ("target") with promoter/exon/intron
sink("Diff1p.hyper.anot.txt")
Diff1.25p.hyper.anot <- annotateWithGeneParts(target = as(Diff1p.hyper,"GRanges"),
                                       feature = refseq_anot)
sink()

# View the distance to the nearest Transcription Start Site; the target.row column in the output indicates the row number in the initial target set
dist_tss <- getAssociationWithTSS(Diff1.25p.hyper.anot)
head(dist_tss)

# See whether the differentially methylated CpGs are within promoters,introns or exons; the order is the same as the target set
getMembers(Diff1.25p.hyper.anot)

# This can also be summarized for all differentially methylated CpGs
plotTargetAnnotation(Diff1.25p.hyper.anot, main = "Differential Methylation Annotation")
```
```R
library(methylKit)

#################################################################################################################################################################
#Read in the methylation ratio files
week0_file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR0_E2_1/bsmap/CpG_BR0_E2_1_bsmap_ratios_filtered.txt")
week0=methRead(week0_file.list, 
    sample.id=list("BR0_E2_1"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
week0
head(week0[[1]])
getMethylationStats(week0[[1]], plot=TRUE, both.strands=FALSE)

#Normalisation and filtering
week0.filt <- filterByCoverage(week0,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

week0.filt.norm <- normalizeCoverage(week0.filt, method = "median")
week0.meth <- unite(week0.filt.norm, destrand=FALSE)
week0.meth
# get percent methylation matrix
pm=percMethylation(week0.meth)
# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("histogram.png", width = 800, height = 600, units = "px", res = 300)
hist(sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
week0.meth <- week0.meth[sds > 2]
# This leaves us with this number of CpG sites
nrow(week0.meth)

#Plot data structure
png("correlation.png", width = 800, height = 600, units = "px", res = 300)
getCorrelation(week0.meth,plot=TRUE)
dev.off()
png("dendogram.png", width = 800, height = 600, units = "px", res = 300)
clusterSamples(week0.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("PCA.png", width = 800, height = 600, units = "px", res = 300)
PCASamples(week0.meth)
dev.off()

#################################################################################################################################################################
#Read in the methylation ratio files
week1_file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_1/bsmap/CpG_BR1_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_2/bsmap/CpG_BR1_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_3/bsmap/CpG_BR1_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_1/bsmap/CpG_NB1_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_2/bsmap/CpG_NB1_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_3/bsmap/CpG_NB1_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT1_E2_1/bsmap/CpG_AT1_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT1_E2_2/bsmap/CpG_AT1_E2_2_bsmap_ratios_filtered.txt")
week1=methRead(week1_file.list, 
    sample.id=list("BR1_E2_1","BR1_E2_2","BR1_E2_3","NB1_E2_1","NB1_E2_2","NB1_E2_3","AT1_E2_1","AT1_E2_2"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1,2,2),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
week1
head(week1[[1]])
getMethylationStats(week1[[1]], plot=TRUE, both.strands=FALSE)

#Normalisation and filtering
week1.filt <- filterByCoverage(week1,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

week1.filt.norm <- normalizeCoverage(week1.filt, method = "median")
week1.meth <- unite(week1.filt.norm, destrand=FALSE)
week1.meth
# get percent methylation matrix
pm=percMethylation(week1.meth)
# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("histogram1.png", res = 1080)
hist(sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
week1.meth <- week1.meth[sds > 2]
# This leaves us with this number of CpG sites
nrow(week1.meth)

#Plot data structure
png("correlation1.png", res = 1080)
getCorrelation(week1.meth,plot=TRUE)
dev.off()
cat("bump")
png("dendogram1.png", res = 1080)
clusterSamples(week1.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("PCA1.png", res = 1080)
PCASamples(week1.meth)
dev.off()

#Identify differential methylation
Diff1 <- calculateDiffMeth(week1.meth,
                            treatment=c(0,0,0,1,1,1,2,2,2),
                            overdispersion = "MN",
                            adjust="BH")
Diff1
png("Volcano1.png", res = 1080)
plot(Diff1$meth.diff, -log10(Diff1$qvalue))
abline(v=0)
dev.off()

diffMethPerChr(Diff1)

# get hyper methylated bases and order by qvalue
Diff1.25p.hyper <- getMethylDiff(Diff1,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
Diff1.25p.hyper <- Diff1.25p.hyper[order(Diff1.25p.hyper$qvalue),]

# get hypo methylated bases and order by qvalue
Diff1.25p.hypo <- getMethylDiff(Diff1,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
Diff1.25p.hypo <- Diff1.25p.hypo[order(Diff1.25p.hypo$qvalue),]

# get all differentially methylated bases and order by qvalue
Diff1.25p <- getMethylDiff(Diff1,
                        difference=25,
                        qvalue=0.01)
Diff1.25p <- Diff25.25p[order(Diff1.25p$qvalue),]

#Annotate CpGs in genic regions
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed")
#Annotate hypermethylated CpGs ("target") with promoter/exon/intron
sink("Diff1p.hyper.anot.txt")
Diff1.25p.hyper.anot <- annotateWithGeneParts(target = as(Diff1p.hyper,"GRanges"),
                                       feature = refseq_anot)
sink()

# View the distance to the nearest Transcription Start Site; the target.row column in the output indicates the row number in the initial target set
dist_tss <- getAssociationWithTSS(Diff1.25p.hyper.anot)
head(dist_tss)

# See whether the differentially methylated CpGs are within promoters,introns or exons; the order is the same as the target set
getMembers(Diff1.25p.hyper.anot)

# This can also be summarized for all differentially methylated CpGs
plotTargetAnnotation(Diff1.25p.hyper.anot, main = "Differential Methylation Annotation")
#################################################################################################################################################################
#Read in the methylation ratio files
week3_file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_1/bsmap/CpG_BR3_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_2/bsmap/CpG_BR3_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR3_E2_3/bsmap/CpG_BR3_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB3_E2_1/bsmap/CpG_NB3_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB3_E2_2/bsmap/CpG_NB3_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB3_E2_3/bsmap/CpG_NB3_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT3_E2_1/bsmap/CpG_AT3_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT3_E2_2/bsmap/CpG_AT3_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT3_E2_3/bsmap/CpG_AT3_E2_3_bsmap_ratios_filtered.txt")
week3=methRead(week3_file.list, 
    sample.id=list("BR3_E2_1","BR3_E2_2","BR3_E2_3","NB3_E2_1","NB3_E2_2","NB3_E2_3","AT3_E2_1","AT3_E2_2","AT3_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1,2,2,2),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
week3
head(week3[[1]])
getMethylationStats(week3[[1]], plot=TRUE, both.strands=FALSE)

#Normalisation and filtering
week3.filt <- filterByCoverage(week3,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

week3.filt.norm <- normalizeCoverage(week3.filt, method = "median")
week3.meth <- unite(week3.filt.norm, destrand=FALSE)
week3.meth
# get percent methylation matrix
pm=percMethylation(week3.meth)
# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("histogram3.png", width = 800, height = 600, units = "px", res = 300)
hist(sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
week3.meth <- week3.meth[sds > 2]
# This leaves us with this number of CpG sites
nrow(week3.meth)

#Plot data structure
png("correlation3.png", width = 800, height = 600, units = "px", res = 300)
getCorrelation(week3.meth,plot=TRUE)
dev.off()
png("dendogram3.png", width = 800, height = 600, units = "px", res = 300)
clusterSamples(week3.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("PCA3.png", width = 800, height = 600, units = "px", res = 300)
PCASamples(week3.meth)
dev.off()

#################################################################################################################################################################
#Read in the methylation ratio files
week6_file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_1/bsmap/CpG_BR6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_2/bsmap/CpG_BR6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR6_E2_3/bsmap/CpG_BR6_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_1/bsmap/CpG_NB6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_2/bsmap/CpG_NB6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB6_E2_3/bsmap/CpG_NB6_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_1/bsmap/CpG_AT6_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_2/bsmap/CpG_AT6_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT6_E2_3/bsmap/CpG_AT6_E2_3_bsmap_ratios_filtered.txt")
week6=methRead(week6_file.list, 
    sample.id=list("BR6_E2_1","BR6_E2_2","BR6_E2_3","NB6_E2_1","NB6_E2_2","NB6_E2_3","AT6_E2_1","AT6_E2_2","AT6_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1,2,2,2),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
week6
head(week6[[1]])
getMethylationStats(week6[[1]], plot=TRUE, both.strands=FALSE)

#Normalisation and filtering
week6.filt <- filterByCoverage(week6,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

week6.filt.norm <- normalizeCoverage(week6.filt, method = "median")
week6.meth <- unite(week6.filt.norm, destrand=FALSE)
week6.meth
# get percent methylation matrix
pm=percMethylation(week6.meth)
# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("histogram6.png", width = 800, height = 600, units = "px", res = 300)
hist(sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
week6.meth <- week6.meth[sds > 2]
# This leaves us with this number of CpG sites
nrow(week6.meth)

#Plot data structure
png("correlation6.png", width = 800, height = 600, units = "px", res = 300)
getCorrelation(week6.meth,plot=TRUE)
dev.off()
png("dendogram6.png", width = 800, height = 600, units = "px", res = 300)
clusterSamples(week6.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("PCA6.png", width = 800, height = 600, units = "px", res = 300)
PCASamples(week6.meth)
dev.off()

#################################################################################################################################################################
#Read in the methylation ratio files
week9_file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_1/bsmap/CpG_BR9_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_2/bsmap/CpG_BR9_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR9_E2_3/bsmap/CpG_BR9_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB9_E2_1/bsmap/CpG_NB9_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB9_E2_2/bsmap/CpG_NB9_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT9_E2_1/bsmap/CpG_AT9_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT9_E2_2/bsmap/CpG_AT9_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT9_E2_3/bsmap/CpG_AT9_E2_3_bsmap_ratios_filtered.txt")
week9=methRead(week9_file.list, 
    sample.id=list("BR9_E2_1","BR9_E2_2","BR9_E2_3","NB9_E2_1","NB9_E2_2","AT9_E2_1","AT9_E2_2","AT9_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,2,2,2),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
week9
head(week9[[1]])
getMethylationStats(week9[[1]], plot=TRUE, both.strands=FALSE)

#Normalisation and filtering
week9.filt <- filterByCoverage(week9,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

week9.filt.norm <- normalizeCoverage(week9.filt, method = "median")
week9.meth <- unite(week9.filt.norm, destrand=FALSE)
week9.meth
# get percent methylation matrix
pm=percMethylation(week9.meth)
# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("histogram9.png", width = 800, height = 600, units = "px", res = 300)
hist(sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
week9.meth <- week9.meth[sds > 2]
# This leaves us with this number of CpG sites
nrow(week9.meth)

#Plot data structure
png("correlation9.png", width = 800, height = 600, units = "px", res = 300)
getCorrelation(week9.meth,plot=TRUE)
dev.off()
png("dendogram9.png", width = 800, height = 600, units = "px", res = 300)
clusterSamples(week9.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("PCA9.png", width = 800, height = 600, units = "px", res = 300)
PCASamples(week9.meth)
dev.off()

#################################################################################################################################################################
#Read in the methylation ratio files
week25_file.list <- list("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_1/bsmap/CpG_BR25_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_2/bsmap/CpG_BR25_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR25_E2_3/bsmap/CpG_BR25_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB25_E2_1/bsmap/CpG_NB25_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB25_E2_2/bsmap/CpG_NB25_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB25_E2_3/bsmap/CpG_NB25_E2_3_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT25_E2_1/bsmap/CpG_AT25_E2_1_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT25_E2_2/bsmap/CpG_AT25_E2_2_bsmap_ratios_filtered.txt",
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT25_E2_3/bsmap/CpG_AT25_E2_3_bsmap_ratios_filtered.txt")
week25=methRead(week25_file.list, 
    sample.id=list("BR25_E2_1","BR25_E2_2","BR25_E2_3","NB25_E2_1","NB25_E2_2","NB25_E2_3","AT25_E2_1","AT25_E2_2","AT25_E2_3"),
    assembly="O_v2",
    header=TRUE,
    treatment=c(0,0,0,1,1,1,2,2,2),
    mincov = 4,
    context="CpG",
    resolution="base",
    pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
week25
head(week25[[1]])
getMethylationStats(week25[[1]], plot=TRUE, both.strands=FALSE)

#Normalisation and filtering
week25.filt <- filterByCoverage(week25,
                      lo.count=4,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

week25.filt.norm <- normalizeCoverage(week25.filt, method = "median")
week25.meth <- unite(week25.filt.norm, destrand=FALSE)
week25.meth
# get percent methylation matrix
pm=percMethylation(week25.meth)
# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)
# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
png("histogram25.png", width = 800, height = 600, units = "px", res = 300)
hist(sds, breaks = 100)
dev.off()
# keep only CpG with standard deviations larger than 2%
week25.meth <- week25.meth[sds > 2]
# This leaves us with this number of CpG sites
nrow(week25.meth)

#Plot data structure
png("correlation25.png", width = 800, height = 600, units = "px", res = 300)
getCorrelation(week25.meth,plot=TRUE)
dev.off()
png("dendogram25.png", width = 800, height = 600, units = "px", res = 300)
clusterSamples(week25.meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
png("PCA25.png", width = 800, height = 600, units = "px", res = 300)
PCASamples(week25.meth)
dev.off()

#Identify differential methylation
Diff25 <- calculateDiffMeth(week25.meth,
                            treatment=c(0,0,0,1,1,1,2,2,2),
                            overdispersion = "MN",
                            adjust="BH")
Diff25
png("Volcano25.png", width = 800, height = 600, units = "px", res = 300)
plot(Diff25$meth.diff, -log10(Diff25$qvalue))
abline(v=0)
dev.off()

diffMethPerChr(Diff25)

# get hyper methylated bases and order by qvalue
Diff25.25p.hyper <- getMethylDiff(Diff25,
                              difference=25,
                              qvalue=0.01,
                              type="hyper")
Diff25.25p.hyper <- Diff25.25p.hyper[order(Diff25.25p.hyper$qvalue),]

# get hypo methylated bases and order by qvalue
Diff25.25p.hypo <- getMethylDiff(Diff25,
                             difference=25,
                             qvalue=0.01,
                             type="hypo")
Diff25.25p.hypo <- Diff25.25p.hypo[order(Diff25.25p.hypo$qvalue),]

# get all differentially methylated bases and order by qvalue
Diff25.25p <- getMethylDiff(Diff25,
                        difference=25,
                        qvalue=0.01)
Diff25.25p <- Diff25.25p[order(Diff25.25p$qvalue),]

#Annotate CpGs in genic regions
refseq_anot <- readTranscriptFeatures("/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.bed")
#Annotate hypermethylated CpGs ("target") with promoter/exon/intron
sink("Diff25p.hyper.anot.txt")
Diff25.25p.hyper.anot <- annotateWithGeneParts(target = as(Diff25p.hyper,"GRanges"),
                                       feature = refseq_anot)
sink()

# View the distance to the nearest Transcription Start Site; the target.row column in the output indicates the row number in the initial target set
dist_tss <- getAssociationWithTSS(Diff25.25p.hyper.anot)
head(dist_tss)

# See whether the differentially methylated CpGs are within promoters,introns or exons; the order is the same as the target set
getMembers(Diff25.25p.hyper.anot)

# This can also be summarized for all differentially methylated CpGs
plotTargetAnnotation(Diff25.25p.hyper.anot, main = "Differential Methylation Annotation")
```


BSsmooth
```bash
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/bsseq1.38.0.sif R
#57758169
```
test:
```bash
for file in \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_1/bismark/BR1_E2_1_bismark.deduplicated.CpG_report.txt" \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_2/bismark/BR1_E2_2_bismark.deduplicated.CpG_report.txt" \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/BR1_E2_3/bismark/BR1_E2_3_bismark.deduplicated.CpG_report.txt" \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_1/bismark/NB1_E2_1_bismark.deduplicated.CpG_report.txt" \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/NB1_E2_3/bismark/NB1_E2_3_bismark.deduplicated.CpG_report.txt" \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT1_E2_1/bismark/AT1_E2_1_bismark.deduplicated.CpG_report.txt" \
"/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/alignment/Myzus/persicae/WGBS/Archana_Mar2021/AT1_E2_2/bismark/AT1_E2_2_bismark.deduplicated.CpG_report.txt"; do
    tail -n 100000 $file > ./temp_$(basename $file)
done
```
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

sapply(assays(bsseq, withDimnames = FALSE), class)
print("Bsseq object:")
bsseq
pData(bsseq)

print("The average coverage of CpGs:")
round(colMeans(getCoverage(bsseq)), 1)
print("The number of CpGs:")
length(bsseq)
print("Number of CpGs which are covered by at least 1 read in all samples:")
sum(rowSums(getCoverage(bsseq) >= 1) == 7)
print("Number of CpGs with 0 coverage in all samples:")
sum(rowSums(getCoverage(bsseq)) == 0)


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

#Remove CpGs with coverage below 4 in all samples:
BS.cov <- getCoverage(bssmooth)
keepLoci.ex <- which(rowSums(BS.cov[, bsseq$treatment == "control"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "NB"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "AT"] >= 4) >= 3)
print("The number of CpGs with coverage >=4 in all samples:")
length(keepLoci.ex)
bssmooth <- bssmooth[keepLoci.ex,]


#Compute t-statistics based on smoothed whole-genome bisulfite sequencing data:
print("Compute t-stats vs AT:")
AT.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("AT1_E2_1", "AT1_E2_2", "AT1_E2_3"),  
                                    group2 = c("BR1_E2_1", "BR1_E2_2", "BR1_E2_3"),
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE)
AT.tstat
print("Compute t-stats vs NB:")
NB.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("NB1_E2_1", "NB1_E2_2", "NB1_E2_3"),  
                                    group2 = c("BR1_E2_1", "BR1_E2_2", "BR1_E2_3"),
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE)
NB.tstat


#Finding DMRs
print("Find DMRs vs NB:")
dmrs0 <- dmrFinder(AT.tstat, qcutoff = c(0.025, 0.975), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation between normal and cancers of at least 0.1.
AT.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(AT.dmrs)
head(AT.dmrs, n = 3)
write(AT.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT1_dmrs.txt")

print("Find DMRs vs AT:")
dmrs0 <- dmrFinder(AT.tstat, qcutoff = c(0.025, 0.975), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation between normal and cancers of at least 0.1.
NB.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(NB.dmrs)
head(NB.dmrs, n = 3)
write(NB.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB1_dmrs.txt")

#Plot the top 1000 DMRs
print("Plot top NB DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT1_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, AT.dmrs[1:1000,], extend = 5000, 
                addRegions = AT.dmrs)
dev.off()
print("Plot top AT DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB1_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, NB.dmrs[1:1000,], extend = 5000, 
                addRegions = NB.dmrs)
dev.off()

##################################################################################################################

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

sapply(assays(bsseq, withDimnames = FALSE), class)
print("Bsseq object:")
bsseq
pData(bsseq)

print("The average coverage of CpGs:")
round(colMeans(getCoverage(bsseq)), 1)
print("The number of CpGs:")
length(bsseq)
print("Number of CpGs which are covered by at least 1 read in all samples:")
sum(rowSums(getCoverage(bsseq) >= 1) == 7)
print("Number of CpGs with 0 coverage in all samples:")
sum(rowSums(getCoverage(bsseq)) == 0)


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

#Remove CpGs with coverage below 4 in all samples:
BS.cov <- getCoverage(bssmooth)
keepLoci.ex <- which(rowSums(BS.cov[, bsseq$treatment == "control"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "NB"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "AT"] >= 4) >= 3)
print("The number of CpGs with coverage >=4 in all samples:")
length(keepLoci.ex)
bssmooth <- bssmooth[keepLoci.ex,]


#Compute t-statistics based on smoothed whole-genome bisulfite sequencing data:
print("Compute t-stats vs AT:")
AT.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("AT3_E2_1", "AT3_E2_2", "AT3_E2_3"),  
                                    group2 = c("BR3_E2_1", "BR3_E2_2", "BR3_E2_3"),
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE)
AT.tstat
print("Compute t-stats vs NB:")
NB.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("NB3_E2_1", "NB3_E2_2", "NB3_E2_3"),  
                                    group2 = c("BR3_E2_1", "BR3_E2_2", "BR3_E2_3"),
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE)
NB.tstat


#Finding DMRs
print("Find DMRs vs NB:")
dmrs0 <- dmrFinder(AT.tstat, qcutoff = c(0.025, 0.975), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation between normal and cancers of at least 0.1.
AT.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(AT.dmrs)
head(AT.dmrs, n = 3)
write(AT.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT3_dmrs.txt")

print("Find DMRs vs AT:")
dmrs0 <- dmrFinder(AT.tstat, qcutoff = c(0.025, 0.975), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation between normal and cancers of at least 0.1.
NB.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(NB.dmrs)
head(NB.dmrs, n = 3)
write(NB.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB3_dmrs.txt")

#Plot the top 1000 DMRs
print("Plot top NB DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT3_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, AT.dmrs[1:1000,], extend = 5000, 
                addRegions = AT.dmrs)
dev.off()
print("Plot top AT DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB3_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, NB.dmrs[1:1000,], extend = 5000, 
                addRegions = NB.dmrs)
dev.off()

##################################################################################################################

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

sapply(assays(bsseq, withDimnames = FALSE), class)
print("Bsseq object:")
bsseq
pData(bsseq)

print("The average coverage of CpGs:")
round(colMeans(getCoverage(bsseq)), 1)
print("The number of CpGs:")
length(bsseq)
print("Number of CpGs which are covered by at least 1 read in all samples:")
sum(rowSums(getCoverage(bsseq) >= 1) == 7)
print("Number of CpGs with 0 coverage in all samples:")
sum(rowSums(getCoverage(bsseq)) == 0)


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

#Remove CpGs with coverage below 4 in all samples:
BS.cov <- getCoverage(bssmooth)
keepLoci.ex <- which(rowSums(BS.cov[, bsseq$treatment == "control"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "NB"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "AT"] >= 4) >= 3)
print("The number of CpGs with coverage >=4 in all samples:")
length(keepLoci.ex)
bssmooth <- bssmooth[keepLoci.ex,]


#Compute t-statistics based on smoothed whole-genome bisulfite sequencing data:
print("Compute t-stats vs AT:")
AT.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("AT6_E2_1", "AT6_E2_2", "AT6_E2_3"),  
                                    group2 = c("BR6_E2_1", "BR6_E2_2", "BR6_E2_3"),
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE)
AT.tstat
print("Compute t-stats vs NB:")
NB.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("NB6_E2_1", "NB6_E2_2", "NB6_E2_3"),  
                                    group2 = c("BR6_E2_1", "BR6_E2_2", "BR6_E2_3"),
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE)
NB.tstat


#Finding DMRs
print("Find DMRs vs NB:")
dmrs0 <- dmrFinder(AT.tstat, qcutoff = c(0.025, 0.975), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation between normal and cancers of at least 0.1.
AT.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(AT.dmrs)
head(AT.dmrs, n = 3)
write(AT.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT6_dmrs.txt")

print("Find DMRs vs AT:")
dmrs0 <- dmrFinder(AT.tstat, qcutoff = c(0.025, 0.975), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation between normal and cancers of at least 0.1.
NB.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(NB.dmrs)
head(NB.dmrs, n = 3)
write(NB.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB6_dmrs.txt")

#Plot the top 1000 DMRs
print("Plot top NB DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT6_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, AT.dmrs[1:1000,], extend = 5000, 
                addRegions = AT.dmrs)
dev.off()
print("Plot top AT DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB6_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, NB.dmrs[1:1000,], extend = 5000, 
                addRegions = NB.dmrs)
dev.off()

##################################################################################################################

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

sapply(assays(bsseq, withDimnames = FALSE), class)
print("Bsseq object:")
bsseq
pData(bsseq)

print("The average coverage of CpGs:")
round(colMeans(getCoverage(bsseq)), 1)
print("The number of CpGs:")
length(bsseq)
print("Number of CpGs which are covered by at least 1 read in all samples:")
sum(rowSums(getCoverage(bsseq) >= 1) == 7)
print("Number of CpGs with 0 coverage in all samples:")
sum(rowSums(getCoverage(bsseq)) == 0)


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

#Remove CpGs with coverage below 4 in all samples:
BS.cov <- getCoverage(bssmooth)
keepLoci.ex <- which(rowSums(BS.cov[, bsseq$treatment == "control"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "NB"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "AT"] >= 4) >= 3)
print("The number of CpGs with coverage >=4 in all samples:")
length(keepLoci.ex)
bssmooth <- bssmooth[keepLoci.ex,]


#Compute t-statistics based on smoothed whole-genome bisulfite sequencing data:
print("Compute t-stats vs AT:")
AT.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("AT9_E2_1", "AT9_E2_2", "AT9_E2_3"),  
                                    group2 = c("BR9_E2_1", "BR9_E2_2", "BR9_E2_3"),
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE)
AT.tstat
print("Compute t-stats vs NB:")
NB.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("NB9_E2_1", "NB9_E2_2", "NB9_E2_3"),  
                                    group2 = c("BR9_E2_1", "BR9_E2_2", "BR9_E2_3"),
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE)
NB.tstat


#Finding DMRs
print("Find DMRs vs NB:")
dmrs0 <- dmrFinder(AT.tstat, qcutoff = c(0.025, 0.975), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation between normal and cancers of at least 0.1.
AT.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(AT.dmrs)
head(AT.dmrs, n = 3)
write(AT.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT9_dmrs.txt")

print("Find DMRs vs AT:")
dmrs0 <- dmrFinder(AT.tstat, qcutoff = c(0.025, 0.975), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation between normal and cancers of at least 0.1.
NB.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(NB.dmrs)
head(NB.dmrs, n = 3)
write(NB.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB9_dmrs.txt")

#Plot the top 1000 DMRs
print("Plot top NB DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT9_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, AT.dmrs[1:1000,], extend = 5000, 
                addRegions = AT.dmrs)
dev.off()
print("Plot top AT DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB9_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, NB.dmrs[1:1000,], extend = 5000, 
                addRegions = NB.dmrs)
dev.off()

##################################################################################################################

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

sapply(assays(bsseq, withDimnames = FALSE), class)
print("Bsseq object:")
bsseq
pData(bsseq)

print("The average coverage of CpGs:")
round(colMeans(getCoverage(bsseq)), 1)
print("The number of CpGs:")
length(bsseq)
print("Number of CpGs which are covered by at least 1 read in all samples:")
sum(rowSums(getCoverage(bsseq) >= 1) == 7)
print("Number of CpGs with 0 coverage in all samples:")
sum(rowSums(getCoverage(bsseq)) == 0)


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

#Remove CpGs with coverage below 4 in all samples:
BS.cov <- getCoverage(bssmooth)
keepLoci.ex <- which(rowSums(BS.cov[, bsseq$treatment == "control"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "NB"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "AT"] >= 4) >= 3)
print("The number of CpGs with coverage >=4 in all samples:")
length(keepLoci.ex)
bssmooth <- bssmooth[keepLoci.ex,]


#Compute t-statistics based on smoothed whole-genome bisulfite sequencing data:
print("Compute t-stats vs AT:")
AT.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("AT25_E2_1", "AT25_E2_2", "AT25_E2_3"),  
                                    group2 = c("BR25_E2_1", "BR25_E2_2", "BR25_E2_3"),
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE)
AT.tstat
print("Compute t-stats vs NB:")
NB.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("NB25_E2_1", "NB25_E2_2", "NB25_E2_3"),  
                                    group2 = c("BR25_E2_1", "BR25_E2_2", "BR25_E2_3"),
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE)
NB.tstat


#Finding DMRs
print("Find DMRs vs NB:")
dmrs0 <- dmrFinder(AT.tstat, qcutoff = c(0.025, 0.975), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation between normal and cancers of at least 0.1.
AT.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(AT.dmrs)
head(AT.dmrs, n = 3)
write(AT.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT25_dmrs.txt")

print("Find DMRs vs AT:")
dmrs0 <- dmrFinder(AT.tstat, qcutoff = c(0.025, 0.975), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation between normal and cancers of at least 0.1.
NB.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(NB.dmrs)
head(NB.dmrs, n = 3)
write(NB.dmrs, file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB25_dmrs.txt")

#Plot the top 1000 DMRs
print("Plot top NB DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/AT25_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, AT.dmrs[1:1000,], extend = 5000, 
                addRegions = AT.dmrs)
dev.off()
print("Plot top AT DMRs:")
pdf(file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/analysis/Myzus/persicae/WGBS/Archana_Mar2021/bsmooth/NB25_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, NB.dmrs[1:1000,], extend = 5000, 
                addRegions = NB.dmrs)
dev.off()

##################################################################################################################




#Test set
col_names <- c("treatment", "replicate", "col")
row_names <- c("BR1_E2_1", "BR1_E2_2", "BR1_E2_3", "NB1_E2_1", "NB1_E2_3", "AT1_E2_1", "AT1_E2_2")
data <- matrix(c("control", "control", "control", "NB", "NB", "AT", "AT", 1, 2, 3, 1, 3, 1, 2, "blue", "blue", "blue", "red", "red", "green", "green"), nrow = length(row_names), ncol = length(col_names))
df <- data.frame(data, row.names = row_names)
colnames(df) <- col_names
print(df)

bsseq <- read.bismark(files = c("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_BR1_E2_1_bismark.deduplicated.CpG_report.txt", "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_BR1_E2_2_bismark.deduplicated.CpG_report.txt", "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_BR1_E2_3_bismark.deduplicated.CpG_report.txt", "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_NB1_E2_1_bismark.deduplicated.CpG_report.txt", "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_NB1_E2_3_bismark.deduplicated.CpG_report.txt", "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_AT1_E2_1_bismark.deduplicated.CpG_report.txt", "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/temp_AT1_E2_2_bismark.deduplicated.CpG_report.txt"),
    colData = df,
    rmZeroCov = FALSE,
    strandCollapse = TRUE,
    verbose = TRUE)
# This is a matrix-backed BSseq object.
sapply(assays(bsseq, withDimnames = FALSE), class)
bsseq
pData(bsseq)

#The average coverage of CpGs:
round(colMeans(getCoverage(bsseq)), 1)
#The number of CpGs:
length(bsseq)
#Number of CpGs which are covered by at least 1 read in all samples
sum(rowSums(getCoverage(bsseq) >= 1) == 7)
#Number of CpGs with 0 coverage in all samples
sum(rowSums(getCoverage(bsseq)) == 0)

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

#Remove CpGs with coverage below 4 in all samples:
BS.cov <- getCoverage(bssmooth)
keepLoci.ex <- which(rowSums(BS.cov[, bsseq$treatment == "control"] >= 4) >= 3 &
                     rowSums(BS.cov[, bsseq$treatment == "NB"] >= 4) >= 2 &
                     rowSums(BS.cov[, bsseq$treatment == "AT"] >= 4) >= 2)
length(keepLoci.ex)
bssmooth <- bssmooth[keepLoci.ex,]

#Compute t-statistics based on smoothed whole-genome bisulfite sequencing data:
AT.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("AT1_E2_1", "AT1_E2_2"),  
                                    group2 = c("BR1_E2_1", "BR1_E2_2", "BR1_E2_3"),
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE)

NB.tstat <- BSmooth.tstat(bssmooth,
                                    group1 = c("NB1_E2_1", "NB1_E2_3"),  
                                    group2 = c("BR1_E2_1", "BR1_E2_2", "BR1_E2_3"),
                                    estimate.var = "group2",
                                    local.correct = TRUE,
                                    verbose = TRUE)

#Finding DMRs
dmrs0 <- dmrFinder(AT.tstat, qcutoff = c(0.025, 0.975), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation between normal and cancers of at least 0.1.
AT.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(AT.dmrs)
head(AT.dmrs, n = 3)

dmrs0 <- dmrFinder(AT.tstat, qcutoff = c(0.025, 0.975), maxGap=300, verbose = TRUE)
#Filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference (across the DMR) in methylation between normal and cancers of at least 0.1.
NB.dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(NB.dmrs)
head(NB.dmrs, n = 3)

#Plot the top 1000 DMRs
pdf(file = "AT1_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, AT.dmrs[1:1000,], extend = 5000, 
                addRegions = AT.dmrs)
dev.off()
pdf(file = "NB1_dmrs_top1000.pdf", width = 10, height = 5)
plotManyRegions(bssmooth, NB.dmrs[1:1000,], extend = 5000, 
                addRegions = NB.dmrs)
dev.off()
```
```R
#Convert from bsmap format to bsseq input format. 

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
```bash

```
## Check for transcriptional changes over the course of the experiments 1 & 2
### QC 
#### fastqc 

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
### Trimming

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

### Mapping
We don't have a reference transcriptome for M.persicae for alignment free approach.

#### Alignment free approach

Could generate de novo transcriptome with trinity. Trinity could find novel transcritps from George's data however discovering novel transcripts is not the focus of this study.

##### Trinity

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

#Run phase 2 step for these 23 with increased heap space:
for file in $(cat temp_trinity_23.txt); do
#Trinity --single $file --output ${file}.out --CPU 1 --max_memory 1G --run_as_paired --seqType fa --trinity_complete --full_cleanup --min_kmer_cov 2 --verbose --bflyHeapSpaceMax 100G
sbatch ~/git_repos/Wrappers/NBI/temp4.sh $file
done
#57795755-77

#Trinity's checkpoints only work if the original command hasn't been modified, --FORCE must therefore be used to force inclusion of both the 401026 and 23 in the final transcriptome:
Trinity --full_cleanup --seqType fq --CPU 64 --max_memory 4030G --min_kmer_cov 2 --normalize_by_read_set --verbose --FORCE\
--left $forward \
--right $reverse \
--output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/M_persicae_trinity_transcriptome 2>&1 >> logs/trinity_run_log2.txt
#57765347 - check for errors after the 99.994% sticking point that will now be forced but have not been run in isolation?

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/M_persicae_trinity_transcriptome.Trinity.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/M_persicae_trinity_transcriptome.Trinity-force1.fasta



#du --max-depth=1 --total -h /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/M_persicae_trinity_transcriptome
#468G    total

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/M_persicae_trinity_transcriptome/read_partitions/Fb_*/CBin_*/*.trinity.reads.fa); do
if  [ ! -e ${file}.out.Trinity.fasta ]; then
echo ${file}
fi
done
Finished.  Final Trinity assemblies are written to /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/M_persicae_trinity_transcriptome/read_partitions/Fb_*/CBin_*/*.trinity.reads.fa.out.Trinity.fasta




Trinity --single "$file" --output "${file}.out" --CPU 1 --max_memory 1G --run_as_paired --seqType fa --trinity_complete --full_cleanup --min_kmer_cov 2 --verbose --bflyHeapSpaceMax 100G
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

#### Alignment approach

##### STAR
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



##### Stringtie

Stringtie generation of transcriptome for alignment free appraoch requires alignments...

```bash

```