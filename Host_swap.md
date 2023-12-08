## Collect data
All files that look like analysis not raw data have been removed/not copied over in order to save space.
```bash
#Organ data raw data:
cp -r /ei/projects/a/a93e4b69-58e3-495b-ac4d-04e978fed5d1/data/organ_data/raw_data/for_archana/* /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/archana_organ_data/.

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

cp -r /ei/projects/b/bda78a5f-7de9-4e0c-813f-1b1105bd6c24/host_adjustment/HostAdaptation/raw_data/X201SC20052390-Z01-F001_0*_20200727.zip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana/.
#All files unzipped
mkdir At
mkdir Br
mkdir Nb
mv AT* At/.
mv BR* Br/.
mv NB* Nb/.

cp -r /ei/projects/b/bda78a5f-7de9-4e0c-813f-1b1105bd6c24/RNA_Mar2021/* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Mar2021/.
#All files unzipped
mkdir BR
mkdir NB
mv BR_* BR/.
mv NB_* NB/.

cp -r /ei/projects/b/bda78a5f-7de9-4e0c-813f-1b1105bd6c24/host_adjustment/Data_Dec2020/* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Archana_Dec2020/.
#All files unzipped
mkdir At
mkdir Br
mkdir Nb
mv AT* At/.
mv BR* Br/.
mv NB* Nb/.


#This is analysis - removed to make room for all raw data
cp -r analysis /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/Archana_hostadaptation_analysis/.
cp -r Effector /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/Archana_hostadaptation_analysis/.
cp -r * /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/Archana_hostadaptation_analysis/alternative_splicing/.
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/Archana_hostadaptation_analysis
```
## Check mutations in the aphid genomes over the course of the experiment
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
# WGBS
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
    if 
    [ ! -e "${OutDir}/${OutFile}.bam.bai" ] || [ ! -s "${OutDir}/${OutFile}.bam.bai" ]; then
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
#57586990-
```