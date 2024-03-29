# Gene prediction
We are interested whether all genes have been successfully annotated in the M. persicae genome; v2.0 annotations were initially performed by T. mathers using braker2, David swarbreck performed v2.1 annotations - his methods for this are not shared - this is our best current annotation, M. Marzo reports success with the helixer annotator.

It may be possible to improve annotations with the addition of George's organ transcriptome data.

It may be possible to add to annotations with new/different methods: trinity assembly free prediction, braker3, EASEL

'Bicycle genes' may be missing from the annotations performed to data - found in other aphid species.

## Collecting data
For M. persicae, we aligned 25 RNA-seq libraries. Specifically, we used a high-coverage (200 million reads), strand-specific, RNA-seq library generated from mixed whole bodies of apterous M. persicae clone O asexual females (Mathers et al. 2017 - PRJEB11304) as well as newly generated (PRJNA613055) and publicly available (Mathers et al. 2019 - PRJNA437622) unstranded RNA-seq data for M. persicae clone O nymphs (derived from apterous asexual females), alate asexual females, apterous asexual females and males (six biological replicates each).

Key for the RNA-Seq data used in Mathers et al. (2020) assembly annotation, 209G raw data in total:
```bash
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/Mathers_2020_RNA_datasets.xlsx
```
### Raw data locations
RNA-seq used for gene annotation of M. persicae by Mathers et al. (2020) were downloaded from NCBI (PRJNA613055). (12/25)
```bash
#Myzus persicae clone O asexual female nymph whole body RNA-seq
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322156_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322156_2.fastq.gz
#Myzus persicae clone O asexual female nymph whole body RNA-seq
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322157_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322157_2.fastq.gz
#Myzus persicae clone O asexual female nymph whole body RNA-seq
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322159_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322159_2.fastq.gz
#Myzus persicae clone O asexual female nymph whole body RNA-seq
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322160_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322160_2.fastq.gz
#Myzus persicae clone O asexual female nymph whole body RNA-seq
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322161_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322161_2.fastq.gz
#Myzus persicae clone O asexual female nymph whole body RNA-seq
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322162_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322162_2.fastq.gz
#Myzus persicae clone O winged female whole body RNA-seq
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322163_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322163_2.fastq.gz
#Myzus persicae clone O winged female whole body RNA-seq
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322164_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322164_2.fastq.gz
#Myzus persicae clone O winged female whole body RNA-seq
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322165_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322165_2.fastq.gz
#Myzus persicae clone O winged female whole body RNA-seq
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322166_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322166_2.fastq.gz
#Myzus persicae clone O winged female whole body RNA-seq
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322167_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322167_2.fastq.gz
#Myzus persicae clone O winged female whole body RNA-seq
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322168_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR11322168_2.fastq.gz
```
RNA-seq used for gene annotation of M. persicae by Mathers et al. (2020) were downloaded from NCBI (PRJNA437622). (12/25)
```bash
#Myzus persicae clone O winged males
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821683_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821683_2.fastq.gz 
#Myzus persicae clone O winged males
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821684_1.fastq.gz 
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821684_2.fastq.gz 
#Myzus persicae clone O asexual unwinged females
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821685_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821685_2.fastq.gz 
#Myzus persicae clone O asexual unwinged females
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821686_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821686_2.fastq.gz 
#Myzus persicae clone O asexual unwinged females
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821687_1.fastq.gz 
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821687_2.fastq.gz
#Myzus persicae clone O asexual unwinged females
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821688_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821688_2.fastq.gz 
#Myzus persicae clone O asexual unwinged females
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821689_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821689_2.fastq.gz 
#Myzus persicae clone O asexual unwinged females
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821690_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821690_2.fastq.gz 
#Myzus persicae clone O winged males
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821691_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821691_2.fastq.gz 
#Myzus persicae clone O winged males
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821692_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821692_2.fastq.gz 
#Myzus persicae clone O winged males
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821693_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821693_2.fastq.gz  
#Myzus persicae clone O winged males
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821694_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR6821694_2.fastq.gz
```
RNA-seq used for gene annotation of M. persicae by Mathers et al. (2020) were downloaded from NCBI (PRJEB11304). (1/25)
```bash
#Myzus persicae clone O, Unwinged asexual female: mixed whole adults and nymphs - Illumina HiSeq 2000 SS PE 
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661483_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661483_2.fastq.gz  

#The following from this study (PRJEB11304) were not used for annotation by Mathers et al. (2020)
#Myzus persicae clone O, adult asexual females reared on B. rapa (1 year) - Cb3
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661479_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661479_2.fastq.gz
#Myzus persicae clone O, adult asexual females reared on N. benthamiana (1 year) - Nb1
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661478_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661478_2.fastq.gz
#Myzus persicae clone O, adult asexual females reared on B. rapa (1 year) - Cb1
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661477_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661477_2.fastq.gz
#Myzus persicae clone O, adult asexual females reared on N. benthamiana (1 year) - Nb4
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661482_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661482_2.fastq.gz
#Myzus persicae clone O, adult asexual females reared on N. benthamiana (1 year) - Nb2
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661480_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661480_2.fastq.gz
#Myzus persicae clone O, adult asexual females reared on B. rapa (1 year) - Cb4
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661481_1.fastq.gz
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661481_2.fastq.gz
#Myzus persicae clone O, nymphs, miRNA-Seq
ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1637389.fastq.gz
```
### Linking raw data into working directory
```bash
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Mathers2020
ln -s /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Mathers2020/.
ln -s /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661483* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Mathers2020/.
```
## Quality control
The raw RNAseq reads were subjected to a quality control check using FastQC:
```bash
for RawData in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Mathers2020/*.fastq.gz); do
ProgDir=~/git_repos/Wrappers/NBI
OutFile=$(basename $RawData .fastq.gz)
OutDir=$(dirname $RawData)
echo $OutFile
sbatch $ProgDir/run_fastqc.sh $RawData $OutDir $OutFile
done
```
The raw RNAseq reads were trimmed to remove adapters an low quality regions:
```bash
#NOTE: to save space the script has been edited to delete input files upon production of outputs. - make sure that you check that it runs correctly first before deleting all input data - in this case only the symlinks should be deleted, these will be replaced afterwards
for Fread in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Mathers2020/*1.fastq.gz); do
    Rread=$(ls $Fread | sed 's@_1.@_2.@g')
    Fread2=$(ls $Fread | sed 's@_1.@_3.@g')
    Rread2=$(ls $Fread | sed 's@_1.@_4.@g')
    OutDir=$(dirname $Fread | sed 's@raw_data@rna_qc@g')trim_galore
    OutFile=$(basename $Fread _1.fastq.gz)
    Quality=20
    Length=20
    ProgDir=~/git_repos/Wrappers/NBI
    Jobs=$(squeue -u did23faz| grep 'trim_g'  | wc -l)
    echo x
    while [ $Jobs -gt 19 ]; do
        sleep 300s
        printf "."
        Jobs=$(squeue -u did23faz| grep 'trim_g'  | wc -l)
    done
    mkdir -p $OutDir
    echo $OutFile >> logs/rna_trim_galore_report.txt
    sbatch $ProgDir/run_trim_galore.sh $OutDir $OutFile $Quality $Length $Fread $Rread $Fread2 $Rread2  2>&1 >> logs/rna_trim_galore_report.txt
done

ln -s /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/SRR* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Mathers2020/.
ln -s /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/CloneORNA/ERR1661483* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/RNA_Seq/Mathers2020/.
```
The trimmed reads were then re-assessed for quality following trimming, using FastQC:
```bash
for QCData in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/rna_qc/Myzus/persicae/RNA_Seq/Mathers2020/trim_galore/*.fq.gz); do
ProgDir=~/git_repos/Wrappers/NBI
OutFile=$(basename $QCData .fastq.gz)
OutDir=$(dirname $QCData)
echo $OutFile
sbatch $ProgDir/run_fastqc.sh $QCData $OutDir $OutFile
done
```
## Transcriptome assembly

#### Trinity
Trinity was run to assembly a transcriptome for M. persicae using the RNASeq data, this wrapper script will do this twice: once with the --genome_guided_bam flag, and once without reference to the genome assembly.
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/rna_qc/Myzus/persicae/RNA_Seq/Mathers2020/trim_galore); do
OutDir=$(echo $ReadDir | sed 's@rna_qc@assembly/transcriptome@g' | sed 's@RNA_Seq/Mathers2020/trim_galore@trinity_2.9.1@g')
mkdir -p $OutDir
OutFile=Mathers2020_25libs
Freads_list=$(ls ${ReadDir}/*.fq.gz | grep '_1.fq.gz' | grep -v 'ERR1661483' | tr '\n' ','| sed 's/,$//')
Rreads_list=$(ls ${ReadDir}/*.fq.gz | grep '_2.fq.gz' | grep -v 'ERR1661483' | tr '\n' ','| sed 's/,$//')
SSFreads_list=$(ls ${ReadDir}/*.fq.gz | grep '_1.fq.gz' | grep 'ERR1661483' | tr '\n' ','| sed 's/,$//')
SSRreads_list=$(ls ${ReadDir}/*.fq.gz | grep '_2.fq.gz' | grep 'ERR1661483' | tr '\n' ','| sed 's/,$//')
Max_intron=10000
Genome=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
ProgDir=~/git_repos/Wrappers/NBI
#sbatch $ProgDir/run_trinity.sh $OutDir $OutFile $Freads_list $Rreads_list $SSFreads_list $SSRreads_list $Max_intron $Genome
done #56025935, 56149160

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/rna_qc/Myzus/persicae/RNA_Seq/Mathers2020/trim_galore); do
OutDir=$(echo $ReadDir | sed 's@rna_qc@assembly/transcriptome@g' | sed 's@RNA_Seq/Mathers2020/trim_galore@trinity_2.9.1@g')
mkdir -p $OutDir
OutFile=Mathers2020_25libs
Freads_list=$(ls ${ReadDir}/*.fq.gz | grep '_1.fq.gz' | grep -v 'ERR1661483' | tr '\n' ','| sed 's/,$//')
Rreads_list=$(ls ${ReadDir}/*.fq.gz | grep '_2.fq.gz' | grep -v 'ERR1661483' | tr '\n' ','| sed 's/,$//')
SSFreads_list=$(ls ${ReadDir}/*.fq.gz | grep '_1.fq.gz' | grep 'ERR1661483' | tr '\n' ','| sed 's/,$//')
SSRreads_list=$(ls ${ReadDir}/*.fq.gz | grep '_2.fq.gz' | grep 'ERR1661483' | tr '\n' ','| sed 's/,$//')
Max_intron=10000
Genome=/jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_trinity.sh $OutDir $OutFile $Freads_list $Rreads_list $SSFreads_list $SSRreads_list $Max_intron $Genome
done #56265401 - run with de novo hashed out and workdir left
```
##### Trinity (from 9_Host_swap.md)

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
