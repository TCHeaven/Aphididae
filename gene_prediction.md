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
Trim Galore v0.4.5, retaining reads where both members of the pair are at least 20-bp long.
```bash

```
```bash
OutDir=
OutFile=
Freads_list=
Rreads_list=
Max_intron=10000
Genome=
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_trinity.sh $OutDir $OutFile $Freads_list $Rreads_list $Max_intron $Genome
```