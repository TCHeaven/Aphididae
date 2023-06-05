```bash
mkdir /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/blast_nt_db
cd /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/blast_nt_db
gzip -cd /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa.gz > Myzus_persicae_O_v2.0.scaffolds.fa
source package 37f0ffda-9f66-4391-87e2-38ccd398861d
makeblastdb -in Myzus_persicae_O_v2.0.scaffolds.fa -input_type fasta -dbtype nucl  -title MP  -parse_seqids -out MP
rm Myzus_persicae_O_v2.0.scaffolds.fa
```

SWeeD
```bash
cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Australia\|Asia\|Africa\|America' | awk -F',' '{print $2}' > temp1.txt
wc -l temp1.txt #59
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp1.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=Europe_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523176, 54541557, 54650460

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Australia\|Asia\|Africa\|Europe' | awk -F',' '{print $2}' > temp2.txt
wc -l temp2.txt #198
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp2.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=America_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523178, 54541563

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Australia\|Asia\|America\|Europe' | awk -F',' '{print $2}' > temp3.txt
wc -l temp3.txt #194
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp3.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=Africa_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523179, 54541566

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Australia\|Africa\|America\|Europe' | awk -F',' '{print $2}' > temp4.txt
wc -l temp4.txt #194
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp4.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=Asia_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523181, 54541571

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Asia\|Africa\|America\|Europe' | awk -F',' '{print $2}' > temp5.txt
wc -l temp5.txt #191
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp5.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=Australia_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523183, 54541575

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'JIC' | awk -F',' '{print $2}' > temp6.txt
wc -l temp6.txt #79
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp6.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=Bass_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523186, 54541577, 54586869, 54851916

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Bass' | awk -F',' '{print $2}' > temp7.txt
wc -l temp7.txt #130
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp7.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=JIC_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523187, 54541579, 54586885

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Brassica\|Prunus\|Solanacae\|Tobacco\|Bass\|Others\|Sugar' | awk -F',' '{print $2}' > temp8.txt
wc -l temp8.txt #197
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp8.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=JIC_Lab_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523189, 54541614

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Brassica\|Prunus\|Solanacae\|Tobacco\|Others\|Sugar' | awk -F',' '{print $2}' > temp9.txt
wc -l temp9.txt #195
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp9.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=Lab_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523190, 54542450

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Brassica\|Prunus\|Solanacae\|Tobacco\|Lab\|Sugar' | awk -F',' '{print $2}' > temp10.txt
wc -l temp10.txt #198
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp10.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=others_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523191, 54542471

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Brassica\|Prunus\|Solanacae\|Tobacco\|Lab\|Others' | awk -F',' '{print $2}' > temp11.txt
wc -l temp11.txt #162
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp11.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=Sugar_beet_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523192, 54542486

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Brassica\|Prunus\|Solanacae\|Lab\|Others\|Sugar' | awk -F',' '{print $2}' > temp12.txt
wc -l temp12.txt #182
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp12.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=Tabacco_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523193, 54542512

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Brassica\|Prunus\|Tobacco\|Lab\|Others\|Sugar' | awk -F',' '{print $2}' > temp13.txt
wc -l temp13.txt #185
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp13.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=Solanacae_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523194, 54544097

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Brassica\|Solanacae\|Tobacco\|Lab\|Others\|Sugar' | awk -F',' '{print $2}' > temp14.txt
wc -l temp14.txt #155
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp14.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=Prunus_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523195, 54544098

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Prunus\|Solanacae\|Tobacco\|Lab\|Others\|Sugar' | awk -F',' '{print $2}' > temp15.txt
wc -l temp15.txt #177
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp15.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=Brassica_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523197, 54544100
rm temp1.txt
rm temp2.txt
rm temp3.txt
rm temp4.txt
rm temp5.txt
rm temp6.txt
rm temp7.txt
rm temp8.txt
rm temp9.txt
rm temp10.txt
rm temp11.txt
rm temp12.txt
rm temp13.txt
rm temp14.txt
rm temp15.txt
echo done
```
```R
library(ggplot2)
library(dplyr)
library(tidyr)

# First 6 CHROMOSOMES ONLY:
#Bass
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Bass")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "Bass samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#JIC
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.JIC")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "JIC samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#Africa
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Africa")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "Africa samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#America
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.America")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "America samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#Asia
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Asia")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "Asia samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#Australia
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Australia")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "Australia samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#Europe
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Europe")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "Europe samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#Lab
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Lab")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "Lab samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#JIC Lab
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.JIC_Lab")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "JIC Lab samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#Brassica
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Brassica")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "Brassica samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#others
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.others")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "others samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#Prunus
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Prunus")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "Prunus samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#Solanacae
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Solanacae")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "Solanacae samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#Sugar_beet
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Sugar_beet")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "Sugar_beet samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

#Tabacco
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Tabacco")
num_searches <- sum(grepl("^//", lines))
data_list <- vector("list", length = num_searches)
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- mutate(search_data, Search = i)
  data_list[[i]] <- search_data
}
data <- bind_rows(data_list[1:min(6, num_searches)])
ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  facet_wrap(~Search, scales = "free_x", nrow = 1) +
  labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = "Tabacco samples 100") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

```