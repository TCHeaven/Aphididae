# Population genomics as part of the M. persicae population genomics project

## Contents

1.[SWeeD](#2)
2.[admixture](#3)
3.[STRUCTURE](#4)
4.[FST](#5)
5.[Gemma](#7)

#### SWeeD <a name="2"></a>

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
wc -l temp6.txt #92
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp6.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=Bass_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523186, 54541577, 54586869, 54851916, 55211023, 55211079

cat snp_calling/Myzus/persicae/biello/PCA_file_host.csv | grep 'Bass' | awk -F',' '{print $2}' > temp7.txt
wc -l temp7.txt #117
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz);do
    Exclusion_list=temp7.txt
    ProgDir=~/git_repos/Wrappers/NBI
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD
    OutFile=JIC_100000
    mkdir $OutDir
    sbatch $ProgDir/run_SweeD.sh $InFile $OutDir $OutFile $Exclusion_list
done #54523187, 54541579, 54586885, 55211032, 55211084

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
Plots of SweeD with 100 regions
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
Plots of SweeD with 100,000 regions
```R
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("C:/Users/did23faz/OneDrive - Norwich Bioscience Institutes/Desktop/R")

# First 6 searches (1-6) plotted separately and saved as PDF files
lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Bass_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("Bass samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("Bass_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("Bass_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.JIC_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("JIC samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("JIC_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("JIC_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Africa_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("Africa samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("Africa_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("Africa_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.America_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("America samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("America_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("America_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Asia_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("Asia samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("Asia_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("Asia_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Australia_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("Australia samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("Australia_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("Australia_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Europe_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("Europe samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("Europe_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("Europe_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Lab_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("Lab samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("Lab_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("Lab_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.others_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("others samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("others_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("others_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Prunus_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("Prunus samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("Prunus_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("Prunus_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Solanacae_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("Solanacae samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("Solanacae_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("Solanacae_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Sugar_beet_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("Suagr beet samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("Sugar_beet_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("Sugar_beet_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Tabacco_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("Tabacco samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("Tabacco_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("Tabacco_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.Brassica_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("Brassica samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("Brassica_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("Brassica_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}

lines <- readLines("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/SweeD/SweeD_Report.JIC_Lab_100000")
num_searches <- sum(grepl("^//", lines))
for (i in 1:min(6, num_searches)) {
  start_index <- grep(paste0("^//", i), lines) + 2
  end_index <- ifelse(i < num_searches, grep(paste0("^//", i + 1), lines) - 1, length(lines))
  search_lines <- lines[start_index:end_index]
  search_data <- read.table(text = search_lines, header = TRUE, stringsAsFactors = FALSE)
  colnames(search_data) <- c("Position", "Likelihood", "Alpha")
  search_data$Position <- search_data$Position / 1e6
  search_data <- dplyr::mutate(search_data, Search = i)
  
  # Plot the data with custom point size
  plot <- ggplot(search_data, aes(x = Position, y = Likelihood)) +
    geom_point(size = 0.5) +  # Adjust the point size here (e.g., size = 3)
    labs(x = "Position (Mb)", y = "Likelihood", color = "Search", title = paste("JIC Lab samples 100,000 - Search", i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0, 110)) +  # Set x-axis limits to 0-120
    scale_y_continuous(limits = c(0, 400))  # Set y-axis limits to 0-400
  
  # Save the plot as PNG image
  png_file <- paste0("JIC_lab_100000_stretch_chromosome_", i, ".png")
  ggsave(png_file, plot, width = 48, height = 4, dpi = 600)    
  # Save the plot with custom dimensions
  pdf_file <- paste0("JIC_lab_100000_stretch_chromosome_", i, ".pdf")
  ggsave(pdf_file, plot, width = 48, height = 6)  # Adjust the width and height values as needed
}
```

#### admixture (INCOMPLETE) <a name="3"></a>
```bash
source package /tgac/software/testing/bin/admixture-1.3.0
source package /nbi/software/testing/bin/bcftools-1.8
source package /nbi/software/testing/bin/plink-1.9 
source package 01ef5a53-c149-4c9e-b07d-0b9a46176cc0

#there are 360 scaffolds in the assembly, however plink cannot process more than 95 chromosome, therefore will use only the large scaffolds 1-6 for admixture analysis.
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz | head -n 5 > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs.vcf
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz | grep -w '##ALT\|##INFO\|##FORMAT\|#CHROM\|##bcftools_callVersion\|##bcftools_callCommand\|##bcftools_concatVersion\|##bcftools_concatCommand\|##bcftools_filterVersion\|##bcftools_filterCommand\|##bcftools_viewVersion\|##bcftools_viewCommand\|scaffold_1\|scaffold_2\|scaffold_3\|scaffold_4\|scaffold_5\|scaffold_6' | sed 's@scaffold_@chr@g' >> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs.vcf
bgzip -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs.vcf > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs.vcf.gz
bcftools sort /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs.vcf.gz -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted.vcf.gz -Oz
plink --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted.vcf.gz --make-bed --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted
plink --bfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted --indep 50 5 2 --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted_pruned
#Pruned 2204765 variants from chromosome 1, leaving 262033.
#Pruned 2626597 variants from chromosome 2, leaving 247859.
#Pruned 2127203 variants from chromosome 3, leaving 212903.
#Pruned 1902203 variants from chromosome 4, leaving 193359.
#Pruned 963227 variants from chromosome 5, leaving 105113.
#Pruned 898708 variants from chromosome 6, leaving 95341.
#Pruning complete.  10722703 of 11839311 variants removed.
plink --bfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted --extract /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted_pruned.prune.in --make-bed --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted_pruned_set

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40; \
do admixture --cv=10 --seed=1234 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted_pruned_set.bed $K | tee /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/admixture/log${K}.out; done

admixture --cv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted_pruned_set.bed 2 #Submitted batch job 56002358

for bedfile in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted_pruned_set.bed); do
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/admixture
Mink=2
Maxk=40
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_admixture_cross_validation.sh $bedfile $Mink $Maxk $OutDir
done #56002475, 56144348 from 12 up

for K in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40; do 
bedfile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted_pruned_set.bed
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/admixture2
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_admixture_cross_validation.sh $bedfile $K $OutDir
done #56178239-56178277, 56258997-56259019, 56359254-56359276, 56607044 (38), 56607045 (39), 56607048 (40), 56607073-6 (27 28 29 30)

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/admixture2/log*); do
CV=$(grep 'CV error' $file | sed 's@CV error (@@g'| sed 's@):@@g')
echo $CV
done

Pruned_vcf=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted_pruned_set
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/admixture
OutFile=193s.M_persicae.genomicSNPs
Mink=21
Maxk=21
Bootstraps=200
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_admixture.sh $OutDir $OutFile $Pruned_vcf $Mink $Maxk $Bootstraps #57059549
```
#### STRUCTURE (INCOMPLETE) <a name="4"></a>
```bash
# Extract individual IDs and genotype data
source package /nbi/software/testing/bin/vcftools-0.1.15
vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs.vcf --plink --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/structure/193s.M_persicae

# Convert PLINK files to STRUCTURE format
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/structure
source package /nbi/software/testing/bin/plink-1.9 
plink --file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/structure/193s.M_persicae --recode structure --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/structure/193s.M_persicae
```
#### FST (INCOMPLETE) <a name="5"></a>
```bash
source package d37013e7-5691-40b6-8885-f029fe5fad54
mkdir snp_calling/Myzus/persicae/biello/gatk/FST
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz --weir-fst-pop population1.txt --weir-fst-pop population2.txt 
#Weir and Cockerham mean Fst estimate: 0.047368
#Weir and Cockerham weighted Fst estimate: 0.10435
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz --weir-fst-pop population1.txt --weir-fst-pop population2.txt --fst-window-size 10000 --fst-window-step 5000 --out ./snp_calling/Myzus/persicae/biello/gatk/FST/output_fst.txt
#Warning: Expected at least 2 parts in INFO entry: ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
#Warning: Expected at least 2 parts in INFO entry: ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
#Warning: Expected at least 2 parts in INFO entry: ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
#Output file is empty

for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=snp_calling/Myzus/persicae/biello/gatk/FST
    OutFile=output_fst_10000-5000.txt
    Populations=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/two.txt
    WindowSize=10000
    WindowSlide=5000
    sbatch $ProgDir/run_FST.sh $vcf $OutDir $OutFile $Populations $WindowSize $WindowSlide
done #55637450, 55675621


Traceback (most recent call last):
  File "/hpc-home/did23faz/git_repos/Scripts/NBI/FST.py", line 13, in <module>
    with open(sample_population_file, 'r') as pop_file:
FileNotFoundError: [Errno 2] No such file or directory: 'population.txt'
cp: cannot stat 'fst.txt': No such file or directory
population.txt
```
## GWAS <a name="6"></a>
### Gemma <a name="7"></a>

GEMMA: Genome-wide Efficient Mixed Model Association
```bash
perl ~/git_repos/Scripts/NBI/bcf2bbgeno.pl -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.bbgeno -p H-W -s -r
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/
gzip 193s.M_persicae.onlySNPs-genic-regions.bbgeno
#57527288
```
```bash
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz | awk 'BEGIN{OFS="\t"} {if ($3 == ".") $3 = $1"_"$2"_"$5; print}' | gzip > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.vcf.gz

zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.vcf.gz | sed 's/scaffold_//g' | gzip > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.tmp.vcf.gz && mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.tmp.vcf.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.mod.vcf.gz

##################################################################################################################
source package 01ef5a53-c149-4c9e-b07d-0b9a46176cc0
bgzip -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf.gz

zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.vcf.gz | awk 'BEGIN{OFS="\t"} {if ($3 == ".") $3 = $1"_"$2"_"$5; print}' | bgzip -c > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf.gz

#57451298

#there are 360 scaffolds in the assembly, however plink cannot process more than 95 chromosome, therefore will use only the large scaffolds 1-6 for admixture analysis.
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf.gz | head -n 5 > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf.gz | grep -w '##ALT\|##INFO\|##FORMAT\|#CHROM\|##bcftools_callVersion\|##bcftools_callCommand\|##bcftools_concatVersion\|##bcftools_concatCommand\|##bcftools_filterVersion\|##bcftools_filterCommand\|##bcftools_viewVersion\|##bcftools_viewCommand\|scaffold_1\|scaffold_2\|scaffold_3\|scaffold_4\|scaffold_5\|scaffold_6' | sed 's@scaffold_@chr@g' >> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/
bgzip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf
#57451608

source package /nbi/software/testing/bin/plink-1.9 
plink --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod.vcf.gz --make-bed --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod
#57452909
plink --bfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod --genome --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod
#57501018
```
```python
import pandas as pd
import numpy as np

ibd_data = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/193s.M_persicae.onlySNPs-mac1.recode.mod.genome", delim_whitespace=True)
df = ibd_data[['FID1', 'IID1', 'FID2', 'IID2', 'PI_HAT']]
kinship_matrix = df.pivot_table(values='PI_HAT', index=['FID1', 'IID1'], columns=['FID2', 'IID2']).values
kinship_matrix = np.nan_to_num(kinship_matrix)
kinship_matrix -= np.mean(kinship_matrix)
np.savetxt("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink2/kinship_matrix.txt", kinship_matrix, fmt='%.6f', delimiter='\t')
```
Graphs were plotted in excel
