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
VCFstats
```bash
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-genic-regions.vcf.gz); do
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity
    OutFile=193_genic
    Exclusion_list=
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_vcfstats.sh $InFile $OutDir $OutFile $Exclusion_list
done #55244504
```
Plot VCFStats
```R
setwd("C:/Users/did23faz/OneDrive - Norwich Bioscience Institutes/Desktop/R")

install.packages("tidyverse")
install.packages("ggplot2")

library(tidyverse)
library(ggplot2)

var_qual <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_site_qual.lqual", delim = "\t",
                       col_names = c("CHROM", "POS", "QUAL"), skip = 1)

var_depth <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_site-mean-depth.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

var_miss <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_missing_sites.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

var_freq <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_allele_freq.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

var_count <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_allele_count.frq.count", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "count", "a2"), skip = 1)

ind_depth <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_depth.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

ind_miss  <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_missing_ind.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

ind_het <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_genic_het.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)


p1 <- ggplot(var_qual, aes(QUAL)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Per Site Quality")
p1 <- p1 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_var_qual.pdf", width=11, height=7)
plot(p1, pdf=T)
dev.off()

p2 <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + labs(title = "Genic SNPs: Mean Depth Per Site")
p2 <- p2 + theme(plot.title = element_text(hjust = 0.5))
summary(var_depth$mean_depth)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  8.249  16.212  19.529  18.710  21.285  30.829 
pdf("193_genic_var_depth.pdf", width=11, height=7)
plot(p2, pdf=T)
dev.off()

p3 <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + xlim(0, 100) + labs(title = "Genic SNPs: Mean Depth Per Site")
p3 <- p3 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_var_depth_x100.pdf", width=11, height=7)
plot(p3, pdf=T)
dev.off()

p4 <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Missingness Per Site")
p4 <- p4 + theme(plot.title = element_text(hjust = 0.5))
summary(var_miss$fmiss)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.005181 0.010363 0.016646 0.025907 0.108808 
pdf("193_genic_var_miss.pdf", width=11, height=7)
plot(p4, pdf=T)
dev.off()

var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
summary(var_freq$maf)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.00000 0.00000 0.01121 0.00000 0.50000 
p5 <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Allele Frequency for Each Site")
p5 <- p5 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_var_freq.pdf", width=11, height=7)
plot(p5, pdf=T)
dev.off()
p51 <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Allele Frequency for Each Site") + xlim(0, 0.1)
p51 <- p51 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_var_freq_0.1.pdf", width=11, height=7)
plot(p51, pdf=T)
dev.off()

p6 <- ggplot(var_count, aes(x = count)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Raw Allele Counts for Each Site")
p6 <- p6 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_allele_count_density.pdf", width = 11, height = 7)
plot(p6, pdf=T)
dev.off()

p7 <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Mean Depth Per Individual")
p7 <- p7 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_ind_depth.pdf", width=11, height=7)
plot(p7, pdf=T)
dev.off()

p8 <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Missingness Per Individual")
p8 <- p8 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_ind_miss.pdf", width=11, height=7)
plot(p8, pdf=T)
dev.off()

p9 <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Genic SNPs: Heterozygosity Per Individual")
p9 <- p9 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_genic_ind_het.pdf", width=11, height=7)
plot(p9, pdf=T)
dev.off()
```
Effector SNPs
This file contains Candidate effectors identified from protein-coding genes from the v2.1 annotation of the M. persicae clone 0 based on the v2.0 assembly (Archana, using Augustus). Transcripts showing high expression in salivary glands  (>25TPM from Jia’s dissected samples), and encoding proteins with predicted secretory signal peptides (SingalP v4.0), OR that were detected in the saliva of aphids by mass spectrometry in previously published studies- this led to an identification of a set of 496 candidate effectors (gene IDs in Effector.list.v2.1.xlsx).
```bash
/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/effector_candidates.txt
```
```bash
awk '{printf "%s\\|", $0}' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/effector_candidates.txt

grep 'MYZPE13164_O_EIv2.1_0002220\|MYZPE13164_O_EIv2.1_0037470\|MYZPE13164_O_EIv2.1_0037650\|MYZPE13164_O_EIv2.1_0037670\|MYZPE13164_O_EIv2.1_0043160\|MYZPE13164_O_EIv2.1_0043750\|MYZPE13164_O_EIv2.1_0045990\|MYZPE13164_O_EIv2.1_0049040\|MYZPE13164_O_EIv2.1_0054380\|MYZPE13164_O_EIv2.1_0055470\|MYZPE13164_O_EIv2.1_0056440\|MYZPE13164_O_EIv2.1_0059180\|MYZPE13164_O_EIv2.1_0061940\|MYZPE13164_O_EIv2.1_0063800\|MYZPE13164_O_EIv2.1_0063860\|MYZPE13164_O_EIv2.1_0063890\|MYZPE13164_O_EIv2.1_0064380\|MYZPE13164_O_EIv2.1_0067290\|MYZPE13164_O_EIv2.1_0075010\|MYZPE13164_O_EIv2.1_0080430\|MYZPE13164_O_EIv2.1_0080440\|MYZPE13164_O_EIv2.1_0080950\|MYZPE13164_O_EIv2.1_0082260\|MYZPE13164_O_EIv2.1_0082270\|MYZPE13164_O_EIv2.1_0082380\|MYZPE13164_O_EIv2.1_0082580\|MYZPE13164_O_EIv2.1_0082750\|MYZPE13164_O_EIv2.1_0083010\|MYZPE13164_O_EIv2.1_0083400\|MYZPE13164_O_EIv2.1_0084880\|MYZPE13164_O_EIv2.1_0087020\|MYZPE13164_O_EIv2.1_0087460\|MYZPE13164_O_EIv2.1_0087920\|MYZPE13164_O_EIv2.1_0087930\|MYZPE13164_O_EIv2.1_0087950\|MYZPE13164_O_EIv2.1_0088410\|MYZPE13164_O_EIv2.1_0092760\|MYZPE13164_O_EIv2.1_0092990\|MYZPE13164_O_EIv2.1_0094080\|MYZPE13164_O_EIv2.1_0097540\|MYZPE13164_O_EIv2.1_0098320\|MYZPE13164_O_EIv2.1_0098420\|MYZPE13164_O_EIv2.1_0098460\|MYZPE13164_O_EIv2.1_0099380\|MYZPE13164_O_EIv2.1_0100510\|MYZPE13164_O_EIv2.1_0100520\|MYZPE13164_O_EIv2.1_0106960\|MYZPE13164_O_EIv2.1_0107140\|MYZPE13164_O_EIv2.1_0107850\|MYZPE13164_O_EIv2.1_0111080\|MYZPE13164_O_EIv2.1_0111090\|MYZPE13164_O_EIv2.1_0111120\|MYZPE13164_O_EIv2.1_0123610\|MYZPE13164_O_EIv2.1_0123730\|MYZPE13164_O_EIv2.1_0124690\|MYZPE13164_O_EIv2.1_0124700\|MYZPE13164_O_EIv2.1_0133690\|MYZPE13164_O_EIv2.1_0134410\|MYZPE13164_O_EIv2.1_0135120\|MYZPE13164_O_EIv2.1_0135180\|MYZPE13164_O_EIv2.1_0135320\|MYZPE13164_O_EIv2.1_0135620\|MYZPE13164_O_EIv2.1_0135720\|MYZPE13164_O_EIv2.1_0135840\|MYZPE13164_O_EIv2.1_0136290\|MYZPE13164_O_EIv2.1_0136330\|MYZPE13164_O_EIv2.1_0136390\|MYZPE13164_O_EIv2.1_0136470\|MYZPE13164_O_EIv2.1_0136520\|MYZPE13164_O_EIv2.1_0136530\|MYZPE13164_O_EIv2.1_0136540\|MYZPE13164_O_EIv2.1_0136550\|MYZPE13164_O_EIv2.1_0137070\|MYZPE13164_O_EIv2.1_0138160\|MYZPE13164_O_EIv2.1_0138430\|MYZPE13164_O_EIv2.1_0138570\|MYZPE13164_O_EIv2.1_0138680\|MYZPE13164_O_EIv2.1_0138800\|MYZPE13164_O_EIv2.1_0138820\|MYZPE13164_O_EIv2.1_0139180\|MYZPE13164_O_EIv2.1_0139910\|MYZPE13164_O_EIv2.1_0139920\|MYZPE13164_O_EIv2.1_0140600\|MYZPE13164_O_EIv2.1_0140720\|MYZPE13164_O_EIv2.1_0140880\|MYZPE13164_O_EIv2.1_0140910\|MYZPE13164_O_EIv2.1_0141310\|MYZPE13164_O_EIv2.1_0141320\|MYZPE13164_O_EIv2.1_0141340\|MYZPE13164_O_EIv2.1_0141350\|MYZPE13164_O_EIv2.1_0141360\|MYZPE13164_O_EIv2.1_0142140\|MYZPE13164_O_EIv2.1_0143240\|MYZPE13164_O_EIv2.1_0144130\|MYZPE13164_O_EIv2.1_0144560\|MYZPE13164_O_EIv2.1_0144610\|MYZPE13164_O_EIv2.1_0145380\|MYZPE13164_O_EIv2.1_0145510\|MYZPE13164_O_EIv2.1_0145520\|MYZPE13164_O_EIv2.1_0145590\|MYZPE13164_O_EIv2.1_0146220\|MYZPE13164_O_EIv2.1_0146520\|MYZPE13164_O_EIv2.1_0146530\|MYZPE13164_O_EIv2.1_0146650\|MYZPE13164_O_EIv2.1_0147220\|MYZPE13164_O_EIv2.1_0147240\|MYZPE13164_O_EIv2.1_0147400\|MYZPE13164_O_EIv2.1_0147410\|MYZPE13164_O_EIv2.1_0148500\|MYZPE13164_O_EIv2.1_0148540\|MYZPE13164_O_EIv2.1_0148810\|MYZPE13164_O_EIv2.1_0149720\|MYZPE13164_O_EIv2.1_0151640\|MYZPE13164_O_EIv2.1_0151650\|MYZPE13164_O_EIv2.1_0151660\|MYZPE13164_O_EIv2.1_0151700\|MYZPE13164_O_EIv2.1_0151720\|MYZPE13164_O_EIv2.1_0151730\|MYZPE13164_O_EIv2.1_0152360\|MYZPE13164_O_EIv2.1_0152530\|MYZPE13164_O_EIv2.1_0152540\|MYZPE13164_O_EIv2.1_0152630\|MYZPE13164_O_EIv2.1_0153210\|MYZPE13164_O_EIv2.1_0153400\|MYZPE13164_O_EIv2.1_0153420\|MYZPE13164_O_EIv2.1_0154480\|MYZPE13164_O_EIv2.1_0154750\|MYZPE13164_O_EIv2.1_0154880\|MYZPE13164_O_EIv2.1_0155220\|MYZPE13164_O_EIv2.1_0155260\|MYZPE13164_O_EIv2.1_0159250\|MYZPE13164_O_EIv2.1_0159280\|MYZPE13164_O_EIv2.1_0160140\|MYZPE13164_O_EIv2.1_0163620\|MYZPE13164_O_EIv2.1_0163630\|MYZPE13164_O_EIv2.1_0163650\|MYZPE13164_O_EIv2.1_0163800\|MYZPE13164_O_EIv2.1_0163810\|MYZPE13164_O_EIv2.1_0164300\|MYZPE13164_O_EIv2.1_0164440\|MYZPE13164_O_EIv2.1_0164630\|MYZPE13164_O_EIv2.1_0164660\|MYZPE13164_O_EIv2.1_0164710\|MYZPE13164_O_EIv2.1_0165190\|MYZPE13164_O_EIv2.1_0165320\|MYZPE13164_O_EIv2.1_0165390\|MYZPE13164_O_EIv2.1_0166420\|MYZPE13164_O_EIv2.1_0166450\|MYZPE13164_O_EIv2.1_0167360\|MYZPE13164_O_EIv2.1_0167600\|MYZPE13164_O_EIv2.1_0169070\|MYZPE13164_O_EIv2.1_0169320\|MYZPE13164_O_EIv2.1_0169370\|MYZPE13164_O_EIv2.1_0169450\|MYZPE13164_O_EIv2.1_0169480\|MYZPE13164_O_EIv2.1_0171150\|MYZPE13164_O_EIv2.1_0171450\|MYZPE13164_O_EIv2.1_0171820\|MYZPE13164_O_EIv2.1_0171940\|MYZPE13164_O_EIv2.1_0172040\|MYZPE13164_O_EIv2.1_0172050\|MYZPE13164_O_EIv2.1_0173370\|MYZPE13164_O_EIv2.1_0173410\|MYZPE13164_O_EIv2.1_0174480\|MYZPE13164_O_EIv2.1_0174560\|MYZPE13164_O_EIv2.1_0175570\|MYZPE13164_O_EIv2.1_0177210\|MYZPE13164_O_EIv2.1_0178420\|MYZPE13164_O_EIv2.1_0178600\|MYZPE13164_O_EIv2.1_0178740\|MYZPE13164_O_EIv2.1_0179580\|MYZPE13164_O_EIv2.1_0179840\|MYZPE13164_O_EIv2.1_0181430\|MYZPE13164_O_EIv2.1_0181660\|MYZPE13164_O_EIv2.1_0181980\|MYZPE13164_O_EIv2.1_0182210\|MYZPE13164_O_EIv2.1_0182260\|MYZPE13164_O_EIv2.1_0182520\|MYZPE13164_O_EIv2.1_0183050\|MYZPE13164_O_EIv2.1_0183710\|MYZPE13164_O_EIv2.1_0183960\|MYZPE13164_O_EIv2.1_0183970\|MYZPE13164_O_EIv2.1_0184720\|MYZPE13164_O_EIv2.1_0185380\|MYZPE13164_O_EIv2.1_0185430\|MYZPE13164_O_EIv2.1_0186140\|MYZPE13164_O_EIv2.1_0186640\|MYZPE13164_O_EIv2.1_0186650\|MYZPE13164_O_EIv2.1_0187090\|MYZPE13164_O_EIv2.1_0187420\|MYZPE13164_O_EIv2.1_0188030\|MYZPE13164_O_EIv2.1_0188050\|MYZPE13164_O_EIv2.1_0188720\|MYZPE13164_O_EIv2.1_0189060\|MYZPE13164_O_EIv2.1_0190050\|MYZPE13164_O_EIv2.1_0190480\|MYZPE13164_O_EIv2.1_0191000\|MYZPE13164_O_EIv2.1_0193570\|MYZPE13164_O_EIv2.1_0195210\|MYZPE13164_O_EIv2.1_0195220\|MYZPE13164_O_EIv2.1_0195600\|MYZPE13164_O_EIv2.1_0196170\|MYZPE13164_O_EIv2.1_0199840\|MYZPE13164_O_EIv2.1_0200020\|MYZPE13164_O_EIv2.1_0200040\|MYZPE13164_O_EIv2.1_0200140\|MYZPE13164_O_EIv2.1_0200400\|MYZPE13164_O_EIv2.1_0200950\|MYZPE13164_O_EIv2.1_0200980\|MYZPE13164_O_EIv2.1_0201680\|MYZPE13164_O_EIv2.1_0204830\|MYZPE13164_O_EIv2.1_0206480\|MYZPE13164_O_EIv2.1_0206630\|MYZPE13164_O_EIv2.1_0206880\|MYZPE13164_O_EIv2.1_0206970\|MYZPE13164_O_EIv2.1_0207040\|MYZPE13164_O_EIv2.1_0207270\|MYZPE13164_O_EIv2.1_0207370\|MYZPE13164_O_EIv2.1_0207560\|MYZPE13164_O_EIv2.1_0207890\|MYZPE13164_O_EIv2.1_0207900\|MYZPE13164_O_EIv2.1_0208060\|MYZPE13164_O_EIv2.1_0208180\|MYZPE13164_O_EIv2.1_0209200\|MYZPE13164_O_EIv2.1_0209240\|MYZPE13164_O_EIv2.1_0209300\|MYZPE13164_O_EIv2.1_0209730\|MYZPE13164_O_EIv2.1_0210180\|MYZPE13164_O_EIv2.1_0210270\|MYZPE13164_O_EIv2.1_0210920\|MYZPE13164_O_EIv2.1_0210940\|MYZPE13164_O_EIv2.1_0210990\|MYZPE13164_O_EIv2.1_0211040\|MYZPE13164_O_EIv2.1_0211490\|MYZPE13164_O_EIv2.1_0211990\|MYZPE13164_O_EIv2.1_0212120\|MYZPE13164_O_EIv2.1_0213140\|MYZPE13164_O_EIv2.1_0213480\|MYZPE13164_O_EIv2.1_0213680\|MYZPE13164_O_EIv2.1_0213980\|MYZPE13164_O_EIv2.1_0216020\|MYZPE13164_O_EIv2.1_0216250\|MYZPE13164_O_EIv2.1_0216310\|MYZPE13164_O_EIv2.1_0216320\|MYZPE13164_O_EIv2.1_0216330\|MYZPE13164_O_EIv2.1_0216740\|MYZPE13164_O_EIv2.1_0217290\|MYZPE13164_O_EIv2.1_0217510\|MYZPE13164_O_EIv2.1_0217890\|MYZPE13164_O_EIv2.1_0218140\|MYZPE13164_O_EIv2.1_0219900\|MYZPE13164_O_EIv2.1_0220560\|MYZPE13164_O_EIv2.1_0220790\|MYZPE13164_O_EIv2.1_0220980\|MYZPE13164_O_EIv2.1_0221360\|MYZPE13164_O_EIv2.1_0221370\|MYZPE13164_O_EIv2.1_0221890\|MYZPE13164_O_EIv2.1_0222540\|MYZPE13164_O_EIv2.1_0223090\|MYZPE13164_O_EIv2.1_0223940\|MYZPE13164_O_EIv2.1_0224920\|MYZPE13164_O_EIv2.1_0225090\|MYZPE13164_O_EIv2.1_0226260\|MYZPE13164_O_EIv2.1_0226530\|MYZPE13164_O_EIv2.1_0227980\|MYZPE13164_O_EIv2.1_0228650\|MYZPE13164_O_EIv2.1_0228810\|MYZPE13164_O_EIv2.1_0228840\|MYZPE13164_O_EIv2.1_0229070\|MYZPE13164_O_EIv2.1_0229570\|MYZPE13164_O_EIv2.1_0230870\|MYZPE13164_O_EIv2.1_0231260\|MYZPE13164_O_EIv2.1_0231810\|MYZPE13164_O_EIv2.1_0232600\|MYZPE13164_O_EIv2.1_0233810\|MYZPE13164_O_EIv2.1_0234500\|MYZPE13164_O_EIv2.1_0234520\|MYZPE13164_O_EIv2.1_0235230\|MYZPE13164_O_EIv2.1_0235660\|MYZPE13164_O_EIv2.1_0236470\|MYZPE13164_O_EIv2.1_0236800\|MYZPE13164_O_EIv2.1_0238080\|MYZPE13164_O_EIv2.1_0241310\|MYZPE13164_O_EIv2.1_0241990\|MYZPE13164_O_EIv2.1_0242890\|MYZPE13164_O_EIv2.1_0242910\|MYZPE13164_O_EIv2.1_0242970\|MYZPE13164_O_EIv2.1_0243080\|MYZPE13164_O_EIv2.1_0243090\|MYZPE13164_O_EIv2.1_0243160\|MYZPE13164_O_EIv2.1_0245060\|MYZPE13164_O_EIv2.1_0246590\|MYZPE13164_O_EIv2.1_0246720\|MYZPE13164_O_EIv2.1_0248170\|MYZPE13164_O_EIv2.1_0248740\|MYZPE13164_O_EIv2.1_0249490\|MYZPE13164_O_EIv2.1_0249680\|MYZPE13164_O_EIv2.1_0250270\|MYZPE13164_O_EIv2.1_0250940\|MYZPE13164_O_EIv2.1_0251130\|MYZPE13164_O_EIv2.1_0253480\|MYZPE13164_O_EIv2.1_0253730\|MYZPE13164_O_EIv2.1_0254200\|MYZPE13164_O_EIv2.1_0254240\|MYZPE13164_O_EIv2.1_0254250\|MYZPE13164_O_EIv2.1_0254760\|MYZPE13164_O_EIv2.1_0255000\|MYZPE13164_O_EIv2.1_0256470\|MYZPE13164_O_EIv2.1_0258410\|MYZPE13164_O_EIv2.1_0258780\|MYZPE13164_O_EIv2.1_0259530\|MYZPE13164_O_EIv2.1_0259540\|MYZPE13164_O_EIv2.1_0259560\|MYZPE13164_O_EIv2.1_0259700\|MYZPE13164_O_EIv2.1_0260190\|MYZPE13164_O_EIv2.1_0260300\|MYZPE13164_O_EIv2.1_0260310\|MYZPE13164_O_EIv2.1_0261410\|MYZPE13164_O_EIv2.1_0261420\|MYZPE13164_O_EIv2.1_0262120\|MYZPE13164_O_EIv2.1_0266490\|MYZPE13164_O_EIv2.1_0266500\|MYZPE13164_O_EIv2.1_0266760\|MYZPE13164_O_EIv2.1_0267430\|MYZPE13164_O_EIv2.1_0267480\|MYZPE13164_O_EIv2.1_0268430\|MYZPE13164_O_EIv2.1_0268970\|MYZPE13164_O_EIv2.1_0269400\|MYZPE13164_O_EIv2.1_0270740\|MYZPE13164_O_EIv2.1_0271340\|MYZPE13164_O_EIv2.1_0272280\|MYZPE13164_O_EIv2.1_0272350\|MYZPE13164_O_EIv2.1_0272650\|MYZPE13164_O_EIv2.1_0272700\|MYZPE13164_O_EIv2.1_0273350\|MYZPE13164_O_EIv2.1_0273540\|MYZPE13164_O_EIv2.1_0275060\|MYZPE13164_O_EIv2.1_0275340\|MYZPE13164_O_EIv2.1_0275790\|MYZPE13164_O_EIv2.1_0276040\|MYZPE13164_O_EIv2.1_0276660\|MYZPE13164_O_EIv2.1_0277210\|MYZPE13164_O_EIv2.1_0277330\|MYZPE13164_O_EIv2.1_0279010\|MYZPE13164_O_EIv2.1_0279090\|MYZPE13164_O_EIv2.1_0279280\|MYZPE13164_O_EIv2.1_0279710\|MYZPE13164_O_EIv2.1_0279720\|MYZPE13164_O_EIv2.1_0280140\|MYZPE13164_O_EIv2.1_0281260\|MYZPE13164_O_EIv2.1_0281280\|MYZPE13164_O_EIv2.1_0281380\|MYZPE13164_O_EIv2.1_0281480\|MYZPE13164_O_EIv2.1_0282080\|MYZPE13164_O_EIv2.1_0282100\|MYZPE13164_O_EIv2.1_0282500\|MYZPE13164_O_EIv2.1_0282560\|MYZPE13164_O_EIv2.1_0282590\|MYZPE13164_O_EIv2.1_0282600\|MYZPE13164_O_EIv2.1_0283700\|MYZPE13164_O_EIv2.1_0283790\|MYZPE13164_O_EIv2.1_0284230\|MYZPE13164_O_EIv2.1_0284250\|MYZPE13164_O_EIv2.1_0284700\|MYZPE13164_O_EIv2.1_0285760\|MYZPE13164_O_EIv2.1_0285960\|MYZPE13164_O_EIv2.1_0287290\|MYZPE13164_O_EIv2.1_0287910\|MYZPE13164_O_EIv2.1_0287930\|MYZPE13164_O_EIv2.1_0289020\|MYZPE13164_O_EIv2.1_0289480\|MYZPE13164_O_EIv2.1_0290630\|MYZPE13164_O_EIv2.1_0292910\|MYZPE13164_O_EIv2.1_0292940\|MYZPE13164_O_EIv2.1_0292960\|MYZPE13164_O_EIv2.1_0293740\|MYZPE13164_O_EIv2.1_0294670\|MYZPE13164_O_EIv2.1_0295400\|MYZPE13164_O_EIv2.1_0295790\|MYZPE13164_O_EIv2.1_0295860\|MYZPE13164_O_EIv2.1_0295880\|MYZPE13164_O_EIv2.1_0296460\|MYZPE13164_O_EIv2.1_0297030\|MYZPE13164_O_EIv2.1_0297270\|MYZPE13164_O_EIv2.1_0298300\|MYZPE13164_O_EIv2.1_0298620\|MYZPE13164_O_EIv2.1_0300870\|MYZPE13164_O_EIv2.1_0300970\|MYZPE13164_O_EIv2.1_0301280\|MYZPE13164_O_EIv2.1_0302000\|MYZPE13164_O_EIv2.1_0302010\|MYZPE13164_O_EIv2.1_0303250\|MYZPE13164_O_EIv2.1_0304300\|MYZPE13164_O_EIv2.1_0304970\|MYZPE13164_O_EIv2.1_0306100\|MYZPE13164_O_EIv2.1_0306680\|MYZPE13164_O_EIv2.1_0309020\|MYZPE13164_O_EIv2.1_0309030\|MYZPE13164_O_EIv2.1_0309330\|MYZPE13164_O_EIv2.1_0310090\|MYZPE13164_O_EIv2.1_0311040\|MYZPE13164_O_EIv2.1_0312120\|MYZPE13164_O_EIv2.1_0317580\|MYZPE13164_O_EIv2.1_0318530\|MYZPE13164_O_EIv2.1_0318960\|MYZPE13164_O_EIv2.1_0320460\|MYZPE13164_O_EIv2.1_0320470\|MYZPE13164_O_EIv2.1_0323070\|MYZPE13164_O_EIv2.1_0325490\|MYZPE13164_O_EIv2.1_0325510\|MYZPE13164_O_EIv2.1_0326160\|MYZPE13164_O_EIv2.1_0326560\|MYZPE13164_O_EIv2.1_0328570\|MYZPE13164_O_EIv2.1_0328670\|MYZPE13164_O_EIv2.1_0329710\|MYZPE13164_O_EIv2.1_0332210\|MYZPE13164_O_EIv2.1_0333430\|MYZPE13164_O_EIv2.1_0333880\|MYZPE13164_O_EIv2.1_0333970\|MYZPE13164_O_EIv2.1_0336300\|MYZPE13164_O_EIv2.1_0336990\|MYZPE13164_O_EIv2.1_0337240\|MYZPE13164_O_EIv2.1_0337250\|MYZPE13164_O_EIv2.1_0337260\|MYZPE13164_O_EIv2.1_0337270\|MYZPE13164_O_EIv2.1_0337580\|MYZPE13164_O_EIv2.1_0337590\|MYZPE13164_O_EIv2.1_0338150\|MYZPE13164_O_EIv2.1_0338940\|MYZPE13164_O_EIv2.1_0339430\|MYZPE13164_O_EIv2.1_0339930\|MYZPE13164_O_EIv2.1_0340140\|MYZPE13164_O_EIv2.1_0341080\|MYZPE13164_O_EIv2.1_0342070\|MYZPE13164_O_EIv2.1_0343870\|MYZPE13164_O_EIv2.1_0343990\|MYZPE13164_O_EIv2.1_0345790\|MYZPE13164_O_EIv2.1_0345820\|MYZPE13164_O_EIv2.1_0345840\|MYZPE13164_O_EIv2.1_0345860\|MYZPE13164_O_EIv2.1_0345880\|MYZPE13164_O_EIv2.1_0346220\|MYZPE13164_O_EIv2.1_0346420\|MYZPE13164_O_EIv2.1_0347830\|MYZPE13164_O_EIv2.1_0347890\|MYZPE13164_O_EIv2.1_0348170\|MYZPE13164_O_EIv2.1_0348220\|MYZPE13164_O_EIv2.1_0348550\|MYZPE13164_O_EIv2.1_0348900\|MYZPE13164_O_EIv2.1_0349310\|MYZPE13164_O_EIv2.1_0349540\|MYZPE13164_O_EIv2.1_0349610\|MYZPE13164_O_EIv2.1_0350040\|MYZPE13164_O_EIv2.1_0350050\|MYZPE13164_O_EIv2.1_0350350\|MYZPE13164_O_EIv2.1_0350610\|MYZPE13164_O_EIv2.1_0350640\|MYZPE13164_O_EIv2.1_0351370\|MYZPE13164_O_EIv2.1_0351790\|MYZPE13164_O_EIv2.1_0353780\|MYZPE13164_O_EIv2.1_0353810\|MYZPE13164_O_EIv2.1_0354030\|MYZPE13164_O_EIv2.1_0357250\|MYZPE13164_O_EIv2.1_0357590\|MYZPE13164_O_EIv2.1_0358240\|MYZPE13164_O_EIv2.1_0358520\|MYZPE13164_O_EIv2.1_0358840\|MYZPE13164_O_EIv2.1_0358890\|MYZPE13164_O_EIv2.1_0359910\|MYZPE13164_O_EIv2.1_0359970\|MYZPE13164_O_EIv2.1_0361760\|MYZPE13164_O_EIv2.1_0362100\|MYZPE13164_O_EIv2.1_0362450\|MYZPE13164_O_EIv2.1_0362830\|MYZPE13164_O_EIv2.1_0362930\|MYZPE13164_O_EIv2.1_0362940\|MYZPE13164_O_EIv2.1_0363870\|MYZPE13164_O_EIv2.1_0364310\|MYZPE13164_O_EIv2.1_0364360\|MYZPE13164_O_EIv2.1_0365150\|MYZPE13164_O_EIv2.1_0365510\|MYZPE13164_O_EIv2.1_0365860\|MYZPE13164_O_EIv2.1_0365920\|MYZPE13164_O_EIv2.1_0366100\|MYZPE13164_O_EIv2.1_0366700\|MYZPE13164_O_EIv2.1_0366730\|MYZPE13164_O_EIv2.1_0366970\|MYZPE13164_O_EIv2.1_0367250\|MYZPE13164_O_EIv2.1_0367640\|MYZPE13164_O_EIv2.1_0367680\|MYZPE13164_O_EIv2.1_0367740\|MYZPE13164_O_EIv2.1_0367750\|MYZPE13164_O_EIv2.1_0367760\|MYZPE13164_O_EIv2.1_0367810\|MYZPE13164_O_EIv2.1_0367980\|MYZPE13164_O_EIv2.1_0369920\|sequence-region\|gff-version' snp_calling/Myzus/persicae/biello/gatk/filtered/MYZPE13164_O_EIv2.1.annotation.gff3 > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/effector_candidates.gff3

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/pybed.simg bedtools intersect \
-a snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz \
-b /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/effector_candidates.gff3 \
-header > snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-effector-genic-regions.vcf
 source package /nbi/software/
bcftools view -v snps snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-effector-genic-regions.vcf | grep -c -v "^#"
#1,112,435
source package 01ef5a53-c149-4c9e-b07d-0b9a46176cc0
bgzip -c snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-effector-genic-regions.vcf > snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-effector-genic-regions.vcf.gz
rm snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-effector-genic-regions.vcf

for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs-effector-genic-regions.vcf.gz); do
    InFile=$vcf
    OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity
    OutFile=193_effector_genic
    Exclusion_list=
    ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_vcfstats.sh $InFile $OutDir $OutFile $Exclusion_list
done #55250981
```
Plot VCFStats
```R
setwd("C:/Users/did23faz/OneDrive - Norwich Bioscience Institutes/Desktop/R")

install.packages("tidyverse")
install.packages("ggplot2")

library(tidyverse)
library(ggplot2)

var_qual <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_effector_genic_site_qual.lqual", delim = "\t",
                       col_names = c("CHROM", "POS", "QUAL"), skip = 1)

var_depth <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_effector_genic_site-mean-depth.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

var_miss <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_effector_genic_missing_sites.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

var_freq <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_effector_genic_allele_freq.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

var_count <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_effector_genic_allele_count.frq.count", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "count", "a2"), skip = 1)

ind_depth <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_effector_genic_depth.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

ind_miss  <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_effector_genic_missing_ind.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

ind_het <- read_delim("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/SNP_diversity/193_effector_genic_het.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)


p1 <- ggplot(var_qual, aes(QUAL)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Effector SNPs: Per Site Quality")
p1 <- p1 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_effector_genic_var_qual.pdf", width=11, height=7)
plot(p1, pdf=T)
dev.off()

p2 <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + labs(title = "Effector SNPs: Mean Depth Per Site")
p2 <- p2 + theme(plot.title = element_text(hjust = 0.5))
summary(var_depth$mean_depth)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  8.249  16.212  19.529  18.710  21.285  30.829 
pdf("193_effector_genic_var_depth.pdf", width=11, height=7)
plot(p2, pdf=T)
dev.off()

p3 <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + xlim(0, 100) + labs(title = "Effector SNPs: Mean Depth Per Site")
p3 <- p3 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_effector_genic_var_depth_x100.pdf", width=11, height=7)
plot(p3, pdf=T)
dev.off()

p4 <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Effector SNPs: Missingness Per Site")
p4 <- p4 + theme(plot.title = element_text(hjust = 0.5))
summary(var_miss$fmiss)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.005181 0.010363 0.016646 0.025907 0.108808 
pdf("193_effector_genic_var_miss.pdf", width=11, height=7)
plot(p4, pdf=T)
dev.off()

var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
summary(var_freq$maf)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.00000 0.00000 0.01121 0.00000 0.50000 
p5 <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Effector SNPs: Allele Frequency for Each Site")
p5 <- p5 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_effector_genic_var_freq.pdf", width=11, height=7)
plot(p5, pdf=T)
dev.off()

p6 <- ggplot(var_count, aes(x = count)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Effector SNPs: Raw Allele Counts for Each Site")
p6 <- p6 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_effector_genic_allele_count_density.pdf", width = 11, height = 7)
plot(p6, pdf=T)
dev.off()

p7 <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Effector SNPs: Mean Depth Per Individual")
p7 <- p7 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_effector_genic_ind_depth.pdf", width=11, height=7)
plot(p7, pdf=T)
dev.off()

p8 <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Effector SNPs: Missingness Per Individual")
p8 <- p8 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_effector_genic_ind_miss.pdf", width=11, height=7)
plot(p8, pdf=T)
dev.off()

p9 <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + labs(title = "Effector SNPs: Heterozygosity Per Individual")
p9 <- p9 + theme(plot.title = element_text(hjust = 0.5))
pdf("193_effector_genic_ind_het.pdf", width=11, height=7)
plot(p9, pdf=T)
dev.off()
```

















Check for SNPs of clone O versus itself to control for diploidy
```bash
for RawData in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/wouters/O/Mp_O_ERR1145176_R*.fastq.gz); do
echo $RawData
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $RawData)
Outfile=$(basename -a -s .fastq.gz $RawData)
sbatch $ProgDir/run_fastqc.sh $RawData $OutDir $Outfile
done #55377790,91

source package /nbi/software/production/bin/python-2.7.11
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/wouters/O);
do
 read1=$(ls $ReadDir/*R1.fastq.gz)
 read2=$(ls $ReadDir/*R2.fastq.gz)
ls $read1
ls $read2
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(echo $ReadDir|sed 's@rawdata@dna_qc@g')
Prefix=$(echo $ReadDir|cut -f9,10,12 -d '/' --output-delimiter '-')
echo $Prefix
mkdir -p $OutDir
sbatch $ProgDir/run_trim_galore.sh $read1 $read2 $OutDir $Prefix
done #55410035

source package /nbi/software/testing/bin/vcftools-0.1.15
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz --indv O --counts
grep -v "^#" out.frq.count| cut -f3 | grep -cv "0"

source package /nbi/software/testing/bin/bcftools-1.8
bcftools view -c1 -s O -Ou /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz | bcftools query -e 'REF=ALT' -f '%CHROM\t%POS\n' | wc -l
#272,953
bcftools view -c1 -s O -Ou /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz | bcftools query -e 'REF=ALT' -f '%CHROM\t%POS\n' > temp.vcf
source package 01ef5a53-c149-4c9e-b07d-0b9a46176cc0
bgzip -cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz > 193s.M_persicae.onlySNPs.vcf
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
# Read the first file and store the values in a set
first_file_lines = set()
with open('temp.vcf', 'r') as f1:
    for line in f1:
        columns = line.strip().split('\t')
        if len(columns) >= 2:
            first_file_lines.add((columns[0], columns[1]))

# Read the second file and write lines that don't match the criteria to a new file
with open('193s.M_persicae.onlySNPs.vcf', 'r') as f2, open('192s.M_persicae.onlySNPs.vcf', 'w') as output:
    for line in f2:
        columns = line.strip().split('\t')
        if len(columns) >= 2 and (columns[0], columns[1]) not in first_file_lines:
            output.write(line)

exit()
grep '#' 193s.M_persicae.onlySNPs.vcf >> 193s-no0.M_persicae.onlySNPs.vcf
grep -v '#' 192s.M_persicae.onlySNPs.vcf >> 193s-no0.M_persicae.onlySNPs.vcf
rm 192s.M_persicae.onlySNPs.vcf 193s.M_persicae.onlySNPs.vcf
bgzip 193s-no0.M_persicae.onlySNPs.vcf
mv 193s-no0.M_persicae.onlySNPs.vcf.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/.


bcftools view -s O -Ou /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz | bcftools query -e 'REF=ALT' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT\t%O[\t]\n' > temp.vcf

bcftools query -f '%CHROM\t%POS\n' -i 'ALT!="."' -s O /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz > snps_to_remove.txt

source package 01ef5a53-c149-4c9e-b07d-0b9a46176cc0
bgzip -cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz > 193s.M_persicae.onlySNPs.vcf
bcftools view -s O -v snps -f 'GT="0/0"' 193s.M_persicae.onlySNPs.vcf > filtered_file.vcf
rm 193s.M_persicae.onlySNPs.vcf
bgzip filtered_file.vcf
tabix -p vcf filtered_file.vcf.gz
echo done








bcftools view -H -s O -v snps /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz | cut -f1,2 > snps_to_remove.txt

bcftools view -H -s O -v snps -u /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz | bcftools query -f '%CHROM\t%POS\n' > snps_to_remove.txt
bcftools view -H -s O -v snps /jic/research-groups/Saskia-Hogenhout/TCHeaven/PopGen/M_persicae_SNP_population/193s.M_persicae.onlySNPs.vcf.gz | bcftools query -f '%CHROM\t%POS\n' > snps_to_remove.txt


bcftools view -T ^snps_to_remove.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/193s.M_persicae.onlySNPs.vcf.gz -Oz -o filtered.vcf.gz


```
```bash
#find gene length:
for line in $(cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_snp_report.txt); do
gene=$(echo $line | cut -d '_' -f2,3,4,5)
length=$(sed -n '2p' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/homo_gene_fastas/hom_${gene}.fa | wc -c)
echo ${line},${length} >> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_snp_report2.txt
done
sed -i 's/dedup_MYZPE13164_O_EIv2.1_//g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_snp_report2.txt
sed -i 's/_snps.vcf//g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_snp_report2.txt
sed -i 's/dedup_//g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/gene_snp_report2.txt
```
#### admixture
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

Pruned_vcf=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/plink/193s.M_persicae.onlySNPs_sorted_pruned_set
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/admixture
OutFile=193s.M_persicae.genomicSNPs
Mink=2
Maxk=3
Bootstraps=200
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_admixture.sh $OutDir $OutFile $Pruned_vcf $Mink $Maxk $Bootstraps #55755934
```
#### FST
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
```bash
/jic/research-groups/Saskia-Hogenhout/reads/genomic/nanopore/Clone_O/cat_08_09_10_11_guppy.min_1kb.fastq
pilon_r3_curated_cov_90.fasta
mkdir Clone_O
cd Clone_O
~/git_repos/Wrappers/NBI/submit-slurm_v1.1.pl -q jic-medium -m 100000 -c 32 -t 1-00:00 -e -j M_perO_wtdbg_p0_k17_L15kb -i "source wtdbg2-2.3;wtdbg2 -x ont -p 0 -k 17 -L 15000 -g 409m -t 32 -i /jic/research-groups/Saskia-Hogenhout/reads/genomic/nanopore/Clone_O/cat_08_09_10_11_guppy.min_1kb.fastq.gz -fo wtdbg2;wtpoa-cns -t 32 -i wtdbg2.ctg.lay.gz -fo wtdbg2_1.ctg.fa" #55681945
~/git_repos/Wrappers/NBI/submit-slurm_v1.1.pl -q jic-medium -m 100000 -c 32 -t 1-00:00 -e -j M_perO_wtdbg_p19_k0_L15kb -i "source wtdbg2-2.3;wtdbg2 -x ont -p 19 -k 0 -L 15000 -g 409m -t 32 -i /jic/research-groups/Saskia-Hogenhout/reads/genomic/nanopore/Clone_O/cat_08_09_10_11_guppy.min_1kb.fastq.gz -fo wtdbg2;wtpoa-cns -t 32 -i wtdbg2.ctg.lay.gz -fo wtdbg2_2.ctg.fa" #55681946
ln -s /tgac/software/production/MUMmer/3.23/x86_64/bin/delta-filter ~/progs/quickmerge/.
ln -s /software/391fb775-22fd-4b97-9e2a-6d213794cc1d/bin/quickmerge ~/progs/quickmerge/.
ln -s /tgac/software/production/bin/core/../..//MUMmer/3.23/x86_64/bin/nucmer ~/progs/quickmerge/.
#~/git_repos/Wrappers/NBI/submit-slurm_v1.1.pl -q jic-long -m 75000 -c 4 -t 7-00:00 -e -j quickmerge_wtdbg2_M_perO -i "export PATH=/hpc-home/did23faz/progs/quickmerge;source package 391fb775-22fd-4b97-9e2a-6d213794cc1d;source package /tgac/software/production/bin/MUMmer-3.23;singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/merge_wrapper.py -pre merged -l 1837291 -ml 10000 wtdbg2_1.ctg.fa wtdbg2_2.ctg.fa" 
~/git_repos/Wrappers/NBI/submit-slurm_v1.1.pl -q jic-long -m 75000 -c 4 -t 7-00:00 -e -j quickmerge_wtdbg2_M_perO -i "export PATH=/hpc-home/did23faz/progs/quickmerge:/hpc-home/did23faz/progs/quickmerge/MUMmer3.23:$PATH;merge_wrapper.py -pre merged -l 1837291 -ml 10000 wtdbg2_1.ctg.fa wtdbg2_2.ctg.fa" 
~/git_repos/Wrappers/NBI/submit-slurm_v1.1.pl -q jic-long -m 200000 -c 16 -t 2-00:00 -e -j M_perO_wtdbg_racon_3r_kat_comp_busco -i "~/git_repos/Wrappers/NBI/racon_3_rounds_pilon_3_rounds_kat_comp_busco.sh merged_merged.fasta /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/SRR11322170.SRA/SRR11322170/SRR11322170_1.fastq.gz /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/SRR11322170.SRA/SRR11322170/SRR11322170_2.fastq.gz /jic/research-groups/Saskia-Hogenhout/reads/genomic/nanopore/Clone_O/cat_08_09_10_11_guppy.min_1kb.fastq.gz map-ont"
Submitted batch job 55755844


source package /nbi/software/testing/bin/racon-1.3.1


cd /usr/users/JIC_a1/tmathers
find . -name "*Mp_057_R2_val_2*" 2>/dev/null

cd /jic/scratch/groups/Saskia-Hogenhout
find . -name "*Mp_057_R2_val_2*" 2>/dev/null

cd ~/../..
find . -name "*Mp_057_R2_val_2*" 2>/dev/null

#PE PCR free assembly
~/BBSRC_aphid_adaptation/Clone_O_assembly_improvement/discovar_lib/reads/Mp_057_R1_val_1.fq ~/BBSRC_aphid_adaptation/Clone_O_assembly_improvement/discovar_lib/reads/Mp_057_R2_val_2.fq

"export PATH=/hpc-home/tmathers/quickmerge-0.3:/hpc-home/tmathers/quickmerge-0.3/MUMmer3.23:~/git_repos/Scripts/NBI/merge_wrapper.py -pre merged -l 1837291 -ml 10000 wtdbg2_1.ctg.fa  wtdbg2_2.ctg.fa" #55687700
delta-filter', and 'merger


submit-slurm_v1.1.pl -q jic-largemem -m 1000000 -c 32 -t 7-00:00 -e -j M_perO_flye_default -i "source anaconda3-5.2.0;source activate flye; flye --nano-raw /jic/research-groups/Saskia-Hogenhout/reads/genomic/nanopore/Clone_O/cat_08_09_10_11_guppy.min_1kb.fastq -g 409m -o out --threads 32"

[tmathers@NBI-HPC out]$ submit-slurm_v1.1.pl -q ei-medium -m 100000 -c 32 -t 2-00:00 -e -j O_flye_pilon_r3_purge_haplotigs -i "source anaconda3-5.2.0; source activate purge_haplotigs; minimap2 -ax map-ont -t 32 pilon_r3.fasta /jic/research-groups/Saskia-Hogenhout/reads/genomic/nanopore/Clone_O/cat_08_09_10_11_guppy.min_1kb.fastq | samtools view -hF 256 - | samtools sort --threads 32 -m 1G -o aligned.bam -T tmp.ali;purge_haplotigs readhist -b aligned.bam -g pilon_r3.fasta -t 32"
Submitted batch job 19555266

[tmathers@NBI-HPC out]$ submit-slurm_v1.1.pl -q ei-medium -m 100000 -c 32 -t 2-00:00 -e -j O_flye_pilon_r3_purge_haplotigs -i "source anaconda3-5.2.0; source activate purge_haplotigs; purge_haplotigs contigcov -i aligned.bam.gencov -l 9 -m 45  -h 92 -o minimnap2_coverage_stats.csv; purge_haplotigs purge -t 32 -g pilon_r3.fasta -c minimnap2_coverage_stats.csv -a 90 -o pilon_r3_curated_cov_90"
Submitted batch job 19555769

1891 97.3%

```