All files that look like analysis not raw data have been removed/not copied over in order to save space.
```bash
#Organ data raw data:
cp -r /ei/projects/a/a93e4b69-58e3-495b-ac4d-04e978fed5d1/data/organ_data/raw_data/for_archana/* /jic/research-groups/Saskia-Hogenhout/TCHeaven/Raw_Data/archana_organ_data/.

#WGS Raw data:
cp -r /ei/projects/b/bda78a5f-7de9-4e0c-813f-1b1105bd6c24/WGS_Feb2021/* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/.
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021
unzip X201SC21011755-Z01-F001.zip

#WGBS Raw data:
cp /ei/projects/b/bda78a5f-7de9-4e0c-813f-1b1105bd6c24/WGBS_Mar2021/X201SC21011757-Z01-F001_*.zip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGBS/Archana_Mar2021/.

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
