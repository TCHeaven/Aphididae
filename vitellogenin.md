# Vitellogenins
Find vitellogenin genes for Reuben

Download vitellogenin protein sequences from NCBI (It seems there is no reference for M. persicae):
```bash
#Acipenser_transmontanus:>AAA87392.1 vitellogenin, partial [Acipenser transmontanus]
#Aedes_aegypti:>AAA18221.1 vitellogenin [Aedes aegypti]
#Anopheles_gambiae:>AAF82131.1 vitellogenin 1 [Anopheles gambiae]
#Anthonomus_grandis:>AAA27740.1 vitellogenin [Anthonomus grandis]
#Aphis_citricidus:>AVP41182.1 vitellogenin, partial [Aphis citricidus]
#Aphis_gossypii:>QFU80926.1 vitellogenin [Aphis gossypii]
#Apis_mellifera:>CAD56944.1 vitellogenin precursor [Apis mellifera]
#Apolygus_lucorum:>AGT39945.1 vitellogenin [Apolygus lucorum]
#Blattella_germanica:>CAA06379.2 vitellogenin [Blattella germanica]
#Bombyx_mandarina:>BAB32642.2 vitellogenin [Bombyx mandarina]
#Bombyx_mori:>BAA02444.1 vitellogenin precursor [Bombyx mori]
#Caenorhabditis_elegans:>NP_509305.1 Vitellogenin-1 [Caenorhabditis elegans]
#Cimex_lectularius:>BAU36889.1 vitellogenin [Cimex lectularius]
#Culex_quinquefasciatus:>AAV31930.1 vitellogenin C1 [Culex quinquefasciatus]
#Diuraphis_noxia:>XP_015366382.1 PREDICTED: vitellogenin isoform X1 [Diuraphis noxia]
#Encarsia_formosa:>AAT48601.1 vitellogenin [Encarsia formosa]
#Gallus_gallus:>CAA31942.1 vitellogenin [Gallus gallus]
#Homalodisca_vitripennis:>AAZ06771.1 vitellogenin [Homalodisca vitripennis]
#Laodelpha_ striatellus:>AGJ26477.1 vitellogenin [Laodelphax striatellus]
#Mythimna_separata:>AHG29547.1 vitellogenin protein [Mythimna separata]
#Nilaparvat_lugens:>BAF75351.1 vitellogenin [Nilaparvata lugens]
#Periplaneta_americana:>BAB32673.1 vitellogenin-2 [Periplaneta americana]
#Periplaneta_americana2:>BAA86656.1 vitellogenin [Periplaneta americana]
#Plautia_stali:>BAA88077.1 vitellogenin-3 [Plautia stali]
#Plautia_stali2:>BAA88075.1 vitellogenin-1 [Plautia stali]
#Plautia_stali3:>BAA88076.1 vitellogenin-2 [Plautia stali]
#Rhyparobia_maderae:>BAB19327.1 vitellogenin [Rhyparobia maderae]
#Solenopsis_invicta:>AAP47155.1 vitellogenin [Solenopsis invicta]
#Solenopsis_invicta2:>AAY22960.1 vitellogenin-2 [Solenopsis invicta]
#Solenopsis_invicta3:>AAY22961.1 vitellogenin-3 [Solenopsis invicta]
#Spodoptera_litura:>ABU68426.1 vitellogenin [Spodoptera litura]
#Tenebrio_molitor:>AAU20328.2 vitellogenin precursor [Tenebrio molitor]
#Tribolium_castaneum:>XP_971398.1 PREDICTED: vitellogenin [Tribolium castaneum]
#Xenopus_laevis:>AAA49982.1 vitellogenin [Xenopus laevis]
#Zootermopsis_nevadensis:>KDR08462.1 Vitellogenin-1, partial [Zootermopsis nevadensis]
#Zootermopsis_nevadensis2:>KDR08463.1 Vitellogenin-2 [Zootermopsis nevadensis]
```
Build a blast database from the downloaded vitellogenin protein sequences and blast search the predicted proteins (annotations v2.1) from myzus persicae against this database:
```bash
mkdir VT
cat * > VT/vitellogenins.faa
cd VT
source package 37f0ffda-9f66-4391-87e2-38ccd398861d
makeblastdb -in vitellogenins.aa.fa -input_type fasta -dbtype prot  -title VT  -parse_seqids -out VT
blastp -query /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.aa.fa -db VT -out results.txt -evalue 1e-5 -outfmt 6 -num_threads 1

#Hits returned:
#MYZPE13164_O_EIv2.1_0031810.1
#MYZPE13164_O_EIv2.1_0059260.1
#MYZPE13164_O_EIv2.1_0133690.1
#MYZPE13164_O_EIv2.1_0133690.2
#MYZPE13164_O_EIv2.1_0133690.3
#MYZPE13164_O_EIv2.1_0133690.4
#MYZPE13164_O_EIv2.1_0209000.2 #most common hit
#MYZPE13164_O_EIv2.1_0213490.1 #2nd most common hit
#MYZPE13164_O_EIv2.1_0213490.2 #2nd most common hit
#MYZPE13164_O_EIv2.1_0248680.1
#MYZPE13164_O_EIv2.1_0248680.2
#MYZPE13164_O_EIv2.1_0248680.3
#MYZPE13164_O_EIv2.1_0248680.4
#MYZPE13164_O_EIv2.1_0248680.5

awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.aa.fa > MYZPE13164_O_EIv2.1.annotation.gff3.aa.fa
source package 46a62eca-4f8f-45aa-8cc2-d4efc99dd9c6
grep -A 1 'MYZPE13164_O_EIv2.1_0031810.1\|MYZPE13164_O_EIv2.1_0059260.1\|MYZPE13164_O_EIv2.1_0133690.1\|MYZPE13164_O_EIv2.1_0133690.2\|MYZPE13164_O_EIv2.1_0133690.3\|MYZPE13164_O_EIv2.1_0133690.4\|MYZPE13164_O_EIv2.1_0209000.2\|MYZPE13164_O_EIv2.1_0213490.1\|MYZPE13164_O_EIv2.1_0213490.2\|MYZPE13164_O_EIv2.1_0248680.1\|MYZPE13164_O_EIv2.1_0248680.2\|MYZPE13164_O_EIv2.1_0248680.3\|MYZPE13164_O_EIv2.1_0248680.4\|MYZPE13164_O_EIv2.1_0248680.5' MYZPE13164_O_EIv2.1.annotation.gff3.aa.fa > hits.aa.fa
```
Blasting these sequences on NCBI web browser:
```bash
#MYZPE13164_O_EIv2.1_0031810.1 #tops hits are hypothetical and unknown proteins; aphis glycine,  Macrosiphum euphorbiae,     Aphis gossypii, Aphis craccivora etc.
#MYZPE13164_O_EIv2.1_0059260.1 #top hits are hypothetical and unknown proteins;  Myzus persicae, Acyrthosiphon pisum,    Aphis gossypii, Melanaphis sacchari etc.
#MYZPE13164_O_EIv2.1_0133690.1 #top hits are apolipophorins-like proteins and unknown and unclassified proteins; Myzus persicae, Acyrthosiphon pisum, Rhopalosiphum maidis, Diuraphis noxia, aphis glycine etc.
#MYZPE13164_O_EIv2.1_0133690.2 #top hits are apolipophorins-like proteins and unknown and unclassified proteins; Myzus persicae, Acyrthosiphon pisum, Rhopalosiphum maidis, Diuraphis noxia, aphis glycine etc.
#MYZPE13164_O_EIv2.1_0133690.3 #top hits are apolipophorins-like proteins and unknown and unclassified proteins; Myzus persicae, Acyrthosiphon pisum, Rhopalosiphum maidis, Diuraphis noxia, aphis glycine etc.
#MYZPE13164_O_EIv2.1_0133690.4 #top hits are apolipophorins-like proteins and unknown and unclassified proteins; Myzus persicae, Acyrthosiphon pisum, Rhopalosiphum maidis, Diuraphis noxia, aphis glycine etc.
#MYZPE13164_O_EIv2.1_0209000.2 #BMP-binding endothelial regulator protein and unknown and unclassified proteins; Myzus persicae, Macrosiphum euphorbiae, Diuraphis noxia, Acyrthosiphon pisum etc.
#MYZPE13164_O_EIv2.1_0213490.1 #top hits are uncharacterised and vitellogenin proteins; Myzus persicae, Diuraphis noxia, Acyrthosiphon pisum, Macrosiphum euphorbiae etc.
#MYZPE13164_O_EIv2.1_0213490.2 #top hits are uncharacterised and vitellogenin proteins; Myzus persicae, Diuraphis noxia, Acyrthosiphon pisum, Macrosiphum euphorbiae etc.
#MYZPE13164_O_EIv2.1_0248680.1 #tops hits are vacuolar-sorting protein, unnamed protein, ELL complex EAP30 subunit-like, ESCRT-2 complex, Snf8,EAP30,Winged helix-turn-helix DNA-binding domain etc.; Macrosiphum euphorbiae, Myzus persicae, Acyrthosiphon pisum, Diuraphis noxia,  Rhopalosiphum maidis etc.
#MYZPE13164_O_EIv2.1_0248680.2 #tops hits are vacuolar-sorting protein, unnamed protein, ELL complex EAP30 subunit-like, ESCRT-2 complex, Snf8,EAP30,Winged helix-turn-helix DNA-binding domain etc.; Macrosiphum euphorbiae, Myzus persicae, Acyrthosiphon pisum, Diuraphis noxia,  Rhopalosiphum maidis etc.
#MYZPE13164_O_EIv2.1_0248680.3 #tops hits are vacuolar-sorting protein, unnamed protein, ELL complex EAP30 subunit-like, ESCRT-2 complex, Snf8,EAP30,Winged helix-turn-helix DNA-binding domain etc.; Macrosiphum euphorbiae, Myzus persicae, Acyrthosiphon pisum, Diuraphis noxia,  Rhopalosiphum maidis etc.
#MYZPE13164_O_EIv2.1_0248680.4 #tops hits are vacuolar-sorting protein, unnamed protein, ELL complex EAP30 subunit-like, ESCRT-2 complex, Snf8,EAP30,Winged helix-turn-helix DNA-binding domain etc.; Macrosiphum euphorbiae, Myzus persicae, Acyrthosiphon pisum, Diuraphis noxia,  Rhopalosiphum maidis etc.
#MYZPE13164_O_EIv2.1_0248680.5 #tops hits are vacuolar-sorting protein, unnamed protein, ELL complex EAP30 subunit-like, ESCRT-2 complex, Snf8,EAP30,Winged helix-turn-helix DNA-binding domain etc.; Macrosiphum euphorbiae, Myzus persicae, Acyrthosiphon pisum, Diuraphis noxia,  Rhopalosiphum maidis etc.
```
Based upon these blast hits MYZPE13164_O_EIv2.1_0213490 is most likely the M.persicae vitellogenin gene

NOTE: there are two versions of this gene predicted: MYZPE13164_O_EIv2.1_0213490.1 and MYZPE13164_O_EIv2.1_0213490.2; version 1 is 1,328 amino acids in length, version 2 is 1,203 amino acids in length. 
```bash
#Convert nucleotide fasta to single line format:
awk ' {if (NR==1) {print $0} else {if ($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.gene.fa > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation_singleline.gff3.nt.gene.fa #21,011 nucleotides in length; blasts to uncharacterised  Myzus persicae and Acyrthosiphon pisum proteins and 'vitellogenin like' proteins in Melanaphis sacchari and Diuraphis noxia, etc.
awk ' {if (NR==1) {print $0} else {if ($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.nt.CDS.fa > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation_singleline.gff3.nt.CDS.fa

#Reference gene sequence:
grep -A 1 'MYZPE13164_O_EIv2.1_0213490' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation_singleline.gff3.nt.gene.fa > MYZPE13164_O_EIv2.1_0213490.fa
#ReferenceCDS sequence:
grep -A 1 'MYZPE13164_O_EIv2.1_0213490.1' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation_singleline.gff3.nt.CDS.fa > MYZPE13164_O_EIv2.1_0213490.1.fa
grep -A 1 'MYZPE13164_O_EIv2.1_0213490.2' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation_singleline.gff3.nt.CDS.fa > MYZPE13164_O_EIv2.1_0213490.2.fa

#Gene sequence multifasta for 193 individuals:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas/hom_MYZPE13164_O_EIv2.1_0213490.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas/het_MYZPE13164_O_EIv2.1_0213490.fa
#CDS sequence multifasta for 193 individuals:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/hom_MYZPE13164_O_EIv2.1_0213490.1_CDS.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/hom_MYZPE13164_O_EIv2.1_0213490.2_CDS.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas/het_MYZPE13164_O_EIv2.1_0213490.1_CDS.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas/het_MYZPE13164_O_EIv2.1_0213490.2_CDS.fa
#VCF SNP file, gene region, for 193 individuals:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/dedup_MYZPE13164_O_EIv2.1_0213490_snps.vcf
#VCF SNP file, CDS region, for 193 individuals:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/dedup_MYZPE13164_O_EIv2.1_CDS_0213490.1_snps.vcf
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/dedup_MYZPE13164_O_EIv2.1_CDS_0213490.2_snps.vcf

#Convert amino acid fasta to single line format:
#Reference protein sequence:
awk ' {if (NR==1) {print $0} else {if ($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/MYZPE13164_O_EIv2.1.annotation.gff3.aa.fa > /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation_singleline.gff3.aa.fa
grep -A 1 'MYZPE13164_O_EIv2.1_0213490.1' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation_singleline.gff3.aa.fa > vitellogenin_1.aa.fa
grep -A 1 'MYZPE13164_O_EIv2.1_0213490.2' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation_singleline.gff3.aa.fa > vitellogenin_2.aa.fa

source package /nbi/software/production/bin/emboss-6.5.7
transeq -sequence /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/hom_MYZPE13164_O_EIv2.1_0213490.1_CDS.fa -outseq hom_MYZPE13164_O_EIv2.1_0213490.1_CDS.aa.fa -frame 1
transeq -sequence /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/hom_MYZPE13164_O_EIv2.1_0213490.2_CDS.fa -outseq hom_MYZPE13164_O_EIv2.1_0213490.2_CDS.aa.fa -frame 1
transeq -sequence /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas/het_MYZPE13164_O_EIv2.1_0213490.1_CDS.fa -outseq het_MYZPE13164_O_EIv2.1_0213490.1_CDS.aa.fa -frame 1
transeq -sequence /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas/het_MYZPE13164_O_EIv2.1_0213490.2_CDS.fa -outseq het_MYZPE13164_O_EIv2.1_0213490.2_CDS.aa.fa -frame 1
```
### White gene
Reuben says that MYZPE13164_O_EIv2.1_0184430 is the white gene
```bash
#Reference gene sequence:
grep -A 1 'MYZPE13164_O_EIv2.1_0184430' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation_singleline.gff3.nt.gene.fa > MYZPE13164_O_EIv2.1_0184430.fa
#ReferenceCDS sequence:
grep -A 1 'MYZPE13164_O_EIv2.1_0184430.1' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation_singleline.gff3.nt.CDS.fa > MYZPE13164_O_EIv2.1_0184430.1.fa
#VCF SNP file, gene region, for 193 individuals:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/dedup_MYZPE13164_O_EIv2.1_0184430_snps.vcf

#Gene sequence multifasta for 193 individuals:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/hom_gene_fastas/hom_MYZPE13164_O_EIv2.1_0184430.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_gene/het_gene_fastas/het_MYZPE13164_O_EIv2.1_0184430.fa
#CDS sequence multifasta for 193 individuals:
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/hom_MYZPE13164_O_EIv2.1_0184430.1_CDS+.fa .
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas/het_MYZPE13164_O_EIv2.1_0184430.1_CDS+.fa .

#Get aa sequences:
source package /nbi/software/production/bin/emboss-6.5.7
grep -A 1 'MYZPE13164_O_EIv2.1_0184430.1' /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Myzus/persicae/O_v2/MYZPE13164_O_EIv2.1.annotation_singleline.gff3.aa.fa > white_1.aa.fa
transeq -sequence /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/hom_CDS_fastas/hom_MYZPE13164_O_EIv2.1_0184430.1_CDS+.fa -outseq hom_MYZPE13164_O_EIv2.1_0184430.1_CDS.aa.fa
transeq -sequence /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/snps_per_CDS/het_CDS_fastas/het_MYZPE13164_O_EIv2.1_0184430.1_CDS+.fa -outseq het_MYZPE13164_O_EIv2.1_0184430.1_CDS.aa.fa

echo scaffold_2 > temp.txt
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file ./temp.txt --input /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Myzus_persicae/O_v2/Myzus_persicae_O_v2.0.scaffolds.fa --output scaffold2.fa
```
```python
def extract_characters_from_fasta(fasta_file, start, end):
	sequence = ""
	with open(fasta_file, 'r') as f:
		header = None
		for line in f:
			if line.startswith('>'): 
				if sequence:  
					return sequence[start - 1:end]  
				header = line.rstrip()
				sequence = ""
			else:
				sequence += line.strip()
		if sequence:
			return sequence[start - 1:end]

# Example usage:
fasta_file = "scaffold2.fa"
start_pos = 63553919
end_pos = 63554147
extracted_sequence = extract_characters_from_fasta(fasta_file, start_pos, end_pos)
print(extracted_sequence)

```
scaffold_2      MYZPE13164_O_EIv2.1     gene    63534977        63556269        1.05e+03        +       .       ID=MYZPE13164_O_EIv2.1_0184430;Name=MYZPE13164_O_EIv2.1_0184430;biotype=pro$
scaffold_2      MYZPE13164_O_EIv2.1     mRNA    63534977        63556269        1.05e+03        +       .       ID=MYZPE13164_O_EIv2.1_0184430.1;Parent=MYZPE13164_O_EIv2.1_0184430;Name=MY$
scaffold_2      MYZPE13164_O_EIv2.1     five_prime_UTR  63534977        63535241        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.five_prime_UTR1;Parent=MYZPE13164_O_EIv2.1$
scaffold_2      MYZPE13164_O_EIv2.1     exon    63534977        63535340        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.exon1;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     CDS     63535242        63535340        .       +       0       ID=MYZPE13164_O_EIv2.1_0184430.1.CDS1;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     CDS     63547435        63547612        .       +       0       ID=MYZPE13164_O_EIv2.1_0184430.1.CDS2;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     exon    63547435        63547612        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.exon2;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     CDS     63548087        63548325        .       +       2       ID=MYZPE13164_O_EIv2.1_0184430.1.CDS3;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     exon    63548087        63548325        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.exon3;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     CDS     63549869        63549940        .       +       0       ID=MYZPE13164_O_EIv2.1_0184430.1.CDS4;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     exon    63549869        63549940        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.exon4;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     CDS     63550024        63550122        .       +       0       ID=MYZPE13164_O_EIv2.1_0184430.1.CDS5;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     exon    63550024        63550122        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.exon5;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     CDS     63550361        63550599        .       +       0       ID=MYZPE13164_O_EIv2.1_0184430.1.CDS6;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     exon    63550361        63550599        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.exon6;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     CDS     63552489        63552672        .       +       1       ID=MYZPE13164_O_EIv2.1_0184430.1.CDS7;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     exon    63552489        63552672        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.exon7;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     CDS     63553055        63553210        .       +       0       ID=MYZPE13164_O_EIv2.1_0184430.1.CDS8;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     exon    63553055        63553210        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.exon8;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     CDS     63553278        63553409        .       +       0       ID=MYZPE13164_O_EIv2.1_0184430.1.CDS9;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     exon    63553278        63553409        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.exon9;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     CDS     63553919        63554147        .       +       0       ID=MYZPE13164_O_EIv2.1_0184430.1.CDS10;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     exon    63553919        63554147        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.exon10;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     CDS     63554994        63555099        .       +       2       ID=MYZPE13164_O_EIv2.1_0184430.1.CDS11;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     exon    63554994        63555099        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.exon11;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     CDS     63555299        63555470        .       +       1       ID=MYZPE13164_O_EIv2.1_0184430.1.CDS12;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     exon    63555299        63555470        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.exon12;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     CDS     63555653        63555751        .       +       0       ID=MYZPE13164_O_EIv2.1_0184430.1.CDS13;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     exon    63555653        63556269        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.exon13;Parent=MYZPE13164_O_EIv2.1_0184430.1
scaffold_2      MYZPE13164_O_EIv2.1     three_prime_UTR 63555752        63556269        .       +       .       ID=MYZPE13164_O_EIv2.1_0184430.1.three_prime_UTR1;Parent=MYZPE13164_O_EIv2.$