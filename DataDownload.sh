###
#### annotate giant clam genome
nohup python2.7 emapper.py -m diamond --output_dir /home/ruiqi/ruiqi_data/GiantClamGenome/ -o Tma_Annotation --cpu 28 -i ./Tma.faa &


### Space clean up
### Move to hard drive
/home/ruiqi/ruiqi_data/Transcriptome
/home/ruiqi/ruiqi_data/
/home/ruiqi/ruiqi_data/PanamaData/Trimmed_Reads

#### cd /home/ruiqi/ruiqi_data/GiantClamGenome/OtherGenomes
### Download other genomes (protein sequence)
# scallop: Mizuhopecten yessoensis
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/113/885/GCF_002113885.1_ASM211388v2/GCF_002113885.1_ASM211388v2_protein.faa.gz
mv GCF_002113885.1_ASM211388v2_protein.faa.gz Mye.faa.gz
# pacific oyster: Crassostrea gigas
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/806/645/GCF_902806645.1_cgigas_uk_roslin_v1/GCF_902806645.1_cgigas_uk_roslin_v1_protein.faa.gz
mv GCF_902806645.1_cgigas_uk_roslin_v1_protein.faa.gz Cgi.faa.gz
# Only CDS available, translate them with TransDecoder

# Sun, J., Zhang, Y., Xu, T. et al. Adaptation to deep-sea chemosynthetic environments as revealed by mussel genomes. Nat Ecol Evol 1, 0121 (2017)
# https://datadryad.org/stash/dataset/doi:10.5061%2Fdryad.h9942
# Gpl, Gigantidas platifrons (mussel);
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/080/005/GCA_002080005.1_Bpl_v1.0/GCA_002080005.1_Bpl_v1.0_genomic.fna.gz
mv GCA_002080005.1_Bpl_v1.0_genomic.fna.gz Gpl.fna.gz
# Mph, Modiolus philippinarum (shallow imparidentian); Philippine horse mussel
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/080/025/GCA_002080025.1_Mph_v1.0/GCA_002080025.1_Mph_v1.0_genomic.fna.gz
mv GCA_002080025.1_Mph_v1.0_genomic.fna.gz Mph.fna.gz


# Rph, Ruditapes philippinarum (shallow imparidentian) :
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/026/015/GCA_009026015.1_ASM902601v1/GCA_009026015.1_ASM902601v1_genomic.fna.gz
mv GCA_009026015.1_ASM902601v1_genomic.fna.gz Rph.fna.gz

# Sco, Sinonovacula constricta (shallow imparidentian); Chinese razor clam
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.qv9s4mwc6
#Ran, Zhaoshou, et al. "Chromosome‐level genome assembly of the razor clam Sinonovacula constricta (Lamarck, 1818)." Molecular ecology resources 19.6 (2019): 1647-1658.
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/844/125/GCA_007844125.1_ASM784412v1/GCA_007844125.1_ASM784412v1_genomic.fna.gz
mv GCA_007844125.1_ASM784412v1_genomic.fna.gz Sco.fna.gz

####### probably should exclude this because it is symbiotic
# Ama, Archivesica marissinica (pliocardiine).deep-sea symbiotic clam
# Jack Chi-Ho Ip, Host–Endosymbiont Genome Integration in a Deep-Sea Chemosymbiotic Clam, Molecular Biology and Evolution, Volume 38, Issue 2, February 2021, Pages 502–518, https://doi.org/10.1093/molbev/msaa241
# https://figshare.com/articles/dataset/Host-Endosymbiont_Genome_Integration_in_a_Deep-Sea_Chemosymbiotic_Clam/12198987
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/843/695/GCA_014843695.1_ASM1484369v1/GCA_014843695.1_ASM1484369v1_genomic.fna.gz
mv GCA_014843695.1_ASM1484369v1_genomic.fna.gz Ama.fna.gz

### Download other genomes from MolluscaDB
# Bathymodiolus platifrons: deep sea mussel



#outgroup: capitella teleta, annelida
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/328/365/GCA_000328365.1_Capca1/GCA_000328365.1_Capca1_protein.faa.gz
mv GCA_000328365.1_Capca1_protein.faa.gz Cte.faa.gz

# outgroup: lingula anatina, Brachiopoda
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/039/355/GCF_001039355.2_LinAna2.0/GCF_001039355.2_LinAna2.0_protein.faa.gz
mv GCF_001039355.2_LinAna2.0_protein.faa.gz Lbr.faa.gz

# outgroup: Lottia gigantea (Owl limpet)
# Outgroup: Pomacea_canaliculata (Apple snail)



gunzip *.gz

#### Don't do these in the future, find protein sequence from original paper, or use protein prediction tools like maker
# CDS to protein
# conda create -n transdecoder
# conda activate transdecoder
# conda install -c bioconda transdecoder
nohup TransDecoder.LongOrfs -t Ama.fna &
nohup TransDecoder.LongOrfs -t Mph.fna &
nohup TransDecoder.LongOrfs -t Rph.fna &
nohup TransDecoder.LongOrfs -t Sco.fna &
nohup TransDecoder.LongOrfs -t Gpl.fna &
# Phylogenetic analyses
## Orthogroups with OrthoFinder
conda create -n OrthoFinder
conda activate OrthoFinder
conda install -c bioconda orthofinder

# copy the giant clam genome to the folder, then run
nohup orthofinder -S diamond_ultra_sens -f subset_test &

## Will need to add species name to each sequence before run the orthofinder
awk '/>/{sub(">","&"FILENAME"_");sub(/\.faa/,x)}1' test.faa > test2.faa
## Run OrthoFinder on the fullset
nohup orthofinder -S diamond_ultra_sens -f fullset &

# Continue testing
# select genes that are found in all 4 species\

# let's use mafft in case of forward and reverse sequence directions
# https://onestopdataanalysis.com/multiple-sequence-alignment-msa-reverse-complement/
#conda create -n muscle
#conda activate muscle
#conda install -c bioconda muscle

conda create -n mafft
conda activate mafft
conda install -c bioconda mafft
conda install -c bioconda trimal

# use mafft to align each orthogroup, then use trimAL to trim them

for i in *.fa
do
    mafft --adjustdirection $i > ${i%.fa}.aligned.fasta
    trimal -in ${i%.fa}.aligned.fasta -out ${i%.fa}.trimmed.fasta -automated1
done
### change the name in the alignment file, before: >Spe_xxxGene; After: >Spe
for i in *.trimmed.fasta
do
  cut -f 1 -d "_" $i > ${i%.trimmed.fasta}_clean.fasta
done
### To concatenate the alignments: catfasta2phyml
#From IQ-Tree version 2 you can input a directory of alignment files. IQ-TREE 2 will load and concatenate all alignments within the directory, eliminating the need for users to manually perform this step
# --concatenate : Concatenate files even when number of taxa differ among alignments. Missing data will be filled with all gap (-) sequences.
catfasta2phyml/catfasta2phyml.pl --concatenate --verbose --fasta *_clean.fasta > Subset.fasta

#Phylogeny with RAxML
conda create -n raxml
conda activate raxml
conda install -c bioconda raxml


# -m PROTGAMMAAUTO tells the program that these are protein sequences, and tells the program how to model protein evolution

run1=all1
run2=all2
run=all_1
raxmlHPC-PTHREADS -T 36 -p 12345 -N 20 -s Subset.fasta -m PROTGAMMAAUTO -n $run1
raxmlHPC-PTHREADS -T 36 -p 12345 -b 12345 -s Subset.fasta -N 100 -m PROTGAMMAAUTO -n $run2
raxmlHPC-PTHREADS -T 36 -f b -t RAxML_bestTree.${run1}  -z RAxML_bootstrap.${run2} -m PROTGAMMAAUTO -n $run
ls


# IQ Tree
conda create -n iqtree
conda activate iqtree
conda install -c bioconda iqtree

# Run IQ Tree
# -m is the option to specify the model name to use during the analysis. The special MFP key word stands for ModelFinder Plus, which tells IQ-TREE to perform ModelFinder and the remaining analysis using the selected model
# example.phy.model: log-likelihoods for all models tested. It serves as a checkpoint file to recover an interrupted model selection
# -T allows specifying the number of CPU cores to use
# -s input sequences
# -B specifies the number of bootstrap replicates where 1000 is the minimum number recommended
nohup iqtree -s Subset.fasta -m MFP -T 34 -B 1000 &

## prepare input files for CAFE
# /home/ruiqi/ruiqi_data/GiantClamGenome/OtherGenomes/fullset/OrthoFinder/Results_Jul19/CAFE
## Tutorial: https://github.com/harish0201/Analyses_Pipelines/blob/main/7.CAFE.sh
conda activate OrthoFinder

# Convert SpeciesTree_rooted.txt to ultrametric
# Find time of outgroup or root node from TimeTree.org: age-> of root node
# Capitella teleta vs. Mizuhopecten yessoensis Median Time: 626 MYA CI: (545.0 - 681.5 MYA) Adjusted Time: 574 MYA
/home/ruiqi/miniconda3/pkgs/orthofinder-2.5.1-0/bin/make_ultrametric.py -r 574 SpeciesTree_rooted.txt

# Add a column
awk -F'\t' '{print "(null)\t"$0}' Orthogroups.GeneCount.tsv > tmp.tsv
# Change the header (null) to Desc and save
vim tmp.tsv
#remove the total column from above, without needed to figure out column numbers.
awk -F'\t' '{$NF=""; print $0}' tmp.tsv | rev | sed 's/^\s*//g' | rev | tr ' ' '\t' > mod.tsv

# Get the python_scripts
git clone https://github.com/hahnlab/cafe_tutorial.git

#run size filter finally
#filter the Orthogroups.GeneCount.tsv file to remove OG that have more than 100 proteins in a particular species:
python2.7 ./cafe_tutorial/python_scripts/cafetutorial_clade_and_size_filter.py -i mod.tsv -s -o cafe.input.tsv

# Prepare lamda tree structure: will look like this:
#(6,(5,(4,(3,((2,2)2,((1,1)1,((,),1)1)2)3)4)5)6)
cafe5 -i cafe.input.tsv -t SpeciesTree_rooted.txt.ultrametric.tre
# run cafe: inputs the tsv output from above; and using the rounded off tree.nwk

#
conda create -n cafe
conda activate cafe
conda install -c bioconda cafe

## results
# https://github.com/hahnlab/CAFE5/blob/master/docs/tutorial/tutorial.md


#### all Genomes
## Will need to add species name to each sequence before run the orthofinder
for i in *.faa
do
  awk '/>/{sub(">","&"FILENAME"_");sub(/\.faa/,x)}1' $i > ${i%.faa}.pep
done

## run orthofinder
nohup orthofinder -S diamond_ultra_sens -f fullset &

###Download the gene count files, select genes that are presented in at least 8/11 species in R
scp lilab:/home/ruiqi/ruiqi_data/GiantClamGenome/OtherGenomes/fullset/OrthoFinder/Results_Jul19/Orthogroups/Orthogroups.GeneCount.tsv

## in Rstudio
# /Users/ruiqi/Downloads/GiantClamGenome_Orthofinder/80PercentSingleCopyOrtholog.R
# Output: a list of orthogroups that cover at least 8/11 species

# cut the first column (Orthorgroup names)
cut -f1 FilteredSCO.txt > SCO_names.txt
# Add .fa to each name, $ stand for the end of each string
sed -i 's/$/.fa/' SCO_names.txt
# Upload the file
scp SCO_names.txt lilab:/home/ruiqi/ruiqi_data/GiantClamGenome/OtherGenomes/fullset/OrthoFinder/Results_Jul19/Orthogroup_Sequences
# Cope the filtered SCO to the directory
mkdir FilteredSCO
cp $(cat SCO_names.txt) ../FilteredSCO/

# let's use mafft in case of forward and reverse sequence directions
## align, trim
conda activate mafft
for i in *.fa
do
    mafft --adjustdirection $i > ${i%.fa}.aligned.fasta
    trimal -in ${i%.fa}.aligned.fasta -out ${i%.fa}.trimmed.fasta -automated1
done
### change the name in the alignment file, before: >Spe_xxxGene; After: >Spe
for i in *.trimmed.fasta
do
  cut -f 1 -d "_" $i > ${i%.trimmed.fasta}_clean.fasta
done

### To concatenate the alignments: catfasta2phyml
#From IQ-Tree version 2 you can input a directory of alignment files. IQ-TREE 2 will load and concatenate all alignments within the directory, eliminating the need for users to manually perform this step
# --concatenate : Concatenate files even when number of taxa differ among alignments. Missing data will be filled with all gap (-) sequences
git clone https://github.com/nylander/catfasta2phyml.git
catfasta2phyml/catfasta2phyml.pl --concatenate --verbose --fasta *_clean.fasta > fullset.fasta

# IQ Tree
conda activate iqtree

# Run IQ Tree
# -m is the option to specify the model name to use during the analysis. The special MFP key word stands for ModelFinder Plus, which tells IQ-TREE to perform ModelFinder and the remaining analysis using the selected model
# example.phy.model: log-likelihoods for all models tested. It serves as a checkpoint file to recover an interrupted model selection
# -T allows specifying the number of CPU cores to use
# -s input sequences
# -B specifies the number of bootstrap replicates where 1000 is the minimum number recommended
# /home/ruiqi/ruiqi_data/GiantClamGenome/OtherGenomes/fullset/iqtree_fullset
nohup iqtree -s fullset.fasta -m MFP -T 34 -B 1000 &

# Visualise the tree XXX.fa.treefile




# Gene Family Analyses
## prepare input files for CAFE
# pwd: /home/ruiqi/ruiqi_data/GiantClamGenome/OtherGenomes/fullset/OrthoFinder/Results_Jul19/CAFE
## Tutorial: https://github.com/harish0201/Analyses_Pipelines/blob/main/7.CAFE.sh
conda activate OrthoFinder

# Convert SpeciesTree_rooted.txt to ultrametric
# Find time of outgroup or root node from TimeTree.org: age-> of root node
# Capitella teleta vs. Modiolus philippinarum Median Time: 626 MYA CI: (545.0 - 681.5 MYA) Adjusted Time: 574 MYA
/home/ruiqi/miniconda3/pkgs/orthofinder-2.5.1-0/bin/make_ultrametric.py -r 574 SpeciesTree_rooted.txt

# Add a column
awk -F'\t' '{print "(null)\t"$0}' Orthogroups.GeneCount.tsv > tmp.tsv
# Change the header (null) to Desc and save
vim tmp.tsv
#remove the total column from above, without needed to figure out column numbers.
awk -F'\t' '{$NF=""; print $0}' tmp.tsv | rev | sed 's/^\s*//g' | rev | tr ' ' '\t' > mod.tsv

# Get the python_scripts
git clone https://github.com/hahnlab/cafe_tutorial.git

#run size filter finally
#filter the Orthogroups.GeneCount.tsv file to remove OG that have more than 100 proteins in a particular species:
python2.7 ./cafe_tutorial/python_scripts/cafetutorial_clade_and_size_filter.py -i mod.tsv -s -o cafe.input.tsv

# Prepare lamda tree structure: will look like this:
#(6,(5,(4,(3,((2,2)2,((1,1)1,((,),1)1)2)3)4)5)6)
conda activate cafe
cafe5 -i cafe.input.tsv -t SpeciesTree_rooted.txt.ultrametric.tre
# run cafe: inputs the tsv output from above; and using the rounded off tree.nwk

## results
# https://github.com/hahnlab/CAFE5/blob/master/docs/tutorial/tutorial.md


###
### Is PAML too old to install with Conda?
conda create -n paml
conda activate paml
conda install -c bioconda paml






## annotation with eggnog mapper
 ./emapper.py  -m diamond --output_dir /home/ruiqi/ruiqi_data/GiantClamGenome/ -o Tma_Annotation --cpu 28 -i ./Tma.faa

 ###subset gene function and seuqnece number
 cut -f1,22  Tma_Annotation.emapper.annotations > Tma_Annotation.subset

 ### Annotate other genomes
 for genome in ../../GiantClamGenome/OtherGenomes/fullset/*.pep
 do
   echo $i
   nohup ./emapper.py  -m diamond --output_dir /home/ruiqi/ruiqi_data/GiantClamGenome/eggnog_annotation/ -o ${genome%.pep}_Annotation --cpu 28 -i $genome
 done

 ### Get gene functions and gene names
 ### delete lines starting with #
 for i in  *_Annotation.emapper.annotations
 do
  echo $i
  sed '/^#/d' $i | cut -f1,22 | tr '\t' ',' > ${i%_Annotation.emapper.annotations}_Annotation.subset.csv
 done

 # concatenate csv files
 cat *_Annotation.subset.csv > othergenomes_Annotation.subset.csv


 ### annotation results to orthogroups
 ## pwd /home/ruiqi/ruiqi_data/GiantClamGenome/eggnog_annotation
cp ../OtherGenomes/fullset/OrthoFinder/Results_Jul19/Orthogroups/Orthogroups.tsv ./
## pwd /Users/ruiqi/Downloads/GiantClamGenome_Orthofinder
scp lilab:/home/ruiqi/ruiqi_data/GiantClamGenome/eggnog_annotation/Orthogroups.tsv ./
scp lilab:/home/ruiqi/ruiqi_data/GiantClamGenome/eggnog_annotation/Tma_Annotation.subset ./
scp lilab:/home/ruiqi/ruiqi_data/GiantClamGenome/eggnog_annotation/OtherGenomes_annotation/othergenomes_Annotation.subset.csv ./
### read.table does not read the rows with "()", convert it to .csv
tr '\t' ',' < Tma_Annotation.subset >  Tma_Annotation_subset.csv
 #### Rstudio
 # /Users/ruiqi/Downloads/GiantClamGenome_Orthofinder/Orthogroup_Annotation.R
 # Output: orthogroup-gene fucntion- counts from each species
 ###




 ## Venn Diagram
#pwd /home/ruiqi/ruiqi_data/GiantClamGenome/OtherGenomes/fullset/VennDiagram
conda activate OrthoFinder
nohup orthofinder -S diamond_ultra_sens -f VennDiagram &

### /Users/ruiqi/Documents/GiantClamGenome_Orthofinder/OrthoFinder2Venn
scp lilab:/home/ruiqi/ruiqi_data/GiantClamGenome/OtherGenomes/fullset/VennDiagram/OrthoFinder/Results_Jul21/Orthogroups/Orthogroups.tsv ./
sp_list=("Ama" "Cgi" "Mph" "Pca" "Tma")
for sp in $sp_list
do
  grep $sp Orthogroups.tsv|sed '1d' | cut -f1 > ${sp}_ortho.tsv
done

### Process venn diagram in Rstudio
#pwd /Users/ruiqi/Documents/GiantClamGenome_Orthofinder/OrthoFinder2Venn
#script : OrthoFinder2Venn.R
# output: Venn diagram
# Add percentage manually


### bar graph for gene family expansion and extraction
scp lilab:/home/ruiqi/ruiqi_data/GiantClamGenome/OtherGenomes/fullset/OrthoFinder/Results_Jul19/CAFE/results/Base_clade_results.txt ./
# /Users/ruiqi/Documents/GiantClamGenome_Orthofinder/bargraph_genefamily
### Script:bargraph_genefamily.R
# Rstudio
# Output: bargraph
