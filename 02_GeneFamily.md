# Gene Family Analysis
Pipeline overview

1. OrthoDB (orthologer) mapping
2. OrthoDB database tidying
3. Build a tree for CAFE
3. CAFE

## 1. OrthoDB (orthologer) mapping

### Rename the fasta file headers

`awk '/^>/{print ">Tmax_REANN" ++i; next}{print}' Tmax_ReAnn.fasta > Tmax_ReAnn_renamed.fasta`

`mv Tmax_ReAnn_renamed.fasta Tmax_ReAnn.fasta`

### orthologer setup and mapping
```bash
# set up docker
# pull orthologer image
docker pull ezlabgva/orthologer:v3.0.0
# setup
mkdir odb
docker run -u $(id -u) -v $(pwd)/odb:/odbwork ezlabgva/orthologer:v3.0.0 setup_odb.sh map
docker run -u $(id -u) -v $(pwd)/odb:/odbwork ezlabgva/orthologer:v3.0.0 ./orthologer.sh
# That should print out orthologer.sh help.
 docker run -u $(id -u) -v $(pwd)/odb:/odbwork ezlabgva/orthologer:v3.0.0 orthomapper -c run -p Tmax -f Tmax_ReAnn.fasta -n 6447
```
### Results
Tmax.og.annotations: The first two columns *query* and *ODB_OG*
Combine this file with the file 3 in the next step
```bash
awk 'NR > 1 { print }'  Tmax.og.hits | cut -f1,2 > Tmax_OG2gene.tab
# Check if one gene is assigned to multiple orthogroup
cut -f2 Tmax_OG2gene.tab | uniq | wc -l #21883
wc -l Tmax_OG2gene.tab #21883
```

## 2. OrthoDB database tidying
We will need 3 files in the end:
1. speciesID    speciesName
2. ODB_OG id    OG annotation (descriptions, GO terms, etc.)
3. ODB_OG id    Gene_ID (distinguishable from each species)

### file 1: speciesID  speciesName
*Mollusca_speciesID2name.tab*

```bash
# get mollusca level ID
grep "Mollusca" odb11v0_levels.tab #6447
# get all mollusca species IDs
grep -w "6447" odb11v0_level2species.tab | cut -f2 > Mollusca_species_IDs
# get species ID to species name file
grep -Fwf Mollusca_species_IDs odb11v0_species.tab | cut -f2,3 > Mollusca_speciesID2name.tab
```



### file 2: ODB_OG id  OG annotation (descriptions, GO terms, etc.)
*Mollusca_OG2Ann.tab*

We can get the all the OG *$(pwd)/odb//data/og_descriptions* file *6447.tsv* for all the mollusca OG annotation.

```bash
cp /home/ruiqi/ruiqi_data/GiantClamGenome/OrthoDB_mollusca/odb/data/og_descriptions/6447.tsv /home/ruiqi/ruiqi_data/GiantClamGenome/OrthoDB/Mollusca_OG2Ann.tab
# count how many OGs
wc -l Mollusca_OG2Ann.tab #29695
# Count how many OGs Tmax are in
cut -f1 Tmax_OG2gene.tab | uniq | wc -l #11566
```


### file 3: ODB_OG id Gene_ID (distinguishable from each species)



Method: using species ids grep the OG2gene mapping
```bash
# use species ids grep the OG2gene mapping
nohup grep -Fwf Mollusca_species_IDs odb11v0_OG2genes.tab > Mollusca_OG2gene.tab0&
# Get mollusca specific OG IDs
grep "at6447" Mollusca_OG2gene.tab0 > Mollusca_OG2gene.tab
# check number
cut -f1 Mollusca_OG2gene.tab | uniq | wc -l #29695
# Concatenate Tmax data with other mollusca
# Mollusca geneID format is: 111167_1:000002 (speciesName:#)
# Giant clam format: Tmax_REANN7304
# Replace "_" with ":" in Giant clams first [to seperate species id, then do the count]
sed 's/_/:Tmax_/g' Tmax_OG2gene.tab > Tmax_OG2gene.tab0
# Concatenate the two files
cat Tmax_OG2gene.tab0 Mollusca_OG2gene.tab > Tmax_Mollusca_OG2gene.tab
# splitting by ":", then get OG_id  SpeciesID
awk -F '[:]' '{ print $1 }' Tmax_Mollusca_OG2gene.tab > Tmax_Mollusca_OG2speciesID.tab
# check number of OGs
cut -f1 Tmax_Mollusca_OG2speciesID.tab | sort -k 1,1 -u | wc -l # 29695
```
Then Process in R.
pwd:/Users/ruiqi/OneDrive - UCB-O365/GiantClamGenome/GeneFamily/GeneCount/OG2speciesID2count.R
output: Mollusca_GeneCount.tsv

## 3. Phylogeny
1. Get single copy OG family ID

Process in R

pwd:/Users/ruiqi/OneDrive - UCB-O365/GiantClamGenome/GeneFamily/GeneCount/OG2speciesID2count.R

I select single copy orthologs in at least 14/16 species.

Output: Mollusca_SingCopyOG_14_OG (2059 OGs)

2. Get the sequenceID in that OG, save to separate files
Output genelist in SingCopyOG_14_genelist

```bash
mkdir SingCopyOG_14_genelist
# grep each OG and save to gene Id to a file
while read line
 do
 grep -w $line Tmax_Mollusca_OG2gene.tab | cut -f2 > ./SingCopyOG_14_genelist/${line}_genelist.txt
 done < Mollusca_SingCopyOG_14_OG
# all mollusca Single copy OG sequences
grep -Fwf Mollusca_SingCopyOG_14_OG Tmax_Mollusca_OG2gene.tab | cut -f2 > Tmax_Mollusca_all_geneID
```

3. Retrieve sequences
```bash
# Environment cdbtools
conda activate cdbtools
# just keep the orthoDB gene id on header
awk '/^>/ {$0=$1} 1' odb11v0_all_og_fasta.tab > odb11v_cleanheader.fasta
# Build a library with all mollusca OG sequences
# since the database is huge (>4GB), will need to split them into 11 fasta files
fasta-splitter.pl --n-parts 11 odb11v_cleanheader.fasta
# index and retrieve sequences
for file in odb11v_cleanheader.part*
  do
   cdbfasta $file
 done
#get mollusca genes
cut -f2 ../Mollusca_OG2gene.tab > Mollusca_gene_ID
# check number of genes
wc -l Mollusca_gene_ID # 308488

# Retrieve mollusca specific genes
for file in odb11v_cleanheader.part*.fasta.cidx
do
  cat Mollusca_gene_ID | cdbyank $file > ${file}_Mollusca_gene_seq.fasta
done
cat *_Mollusca_gene_seq.fasta > Mollusca_gene_OrthoDB.fasta
# check num of sequences
grep ">" Mollusca_gene_OrthoDB.fasta | wc -l #308488
# Concatenate mollusca OG database + Tmax
# Chang header + Tmax:Tmax_REANNxxx
sed -i 's/>Tmax/>Tmax:Tmax/g' Tmax_ReAnn.fasta
# Concatenate Tmax and Mollusca Database
cat Mollusca_Gene_OrthoDB/Mollusca_gene_OrthoDB.fasta Tmax_ReAnn.fasta > Tmax_Mollusca.fasta #354957=308488+46469

### Retrieve sequences
cdbfasta Tmax_Mollusca.fasta
mkdir SingCopyOG_14_seqs
# Retrieve seqs from each OG and save to separate fasta files
for file in ./SingCopyOG_14_genelist/*_genelist.txt
do
  filename=`basename $file`
  name=${filename%_*}
  echo $name
  cat $file | cdbyank Tmax_Mollusca.fasta.cidx > ./SingCopyOG_14_seqs/${name}.fasta
done
```
Output sequences in OrthoDB/SingCopyOG_14_seqs

4. Alignment, Trim, Phylogeny Inference

conda enviroment: exon
pwd: /home/ruiqi/ruiqi_data/GiantClamGenome/OrthoDB/SingCopyOG_14_seqs
script: *Phylogeny_Tmax_Mollusca_SingCopyOG_14.sh*
```bash
# OG_list
cp ../Mollusca_SingCopyOG_14_OG ./
## script began here
### Alignment
# align each gene separately,
while read name;
do mafft ${name}.fasta > $name.output
  trimal -in $name.output -out ${name}_trimmed.output -automated1
done < Mollusca_SingCopyOG_14_OG
# Trimal also automatically only keeps the species name (before :)
# To concatenate the alignments: catfasta2phyml
# --concatenate : Concatenate files even when number of taxa differ among alignments. Missing data will be filled with all gap (-) sequences
./catfasta2phyml/catfasta2phyml.pl --concatenate --verbose --fasta *_trimmed.output > Tmax_Mollusca_SingCopyOG_14.fasta
echo "alignment done!"

mkdir Phylogeny
mv Tmax_Mollusca_SingCopyOG_14.fasta ./Phylogeny/
cd Phylogeny
# raxml
run1=all1
run2=all2
run=Tmax_Mollusca_SingCopyOG_14
raxmlHPC-PTHREADS -T 36 -p 12345 -N 20 -s Tmax_Mollusca_SingCopyOG_14.fasta -m PROTGAMMAWAG -n $run1
raxmlHPC-PTHREADS -T 36 -p 12345 -b 12345 -s Tmax_Mollusca_SingCopyOG_14.fasta -N 100 -m PROTGAMMAWAG -n $run2
raxmlHPC-PTHREADS -T 36 -f b -t RAxML_bestTree.${run1}  -z RAxML_bootstrap.${run2} -m PROTGAMMAWAG -n $run

echo "phylogeny done!"
```

Root the tree in Figtree and Export to newick format. Note to get rid of the "%" in the end.
*Tmax_Mollusca_SCOG14_bestTree.rooted*

```bash
# Make it ultrametric
conda activate OrthoFinder
#Find time of outgroup or root node from TimeTree.org: age-> of root node
# python2.7 OrthoFinder/tools/make_ultrametric.py -r <age> SpeciesTree_rooted.txt
python /home/ruiqi/miniconda3/envs/OrthoFinder/bin/make_ultrametric.py -r 543 Tmax_Mollusca_SCOG14_bestTree.rooted

### species ID to name replacement
sed 's/6454_0/Haliotis_rufescens/g; s/6500_0/Aplysia_californica/g; s/6526_0/Biomphalaria_glabrata/g; s/6565_0/Crassostrea_virginica/g; s/6573_0/Mizuhopecten_yessoensis/g; s/6579_0/Pecten_maximus/g; s/6596_0/Mercenaria_mercenaria/g; s/29159_0/Crassostrea_gigas/g; s/36100_0/Haliotis_rubra/g; s/37623_0/Ostrea_edulis/g; s/37653_0/Octopus_bimaculoides/g; s/225164_0/Lottia_gigantea/g; s/400727_0/Pomacea_canaliculata/g; s/1735272_0/Gigantopelta_aegis/g; s/2607531_0/Octopus sinensis/g; s/Tmax/Tridacna_maxima/g' Tmax_Mollusca_SCOG14_bestTree.rooted.ultrametric.tre
```










#### CAFE

[CAFE Tutorial](https://evomics.org/wp-content/uploads/2016/06/cafe_tutorial-1.pdf)
[CAFE 5 Turtorial](https://github.com/hahnlab/CAFE5)
[OrthoFinder Species Tree to tree for CAFE](https://github.com/davidemms/OrthoFinder/blob/master/tools/make_ultrametric.py)
[OrthoFinder CAFE pipeline](https://github.com/harish0201/Analyses_Pipelines/blob/main/7.CAFE.sh)
[OrthoFinder-CAFE pipeline in Chinese](https://www.yuque.com/shawn-yepfo/oxua15/ftrlin)
