# Gene Family Analysis
Pipeline overview

1. OrthoDB (orthologer) mapping
2. OrthoDB database tidying
3. Build a tree for CAFE
4. CAFE
5. Enrichment

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










#### 4. CAFE

Install CAFE from conda in the enviroment cafe

pwd:*/home/ruiqi/ruiqi_data/GiantClamGenome/OrthoDB/CAFE_Tmax_mollusca*

Run CAFE 5
```bash
# -t parameter specifies a file containing the tree that CAFE uses
# -i parameter specifies a list of gene families
# -c Number of processing cores to use
# -k 3 To incorporate among family rate variation with both lambda and alpha estimated and three discrete gamma rate categories
# FAILED: Failed to initialize any reasonable values
#cafe5 -i Mollusca_GeneCount.tsv -t Tmax_Mollusca_SCOG14_bestTree.rooted.ultrametric.tre -k 3 -c 6 -l 0.0001 -o k3
# -l 0.0001 --fixed_lambda, -l :Value (between 0 and 1) for a single user provided lambda value, otherwise lambda is estimated
cafe5 -i Mollusca_GeneCount.tsv -t Tmax_Mollusca_SCOG14_bestTree.rooted.ultrametric.tre -k 3 -c 6 -l 0.0001 -o k3_l1e-4


#-p Use a Poisson distribution for the root frequency distribution. If no -p flag is given, a uniform distribution will be used
# the best (lowest) likelihood score???
nohup cafe5 -i Mollusca_GeneCount.tsv -t Tmax_Mollusca_SCOG14_bestTree.rooted.ultrametric.tre -k 1 -p -c 6 -o k1p&
nohup cafe5 -i Mollusca_GeneCount.tsv -t Tmax_Mollusca_SCOG14_bestTree.rooted.ultrametric.tre -k 3 -p -c 6 -o k3p&
nohup cafe5 -i Mollusca_GeneCount.tsv -t Tmax_Mollusca_SCOG14_bestTree.rooted.ultrametric.tre -k 2 -p -c 2 -o k2p&
nohup cafe5 -i Mollusca_GeneCount.tsv -t Tmax_Mollusca_SCOG14_bestTree.rooted.ultrametric.tre -k 4 -p -c 2 -o k4p&
nohup cafe5 -i Mollusca_GeneCount.tsv -t Tmax_Mollusca_SCOG14_bestTree.rooted.ultrametric.tre -k 5 -p -c 2 -o k5p& # failed
nohup cafe5 -i Mollusca_GeneCount.tsv -t Tmax_Mollusca_SCOG14_bestTree.rooted.ultrametric.tre -k 6 -p -c 2 -o k6p& # failed

```

Results Mining
```bash
# write the just significant families to a file
grep "y" Gamma_family_results.txt > Significant_families.txt # 6794
# filter gene families to a lower p-value (0.01 in the example)
awk '$2 < .01 {print $0}' Gamma_family_results.txt > Sig_at_p.01.txt
wc -l Sig_at_p.01.txt #4959
# From
#Taxon  Increased Decreased
#Tmax<17>        1969    1925
# Only viewing the trees significant changes
echo $'#nexus\nbegin trees;'>Significant_trees.tre
grep "*" Gamma_asr.tre >>Significant_trees.tre
echo "end;">>Significant_trees.tre # 6794 trees

## Only viewing the trees significant changes in Tmax
echo $'#nexus\nbegin trees;'>Tmax_Significant_trees.tre
grep -F "<17>*" Gamma_asr.tre >>Tmax_Significant_trees.tre
echo "end;">>Tmax_Significant_trees.tre # 2994 trees
# Get the gene family IDs that have significant changes in Tmax
grep -o "\w*at6447\w*" Tmax_Significant_trees.tre > Tmax_Significant_trees_IDs

# Get the gene family expaned in Tmax
cat Gamma_change.tab |cut -f1,19|grep "+[1-9]" >Tmax.expanded # 1969
# Get the gene family extracted in Tmax
cat Gamma_change.tab |cut -f1,19|grep "-" >Tmax.contracted # 1925

############ ############ This method has a caveat ############ ############
############ If the gene gamily is significantly changed but not in Tmax, it will still count ############
# Get the list of significant changed gene family IDs  
#grep "at6447" Significant_families.txt | cut -f1 > Sig_at_p.05_IDs #6794
#  Get the gene family SIGNIFICANTLY expaned in Tmax
#grep -Fwf Sig_at_p.05_IDs Tmax.expanded | cut -f1 > Tmax.expanded.Sig_at_p.05 #1418
#  Get the gene family SIGNIFICANTLY contracted in Tmax
#grep -Fwf Sig_at_p.05_IDs Tmax.contracted | cut -f1 > Tmax.contracted.Sig_at_p.05  #1614

############ Significant Branch on the Tree ########################
#  Get the gene family SIGNIFICANTLY expaned in Tmax
grep -Fwf Tmax_Significant_trees_IDs Tmax.expanded | cut -f1 > Tmax.expanded.Sig.treebranch #1417
#  Get the gene family SIGNIFICANTLY contracted in Tmax
grep -Fwf Tmax_Significant_trees_IDs Tmax.contracted | cut -f1 > Tmax.contracted.Sig.treebranch  #1577

############# ############ Valid the method, USE significant branches ############ ############
#grep -vFwf Tmax.expanded.Sig.treebranch Tmax.expanded.Sig_at_p.05 #3453at6447
# Visualize the tree on http://etetoolkit.org/treeview/
#grep -w "3453at6447" Gamma_asr.tre # obviously not expanded in Tmax

```


Setting λ to a previously estimated value to deal with families with large numbers of gene copies

check Model Gamma Final Likelihood (-lnL) in Gamma_results.txt from different runs (k2p，k3p，k5p，k6p), select the largest??? as the best
```bash
#把每个节点收缩扩张的基因数量画在树上
# R package ggtree
https://yanzhongsino.github.io/2022/01/24/bioinfo_evolutionary.tree_ggtree/
```

## 5. Functional Enrichment
Gene Family Functional Description Metadata:
*Mollusca_GeneFamily_metadata.tsv* obtained from */Users/ruiqi/OneDrive - UCB-O365/GiantClamGenome/GeneFamily/GeneCount/OG2speciesID2count.R*

Get the contracted gene family data:

`head -1 ../Mollusca_GeneFamily_metadata.tsv > Tmax.contracted.genefamily.annotation.tsv`
`grep -Fwf Tmax.contracted.Sig.treebranch ../Mollusca_GeneFamily_metadata.tsv >> Tmax.contracted.genefamily.annotation.tsv`

Get the expanded gene family data:
`head -1 ../Mollusca_GeneFamily_metadata.tsv > Tmax.expanded.genefamily.annotation.tsv`
`grep -Fwf Tmax.expanded.Sig.treebranch ../Mollusca_GeneFamily_metadata.tsv >> Tmax.expanded.genefamily.annotation.tsv`

### Annotation of Tmax with emapper
```bash
conda activate emapper
outputdir=/home/ruiqi/ruiqi_data/GiantClamGenome/TmaximaGenomeReannotation/Tmax_ReAnn_emapper
input=/home/ruiqi/ruiqi_data/GiantClamGenome/TmaximaGenomeReannotation/Tmax_ReAnn.fasta
db=/home/ruiqi/ruiqi_data/emapper/databases
nohup emapper.py -m diamond --output_dir $outputdir -o Tmax_ReAnn --cpu 40 -i $input --data_dir $db&
cat Tmax_ReAnn.emapper.annotations | sed '/^#/d' > Tmax_ReAnn.clean.annotations
```
### Retrieve gene sequences for all gene families, then construct gene trees

Retrieve sequences

```bash
# get the mollusca OG list  
cat Mollusca_OG2Ann.tab | cut -f1 > Mollusca_OG_IDs
mkdir Mollusca_Tmax_genelist
# grep each OG and save to gene Id to a file
while read line
 do
 grep -w $line Tmax_Mollusca_OG2gene.tab | cut -f2 > ./Mollusca_Tmax_genelist/${line}_genelist.txt
 done < Mollusca_OG_IDs

### Retrieve sequences (Use Tmax_Mollusca.fasta from the step 3)
# Environment cdbtools
conda activate cdbtools
mkdir Mollusca_Tmax_seqs
# Retrieve seqs from each OG and save to separate fasta files
for file in ./Mollusca_Tmax_genelist/*_genelist.txt
do
  filename=`basename $file`
  name=${filename%_*}
  echo $name
  cat $file | cdbyank Tmax_Mollusca.fasta.cidx > ./Mollusca_Tmax_seqs/${name}.fasta
done
```

Construct gene trees
conda enviroment: exon
pwd: /home/ruiqi/ruiqi_data/GiantClamGenome/OrthoDB/Mollusca_Tmax_seqs
script: *Phylogeny_Tmax_genetree.sh*


```bash
#!/bin/bash
### Alignment
# align each gene separately,
while read name;
do mafft ${name}.fasta > $name.output
  sed -i 's/:/_/g' $name.output
  trimal -in $name.output -out ${name}_trimmed.output -automated1
done < Mollusca_OG_IDs

mkdir Phylogeny_GeneTree_IntermediateFiles
mkdir GeneTrees
cd Phylogeny_GeneTree_IntermediateFiles
# Phylogeny
while read name;
do
  run1=${name}_run1
  run2=${name}_run2
  run=${name}
  raxmlHPC-PTHREADS -T 40 -p 12345 -N 20 -s ../${name}_trimmed.output -m PROTGAMMAWAG -n $run1
  raxmlHPC-PTHREADS -T 40 -p 12345 -b 12345 -s ../${name}_trimmed.output -N 100 -m PROTGAMMAWAG -n $run2
  raxmlHPC-PTHREADS -T 40 -f b -t RAxML_bestTree.${run1}  -z RAxML_bootstrap.${run2} -m PROTGAMMAWAG -n $run
  cp RAxML_bipartitions.${name} ../GeneTrees/${name}_GeneTree.tre
done < ../Mollusca_OG_IDs

mkdir TrimmedAlignment_GeneTree
mv *_trimmed.output TrimmedAlignment_GeneTree/
mkdir Alignment_GeneTree
mv *.output Alignment_GeneTree/
mkdir Mollusca_Tmax_fasta
mv *.fasta Mollusca_Tmax_fasta/
```



### GO/KEGG enrichment
1. Get Tmax gene IDs that are clustered to mollusca gene families  (Universe)
2. Get Tmax gene IDs that below to expanded/contracted gene families


[CAFE Tutorial](https://evomics.org/wp-content/uploads/2016/06/cafe_tutorial-1.pdf)

[CAFE 5 Turtorial](https://github.com/hahnlab/CAFE5)

[OrthoFinder Species Tree to tree for CAFE](https://github.com/davidemms/OrthoFinder/blob/master/tools/make_ultrametric.py)

[OrthoFinder CAFE pipeline](https://github.com/harish0201/Analyses_Pipelines/blob/main/7.CAFE.sh)

[OrthoFinder-CAFE pipeline in Chinese](https://www.yuque.com/shawn-yepfo/oxua15/ftrlin)

[CAFE in Chinese](https://yanzhongsino.github.io/2021/10/29/bioinfo_gene.family_CAFE5/)

Name Code Replacement
### species ID to name replacement
```bash
sed -i 's/6454_0/Haliotis_rufescens/g; s/6500_0/Aplysia_californica/g; s/6526_0/Biomphalaria_glabrata/g; s/6565_0/Crassostrea_virginica/g; s/6573_0/Mizuhopecten_yessoensis/g; s/6579_0/Pecten_maximus/g; s/6596_0/Mercenaria_mercenaria/g; s/29159_0/Crassostrea_gigas/g; s/36100_0/Haliotis_rubra/g; s/37623_0/Ostrea_edulis/g; s/37653_0/Octopus_bimaculoides/g; s/225164_0/Lottia_gigantea/g; s/400727_0/Pomacea_canaliculata/g; s/1735272_0/Gigantopelta_aegis/g; s/2607531_0/Octopus sinensis/g; s/Tmax/Tridacna_maxima/g' file
```
