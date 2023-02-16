### Data Preparation

Download all the bivalve genome from MolluscaDB then upload to the serve, first we are going to select high quality genome based on the BUSCO scores, symbiotic status, and representation of species

Bivalvia: key species: Archivesica marissinica (Chemosymbiotic clam)

pwd:/home/ruiqi/ruiqi_data/GiantClamGenome/Bivalvia
conda activate BUSCO
```bash
#! usr/bin/bash
mkdir ../BUSCO_results
for file in *.fa
do
  name=$(echo ${file} | awk -F'[_.]' '{print $1"_"$2}')
  echo $name
  busco -m protein -i $file -o ../BUSCO_results/${name}_mollusca_BUSCO  -l mollusca_odb10 -c 44
  busco -m protein -i $file -o ../BUSCO_results/${name}_metazoa_BUSCO  -l metazoa_odb10 -c 44
done
```

Other molluscas:
Key species: Elysia_chlorotica (Symbiotic nudibranch); Chrysomallon squamiferum (The Scaly-foot Snail)
For the symbiotic squid, extract CDS from gff file: Euprymna scolopes (bobtail squid)



### get the BUSCO summary figure
```bash
mkdir BUSCO_summaries
cp */short_summary*txt BUSCO_summaries/
cd BUSCO_summaries

db=metazoa_odb10
db=mollusca_odb10

taxon=Other
taxon=Bivalvia

for file in short_summary.specific.*.txt
do
  echo $file
  name=$(echo ${file} | awk -F'[_.]' '{print $6"_"$7}')
  echo $name
  mv $file short_summary.specific.${db}.${name}.txt
done
cd ..
python3 /home/ruiqi/miniconda3/envs/BUSCO/bin/generate_plot.py -wd BUSCO_summaries
mv BUSCO_summaries/busco_figure.png BUSCO_summaries/busco_${db}_${taxon}.png


# move all the figures to one file
mv *_BUSCO_results/*/BUSCO_summaries/*png BUSCO_figures/
```

### List of Bivalvia
Archivesica marissinica (Chemosymbiotic clam)
Amusium japonicum (scallop)
Crassostrea virginica (oyster)
Scapharca_broughtonii (blood clam)
Sinonovacula constricta (razor clam, relatively low quality but closed related)

#### Add other mollusca
Elysia_chlorotica (Symbiotic nudibranch);
Chrysomallon squamiferum (The Scaly-foot Snail);


##### Rename them to G(Genus)sp(species).pep
(OrthoFinder only takes following extensions: fas, fa, faa, fasta, pep)
##### Rename sequence name
```bash
## Will need to add species name to each sequence before run the orthofinder
## Giant clam genome header cleaning
awk -F" " '/^>/ { print $1; next } 1' Tma.fna >Tma2.fna
#
for file in *.fna
do
  awk '/>/{sub(">","&"FILENAME"_");sub(/\.faa/,x)}1' $file > ${file}_1
  mv ${file}_1 $file
done
## Check
head *.fna | grep ">"
#### Orthofinder
## Need to edit the script on the new server: https://github.com/davidemms/OrthoFinder/issues/760
###path: miniconda3/envs/OrthoFinder/bin/scripts_of/newick.py
### env: Orthofinder
#! usr/bin/bash
## Run OrthoFinder on the them
orthofinder -f Biv_Ortho
orthofinder -f Mol_Ortho
```
for file in *.fna
do
  filename="${file%.*}"
  echo $filename
  mv $file ${filename}.pep
done

#### alignment
# use mafft to align each orthogroup, then use trimAL to trim them
```bash
conda activate mafft
# use mafft to align each orthogroup, then use trimAL to trim them
# directory Single_Copy_Orthologue_Sequences
for i in *.fa
do
    mafft --adjustdirection $i > ${i%.fa}.aligned.fasta
    trimal -in ${i%.fa}.aligned.fasta -out ${i%.fa}.trimmed.fasta -automated1
done
### change the name in the alignment file, before: >Spe_xxxGene; After: >Spe
for i in *.trimmed.fasta
do
  cut -f 1 -d "." $i > ${i%.trimmed.fasta}_clean.fasta
done
```



```



#### other giant clams
```bash
nohup busco -m genome -i ./Tcrocea.fasta -o ./BUSCO_results/Tcrocea_metozoa_BUSCO  -l metazoa_odb10 -c 20 &
nohup busco -m genome -i ./Tgigas.fasta -o ./BUSCO_results/Tgigas_metozoa_BUSCO  -l metazoa_odb10 -c 20 &
nohup busco -m genome -i ./Tmaxima_hirise.fasta -o ./BUSCO_results/Tmaxima_metozoa_BUSCO  -l metazoa_odb10 -c 20 &

nohup busco -m genome -i ./Tcrocea.fasta -o ./BUSCO_results/Tridacna_crocea_mollusca_BUSCO  -l mollusca_odb10 -c 20 &
nohup busco -m genome -i ./Tgigas.fasta -o ./BUSCO_results/Tridacna_gigas_mollusca_BUSCO  -l mollusca_odb10 -c 20 &
nohup busco -m genome -i ./Tmaxima_hirise.fasta -o ./BUSCO_results/Tridacna_maxima_mollusca_BUSCO  -l mollusca_odb10 -c 20 &
```
