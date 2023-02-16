### Use Carbonic anhydrase to test OrthoFisher workflow
### Carbonic anhydrase 2-like in the giant clam,Tridacna squamosa
## Reference: Ip, Yuen K., et al. "Carbonic anhydrase 2‐like in the giant clam, Tridacna squamosa: characterization, localization, response to light, and possible role in the transport of inorganic carbon from the host to its symbionts." Physiological reports 5.23 (2017): e13494.
## Translate it on Transeq (EMBOSS)
## Blastp the sequence T_squamosa_CA.fasta on NCBI TBLASTX, download aligned sequence, align them on Clustal Omega
## Download
wget https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/clustalo-I20220714-043729-0955-71534832-p2m/aln-fasta
mv aln-fasta T_squamosa_CA_aln.fasta
## Install OrthoFisher and all the dependecies
conda create -n orthofisher
conda activate orthofisher
conda install -c bioconda clipkit phykit biokit
conda install -c bioconda hmmer iqtree mafft
conda install -c jlsteenwyk orthofisher
# build HMMs
hmmbuild T_squamosa_CA.hmm T_squamosa_CA_aln.fasta
# create a list of hmms
ls *hmm > hmms_list.txt
# create a list of genomes
paste <(ls *faa) <(ls *faa | sed 's/.faa//g') > prot_and_names_list.txt
# run orthofisher
nohup orthofisher -m hmms_list.txt -f prot_and_names_list.txt &
