########1. BUSCO #########
conda activate BUSCO
mkdir BUSCO_results
cd BUSCO_results
#pwd: /home/ruiqi/ruiqi_data/GiantClamGenome/BUSCO_results
nohup busco -m genome -i ../Tmaxima_hirise.fasta -o Tmaxima_mollusca_BUSCO  -l mollusca_odb10 -c 10 &


#######2. GenomeQC ########
# citation: GenomeQC: A quality assessment tool for genome assemblies and gene structure annotations
#Nancy Manchanda, John L. Portwood II, Margaret R. Woodhouse, Arun S. Seetharam,
#Carolyn J. Lawrence-Dill, Carson M. Andorf, Matthew B. Hufford

# https://genomeqc.maizegdb.org
# Estimated genome size (Mb):
# AUGUSTUS Species: fly



####### 3. PDR #####
### Luyu Xie, Limsoon Wong, PDR: a new genome assembly evaluation metric based on genetics concerns
wget https://github.com/XLuyu/PDR/releases/download/v0.6.4/PDRi.jar

conda create -n PDR openjdk=11
conda activae PDR
conda install -c anaconda java-1.8.0-openjdk-devel-cos7-s390x
java -version
conda install -c bioconda bwa

nohup java -jar PDRi.jar GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna ../Tmaxima_hirise.fasta --threads 20 &


######## 4. Quast ######
conda create --name Quast
conda activate Quast
conda install -c bioconda quast

nohup ~/miniconda3/envs/Quast/opt/quast-5.0.2/quast.py ../Tmaxima_hirise.fasta \
        -o quast_output \
        -t 20 \
        --eukaryote \
        --large \
        --k-mer-stats \
        --k-mer-size 19 &

### Error when run with -g
nohup ~/miniconda3/envs/Quast/opt/quast-5.0.2/quast.py ../Tmaxima_hirise.fasta \
        -o quast_output \
        -t 20 \
        --eukaryote \
        --large \
        --k-mer-stats \
        --k-mer-size 19 \
        -g ../Annotation/PO1429_Tridacna_maxima.annotation.gff \
        -r ../Tmaxima_hirise.fasta &


### Kmer analysis
### KAT
conda create --name Kat
conda activate Kat
conda install -c bioconda kat

# Histogram
nohup kat hist -t 20 ../PacBioReads/XDOVE*.fastq.gz &

# GC content and kmer coverage
nohup kat gcp -t 20 ../PacBioReads/XDOVE*.fastq.gz &
