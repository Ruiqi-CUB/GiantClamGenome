########1. BUSCO #########
conda activate BUSCO
mkdir BUSCO_results
cd BUSCO_results
#pwd: /home/ruiqi/ruiqi_data/GiantClamGenome/BUSCO_results
nohup busco -m genome -i ./GCA_943736015.1_xbTriCroc3.1_genomic.fna -o Tcrocea_mollusca_BUSCO  -l mollusca_odb10 -c 20 &
nohup busco -m genome -i ./GCA_943736015.1_xbTriCroc3.1_genomic.fna -o Tcrocea_metazoa_BUSCO  -l metazoa_odb10 -c 20 &

####QUAST
conda activate Quast

nohup ~/miniconda3/envs/Quast/opt/quast-5.0.2/quast.py ./GCA_943736015.1_xbTriCroc3.1_genomic.fna \
        -o quast_output \
        -t 20 \
        --eukaryote \
        --large \
        --k-mer-stats \
        --k-mer-size 19 &


## Reannotation
nohup busco -m protein -i ../Tmax_merged_longestisoform_prot.fasta -o Tmaxima_reannotation_mollusca_BUSCO  -l mollusca_odb10 -c 36 &
nohup busco -m protein -i ../Tmax_merged_longestisoform_prot.fasta -o Tmaxima_reannotation_metazoa_BUSCO  -l metazoa_odb10 -c 44 &
