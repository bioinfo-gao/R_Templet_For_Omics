# projecton
# PCA TSNE UMAP
# RPKM TPM
# featureCount = rsem in function, but rsem better for reads span multiple exons
# conda install -c bioconda cutadapt
# conda install -y bioconda::cutadapt 
# conda install bioconda/label/cf201901::cutadapt

conda -n rna python=3

conda activate rna 

conda install bioconda::sra-tools 
conda install -y bioconda::hisat2

conda install -y bioconda::multiqc 
conda install -y bioconda::trim-galore 


conda config --add channels bioconda
conda config --add channels conda-forge
# conda config --set channel_priority strict



# conda search sra-tools
# conda install bioconda/label/cf201901::sra-tools
# /home/gao/anaconda3/envs/rna/bin/prefetch : 2.10.0
# sudo find . -type f -exec mv {} destination/ \;
# ls ./scRNA/*/*.sra 


prefetch --help
prefetch SRP133642 -O /home/gao/Desktop/Code/Bioinfo/scRNA/



####================
conda env list
conda activate rna
cat SRR_Acc_List.txt | while read id; do (nohup prefetch $id -X 100G  -O ./scRNA &); done     

###==================
mv ./scRNA/*/*.sra ./sra/ # move all file in subdiretory to another folder


###==================
cutadapt=cutadapt-venv/bin/cutadapt
###==================
ls /home/gao/Desktop/Code/Bioinfo/sra/ |
while read id
do 
nohup fastq-dump -gzip --split-3 -O ./fq ${id} &
done     



sudo apt install python3-virtualenv

virtualenv cutadapt-venv
cutadapt-venv/bin/pip --upgrade pip
cutadapt-venv/bin/pip install cutadapt

cutadapt-venv/bin/cutadapt --version