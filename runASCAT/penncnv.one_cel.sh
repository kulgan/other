#!/bin/bash
# take options
uuid="$1"
file="$2"

# SNP6 probe definitions from Affy
cdf='/home/ubuntu/SCRATCH/snp6cnv/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.cdf'
model='/home/ubuntu/SCRATCH/snp6cnv/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.birdseed.models'
special='/home/ubuntu/SCRATCH/snp6cnv/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.specialSNPs'

          # Additional SNP6 supporting files from PennCNV
target='/home/ubuntu/SCRATCH/snp6cnv/PennCNV-1.0.5/affy/libgw6/hapmap.quant-norm.normalization-target.txt'
pfb='/home/ubuntu/SCRATCH/snp6cnv/PennCNV-1.0.5/affy/libgw6/affygw6.hg38.pfb'
genocluster='/home/ubuntu/SCRATCH/snp6cnv/PennCNV-1.0.5/affy/libgw6/hapmap.genocluster'

# real command
cd /home/ubuntu/SCRATCH/snp6cnv/
/home/ubuntu/SCRATCH/snp6cnv/gdc-client download $uuid -t ~/token
cd /home/ubuntu/SCRATCH/snp6cnv/$uuid
/home/ubuntu/SCRATCH/snp6cnv/apt-2.10.2.1-x86_64-intel-linux/bin/apt-probeset-genotype --cdf-file  $cdf --analysis birdseed --read-models-birdseed $model --special-snps $special --out-dir . --cels $file
/home/ubuntu/SCRATCH/snp6cnv/apt-2.10.2.1-x86_64-intel-linux/bin/apt-probeset-summarize --cdf-file $cdf --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch $target --out-dir . --cels $file
/home/ubuntu/SCRATCH/snp6cnv/PennCNV-1.0.5/affy/bin/normalize_affy_geno_cluster.pl $genocluster quant-norm.pm-only.med-polish.expr.summary.txt -locfile $pfb -out ../output/$uuid.lrr_baf.txt
cd ../output/
gzip $uuid.lrr_baf.txt
cd /home/ubuntu/SCRATCH/snp6cnv/
rm -rf /home/ubuntu/SCRATCH/snp6cnv/$uuid
