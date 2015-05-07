#!/bin/bash
#SBATCH --workdir=/home/ubuntu/
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nodelist=compute004

cd ~/SCRATCH
s3cmd get s3://bioinformatics_scratch/phlat-release.tar.gz
tar vxfz phlat-release.tar.gz
mkdir -p aly/input
mv phlat-release aly
rm phlat-release.tar.gz

