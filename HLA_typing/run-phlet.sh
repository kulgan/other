#!/bin/bash
#This is Zhenyu's revised code to run on elasticluster

#SBATCH --workdir=/home/ubuntu/
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nodelist=compute00##XXX##

export PYTHONPATH=/home/ubuntu/lib/python2.7/dist-packages/

phlatdir=/home/ubuntu/SCRATCH/aly/phlat-release   #PATH TO PHLAT DIR
datadir=/run/shm #PATH TO RAM DISK
indexdir=/home/ubuntu/SCRATCH/aly/phlat-release/b2folder #PATH TO HLA INDEX
rsdir=/home/ubuntu/SCRATCH/aly  #RESULTS FOLDER
b2url=/home/ubuntu/bin/bowtie2-2.2.3/bowtie2  #PATH TO BOWTIE2 MAIN EXEC
inputdir=/home/ubuntu/SCRATCH/aly/input

while read tag id tag2 bam; do

	############ FOR EACH PATIENT OR BAM FILE:
	rm /run/shm/*  # CLEAR RAM DISK
	rm $inputdir/*
	cd $inputdir
	s3cmd get s3://tcga_cghub_protected/$id/*

	#WE DO SAMTOOLS EXTRACT chr6:29000000-34000000
	samtools view $bam chr6:29000000-34000000  | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' | pigz --fast --processes 30 >> ${datadir}/tmp.fastq.gz

	#WE DO SAMTOOLS EXTRACT 6:29000000-34000000
	samtools view $bam 6:29000000-34000000 | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' | pigz --fast --processes 30 >> ${datadir}/tmp.fastq.gz

	#RUN HLA TYPING TOOL
	python2.7 -O ${phlatdir}/dist/PHLAT.py -1 ${datadir}/tmp.fastq.gz -index $indexdir -b2url $b2url -tag $tag -e $phlatdir -o $rsdir -p 29 -pe 0

done < ~/aly/xaa
