#!/bin/bash
#Aly Khan's code for HLA typing
#RUN JUST ONCE
phlatdir=/home/ubuntu/SCRATCH/aly/phlat-release   #PATH TO PHLAT DIR
datadir=/run/shm #PATH TO RAM DISK
indexdir=/home/ubuntu/SCRATCH/aly/phlat-release/b2folder #PATH TO HLA INDEX
rsdir=/home/ubuntu/SCRATCH/aly  #RESULTS FOLDER
b2url=/home/ubuntu/bin/bowtie2-2.2.4  #PATH TO BOWTIE2 MAIN EXEC
################


############ FOR EACH PATIENT OR BAM FILE:
rm /run/shm/*  # CLEAR RAM DISK

#WE DO SAMTOOLS EXTRACT chr6:29000000-34000000 
samtools view /home/ubuntu/SCRATCH/0eaa3eb7-e2ed-47bd-9ac6-93441dad3b9a/TCGA-AB-2969-11A-01D-0739-09_whole.bam chr6:29000000-34000000  | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' | pigz --fast --processes 30 >> ${datadir}/tmp.fastq.gz

#WE DO SAMTOOLS EXTRACT 6:29000000-34000000
samtools view /mnt/cinder/ak_project/1ebc27b3-f55f-4ba3-a720-8b4920e1c102/C484.TCGA-02-0055-01A-01D-1490-08.4.bam 6:29000000-34000000 | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' | pigz --fast --processes 30 >> ${datadir}/tmp.fastq.gz

#WE DECLARE PATIENT ID
tag="0eaa3eb7-e2ed-47bd-9ac6-93441dad3b9a”

#RUN HLA TYPING TOOL
python2.7 -O ${phlatdir}/dist/PHLAT.py -1 ${datadir}/tmp.fastq.gz -index $indexdir -b2url $b2url -tag $tag -e $phlatdir -o $rsdir -p 29 -pe 0


