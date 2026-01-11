#!/bin/bash

## Unpack the data
for f in $(ls Data/$1); do
	cd Data/$1/$f;
	gunzip *;
	cd ../../..;
done

## Generate QC reports
mkdir QC-reports/$1;

for f in $(ls Data/$1); do
	mkdir QC-reports/$1/$f;
	fastqc Data/$1/$f/*R1* Data/$1/$f/*R2* -o QC-reports/$1/$f -q
done

## Trim the reads
mkdir Processed-reads/$1;

for f in $(ls Data/$1); do
	readlength.sh in=$(ls Data/$1/$f/*R1*) in2=$(ls Data/$1/$f/*R2*) out=$(echo Processed-reads/${f}/$f.length-preTrimmed.txt);
	mkdir Processed-reads/$1/$f;
	bbduk.sh in1=$(ls Data/$1/$f/*R1*) in2=$(ls Data/$1/$f/*R2*) out1=$(echo Processed-reads/${1}/${f}/${f}_1.trimmed.fastq) \
	out2=$(echo Processed-reads/${1}/${f}/${f}_2.trimmed.fastq) qtrim=rl trimq=15 ref=adaptersPhiX.fa stats=$(echo Processed-reads/$1/$f/$f.stats.log);
done

## Map the reads
mkdir Read-mapping/$1;

for f in $(ls Data/$1); do
	mkdir Read-mapping/$1/$f;
	bwa mem -t 4 SARS-CoV-2-reference.fasta Processed-reads/${1}/${f}/${f}_1.trimmed.fastq Processed-reads/${1}/${f}/${f}_2.trimmed.fastq > Read-mapping/$1/$f/$f.trimmed.sam;
	samtools view -S -b Read-mapping/$1/$f/$f.trimmed.sam > Read-mapping/$1/$f/$f.trimmed.bam;
	rm Read-mapping/$1/$f/$f.trimmed.sam;
	samtools sort Read-mapping/$1/$f/$f.trimmed.bam -o Read-mapping/$1/$f/$f.trimmed.sorted.bam;
done

## Soft clip primer locations
for f in $(ls Data/$1); do
	ivar trim -b nCoV-2019.primer.bed -p Read-mapping/$1/$f/$f.clean.sorted -i Read-mapping/$1/$f/$f.trimmed.sorted.bam -q 15 -e;
done 

for f in $(ls Data/$1); do
	mv Read-mapping/$1/$f/$f.clean.sorted.bam Read-mapping/$1/$f/$f.clean.bam;
	samtools sort Read-mapping/$1/$f/$f.clean.bam -o Read-mapping/$1/$f/$f.clean.sorted.bam;
done 

## Get depth of coverage
mkdir Coverage/Coverage-${1:3}
for f in $(ls Data/$1); do
	mkdir Coverage/Coverage-${1:3}/$f
	samtools depth Read-mapping/$1/$f/$f.clean.sorted.bam > Coverage/Coverage-${1:3}/$f/$f.clean.sorted.coverage.txt;
done

## Call varaints: LoFreq
mkdir Variant-calling-LoFreq/$1;

for f in $(ls Data/$1); do
	lofreq indelqual --dindel -f SARS-CoV-2-reference.fasta Read-mapping/$1/$f/$f.clean.sorted.bam -o Variant-calling-LoFreq/$1/$f.clean.indelqual.bam;
	lofreq call -f SARS-CoV-2-reference.fasta --call-indels -o Variant-calling-LoFreq/$1/$f.clean.vcf Variant-calling-LoFreq/$1/$f.clean.indelqual.bam;
done

## Call variants: iVar
mkdir Variant-calling-iVar/$1;

for f in $(ls Data/$1); do
	samtools mpileup -aa -A -d 600000 -B -Q 0 Read-mapping/$1/$f/$f.clean.sorted.bam | ivar variants -p Variant-calling-iVar/$1/$f.clean \
	-q 20 -t 0 -r SARS-CoV-2-reference.fasta -g SARS-CoV-2-reference.gff3;
done

## Merge variant calls
mkdir Variant-calling-combined/$1;

for f in $(ls Data/$1); do
	python MergeCalls.py -m 0.02 -g SARS-CoV-2-reference.gff3 -o Variant-calling-combined/$1/$f-merged-indel.clean.tsv \
	Variant-calling-iVar/$1/$f.clean.tsv Variant-calling-LoFreq/$1/$f.clean.vcf
done
