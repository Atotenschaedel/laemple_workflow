import os
import subprocess
import logging

def unpack_data(path_to_archive, date, prefix="HHD"):
    pass

def run_qc_check(date, samples, prefix="HHD"):
    for sample in samples:
        try:
            os.mkdir(f"QC-reports/{prefix}{date}/{sample}")
        except:
            pass

        if os.path.isdir(f"QC-reports/{prefix}{date}/{sample}"):
            rc = subprocess.run(["fastqc", f"Data/{prefix}{date}/{sample}/*R1*", f"Data/{prefix}{date}/{sample}/*R1*",
                                 "-o", f"QC-reports/{prefix}{date}/{sample}", "-q"])
            if rc != 0:
                logging.warning("FastQC has returned non-zero exit code.")


def run_read_trimming(date, samples, prefix="HHD"):
    for sample in samples:
        try:
            os.mkdir(f"Processed-reads/{prefix}{date}/{sample}")
        except:
            pass

        if os.path.isdir(f"Processed-reads/{prefix}{date}/{sample}"):
            rc = subprocess.run(["readlength.sh", f"in=$(ls Data/{prefix}{date}/{sample}/*R1*)", f"in2=$(ls Data/{prefix}{date}/{sample}/*R2*)",
                                 f"out=Processed-reads/{prefix}{date}/{sample}/{sample}.length-preTrimmed.txt"], shell=True)
            rc = subprocess.run(["bbduk.sh", f"in1=$(ls Data/{prefix}{date}/{sample}/*R1*)", f"in2=$(ls Data/{prefix}{date}/{sample}/*R2*)",
                                 f"out1=Processed-reads/{prefix}{date}/{sample}/{sample}_1.trimmed.fastq", 
                                 f"out2=Processed-reads/{prefix}{date}/{sample}/{sample}_2.trimmed.fastq",
                                 "qtrim=rl", "trimq=15", "ref=data/adaptersPhiX.fa",
                                 f"stats=Processed-reads/{prefix}{date}/{sample}/{sample}.stats.log"], shell=True)


def run_read_mapping(date, samples, prefix="HHD"):
    for sample in samples:
        try:
            os.mkdir(f"Read-mapping/{prefix}{date}/{sample}")
        except:
            pass

        if os.path.isdir(f"Read-mapping/{prefix}{date}/{sample}"):
            rc = subprocess.run(["bwa", "mem", "-t", "4", SARS-CoV-2-reference.fasta Processed-reads/${1}/${f}/${f}_1.trimmed.fastq Processed-reads/${1}/${f}/${f}_2.trimmed.fastq > Read-mapping/$1/$f/$f.trimmed.sam])

            samtools view -S -b Read-mapping/$1/$f/$f.trimmed.sam > Read-mapping/$1/$f/$f.trimmed.bam;
            rm Read-mapping/$1/$f/$f.trimmed.sam;
            samtools sort Read-mapping/$1/$f/$f.trimmed.bam -o Read-mapping/$1/$f/$f.trimmed.sorted.bam;