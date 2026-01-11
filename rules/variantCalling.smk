# -*- coding: utf-8 -*-
# snakemake workflow
"""
Created on 21-05-2023
@author: aschedl
Description: This workflow uses raw sequencing data and filters, maps and calls variants as well as checks for quality parameter.

Example usage:
snakemake --snakefile variantCalling.smk --use-conda -c1
"""
# get parameters from config file
import yaml
from pathlib import Path

configfile: "config/workflow_config.yaml"

#get parameters from config file
EXPERIMENT_NAME = config["EXPERIMENT_NAME"]

p = Path.cwd() / "experiments" / EXPERIMENT_NAME / "data"
SAMPLES = list(set([x.stem[:-3] for x in p.glob("*") if x.is_file()]))

rule all:
     input:
        expand("experiments/{exp}/results/variantCall/12_phasedLofreqVcf/{sample}_samp.lofreq_filtered_phased.vcf.gz",
                 exp=EXPERIMENT_NAME, sample=SAMPLES), 
        expand("experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}{ext}",
                exp=EXPERIMENT_NAME, sample=SAMPLES, ext=[".stats", "_virus.stats","_cov.png", ".bw", "_clip_viral_stats.txt"])
                # always create stat files as well
         
rule merge_sample:
    # merge different runs of the same sample
    input:  
        # unmapped, paired-end bam file
        bam1="experiments/{exp}/data/{sample}_1.bam",
        bam2="experiments/{exp}/data/{sample}_2.bam"
    output:
        fastq1= "experiments/{exp}/results/variantCall/01_merged/{sample}_npa_1.fastq.gz",
        fastq2= "experiments/{exp}/results/variantCall/01_merged/{sample}_npa_2.fastq.gz"
    conda:
        "../envs/picard.yaml"
    log: 
        # picard sends logging output to stderr
        out_merge= "logs/{exp}/variantCall/{sample}/01_{sample}.merge.out",
        out_sam= "logs/{exp}/variantCall/{sample}/01_{sample}.samfastq.out"
    shell:
        "picard MergeSamFiles -SO unsorted "
        "-I {input.bam1} -I {input.bam2} -O /dev/stdout 2>{log.out_merge} | "
        "picard SamToFastq -INPUT /dev/stdin "
        "-FASTQ {output.fastq1} -SECOND_END_FASTQ {output.fastq2} "
        "2>{log.out_sam}"

rule fastqInput:
    # in case input files are already in fastq format and need no merge
    input:
        # paired - end, forward and reverese reads 
        fastq1="experiments/{exp}/data/{sample}_R1.fastq",
        fastq2="experiments/{exp}/data/{sample}_R2.fastq"
    output:
        fastq1="experiments/{exp}/results/variantCall/01_merged/{sample}_npa_1.fastq.gz",
        fastq2="experiments/{exp}/results/variantCall/01_merged/{sample}_npa_2.fastq.gz"
    shell:
        "gzip -c {input.fastq1} > {output.fastq1} && "
        "gzip -c {input.fastq2} > {output.fastq2} "

rule correct_overlapPE:
    # correct overlapping Pair end
    input:
        fastq1= "experiments/{exp}/results/variantCall/01_merged/{sample}_npa_1.fastq.gz",
        fastq2= "experiments/{exp}/results/variantCall/01_merged/{sample}_npa_2.fastq.gz"
    output:
        cor1= "experiments/{exp}/results/variantCall/02_correctOverlapPE/{sample}_ecco_1.fq.gz",
        cor2= "experiments/{exp}/results/variantCall/02_correctOverlapPE/{sample}_ecco_2.fq.gz"
    params:
        minReadLength=40,
        mininsert=25
    conda:
        "../envs/bbmap.yaml"
    log:
        # bbmerge sends logging output to stderr
        out= "logs/{exp}/variantCall/{sample}/02_{sample}.correctOverlapPE.out"
    shell:
        "bbmerge.sh -Xmx4g in1={input.fastq1} in2={input.fastq2} "
        "out={output.cor1} out2={output.cor2} minlength={params.minReadLength} "
        "mininsert={params.mininsert} tbo ecco mix ordered=t qtrim2=r t=10 "
        "2>{log.out}"

rule generate_comb_genome:
    input:
        ref1=config["REFERENCES"]["REF_GENOME"],
        ref2=config["REFERENCES"]["REF_HUMAN"] 
    output:
        "reference/indices_for_BWA/refSeq_comb.fa"
    shell:
        # conc all necessary genomes
        "sed '/^$/d' {input.ref1} {input.ref2} > {output} "


rule index_genome:
    input:
        "reference/indices_for_BWA/refSeq_comb.fa" # index reference genome for BWA in fasta format
    output:
        "reference/indices_for_BWA/refSeq_comb.fa.amb"
    conda:
        "../envs/bwa.yaml"
    log:
        # bwa tool sends logging output to stderr
        out="logs/reference_index_Refgenmoe.out"
    shell: 
        "bwa index {input} 2>{log.out}"

rule read_mapping:
    # alignment via Burrows-Wheeler transformation (bwa)
    input:
        # needs indexed reference genome and corrected bam files
        ref_genome_amb="reference/indices_for_BWA/refSeq_comb.fa.amb", # viral and human reference sequence combined in fasta format
        ref_genome="reference/indices_for_BWA/refSeq_comb.fa",
        read1= "experiments/{exp}/results/variantCall/02_correctOverlapPE/{sample}_ecco_1.fq.gz",
        read2= "experiments/{exp}/results/variantCall/02_correctOverlapPE/{sample}_ecco_2.fq.gz"
    output: 
        # sorted bam file
        "experiments/{exp}/results/variantCall/03_readMapping/{sample}_sorted.bam"
    conda:
        "../envs/bwa.yaml"
    log: 
        # bwa tool sends logging output to stderr
        bwa_out= "logs/{exp}/variantCall/{sample}/03_{sample}.readMappingBWA.out"
    shell:
        # Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. 
        # The read group ID will be attached to every read in the output.
        "bwa mem -k 17 -r 1.25 -M -t 12 -R "
        "'@RG\\tID:{wildcards.sample}_npa_npp"
        "\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:ILLUMINA"
        "\\tPU:HL37TBBXX.5' " ### what is this flowcell no?
        "{input.ref_genome} {input.read1} {input.read2} "
        "2>{log.bwa_out} "
        "| samtools view -Shb "
        "| samtools sort -T "
        "experiments/{wildcards.exp}/results/variantCall/03_readMapping/{wildcards.sample}tmpBWA "
        " > {output} "

rule mapping_stats_1:
    input:
        # sorted bam file
        "experiments/{exp}/results/variantCall/03_readMapping/{sample}_sorted.bam"
    output:
        # bam index, different mapping stats
        idx="experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}.idxstats",
        stats="experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}.stats"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input} && "
        "samtools idxstats {input} > {output.idx} && "
        "samtools stats {input} > {output.stats}"

rule create_seqDictonary:
    input: 
        config["REFERENCES"]["REF_GENOME"] # reference fasta
    output:
        "reference/SequenceDictionaries/RefSeq_sequence.dict" # ref sequence dictionary format
    conda:
        "../envs/picard.yaml"
    log:
        # picard sends logging output to stderr
        out="logs/seqDict_refGenome.out"
    shell:
        "picard CreateSequenceDictionary -R {input} -O {output} "
        "2> {log.out}"

rule remove_nonPrimAndHuman:
    input:
        # sorted bam file
        bam="experiments/{exp}/results/variantCall/03_readMapping/{sample}_sorted.bam",
        ref_genome_dict="reference/SequenceDictionaries/RefSeq_sequence.dict" # ref sequence dictionary format
    output:
        bam_file="experiments/{exp}/results/variantCall/04_removeNoPrimHuman/{sample}_sorted_viral_rh.bam",
        bam_index="experiments/{exp}/results/variantCall/04_removeNoPrimHuman/{sample}_sorted_viral_rh.bam.bai"
    conda:
        "../envs/picard.yaml"
    log:
        # picard sends logging output to stderr
        pic_out="logs/{exp}/variantCall/{sample}/04_{sample}.removePicard.out"
    shell: # reference sequence dictionary is required first
        "picard ReorderSam -I {input.bam} -SD {input.ref_genome_dict}  "
        "-O /dev/stdout -S true 2>{log.pic_out} | samtools view -b -f 2 -F 256 "
        ">{output.bam_file} && samtools index {output.bam_file}"

rule softclip_primer:
    input:
        primer_bed = config["REFERENCES"]["PRIMER_SET"], # primer set in bed format
        bam_file="experiments/{exp}/results/variantCall/04_removeNoPrimHuman/{sample}_sorted_viral_rh.bam"
    output:
        "experiments/{exp}/results/variantCall/05_softclipPrimer/{sample}_sorted_viral_rh_trim.bam"
    conda:
        "../envs/ivar.yaml"
    log:
        # ivar sends logging output to stderr
        out_ivar="logs/{exp}/variantCall/{sample}/05_{sample}.softclipIvar.out"
    shell:
        "ivar trim -b {input.primer_bed} -s 4 -q 0 -m 20 -e -i "
        "{input.bam_file} -p {output}_tmp 2>{log.out_ivar} && "
        "samtools sort -o {output} {output}_tmp.bam && "
        "samtools index {output} && "
        "rm {output}_tmp.bam "

rule mapping_stats_2:
    input:
        # sorted bam file
        "experiments/{exp}/results/variantCall/05_softclipPrimer/{sample}_sorted_viral_rh_trim.bam"
    output:
        # bam index, different mapping stats
        index="experiments/{exp}/results/variantCall/05_softclipPrimer/{sample}_sorted_viral_rh_trim.bam.bai",
        flagstat="experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}_virus.flagstat",
        idx="experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}_virus.idxstats",
        stats="experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}_virus.stats"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input} && "
        "samtools flagstat {input} > {output.flagstat} && "
        "samtools idxstats {input} > {output.idx} && "
        """samtools view {input} | awk '{{sum+=$5}} END {{print "# Mean MAPQ =",sum/NR}}' >{output.stats} && """
        "samtools stats {input} >>{output.stats}"

rule calc_coverage:
    input:
        # viral only bam files
        bam="experiments/{exp}/results/variantCall/05_softclipPrimer/{sample}_sorted_viral_rh_trim.bam",
        index="experiments/{exp}/results/variantCall/05_softclipPrimer/{sample}_sorted_viral_rh_trim.bam.bai"
    output:
        # only coverage of viral?
        bigwig="experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}.bw",
        bedgraph="experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}.bg"
    conda:
        "../envs/deeptools.yaml"
    log:
        out_bw="logs/{exp}/variantCall/{sample}/00_{sample}.calc_cov_bigwig.out",
        out_bg="logs/{exp}/variantCall/{sample}/00_{sample}.calc_cov_bedgraph.out"
    shell: 
        # bamCoverage sends logging output to stderr
        "bamCoverage --bam {input.bam} -o {output.bigwig} --binSize 10 "
        "2>{log.out_bw} && "
        "bamCoverage --bam {input.bam} -o {output.bedgraph} "
        "--outFileFormat bedgraph --binSize 10 "
        "2>{log.out_bg}"

rule calc_coverage_COV:
    input:
        # viral only bam file
        all_cov="experiments/{exp}/results/variantCall/05_softclipPrimer/{sample}_sorted_viral_rh_trim.bam",
        max_cov="experiments/{exp}/results/variantCall/06_lofreqViterbi/{sample}_real_viterbi_indelQual_max_cov_10000.bam"
    output:
        # coverage file (chrom, pos, coverage) with header
        all_cov = "experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}_cov.tsv",
        max_cov = "experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}_max_cov_10000.tsv"
    conda:
        "../envs/samtools.yaml"
    shell:
        "echo -e '#CHROM\tPOS\t{wildcards.sample}' > {output.all_cov} && "
        "samtools depth -d 0 -Q 20 -a {input.all_cov} >> {output.all_cov} && "
        "echo -e '#CHROM\tPOS\t{wildcards.sample}' > {output.max_cov} && "
        "samtools depth -d 0 -Q 20 -a {input.max_cov} >> {output.max_cov} "

rule plot_coverage_COV:
    input:
        all_cov="experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}_cov.tsv",
        max_cov="experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}_max_cov_10000.tsv"
    output:
        all_cov="experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}_cov.png",
        max_cov="experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}_max_cov_10000.png"
    conda:
        "../envs/python3.yaml"
    log:
        err="logs/{exp}/variantCall/{sample}/00_{sample}.plot_cov.err",
    shell:
        "python3 bin/custom_scripts/plot_coverage.py "
        "-i {input.all_cov} -o {output.all_cov} -s {wildcards.sample} 2>{log} && "
        "python3 bin/custom_scripts/plot_coverage.py "
        "-i {input.max_cov} -o {output.max_cov} -s {wildcards.sample} 2>>{log}"

rule lofreqViterbi:
    # low frequency viterbi realignment
    input:
        ref=config["REFERENCES"]["REF_GENOME"], # reference sequence in fasta format
        # viral only, trimmed bam file
        bam="experiments/{exp}/results/variantCall/05_softclipPrimer/{sample}_sorted_viral_rh_trim.bam"
    output:
        "experiments/{exp}/results/variantCall/06_lofreqViterbi/{sample}_real_viterbi.bam"
    conda:
        "../envs/lofreq.yaml"
    log:
        out = "logs/{exp}/variantCall/{sample}/06_{sample}.viterbi.out",
        err = "logs/{exp}/variantCall/{sample}/06_{sample}.viterbi.err"
    shell:
        "lofreq viterbi -f {input.ref} {input.bam} 2>{log.err} | "
        "samtools sort -T {output}_tmp -o {output} && "
        "samtools index {output} " 

rule lofreqIndelqual:
    # add low frequency indel qualities
    input:
        ref=config["REFERENCES"]["REF_GENOME"], # reference sequence in fasta format
        # low frequency bam file
        bam="experiments/{exp}/results/variantCall/06_lofreqViterbi/{sample}_real_viterbi.bam"
    output:
        bam="experiments/{exp}/results/variantCall/06_lofreqViterbi/{sample}_real_viterbi_indelQual.bam",
        index="experiments/{exp}/results/variantCall/06_lofreqViterbi/{sample}_real_viterbi_indelQual.bam.bai",
        idx="experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}_real_viterbi_indelQual.idxstats",
    conda:
        "../envs/lofreq.yaml"
    log:
        out="logs/{exp}/variantCall/{sample}/06_{sample}.indelqual.out",
        err="logs/{exp}/variantCall/{sample}/06_{sample}.indelqual.err",
    shell:
        "lofreq indelqual --dindel -f {input.ref} -o {output.bam} {input.bam} "
        ">{log.out} 2>{log.err} && "
        "samtools index {output.bam} && "
        "samtools idxstats {output.bam} > {output.idx}"

rule downsample:
    # downsample to 10k reads per region
    input:
        "experiments/{exp}/results/variantCall/06_lofreqViterbi/{sample}_real_viterbi_indelQual.bam",
    output:
        "experiments/{exp}/results/variantCall/06_lofreqViterbi/{sample}_real_viterbi_indelQual_max_cov_10000.bam"
    conda:
        "../envs/downsample.yaml"
    log:
        # python script sends logging output to stderr
        out="logs/{exp}/variantCall/{sample}/06_{sample}.downsample.out"
    shell: 
        "python3 bin/custom_scripts/downsample_to_cov.py --start "
        "--conda-env -b {input} -c 10000 2>{log.out}"

rule clipOverlap:
    input:
        "experiments/{exp}/results/variantCall/06_lofreqViterbi/{sample}_real_viterbi_indelQual_max_cov_10000.bam"
    output:
        bam="experiments/{exp}/results/variantCall/07_clipOverlap/{sample}_ds_clip_viral.bam",
        stats="experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}_clip_viral_stats.txt"
    conda:
        "../envs/bamutil.yaml"
    # no log file created, all logging information written to output.stats file
    shell:
        "bam clipOverlap --poolSize 5000000 --in {input} "
        "--out {output.bam} --stats 2>{output.stats} && "
        "samtools index {output.bam}"

rule lofreqCall:
    input:
        ref=config["REFERENCES"]["REF_GENOME"], # reference sequence in fasta format
        bam="experiments/{exp}/results/variantCall/07_clipOverlap/{sample}_ds_clip_viral.bam"
    output:
        "experiments/{exp}/results/variantCall/08_lofreqCall/{sample}_lofreq.vcf"
    conda:
        "../envs/lofreq.yaml"
    log:
        out="logs/{exp}/variantCall/{sample}/08_{sample}.lofreqCall.out",
        err="logs/{exp}/variantCall/{sample}/08_{sample}.lofreqCall.err"
    shell:
        "echo starting... > {log.out} && "
        "lofreq call -f {input.ref} -C 75 --no-default-filter --call-indels "
        "-o {output} {input.bam} >>{log.out} 2>{log.err} "
        "&& echo finished >>{log.out}" 

rule lofreqFilt:
    input:
        "experiments/{exp}/results/variantCall/08_lofreqCall/{sample}_lofreq.vcf"
    output:
        "experiments/{exp}/results/variantCall/09_lofreqFilt/{sample}_lofreq_filter.vcf"
    conda:
        "../envs/lofreq.yaml"
    log:
        out="logs/{exp}/variantCall/{sample}/09_{sample}.lofreqFilt.out"
    shell:
        "lofreq filter -i {input} -v 75 -a 0.01 --no-defaults -Q 90 --print-all | "
        "bcftools filter -m + -s 'HRUN_gt_3' -e 'HRUN >= 4' > {output} "
        " 2>{log.out}"

rule reformLofreq:
    input:
        "experiments/{exp}/results/variantCall/09_lofreqFilt/{sample}_lofreq_filter.vcf"
    output:
        vcf="experiments/{exp}/results/variantCall/10_lofreqRefor/{sample}_samp.lofreq.vcf",
        gz="experiments/{exp}/results/variantCall/10_lofreqRefor/{sample}_samp.lofreq.vcf.gz",
        tbi="experiments/{exp}/results/variantCall/10_lofreqRefor/{sample}_samp.lofreq.vcf.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        "sh bin/shell_scripts/reformLofreq.sh -i {input} -o {output.vcf} -s {wildcards.sample} && "
        "bgzip -f -k {output.vcf} && "
        "tabix -f -p vcf {output.gz}"


rule refiltLoFreq:
    input:
        "experiments/{exp}/results/variantCall/10_lofreqRefor/{sample}_samp.lofreq.vcf.gz"
    output:
        gz="experiments/{exp}/results/variantCall/11_lofreqrefilt/{sample}_samp.lofreq_filtered.vcf.gz",
        tbi="experiments/{exp}/results/variantCall/11_lofreqrefilt/{sample}_samp.lofreq_filtered.vcf.gz.tbi"
    conda:
        "../envs/lofreq.yaml"
    shell:
        'bcftools view -i "FORMAT/AF >= 0.01 && QUAL > 90" {input} | '
        "bgzip -c >{output.gz} && "
        "tabix -f -p vcf {output.gz} "

rule phasedLofreqVcf:
    input:
        ref=config["REFERENCES"]["REF_GENOME"], # reference sequence in fasta format
        gz="experiments/{exp}/results/variantCall/11_lofreqrefilt/{sample}_samp.lofreq_filtered.vcf.gz",
        bam="experiments/{exp}/results/variantCall/06_lofreqViterbi/{sample}_real_viterbi_indelQual_max_cov_10000.bam"
    output:
        temp("experiments/{exp}/results/variantCall/12_phasedLofreqVcf/{sample}_samp.lofreq_filtered_phased_temp.vcf.gz")
    conda:
        "../envs/freebayes.yaml"
    shell:
        "freebayes -f {input.ref} -F 0.005 -C 1 --pooled-continuous --haplotype-length 3 "
        "--haplotype-basis-alleles {input.gz} {input.bam} >{output} "
        
rule add_afs_freebayes:
    input:
        ref=config["REFERENCES"]["REF_GENOME"], # reference sequence in fasta format
        vcf="experiments/{exp}/results/variantCall/12_phasedLofreqVcf/{sample}_samp.lofreq_filtered_phased_temp.vcf.gz"
    output:
        vcf="experiments/{exp}/results/variantCall/12_phasedLofreqVcf/{sample}_samp.lofreq_filtered_phased.vcf.gz",
        tbi="experiments/{exp}/results/variantCall/12_phasedLofreqVcf/{sample}_samp.lofreq_filtered_phased.vcf.gz.tbi"
    conda:
        "../envs/pysam.yaml"
    shell:
        "python3 bin/custom_scripts/add_afs_freebayes.py -i {input.vcf} | "
        "bcftools norm -m -any -f {input.ref} -O z > {output.vcf} && "
        "tabix -f -p vcf {output.vcf}"