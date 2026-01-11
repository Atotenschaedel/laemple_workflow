# -*- coding: utf-8 -*-
# snakemake workflow
"""
Created on 27-06-2023
@author: aschedl
Description: This workflow generates simulated wastewater sequencing data for testing, 
    either based on a logistic growth model or real timecourse data.

Example usage:
snakemake --snakefile simulation.smk --use-conda -c1
"""
configfile: "./config/workflow_config.yaml"
import os

#get parameters from config file
EXPERIMENT_NAME = config["EXPERIMENT_NAME"]
SET_NAME=config["CURRENT_SET"]

# get newly generated file names
def get_abundances_file_names(wildcards):
    check_output=checkpoints.generateAbundances.get(**wildcards).output[0]
    s_name, = glob_wildcards(os.path.join(check_output, "{sample}.tsv"))
    s_name = [s for s in s_name if not s.endswith("_summary")]
    return expand(os.path.join(check_output, "{SAMPLE}.tsv"), SAMPLE=s_name)

def get_sequences_file_names(wildcards):
    seq_output=checkpoints.generateAbundances.get(**wildcards).output[0]
    s_name, = glob_wildcards(os.path.join(seq_output, "{sample}.tsv"))
    return expand(os.path.join("experiments/"+EXPERIMENT_NAME+"/data", "{SM}_{ext}"), SM=s_name, ext=["R1.fastq","R2.fastq"])

## check qc for simulated files - adapter sequences ect
rule all:
    input:
        expand("experiments/{exp}/simulation/QualityControl/{exp}_multiqc_report.html", exp=EXPERIMENT_NAME),
        expand("experiments/{exp}/simulation/{exp}_metadata.tsv", exp=EXPERIMENT_NAME)

checkpoint generateAbundances:
    output: 
        dir=directory("experiments/{exp}/simulation/abundances/"),
        metafile="experiments/{exp}/simulation/{exp}_metadata.tsv",
        data="experiments/{exp}/simulation/{exp}_data.csv"
    conda:
        "../envs/python3.yaml"
    params:
        seed=config["APP"]["SEED"],
        sample_mode=config["SIMULATION"]["SAMPLE_COLLECTION"]["SAMPLE_MODE"],
        sample_number=config["SIMULATION"]["SAMPLE_COLLECTION"]["NUMBER"],
        real_timecourse=config["SIMULATION"]["DATASET"][SET_NAME]["USE_REAL_TIMECOURSE"],
        real_timecourse_data=config["SIMULATION"]["DATASET"][SET_NAME]["TIMECOURSE_DATA"],
        max_timepoints=config["SIMULATION"]["DATASET"][SET_NAME]["MAX_TIMEPOINTS"],
        variants=config["SIMULATION"]["DATASET"][SET_NAME]["VARIANTS"],
        form=config["SIMULATION"]["DATASET"][SET_NAME]["FORM_PARA"],
        constant_abundance=config["SIMULATION"]["DATASET"][SET_NAME]["CONSTANT_VARIANT_ABUNDANCE"],
        start_date=config["SIMULATION"]["DATASET"][SET_NAME]["START_DATE"]
    log:
        "logs/{exp}/simulation/genSamples.out.log"
    shell: 
        "python3 bin/custom_scripts/generateSamples.py --experiment_name {wildcards.exp} "
        "--seed {params.seed} --sample_mode {params.sample_mode} "
        "--sample_number {params.sample_number} --real_timecourse {params.real_timecourse} "
        "--real_timecourse_data '{params.real_timecourse_data}' "
        '--max_timepoints {params.max_timepoints} -v "{params.variants}" '
        '--form "{params.form}" '
        "--constant_abundance {params.constant_abundance} --start_date {params.start_date} "
        "--output_dir {output.dir} &>{log}"

rule generate_cons_genome:
    input:
        config["SIMULATION"]["REFERENCES_SEQUENCES"]
    output:
        temp("experiments/{exp}/simulation/genomes.fasta")
    shell:
        # conc all necessary genomes. make sure that reads do not contain gaps
        "cat {input}| sed  '/^>/!s/-//g' >{output} "

checkpoint init_swampy:
    input:
        tsv=get_abundances_file_names,
        genome=expand("experiments/{exp}/simulation/genomes.fasta", exp=EXPERIMENT_NAME)
    output: 
        r1="experiments/{exp}/data/{exp}_{sample}_R1.fastq",
        r2="experiments/{exp}/data/{exp}_{sample}_R2.fastq"
    conda: 
        "../envs/swampy.yaml"
    params:
        exp=EXPERIMENT_NAME,
        primer_set=config["SIMULATION"]["PRIMER_SET"],
        n_reads=config["SIMULATION"]["N_READS"],
        seqSys=config["SIMULATION"]["SEQ_SYS"],
        read_length=config["SIMULATION"]["READ_LENGTH"],
        amplicon_distribution=config["SIMULATION"]["AMPLICON_DISTRIBUTION"],
        amplicon_pseudocounts=config["SIMULATION"]["AMPLICON_PSEUDOCOUNTS"]
    log:    
        "logs/{exp}/simulation/SWAMPy_{sample}.out.log"
    shell:
        "python bin/SWAMPy/src/simulate_metagenome.py --genomes_file {input.genome} "
        "--genome_abundances experiments/{wildcards.exp}/simulation/abundances/{wildcards.exp}_{wildcards.sample}.tsv " 
        "--output_folder experiments/{params.exp}/data --output_filename_prefix {wildcards.exp}_{wildcards.sample} "
        "--n_reads {params.n_reads} --seqSys {params.seqSys} --read_length {params.read_length} --primer_set {params.primer_set} "
        "--amplicon_distribution {params.amplicon_distribution} --amplicon_pseudocounts {params.amplicon_pseudocounts} "
        "--temp_folder bin/SWAMPy/temp/{wildcards.sample} --autoremove &>{log} "
        "&& rm experiments/{wildcards.exp}/data/{wildcards.exp}_{wildcards.sample}_amplicon_abundances_summary.tsv "
        "experiments/{wildcards.exp}/data/{wildcards.exp}_{wildcards.sample}_PCR_errors.vcf "
        "experiments/{wildcards.exp}/data/{wildcards.exp}_{wildcards.sample}.log"

rule qualityControl:
    input:
        seq_files=get_sequences_file_names
    output:
        "experiments/{exp}/simulation/QualityControl/{exp}_multiqc_report.html"
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input.seq_files} --quiet -o experiments/{wildcards.exp}/simulation/QualityControl && "
        "multiqc experiments/{wildcards.exp}/simulation/QualityControl/. --outdir experiments/{wildcards.exp}/simulation/QualityControl/ "
        "--quiet --no-data-dir --filename {wildcards.exp}_multiqc_report.html"