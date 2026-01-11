# -*- coding: utf-8 -*-
# snakemake workflow
"""
Created on 30-06-2023
@author: aschedl
Description: This workflow uses output files from variant calling pipeline
    and calls different lineage deconvolution tools based on config file.

Example usage:
snakemake --snakefile lineage_deconvultion.smk --use-conda -c1
"""
# import modules
from pathlib import Path

# get parameters from config file
configfile: "config/workflow_config.yaml"

# get parameters from config file
EXPERIMENT_NAME = config["EXPERIMENT_NAME"]
SET_NAME=config["CURRENT_SET"]

p = Path.cwd() / "experiments" / EXPERIMENT_NAME / "data"
SAMPLES = list(set([x.stem[:-3] for x in p.glob("*") if x.is_file()]))
config["SAMPLES"] = SAMPLES
 
MODULE_NAMES = []
for tool in config["TOOLS"]:
    if config["TOOLS"][tool]["INCLUDE_IN_ANALYSIS"]:
        MODULE_NAMES.append(config["TOOLS"][tool]["TOOL_NAME"])

submodule_all_outputs = [f"experiments/{EXPERIMENT_NAME}/results/postPrediction/simulation_summary.csv"]

for module_name in MODULE_NAMES:
    module:
        name: module_name
        snakefile: f"subworkflow_{module_name}/Snakefile.smk"
        config: config
    
    use rule * from module_name as module_name*
    submodule_all_outputs.append(f"experiments/{EXPERIMENT_NAME}/results/postPrediction/{module_name}_summary.csv")

for rule in workflow.rules:
    print(rule)

rule all:
    input: 
        submodule_all_outputs
    default_target: True

rule pango_sequences_refSet:
    output:
        "reference/pango-sequences/data/pango-consensus-sequences_summary.json"
    params:
        git_repo = "https://github.com/corneliusroemer/pango-sequences.git"
    shell:
        "cd references && git clone {params.git_repo} && cd .. "

rule simulation_summary:
    input: 
        abundances=expand("experiments/{exp}/simulation/abundances/{sample}.tsv", exp=EXPERIMENT_NAME, sample=SAMPLES),
        coverage=expand("experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}_cov.tsv", exp=EXPERIMENT_NAME, sample=SAMPLES),
        stats=expand("experiments/{exp}/results/variantCall/00_stats/{sample}/{sample}_virus.stats", exp=EXPERIMENT_NAME, sample=SAMPLES),
        meta="experiments/{exp}/simulation/{exp}_metadata.tsv"
    output: 
        "experiments/{exp}/results/postPrediction/simulation_summary.csv"
    conda: 
        "../envs/python3.yaml"
    params:
        use_real_timecourse=config["SIMULATION"]["DATASET"][SET_NAME]["USE_REAL_TIMECOURSE"],
        min_read_count=config["POSTPRED"]["MIN_READ_COUNT"],
        lineage_min_threshold = config["POSTPRED"]["LINEAGE_MIN_THRESHOLD"]
    shell:
        "python bin/custom_scripts/SequencingSummary.py "
        "--abundances_files {input.abundances} "
        "--coverage_files {input.coverage} "
        "--stat_files {input.stats} "
        "--meta_file {input.meta} "
        "--real_timecourse {params.use_real_timecourse} "
        "--min_read_count {params.min_read_count} "
        "--lineage_min_threshold {params.lineage_min_threshold} "
        "-o {output}"