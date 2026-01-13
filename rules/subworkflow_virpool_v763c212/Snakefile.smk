"""
Deprecated: in order for this script to run, additional code and configs need to be added for tool: virpool_v763c212. see README for additional information.
virpool can be installed in bin folder following instructions by the orginal authors: https://github.com/fmfi-compbio/virpool, script is based on commit ID: 763c212.
"""
# snakemake workflow
include: "common.smk"

# get parameters from config or common.smk 
merged_config = {**default, **config}

TOOL_NAME="virpool_v763c212"

# Run pipeline solo:
# snakemake --snakefile rules/subworkflow_virpool_v763c212/Snakefile.smk --use-conda -c1"

rule _all:
    input:
        expand("experiments/{exp}/results/postPrediction/virpool_v763c212_summary.csv", 
        exp=merged_config["EXPERIMENT_NAME"])

rule _main:
    input:
        ref=merged_config["REFERENCES"]["REF_GENOME"],
        bam="experiments/{exp}/results/variantCall/07_clipOverlap/{sample}_ds_clip_viral.bam"
    output:
        yaml="experiments/{exp}/results/virpool_v763c212/{sample}.yaml"
    params:
        variants="bin/virpool_v763c212/virpool/profiles/"+ merged_config["VARIANT_PROFILE"] + ".tsv",
        outdir="experiments/{exp}/results/virpool_v763c212/{sample}"
    conda:
        "../../envs/virpool_v763c212.yaml"
    log:
        "logs/{exp}/virpool_v763c212/{sample}.out"
    shell:
        "python3 bin/virpool_v763c212/virpool/src/virpool --variants {params.variants} "
        "--genome {input.ref} --output-dir {params.outdir} {input.bam} && "
        "mv experiments/{wildcards.exp}/results/virpool_v763c212/{wildcards.sample}/estimated_weights.yaml {output.yaml} && "
        "rm -r {params.outdir}"

checkpoint _prepareSummary:
    input:
        expand("experiments/{exp}/results/virpool_v763c212/{sample}.yaml",
        exp=merged_config["EXPERIMENT_NAME"], sample=merged_config["SAMPLES"])
    output:
        "experiments/{exp}/results/postPrediction/virpool_v763c212_summary.csv"
    conda:
        "../../envs/python3.yaml"
    params:
        min_threshold = merged_config["POSTPRED"]["LINEAGE_MIN_THRESHOLD"]
    log:
        "logs/{exp}/virpool_v763c212/PostPred_summary.csv"
    shell:
        "python rules/subworkflow_virpool_v763c212/prepareSummary.py -i {input} "
        "--name virpool_v763c212 --minThreshold {params.min_threshold} -o {output}"