# snakemake workflow
include: "common.smk"

# get parameters from config or common.smk 
merged_config = {**default, **config}

# Run pipeline solo:
# snakemake --snakefile rules/subworkflow_freyja_v1.4.3/Snakefile.smk --use-conda -c1

rule _all:
    input:
        expand("experiments/{exp}/results/postPrediction/freyja_v1.4.3_summary.csv", 
        exp=merged_config["EXPERIMENT_NAME"])

rule _generateInput: # generate depth file, unzip vcf, variants and demix
    input:
        bam="experiments/{exp}/results/variantCall/07_clipOverlap/{sample}_ds_clip_viral.bam",
        ref=merged_config["REFERENCES"]["REF_GENOME"]
    output:
        depth="experiments/{exp}/results/freyja_v1.4.3/{sample}.depth"
    conda: 
        f"../../envs/freyja_v1.4.3.yaml"
    log:
        "logs/{exp}/freyja_v1.4.3/{sample}.input.out"
    shell:  
        "samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f {input.ref} {input.bam} | cut -f1-4 >{output.depth} &>{log}"

rule _main: # start main freyja step
    input:
        vcf="experiments/{exp}/results/variantCall/11_lofreqrefilt/{sample}_samp.lofreq_filtered.vcf.gz",
        bam="experiments/{exp}/results/variantCall/07_clipOverlap/{sample}_ds_clip_viral.bam",
        depth="experiments/{exp}/results/freyja_v1.4.3/{sample}.depth",
        ref=merged_config["REFERENCES"]["REF_GENOME"]
    output:
        tsv="experiments/{exp}/results/freyja_v1.4.3/{sample}.tsv"
    conda: 
        f"../../envs/freyja_v1.4.3.yaml"
    log:
        "logs/{exp}/freyja_v1.4.3/{sample}.demix.out"
    shell: 
        "gunzip -kf {input.vcf} && "
        "freyja variants --ref {input.ref} "
        "--variants experiments/{wildcards.exp}/results/variantCall/11_lofreqrefilt/{wildcards.sample}_samp.lofreq_filtered.vcf "
        "--depths {input.depth} {input.bam} &>{log} && " 
        "freyja demix experiments/{wildcards.exp}/results/variantCall/11_lofreqrefilt/{wildcards.sample}_samp.lofreq_filtered.vcf " 
        "{input.depth} --output {output.tsv} >>{log} 2>&1"

checkpoint _prepareSummary:
    input:
        expand("experiments/{exp}/results/freyja_v1.4.3/{sample}.tsv", 
        exp=merged_config["EXPERIMENT_NAME"], 
        sample=merged_config["SAMPLES"])
    output:
        "experiments/{exp}/results/postPrediction/freyja_v1.4.3_summary.csv"
    conda:
        "../../envs/freyja_v1.4.3.yaml"
    params:
        min_threshold = merged_config["POSTPRED"]["LINEAGE_MIN_THRESHOLD"]
    log:
        "logs/{exp}/freyja_v1.4.3/PostPred_summary.csv"
    shell:
        "python rules/subworkflow_freyja_v1.4.3/prepareSummary.py -i {input} "
        "-n freyja_v1.4.3 -m {params.min_threshold} -o {output}"