# snakemake workflow
include: "common.smk"

# get parameters from config or common.smk 
merged_config = {**default, **config}

TOOL_NAME="VLQ_vd04ad36"

# Run pipeline solo:
# snakemake --snakefile rules/subworkflow_VLQ_vd04ad36/Snakefile.smk --use-conda -c1"

rule _all:
    input:
        expand("experiments/{exp}/results/postPrediction/VLQ_vd04ad36_summary.csv", 
        exp=merged_config["EXPERIMENT_NAME"])

if merged_config["CREATE_REFSET"]:
    rule _generateRefset:
        input:
            metadata=merged_config["GISAID_METADATA"],
            fasta=merged_config["GISAID_SEQUENCES"]
        output:
            metadata="reference/VLQ_reference_set/metadata.tsv",
            fasta="reference/VLQ_reference_set/sequences.fasta",
            kallisto_idx="reference/VLQ_reference_set/sequences.kallisto_idx"
        params:
            ref=merged_config["REFERENCES"]["REF_GENOME"],
            seed = merged_config["APP"]["SEED"],
            continent="Europe",
            ref_set="reference/VLQ_reference_set"
        conda:
            "../../envs/VLQ_vd04ad36.yaml"
        log:
            expand("logs/{exp}/VLQ/generateRefSet.out", exp=merged_config["EXPERIMENT_NAME"])
        shell: # build a reference set --results in sequences.fasta and metadata.tsv
            """
            # cleanUp GISAID download, create UTF-8 encoded text-files
    
            unxz -keep {input.metadata}.tar.xz && mv {input.metadata}.tar {input.metadata}.tsv && 
            cat {input.metadata}.tsv | perl -F"\t" -lane 'if($.==1 || $F[5] =~ m/20\d\d-\d\d-\d\d/) {{ print }}' > {input.metadata}.tsv && 

            # create symbolic link & rename GISAID download files

            mkdir {params.ref_set}/GISAID_dump {params.ref_set}/lineages &&
            ln -s {input.metadata} {params.ref_set}/GISAID_dump/metadata.tsv &&
            ln -s {input.fasta} {params.ref_set}/GISAID_dump/sequences.fasta &&

            # preprocessing reference sequencing data

            python bin/VLQ/wastewater_analysis/pipeline/preprocess_references.py -m {params.ref_set}/GISAID_dump/metadata.tsv 
            -f {params.ref_set}/GISAID_dump/sequences.fasta -k 1000 --seed {params.seed} --continent {params.continent} -o {params.ref_set}/lineages/ >>{log} && 

            # call variants compared to orginal sars cov-2 and compute allele frequencies per lineage
            # bash script does not recognize symbolic links, need orginal ABSOLUTE path
            bash bin/VLQ/wastewater_analysis/pipeline/call_variants.sh {params.ref_set}/lineages {params.ref} >>{log} &&

            # select sequences per lineage with all mutation with al freq >50% were captured at least once
 
            python bin/VLQ/wastewater_analysis/pipeline/select_samples.py -m {params.ref_set}/GISAID_dump/metadata.tsv -f {params.ref_set}/GISAID_dump/sequences.fasta 
            -o {params.ref_set}/ --vcf {params.ref_set}/lineages/*_merged.vcf.gz --freq {params.ref_set}/lineages/*_merged.frq >>{log} &&

            # index reference set
            kallisto index -i {output.kallisto_idx} {output.fasta} >>{log} && 

            # clean-up

            rm {params.ref_set}/lineages/*_merged.vcf.gz {params.ref_set}/lineages/*_merged.frq {params.ref_set}/lineages/*_merged.frq {params.ref_set}/lineages/*_merged.log {params.ref_set}/lineages/*_merged.sites.pi &&
            rm -rf {params.ref_set}/lineages/*
            """

else:
    rule _dummy:
        output:
            metadata="reference/VLQ_reference_set/metadata.tsv",
            kallisto_idx="reference/VLQ_reference_set/sequences.kallisto_idx"

rule _main:
    input:
        kallisto_idx="reference/VLQ_reference_set/sequences.kallisto_idx",
        metadata="reference/VLQ_reference_set/metadata.tsv",
        fastq1="experiments/{exp}/results/variantCall/02_correctOverlapPE/{sample}_ecco_1.fq.gz",
        fastq2="experiments/{exp}/results/variantCall/02_correctOverlapPE/{sample}_ecco_2.fq.gz"
    output:
        folder = temp(directory("experiments/{exp}/results/VLQ_vd04ad36/{sample}/")),
        prediction="experiments/{exp}/results/VLQ_vd04ad36/{sample}_predictions.tsv"
    params:
        threads = 20
    conda:
        "../../envs/VLQ_vd04ad36.yaml"
    log:
        "logs/{exp}/VLQ_vd04ad36/{sample}.out"
    shell:  # calculate abundance per reference sequence and obtain abundances
        "kallisto quant -t {params.threads} -i {input.kallisto_idx} "
        "-o experiments/{wildcards.exp}/results/VLQ_vd04ad36/{wildcards.sample} "
        " {input.fastq1} {input.fastq2} &>{log} && "

        "python bin/VLQ_vd04ad36/VLQ/wastewater_analysis/pipeline/output_abundances.py "
        "-o {output.prediction} --metadata {input.metadata} "
        "experiments/{wildcards.exp}/results/VLQ_vd04ad36/{wildcards.sample}/abundance.tsv >>{log} 2>&1"

checkpoint _prepareSummary:
    input:
        expand("experiments/{exp}/results/VLQ_vd04ad36/{sample}_predictions.tsv",
        exp=merged_config["EXPERIMENT_NAME"], sample=merged_config["SAMPLES"])
    output:
        "experiments/{exp}/results/postPrediction/VLQ_vd04ad36_summary.csv"
    conda:
        "../../envs/python3.yaml"
    params:
        min_threshold = merged_config["POSTPRED"]["LINEAGE_MIN_THRESHOLD"]
    log:
        "logs/{exp}/VLQ_vd04ad36/PostPred_summary.csv"
    shell:
        "python rules/subworkflow_VLQ_vd04ad36/prepareSummary.py -i {input} "
        "--name VLQ_vd04ad36 --minThreshold {params.min_threshold} -o {output}"