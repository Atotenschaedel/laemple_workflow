# snakemake workflow
include: "common.smk"

# get parameters from config or common.smk 
merged_config = {**default, **config}

TOOL_NAME="lollipop_v0.3.0"

# Run pipeline solo:
# snakemake --snakefile rules/subworkflow_lollipop_v0.3.0/Snakefile.smk --use-conda -c1"

# get sample metadata
def get_sample_date(EXPERIMENT_NAME, sample):
    with open("experiments/"+ EXPERIMENT_NAME + "/simulation/" + EXPERIMENT_NAME + "_metadata.tsv") as f:
        lines = f.readlines()
        for l in lines:
            if l.split("\t")[0] == sample:
                return l.split("\t")[1].strip()

rule _all:
    input:
        expand("experiments/{exp}/results/postPrediction/lollipop_v0.3.0_summary.csv", 
        exp=merged_config["EXPERIMENT_NAME"])

if merged_config["CREATE_REFSET"]:
    rule _mutationList:
        input:
            gff3="reference/gffs/Genes_NC_045512.2.GFF3"
        output: 
            mutlist="reference/Lollipop_reference_set/phe-genomics/phe2cojac/mutlist.tsv",
            pangovars="reference/Lollipop_reference_set/phe-genomics/phe2cojac/variants_pangolin.yaml"
        conda:
            "../../envs/lollipop_v0.3.0.yaml"
        params:
            ref_dir="reference/Lollipop_reference_set/phe-genomics",
            git_repo="https://github.com/phe-genomics/variant_definitions.git"
        shell:
            # fetch the repository of standardised variant definitions
            'mkdir {params.ref_dir} && cd {params.ref_dir} && git clone {params.git_repo} && cd ../.. && '

            # generate a YAML for lineages using the corresponding standardised variant definitions in cojac format
            'for filename in {params.ref_dir}/variant_definitions/variant_yaml/*.yml; do '
            'python bin/lollipop/phe2cojac.py {params.ref_dir}/variant_definitions/variant_yaml/"${filename##*/}" '
            '-y {params.ref_dir}/phe2cojac/"${filename##*/}"; done && ' 

            # create mutation list
            'lollipop generate-mutlist --genes {input.gff3} --output {output.mutlist} --out-pangovars {output.pangovars} -- {params.ref_dir}/phe2cojac/*.yml'

else:
    rule _dummy:
        output: 
            mutlist="reference/Lollipop_reference_set/phe-genomics/phe2cojac/mutlist.tsv",
            pangovars="reference/Lollipop_reference_set/phe-genomics/phe2cojac/variants_pangolin.yaml"

rule _mainSingleSample:
    input:
        bam="experiments/{exp}/results/variantCall/05_softclipPrimer/{sample}_sorted_viral_rh_trim.bam",
        mutlist="reference/Lollipop_reference_set/phe-genomics/phe2cojac/mutlist.tsv",
    output:
        basecnt=temp("experiments/{exp}/results/lollipop_v0.3.0/{sample}.basecnt.tsv.gz"),
        cov=temp("experiments/{exp}/results/lollipop_v0.3.0/{sample}.coverage.tsv.gz"),
        mut="experiments/{exp}/results/lollipop_v0.3.0/{sample}.mut.tsv"
    conda:
        "../../envs/lollipop_v0.3.0.yaml"
    params:
        sample_date=lambda wildcards: get_sample_date(merged_config["EXPERIMENT_NAME"], f"{wildcards.sample}"),
        location=merged_config["LOCATION_NAME"]
    log:
        "logs/{exp}/lollipop_v0.3.0/{sample}.out"
    shell:
        # search mutation in a single sample - create basecount table
        "aln2basecnt --basecnt {output.basecnt} --first 1 --coverage {output.cov} --name \"{wildcards.sample}\" {input.bam} &>{log} && "
        "lollipop getmutations from-basecount --based 1 --output {output.mut} "
        "--location \"{params.location}\" --date \"{params.sample_date}\" -m {input.mutlist} -- {output.basecnt} >>{log} 2>&1"

rule _mainCombinedSeries:
    input:
        pangovars="reference/Lollipop_reference_set/phe-genomics/phe2cojac/variants_pangolin.yaml",
        decon_config="bin/lollipop/" + merged_config["KERNEL_DECONVOLUTION_CONFIG"],
        s_mut= lambda wildcards: expand("experiments/{exp}/results/lollipop_v0.3.0/{sample}.mut.tsv", 
        exp=merged_config["EXPERIMENT_NAME"], sample=merged_config["SAMPLES"])
    output:
        mutgz=temp("experiments/{exp}/results/lollipop_v0.3.0/tallymut.tsv.gz"),
        decon="experiments/{exp}/results/lollipop_v0.3.0/deconvoluted.tsv"
    params: 
        seed=merged_config["APP"]["SEED"]
    conda:
        "../../envs/lollipop_v0.3.0.yaml"
    log:
        "logs/{exp}/lollipop_v0.3.0/deconvolution.out"
    shell:
        # combine time series and run deconvolution
        "xsv cat rows {input.s_mut} | xsv fmt --out-delimiter '\\t' | gzip -f >{output.mutgz} && "
        "lollipop deconvolute --output {output.decon} --var {input.pangovars} --seed={params.seed} --dec {input.decon_config} -- {output.mutgz} &>{log}"

checkpoint _prepareSummary:
    input:
        decon="experiments/{exp}/results/lollipop_v0.3.0/deconvoluted.tsv",
        meta="experiments/{exp}/simulation/{exp}_metadata.tsv"
    output:
        "experiments/{exp}/results/postPrediction/lollipop_v0.3.0_summary.csv"
    conda:
        "../../envs/python3.yaml"
    params:
        min_threshold = merged_config["POSTPRED"]["LINEAGE_MIN_THRESHOLD"]
    log:
        "logs/{exp}/lollipop_v0.3.0/PostPred_summary.csv"
    shell:
        "python rules/subworkflow_lollipop_v0.3.0/prepareSummary.py -i {input.decon} "
        "--metafile {input.meta} --name lollipop_v0.3.0 "
        "--minThreshold {params.min_threshold} -o {output}"