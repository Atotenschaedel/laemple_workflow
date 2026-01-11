# snakemake workflow
include: "common.smk"

# get parameters from config or common.smk 
merged_config = {**default, **config}

TOOL_NAME="LCS_vcfc2f02"

# Run pipeline solo:
# snakemake --snakefile rules/subworkflow_LCS_vcfc2f02/Snakefile.smk --use-conda -c1"

rule _all:
    input:
        expand("experiments/{exp}/results/postPrediction/LCS_vcfc2f02_summary.csv", 
        exp=merged_config["EXPERIMENT_NAME"])

rule _folderPrep:
    output:
        marker_table=temp("experiments/{exp}/results/LCS_vcfc2f02/outputs/variants_table/pango-markers-table.tsv"),
        config=temp("experiments/{exp}/results/LCS_vcfc2f02/rules/config.py")
    params:
        pango_designation = "experiments/{exp}/results/LCS_vcfc2f02/data/pre-generated-marker-tables/pango-designation-markers-v1.2.60.tsv.gz",
        samples = merged_config["SAMPLES"]
    conda:
        "../../envs/lcs_vcfc2f02.yaml"
    shell:
        """
        mkdir -p experiments/{wildcards.exp}/results/LCS_vcfc2f02/data/fastq/ && 
        mkdir -p experiments/{wildcards.exp}/results/LCS_vcfc2f02/outputs/variants_table && 
        cp -rp bin/LCS_vcfc2f02/LCS/* experiments/{wildcards.exp}/results/LCS_vcfc2f02/ && 
        cp -p experiments/{wildcards.exp}/data/*.fastq experiments/{wildcards.exp}/results/LCS_vcfc2f02/data/ && 
        for f in experiments/{wildcards.exp}/results/LCS_vcfc2f02/data/*.fastq; do mv "$f" "$(echo "$f" | sed s/_R/_/)"; done && 
        gzip -f experiments/{wildcards.exp}/results/LCS_vcfc2f02/data/*.fastq && 
        mv experiments/{wildcards.exp}/results/LCS_vcfc2f02/data/*.gz experiments/{wildcards.exp}/results/LCS_vcfc2f02/data/fastq/ && 
        zcat {params.pango_designation} > {output.marker_table} && 
        touch  experiments/{wildcards.exp}/results/LCS_vcfc2f02/data/tags_pool_mypool && 
        for s in {params.samples}; do echo "$s" >> experiments/{wildcards.exp}/results/LCS_vcfc2f02/data/tags_pool_mypool; done
        """

rule _setPrimer:
    input:
        "experiments/{exp}/results/LCS_vcfc2f02/rules/config.py"
    output:
        temp("experiments/{exp}/results/LCS_vcfc2f02/{exp}_primer_set.txt")
    params:
        primer = merged_config["PRIMER"][0],
        primer_path =  merged_config["PRIMER"][1]
    run:   
        if params.primer == "articV3":
            shell("""
            cp {params.primer_path} experiments/{wildcards.exp}/results/LCS_vcfc2f02/data/ && 
            sed -i 's/artic_v3_primers.fa/artic_v3_primers.fasta/g' {input} && 
            touch {output}
            """)
        elif params.primer == "articV4":
            shell("""
            cp {params.primer_path} experiments/{wildcards.exp}/results/LCS_vcfc2f02/data/ && 
            sed -i 's/artic_v3_primers.fa/artic_v4_primers.fasta/g' {input} && 
            touch {output}
            """)
        elif params.primer == "articV5":
            shell("""
            cp {params.primer_path} experiments/{wildcards.exp}/results/LCS_vcfc2f02/data/ && 
            sed -i 's/artic_v3_primers.fa/artic_v5.3_primers.fasta/g' {input} && 
            touch {output}
            """)
        else:
            print(params.primer)

rule _main:
    input:
        primer="experiments/{exp}/results/LCS_vcfc2f02/{exp}_primer_set.txt",
        marker_tsv="experiments/{exp}/results/LCS_vcfc2f02/outputs/variants_table/pango-markers-table.tsv"
    output:
        "experiments/{exp}/results/LCS_vcfc2f02/{exp}.out"
    params:
        markers = "pango",
        cores = 2,
        primer_trimming = "True"
    conda:
        "../../envs/lcs_vcfc2f02.yaml"
    log: 
        "logs/{exp}/LCS_vcfc2f02/snakemake.out"
    shell:
        """
        cd experiments/{wildcards.exp}/results/LCS_vcfc2f02/ && 
        snakemake --config markers={params.markers} dataset=mypool primer_trimming={params.primer_trimming} --cores {params.cores} --rerun-incomplete >../../../../{log} 2>&1 && 
        cp outputs/decompose/mypool.out . && mv mypool.out {wildcards.exp}.out && cp outputs/variants_table/pango-markers-table.tsv . && 
        rm -r */ .snakemake/ Snakefile *.yaml *.tsv *.md *.txt && cd ../../../../
        """

checkpoint _prepareSummary:
    input:
        "experiments/{exp}/results/LCS_vcfc2f02/{exp}.out"
    output:
        "experiments/{exp}/results/postPrediction/LCS_vcfc2f02_summary.csv"
    conda:
        "../../envs/python3.yaml"
    params:
        min_threshold = merged_config["POSTPRED"]["LINEAGE_MIN_THRESHOLD"]
    log:
        "logs/{exp}/LCS_vcfc2f02/PostPred_summary.csv"
    shell:
        "python rules/subworkflow_LCS_vcfc2f02/prepareSummary.py -i {input} "
        "--name LCS_vcfc2f02 --minThreshold {params.min_threshold} -o {output}"