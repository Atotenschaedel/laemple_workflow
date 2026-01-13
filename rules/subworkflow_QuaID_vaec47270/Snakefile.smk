"""
Deprecated: in order for this script to run, additional code and configs need to be added for tool: QuaID_vaec47270. see README for additional information.
QuaID can be installed in bin folder following instructions by the orginal authors: https://gitlab.com/treangenlab/quaid, script is based on commit ID: aec47270.
"""
# snakemake workflow
include: "common.smk"

# get parameters from config or common.smk 
merged_config = {**default, **config}

TOOL_NAME="QuaID_vaec47270"

# Run pipeline solo:
# snakemake --snakefile rules/subworkflow_QuaID_vaec47270/Snakefile.smk --use-conda -c1"

# get sample metadata
def get_sample_date(EXPERIMENT_NAME, sample):
    with open("experiments/"+ EXPERIMENT_NAME + "/simulation/" + EXPERIMENT_NAME + "_metadata.tsv") as f:
        lines = f.readlines()
        for l in lines:
            if l.split("\t")[0] == sample:
                return l.split("\t")[1].strip()

rule _all:
    input:
        expand("experiments/{exp}/results/postPrediction/QuaID_vaec47270_summary.csv", 
        exp=merged_config["EXPERIMENT_NAME"])

if merged_config["fix_date_bug"]: # fixes a bug in source code that prevents loading of file for variant call pipeline
    rule _fixDate:
        input:
            expand("experiments/{exp}/simulation/{exp}_metadata.tsv", exp=merged_config["EXPERIMENT_NAME"])
        output:
            "experiments/{exp}/results/QuaID_vaec47270/QuaID_2022/{exp}_metadata_2022.tsv"
        params:
            out_dir = "experiments/{exp}/results/QuaID_vaec47270/QuaID_2022/"
        shell:
            "mkdir -p {params.out_dir} && cp {input} {params.out_dir} && "
            "mv {params.out_dir}{wildcards.exp}_metadata.tsv {output} && "
            "sed -i -e 's/2023/2022/g' -e 's/2024/2023/g' {output}"
        
    rule _main:
        input:
            txt = "experiments/{exp}/results/QuaID_vaec47270/run/Variant-calling-combined/{sample}_finished.txt",
            vdb_ref = "reference/QuaID_reference_set/vdb_lineage_df_week.csv",
            metadata = "experiments/{exp}/results/QuaID_vaec47270/QuaID_2022/{exp}_metadata_2022.tsv"
        output:
            "experiments/{exp}/results/QuaID_vaec47270/QuaID_2022/{sample}.csv"
        params:
            sample_date=lambda wildcards: get_sample_date(merged_config["EXPERIMENT_NAME"], f"{wildcards.sample}"),
            occ_cutoff= merged_config["OCCURENCE_CUTOFF"],
            incl_cutoff= merged_config["INCLUSION_CUTOFF"],
            exc_cutoff= merged_config["EXCLUSION_CUTOFF"],
            time_wind= merged_config["TIME_WINDOW"],
            level= merged_config["LEVEL"]
        conda:
            "../../envs/quaid_vaec47270.yaml"
        log:
            "logs/{exp}/QuaID_vaec47270/quaid_2022/{sample}_quaid.out"
        shell:
            "cd experiments/{wildcards.exp}/results/QuaID_vaec47270/run && mkdir -p output && "
            "python quaid.py ../../../../../{input.vdb_ref} ../../../../../{input.metadata} Variant-calling-combined Coverage "
            "--occurence-cutoff {params.occ_cutoff} --inclusion-cutoff {params.incl_cutoff} --exclusion-cutoff {params.exc_cutoff} "
            "--time-window {params.time_wind} --level {params.level} --verbose >../../../../../{log} 2>&1 && "
            "cd ../../../../../ && "
            "mv experiments/{wildcards.exp}/results/QuaID_vaec47270/run/output/{params.sample_date}/Summary-{params.sample_date}-o{params.occ_cutoff}-i{params.incl_cutoff}-e{params.exc_cutoff}-t{params.time_wind}-l{params.level}-QU-noX-mask.csv "
            "experiments/{wildcards.exp}/results/QuaID_vaec47270/{wildcards.sample}.csv >>../../../../../{log}"

    rule _clean:
        input:
            expand("experiments/{exp}/results/QuaID_vaec47270/QuaID_2022/{sample}.csv", exp=EXPERIMENT_NAME, sample=SAMPLES)
        output:
            temp(expand("experiments/{exp}/results/QuaID_vaec47270/QuaID_2022/cleaned.txt", exp=EXPERIMENT_NAME))
        params:
            exp = EXPERIMENT_NAME
        shell:
            "rm -rf experiments/{params.exp}/results/QuaID_vaec47270/QuaID_2022/run && "
            "touch experiments/{params.exp}/results/QuaID_vaec47270/QuaID_2022/cleaned.txt"

if merged_config["translate_package_list"]:
    rule _envPrep: # translates package list from source repro to conda readable package recipe
        output:
            "../../envs/quaid_vaec47270.yaml"
        params:
            req_file = Path.cwd() /"bin" / "QuaID_vaec47270" / "quaid" / "quaid-env-list.linux64.txt"
        run:
            with open(req_file, 'r') as f:
                lines = f.readlines()
                df = pd.read_csv(f, sep=" ", names=["Name", "Version", "Build", "Channel"], header=None, skiprows=6)
                
            channels = []
            dependencies = []

            for line in lines:
                if line.startswith('#'): pass
                else:
                    line = line.split(" ")
                    line = list(filter(None, line))
                    new_line = line[0] + "=" + line[1] #+ "=" + line[3].replace("\n","")

                    channels.append(line[3].replace("\n",""))
                    dependencies.append(new_line)

            channels=list(set(channels))

            data = {
                'name': "quaid-env",
                'channels': channels,
                'dependencies': dependencies
            }

            with open("../../envs/quaid_vaec47270.yaml", "w+") as yaml_file:
                yaml.dump(data, yaml_file, default_flow_style=False)

rule _folderPrep: # prepares required folder structure 
    output:
        "experiments/{exp}/results/QuaID_vaec47270/run/SARS-CoV-2-reference.fasta.bwt"
    params:
        ref = "SARS-CoV-2-reference.fasta"
    shell:
        "mkdir -p experiments/{wildcards.exp}/results/QuaID_vaec47270/run/QC-reports "
        "experiments/{wildcards.exp}/results/QuaID_vaec47270/run/Coverage "
        "experiments/{wildcards.exp}/results/QuaID_vaec47270/run/Processed-reads "
        "experiments/{wildcards.exp}/results/QuaID_vaec47270/run/Read-mapping "
        "experiments/{wildcards.exp}/results/QuaID_vaec47270/run/Variant-calling-iVar "
        "experiments/{wildcards.exp}/results/QuaID_vaec47270/run/Variant-calling-LoFreq "
        "experiments/{wildcards.exp}/results/QuaID_vaec47270/run/Variant-calling-combined "
        "experiments/{wildcards.exp}/results/QuaID_vaec47270/run/MSA_precomputed && "

        "cp -r bin/QuaID_vaec47270/quaid/* experiments/{wildcards.exp}/results/QuaID_vaec47270/run/ && "

        "cd experiments/{wildcards.exp}/results/QuaID_vaec47270/run/ && "
        "bwa index {params.ref} && cd ../../../../../"

rule _variantCall:
    input:
        fastq1= "experiments/{exp}/data/{sample}_R1.fastq",
        fastq2= "experiments/{exp}/data/{sample}_R2.fastq",
        idx = lambda wildcards: expand("experiments/{exp}/results/QuaID_vaec47270/run/SARS-CoV-2-reference.fasta.bwt", 
        exp=merged_config["EXPERIMENT_NAME"])
    output:
        #tsv = "experiments/{exp}/results/QuaID/run/Variant-calling-combined/{params.sample_date}/{sample}-merged-indel.clean.tsv"
        temp("experiments/{exp}/results/QuaID_vaec47270/run/Variant-calling-combined/{sample}_finished.txt")
    conda:
        "../../envs/quaid_vaec47270.yaml"
    params:
        sample_date=lambda wildcards: get_sample_date(merged_config["EXPERIMENT_NAME"], f"{wildcards.sample}")
    log:
        "logs/{exp}/QuaID_vaec47270/{sample}_varCall.out"
    shell:
        "mkdir -p experiments/{wildcards.exp}/results/QuaID_vaec47270/run/Data/{params.sample_date}/{wildcards.sample} && "
        "cp {input.fastq1} experiments/{wildcards.exp}/results/QuaID_vaec47270/run/Data/{params.sample_date}/{wildcards.sample}/ && "
        "cp {input.fastq2} experiments/{wildcards.exp}/results/QuaID_vaec47270/run/Data/{params.sample_date}/{wildcards.sample}/ && "
        "cd experiments/{wildcards.exp}/results/QuaID_vaec47270/run/ && "
        "./run_analyses.sh {params.sample_date} >../../../../../{log} 2>&1 && "
        "touch Variant-calling-combined/{wildcards.sample}_finished.txt && cd ../../../../../"

if merged_config["CREATE_REFSET"]:
    rule _generateRefset:
        input:
            metadata=merged_config["GISAID_METADATA"],
            msa=merged_config["GISAID_MSA"]
        output: 
            "reference/QuaID_reference_set/vdb_lineage_df_week.csv",
            "reference/GISAID/metadata.tsv"
        params:
            ref_dir = merged_config["REF_DIR"],
            msa_name = merged_config["GISAID_MSA_NAME"],
            metadata_name = merged["GISAID_METADATA_NAME"]
        shell:
            "tar -xOf {input.msa} {params.ref_dir}{params.msa_name}/{params.msa_name}.fasta | bin/quaid/vdbCreate -s -N && "
            "bin/QuaID_vaec47270/quaid/vdb -t vdb_msa_nucl.txt vdb_msa_nucl_trimmed.txt && "
            "python bin/QuaID_vaec47270/quaid/vdb-to-df.py vdb_msa_nucl_trimmed.txt {params.metadata}.tsv"

else:
    rule _dummny:
        output:
            "reference/QuaID_reference_set/vdb_lineage_df_week.csv",
            "reference/GISAID/metadata.tsv"

rule _main:
    input:
        txt = "experiments/{exp}/results/QuaID_vaec47270/run/Variant-calling-combined/{sample}_finished.txt",
        vdb_ref = "reference/QuaID_reference_set/vdb_lineage_df_week.csv",
        metadata = "reference/GISAID/metadata.tsv"
    output:
        "experiments/{exp}/results/QuaID_vaec47270/{sample}.csv"
    params:
        sample_date=lambda wildcards: get_sample_date(merged_config["EXPERIMENT_NAME"], f"{wildcards.sample}"),
        occ_cutoff= merged_config["OCCURENCE_CUTOFF"],
        incl_cutoff= merged_config["INCLUSION_CUTOFF"],
        exc_cutoff= merged_config["EXCLUSION_CUTOFF"],
        time_wind= merged_config["TIME_WINDOW"],
        level= merged_config["LEVEL"]
    conda:
        "../../envs/quaid_vaec47270.yaml"
    log:
        "logs/{exp}/QuaID_vaec47270/{sample}_quaid.out"
    shell:
        "cd experiments/{wildcards.exp}/results/QuaID_vaec47270/run && mkdir -p output && "
        "python quaid.py ../../../../../{input.vdb_ref} ../../../../../{input.metadata} Variant-calling-combined Coverage "
        "--occurence-cutoff {params.occ_cutoff} --inclusion-cutoff {params.incl_cutoff} --exclusion-cutoff {params.exc_cutoff} "
        "--time-window {params.time_wind} --level {params.level} --verbose >../../../../../{log} 2>&1 && "
        "cd ../../../../../ && "
        "mv experiments/{wildcards.exp}/results/QuaID_vaec47270/run/output/{params.sample_date}/Summary-{params.sample_date}-o{params.occ_cutoff}-i{params.incl_cutoff}-e{params.exc_cutoff}-t{params.time_wind}-l{params.level}-QU-noX-mask.csv "
        "experiments/{wildcards.exp}/results/QuaID_vaec47270/{wildcards.sample}.csv"

rule _clean:
    input:
        expand("experiments/{exp}/results/QuaID_vaec47270/{sample}.csv", 
        exp=merged_config["EXPERIMENT_NAME"], sample=merged_config["SAMPLES"])
    output:
        temp(expand("experiments/{exp}/results/QuaID_vaec47270/cleaned.txt", exp=merged_config["EXPERIMENT_NAME"]))
    params:
        exp = merged_config["EXPERIMENT_NAME"]
    shell:
        "rm -rf experiments/{params.exp}/results/QuaID_vaec47270/run && "
        "touch experiments/{params.exp}/results/QuaID_vaec47270/cleaned.txt"

checkpoint _prepareSummary:
    input:
        cleaned_flag=expand("experiments/{exp}/results/QuaID_vaec47270/cleaned.txt", 
        exp=merged_config["EXPERIMENT_NAME"]),
        meta=expand("experiments/{exp}/simulation/{exp}_metadata.tsv", 
        exp=merged_config["EXPERIMENT_NAME"]),
        who_file = merged_config["WHO_VOC_DEFINITION_FILE"],
        files=expand("experiments/{exp}/results/QuaID_vaec47270/{sample}.csv",
        exp=merged_config["EXPERIMENT_NAME"], sample=merged_config["SAMPLES"])
    output:
        "experiments/{exp}/results/postPrediction/QuaID_vaec47270_summary.csv"
    conda:
        "../../envs/python3.yaml"
    params:
        min_threshold = merged_config["POSTPRED"]["LINEAGE_MIN_THRESHOLD"]
    log:
        "logs/{exp}/QuaID_vaec47270/PostPred_summary.log"
    shell:
        "python rules/subworkflow_QuaID_vaec47270/prepareSummary.py -i {input.files} "
        "--name QuaID_vaec47270 --metaFile {input.meta} --whoFile {input.who_file} "
        "--minThreshold {params.min_threshold} -o {output}"