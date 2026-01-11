# snakemake workflow
include: "common.smk"

# get parameters from config or common.smk 
merged_config = {**default, **config}

TOOL_NAME="vaquero_v24d9211"

# Run pipeline solo:
# snakemake --snakefile rules/subworkflow_vaquero_v24d9211/Snakefile.smk --use-conda -c1

rule _all:
    input:
        expand("experiments/{exp}/results/postPrediction/vaquero_v24d9211_summary.csv", 
        exp=merged_config["EXPERIMENT_NAME"])

rule _generateInputAF: #generate tidy table from vcf files 
    input: 
        "experiments/{exp}/results/variantCall/12_phasedLofreqVcf/{sample}_samp.lofreq_filtered_phased.vcf.gz"
    output:
        temp("experiments/{exp}/results/vaquero_v24d9211/{sample}_af.tsv")
    params: 
        sample_suffix = merged_config["SAMPLE_SUFFIX"]
    conda: 
        "../../envs/pysam.yaml"
    shell:
        "python3 bin/vaquero_v24d9211/scripts/vcf2tsv_long.py -i {input} -m 0.1 | "
        "sed 's/{wildcards.sample}/{wildcards.sample}{params.sample_suffix}/g' >{output} "

rule _mergeAFs: # merge tidy tables togehter, remove all but first header
    input: 
        expand("experiments/{exp}/results/vaquero_v24d9211/{sample}_af.tsv", 
        exp=merged_config["EXPERIMENT_NAME"], sample=merged_config["SAMPLES"])
    output: 
        expand("experiments/{exp}/results/vaquero_v24d9211/Input_afs.tsv", 
        exp=merged_config["EXPERIMENT_NAME"])
    shell: 
        "cat {input} | sed -e '2,${{ /^SAMPLEID/d }}' >{output}"

rule _generateInputMetaData: # generate metadata for vaquero based on generic metadata file
    input: 
        "experiments/{exp}/simulation/{exp}_metadata.tsv"
    output: 
        "experiments/{exp}/results/vaquero_v24d9211/Input_metadata.tsv"
    conda: 
        "../../envs/python3.yaml"
    params: 
        bsf_run = merged_config["BSF_RUN"],
        bsf_start_date = merged_config["BSF_START_DATE"],
        locationID = merged_config["LOCATION_ID"],
        locationName=merged_config["LOCATION_NAME"],
        N_in_consensus=merged_config["N_IN_CONSENSUS"],
        adress_town=merged_config["ADDRESS_TOWN"],
        connected_people=merged_config["CONNECTED_PEOPLE"],
        dcpLatitude=merged_config["DCP_LATITUDE"],
        dcpLongitude=merged_config["DCP_LONGTITUTDE"],
        include_in_report=merged_config["INCLUDE_IN_REPORT"],
        report_category=merged_config["REPORT_CATEGORY"],
        status=merged_config["STATUS"],
        sample_suffix=merged_config["SAMPLE_SUFFIX"]
    shell: 
        "python3 rules/subworkflow_vaquero_v24d9211/generateMetadata.py --input {input} --output {output} "
        "--run {params.bsf_run} --start_date {params.bsf_start_date} --locationID {params.locationID} "
        "--locationName {params.locationName} --N_in_consensus {params.N_in_consensus} --adress_town {params.adress_town} "
        "--connected_people {params.connected_people} --latitude {params.dcpLatitude} --longitude {params.dcpLongitude} "
        "--include_in_report {params.include_in_report} --report {params.report_category} --status {params.status} --suffix {params.sample_suffix}"

rule _main: # start main vaquero step
    input:
        afs="experiments/{exp}/results/vaquero_v24d9211/Input_afs.tsv",
        meta="experiments/{exp}/results/vaquero_v24d9211/Input_metadata.tsv"
    output: 
        "experiments/{exp}/results/vaquero_v24d9211/globalFittedData.csv"
    params:
        markerDate = merged_config["MARKER_DATE"],
        markerLocation = merged_config["MARKER_LOCATION"]
    conda: 
        "../../envs/vaquero_v24d9211.yaml"
    log: 
        "logs/{exp}/vaquero_v24d9211/vaquero.out"
    shell:
        "Rscript bin/custom_scripts/install_devtools.R &>{log} && "
        "Rscript bin/vaquero_v24d9211/scripts/VaQuERo_v2.r --metadata {input.meta} --data {input.afs} --inputformat tidy "
        "--marker bin/vaquero_v24d9211/resources/mutations_list_grouped_pango_codonPhased_{params.markerDate}_{params.markerLocation}.csv "
        "--smarker bin/vaquero_v24d9211/resources/mutations_special_2022-12-21.csv "
        "--pmarker bin/vaquero_v24d9211/resources/mutations_problematic_vss1_v3.csv "
        "--dir experiments/{wildcards.exp}/results/vaquero_v24d9211/ >>{log} 2>&1" 

checkpoint _prepareSummary:
    input:
        fitted="experiments/{exp}/results/vaquero_v24d9211/globalFittedData.csv",
        meta="experiments/{exp}/simulation/{exp}_metadata.tsv"
    output:
        "experiments/{exp}/results/postPrediction/vaquero_v24d9211_summary.csv"
    conda:
        "../../envs/python3.yaml"
    params:
        min_threshold = merged_config["POSTPRED"]["LINEAGE_MIN_THRESHOLD"]
    log:
        "logs/{exp}/vaquero_v24d9211/PostPred_summary.csv"
    shell:
        "python rules/subworkflow_vaquero_v24d9211/prepareSummary.py -i {input.fitted} "
        "--metafile {input.meta} --name vaquero_v24d9211 "
        "--minThreshold {params.min_threshold} -o {output}"