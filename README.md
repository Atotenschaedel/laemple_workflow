# Readme

# Laemple Workflow

**a Benchmarking Framework for Virus Lineage Deconvolution Tools for SARS-CoV-2 from Wastewater**

This repository contains a reproducible benchmarking pipeline for evaluating SARS-CoV-2 lineage deconvolution tools from wastewater sequencing data.. It is designed to simulate sequencing data, perform variant calling, apply lineage deconvolution methods, and assess their performance against known ground truth. The pipeline is implemented using **Snakemake**, with configuration-driven execution and downstream analysis in **R**.


## Overview

Wastewater sequencing is increasingly used to monitor viral diversity in populations. Estimating the relative abundance of viral lineages from mixed samples requires specialized deconvolution methods. This project provides a standardized framework to:

- Simulate sequencing data with known lineage compositions
- Perform variant calling on simulated or real sequencing reads
- Run and benchmark lineage deconvolution tools
- Compare predicted lineage abundances to ground truth
- Visualize and summarize benchmarking results

The modular design allows individual components of the pipeline to be adapted or extended as new tools and methods become available.


## Project Structure

The repository is organized into modular workflows, configuration files, and supporting resources:

```text
Project
├── bin                 Tool-specific source code, helper scripts and utilities used by workflows      
├── envs                conda environment definitions for reproducibility
├── references          reference data (e.g. genomes, lineage definitions)
├── config
│   └── workflow_config.yaml            central configuration file for all workflows
├── rules
│   ├── simulation.smk                  Snakemake Workflow for Input Data Simulation
│   ├── variantCalling.smk              Snakemake Workflow for Mutation Calling
│   ├── lineage_deconvolution.smk       Snakemake Workflow to call tool specific subworfklows
│   ├── subworkflow_TOOL NAME           tool specific subworkflow
│   │   ├── Snakefile.smk               tool specific Snakemake workflow
│   │   ├── common.smk                  tool specific config file for subworkflow
│   └── └── prepareSummary.py           Python script to convert tool specfic output to standard output format
└──main.py                 main script to start complete project

# additional folders that will be created after workflow completion:

├── logs
├── experiments
│   ├── EXPERIMENT NAME           
│   │   ├── data                        contains simulated sequencing data in fastq format 
│   │   ├── results                     Python script to convert tool
│   │   │   ├── postPrediction          contains standardized output files from tool specific workflows
│   │   │   ├── variantCall             contains all results from mutation calling pipeline
│   │   │   └──TOOL NAME               tool specific output files 
│   │   ├── simulation                  Python script to convert tool
│   │   │   ├── abundances              abundances per samples in tsv format
│   │   │   ├── QualityControl          QC report per sample
│   │   │   ├── EXPERIMENT_NAME_data.csv        simulated abundances data over complete timecourse
│   │   │   ├── EXPERIMENT_NAME_metadata.tsv    simulated metadata per sample 
└── └── └── └── EXPERIMENT_NAME_plot.png        Plot of simulated timecourse      
```
# Installation

This project consist of two separate installation layers:

1. **Overall pipeline installation** (base setup)
2. **Tool-specific installation** (each lineage deconvolution tool)

The pipeline is designed to be modular: new tools can be added without modifiying the core workflow.

## Overall pipeline installation

1. **Clone the Repository**
2. **Install Snakemake** using `conda`
    ```
    conda env create -f envs/snakemake.yaml
    ```
3. **Configuration**
    All pipeline parameters are defined in `config/worfklow_config.yaml`. 
    This file controls datasets, simulation parameters, seeds, enabled tools and reporting behavior. No code changes are required to modify experimental setups. 

## Minimal Example

The default workflow is configured to include a minimal example of the complete workflow which includes only tools Freyja (https://anaconda.org/channels/bioconda/packages/freyja/overview) and VaQuERo (https://github.com/fabou-uobaf/VaQuERo.git) and uses additional reference sequences from https://github.com/corneliusroemer/pango-sequences.git as well as SWAMPy for simulation of sequence data from wastewater samples (https://github.com/goldman-gp-ebi/SWAMPy). It can be started by activating by:
    
    ```
    conda activate snakemake
    python main.py
    ```

It should end with 6 new experiments with each having three result files in the resultfolder: `freyja_v2.0.0_summary.csv`, `vaquero_v24d9211_summary.csv` and `vaquero_v24d9211_summary.csv`. Final Report can be rendered using `PostPredict_plots_all.Rmd`

It is recommended to update reference sequences when working with newer data.

## Tool-specific installation

Each lineage deconvolution tool is installed independently and integrated into the pipeline in a standardized way. Freyja will be installed 

### General rules for adding a new tool

For a tool named `TOOLNAME`, the following components are required:
1. Conda environment
2. Source code
3. Snakemake subworkflow
4. Output standardization script
5. Configuration entry

### Install the tool itself
Follow the official installation instruction of the tools (e.g. GitHub README, documentation).

### Create a Conda environment

Create a conda yaml file named:  `envs/TOOLNAME.yaml`, which contains the tool, all it dependencies and any required Python/R packages

Example:
```
name: TOOLNAME
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.9
  - toolname
  - numpy
```

### Add tool source code

If required, place the tool's source code in: `bin/TOOLNAME`. this may include any wrapper scripts, helper utilities and configuration templates.

### Create a Snakemake Subworkflow

Each tool needs its own subworfklow directory: 

```
rules/subworkflow_TOOLNAME/
├── Snakefile.smk
├── common.smk
└── prepareSummary.py
```

#### Snakefile.smk
This file should execute the tool, handle all input/output paths, call prepareSummary.py and writes standardized output into: `experiments/<EXPERIMENT_NAME>/postPrediction/TOOLNAME_summary.csv

#### common.smk

This file contains tool-specific configuration defaults and shared variables. Be aware that overlapping parameters from `config/woorkflow_config.yaml`overwrites any tool-specific parameters. 

The configuration must follow this structure: 

```
default = {
  "EXPERIMENT_NAME": "Ex15_03_WideQual",
  "SAMPLES": "Ex15_03_WideQual_simul-1",
  "REFERENCES": {
    "REF_GENOME": ["reference/RefSeq_sequence_Wuhan-Hu1.fa"]
  },
  "POSTPRED": {
    "LINEAGE_MIN_THRESHOLD": 0.01
  }
}
```

#### prepareSummary.py

This script should refprmat the tool's output into a standardized post-prediction format used across the pipeline. 

The output table **must contain the flowwoing columns:**
- timepoint
- lineage 1
- lineage 2
- ...
- others
- sample_name
- tool_name

#### Register the tool in the configuration

Finally, add the tool to `config/workflow_config.yaml` under the `TOOLS` section:

```
TOOLS:
  TOOLNAME_01:
    COLOUR_IN_REPORT: '#123456'
    INCLUDE_IN_ANALYSIS: true
    TOOL_NAME: toolname
    TOOL_LABEL: "ToolName\nv.X.Y.Z"
```
**Configuration Fields**
- `COLOUR_IN_REPORT` Color used in plots and reports
- `INCLUDE_IN_ANALYSIS` Whether the tool is included in benchmarking
- `TOOL_NAME` Internal identifier (used by the pipeline)
- `TOOL_LABEL` Human-readable name for plots and tables

# How to run

`main.py` is the top-level controller script for this project. It automates the execution of multiple Snakemake workflows by iterating over simulation parameters defined in `config/workflow_config.yaml` and launching complete pipeline for each experiment.

Instead of manually running each Snakemake workflow, `main.py`:
- Generates experiment names
- Randomizes seeds
- Updates `config/workflow_config.yaml` between each experiment
- Executes simulation, variant calling and lineage deconvolution workflows sequentially. 

## Experiment Naming Scheme
Experiments are automaiticall named using the format: 
**Ex<*experiment_number*>_<*seed_index*>_<*dataset_name*>**

# Step-by-step guide
## Configure Workflow Parameter

Edit `config/workflow_config.yaml`to define:
- Datasets to be simulated
- Quality score values
- Number of random seeds
- Reference paths and tool settings

## Run the Pipeline

From the project root directory:

```
python main.py
```

**Snakemake Execution Details**

For each experiment the following commands are executed:

```
snakemake --snakefile rules/simulation.smk \
          -c20 --use-conda --rerun-incomplete --rerun-triggers mtime

snakemake --snakefile variantCalling.smk \
          -c20 --use-conda

snakemake --snakefile lineage_deconvolution.smk \
          -c20 --use-conda --keep-going --rerun-incomplete --rerun-triggers mtime

```
- `--use-conda` ensures reproducible environents
- `--rerun-incomplete`restarts failed or interrupted jobs
- `--keep-going`continues independent jobs if one fails.
- `-c` describes number of cores to be used by each worfklow, default: 20

## Analyze Results

Render `PostPredict_plots.Rmd` to visualize and interpret benchmarking outcomes.


# Notes & Best Practices

- Each tool is fully isolated via its Conda environment
- Tools can be enabled or disabled via configuration
- No changes to core workflows are required to add new tools
- Do not run multiple instances of main.py simultaneously

# Credits and third-party software

This repository includes third-party software developed by other authors, which is incorporated to enable the execution of the tool presented here and to support the reproduction of previously published results.

Detailed attribution, licensing information, commit identifiers, and references to the original publications are provided in the [CREDITS](CREDITS) file. Users are encouraged to cite the original works when making use of these components.