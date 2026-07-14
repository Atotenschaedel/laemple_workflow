# Readme

# Læmple Workflow

#### A Benchmarking Framework for Virus Lineage Deconvolution Tools for SARS-CoV-2 from Wastewater

This repository contains a reproducible benchmarking pipeline for evaluating SARS-CoV-2 lineage deconvolution tools from wastewater sequencing data. The pipeline is implemented using **Snakemake** (version ≥9.5.1) and **conda** (version ≥24.7.1), with configuration-driven execution and downstream analysis in **R**. 

It is designed to simulate dynamic lineage abundance timecourse, simulate sequencing data accordingly, perform lineage calling, apply different lineage deconvolution methods, and finally assess their performance against known ground truth, and render a telling report. The pipeline can be run entirely, or only partially. In this sense, simulated time courses can also be provided externally from true time courses, simulated reads can be substituted with real wastewater sequencing reads, accompanied by corresponding ground truth from clinical surveillance.

The purpose of the workflow is to select the most suitable tool for tasks at hand, evaluate the tools in various scenarios (e.g. poor sequencing data due to low virus abundance), and assess the effect of updated reference sets, methodology, or parametrization. To this end, users can adapt the configuration to their needs, and conveniently repeat the task with little hands-on effort.

In the following manual, installation and tailor-made configuration is described in detail.

# Table of Contents

1. [Overview](#overview)
2. [Supported Analysis Modes](#supported-analysis-modes)
3. [Installation](#installation)
4. [Running the Minimal Example](#running-the-minimal-example)
5. [Workflow Output - Report](#workflow-output)
6. [Step-by-step Guide](#step-by-step-guide)
7. [Advanced Customization](#advanced-customization)
8. [Developer Guide: Adding New Deconvolution Tools](#developer-guide-adding-new-deconvolution-tools)
9. [Project Structure](#project-structure)
10. [Troubleshooting](#troubleshooting)
11. [Credits and Third-Party Software](#credits-and-third-party-software)


## Overview

Wastewater sequencing is increasingly used to monitor viral diversity in populations. Estimating the relative abundance of viral lineages from mixed samples requires specialized deconvolution methods. This project provides a standardized framework to:

- Simulate time courses of simple lineage dynamics
- Simulate sequencing data with known lineage compositions
- Perform variant calling on simulated or real sequencing reads
- Run and benchmark lineage deconvolution tools
- Compare predicted lineage abundances to ground truth
- Visualize and summarize benchmarking results

The modular design allows individual components of the pipeline to be adapted or extended as new tools and methods become available. Furthermore, initial steps can also be skipped, by providing own, arbitrarily complex timecourse or providing own sequence data.

## Key Features

### - Reproducible Benchmarking

Entire analyses are configuration-driven and executed through Snakemake workflows and Conda environments.

### - Modular Architecture

New lineage deconvolution tools can be added without modifying the core framework.

### - Flexible Data Sources

Benchmarking can be performed using:

- Simulated or real lineage dynamics
- Simulated sequencing datasets
- Real wastewater sequencing data
- Hybrid analyses combining both

# Supported Analysis Modes

## 1. Simulated Data Benchmarking

Generate synthetic wastewater sequencing datasets with known lineage abundances.

```text
Simulate Lineage Timecourses
    ↓
Read Simulation 
    ↓
Variant Calling
    ↓
Lineage Deconvolution
    ↓
Benchmarking
    ↓
Report Generation
```

Synthetic sequencing reads are generated using SWAMPy (see [credits](credits.md)), which is integrated into the Læmple simulation pipeline. Simulation parameters such as sequencing platform, primer scheme and quality score profiles can be configured through `config/workflow_config.yaml`. The minimal example included with this repository provides reference sequences and primer definitions for the ARTIC v4 amplicon scheme, but alternative primer sets and reference resources can also be used. For detailed information on available simulation options and supported sequencing configurations, please refer to the SWAMPy documentation. Note that while the simulation workflow is flexible with respect to sequencing settings, the downstream variant-calling pipeline is currently designed for amplicon-based sequencing data.

## 2. Real Sequencing Data Analysis

Skip the simulation stage and begin directly from raw sequencing reads.

```text
FASTQ files
    ↓
Variant Calling
    ↓
Lineage Deconvolution
    ↓
Benchmarking
```

See [Using Real Sequencing Data](#using-real-sequencing-data) for more details.


---
# Installation

This project consist of two separate installation layers:

1. **Core pipeline installation** (base setup)
2. **Tool-specific installation** (each lineage deconvolution tool)

The pipeline is designed to be modular: new tools can be added without modifiying the core workflow.

## Overall pipeline installation

#### Clone the Repository

```
git clone https://github.com/Atotenschaedel/Laemple_workflow.git
cd Laemple_workflow/
```

#### Install Snakemake using `conda`

```
conda env create -f envs/snakemake.yaml
```
    
# Running the Minimal Example

For demonstration purpose, the default workflow is configured to include a minimal example of the complete workflow which includes only tools Freyja (https://anaconda.org/channels/bioconda/packages/freyja/overview) and VaQuERo (https://github.com/fabou-uobaf/VaQuERo.git) and uses additional reference sequences from https://github.com/corneliusroemer/pango-sequences.git as well as SWAMPy for simulation of sequence data from wastewater samples (https://github.com/goldman-gp-ebi/SWAMPy). It can be started by using:
    
```
conda activate snakemake
python main.py
Rscript -e "rmarkdown::render('PostPredict_report.Rmd')"
```

It should end with six new simulated experiments with each having two result files in the respective result folders: `freyja_v2.0.0_summary.csv`, and `vaquero_v24d9211_summary.csv`. Final Report can be rendered using `PostPredict_report.Rmd`

Initial expected run time, including generation of required `conda` environments for the minimal example is <2 hours.

# Workflow Output

After post-prediction analysis, with the included R markdown script PostPredict_report.Rmd, and html report in html format is being generated. The report includes the following sections

### Simulation

**Settings & Sample Overview:** A summary of the included tools, their version and any comments, as defined in the config file

**Simulated Lineage Abundance:** A overview of the underlying lineage timecourse, as simulated by the workflow, or as provided by the user in the first place.

**Simulated Coverage:** Distribution of the genomic read coverage for each experiment as controlled by the quality control parameter in the config file. 

### Results

**Comparing replicates:** Foreach simulated experiment and for each assessed tool the consitency between replicas is visualised for each timecourse. For classification results, Venn diagrams are used, for quantiative parameters, density plots are provided. Depicted variables are 

- RMSE, PP, TP, FP, FN, whole genome coverage
- Nr. of detected lineages

**Compare experiments:** Foreach simulated timecourse and for each assessed tool the consitency between experiments is visualised as density plot for the quantiative parameters, namely RMSE, PP, TP, FP, FN, whole genome coverage.

**Predicted timecoures:** Foreach simulated timecourse and for each assessed tool the the detected lineages and there relative deduced abundance are presented. Additionally, the simulated abundance is add for direct comparision on a per lineage basis.

**Lineage identifcation:** Foreach simulated timecourse and for each assessed tool the detected lineages are compared to the ground truth of simulated lineages. Here different ambiguity levels, meaning, how much tolerance for the detection of closely related lineages is given. The classification accuracy is measured with the following metrics: false negative, false positive, predicted positive, true positive, false negative rate, true positive rate, positive predictive value, Jaccard index, F1 score, false discovery rate.

**Library complexity and Sequencing quality"** Foreach simulated timecourse and for each assessed tool the classification accuracy metrics are visualised as a function of library complexity, i.e., number of simualted lineages in the library, and sequencing quality, i.e. observed genome read coverage.

**False Positive - Phylogenetic Relationship:** Foreach simulated timecourse and for each assessed tool the observed false positives are contextualised by plotting the distribution of the closes lineage simulated in the ground truth and its phylogenetic relationship, i.e., ancestor, descendant, or sibling lineage.

### Summary plot

Finally a Summary plot, depciting the **Jaccard index** vs **RMSE**, as a measure for qualitative and quantitative accuracy, respectively, is provided.


# Step-by-step guide

The minimal example provides a quick way to verify that the framework is installed correctly and to familiarize yourself with the workflow. For custom analyses, Læmple can be configured to use different simulation settings, reference resources, datasets, and lineage deconvolution tools. The following guide outlines the required configuration steps and explains how to execute a complete benchmarking analysis.

## Configure Workflow Parameter

Edit `config/workflow_config.yaml` to define:
- Datasets - *each dataset corresponce to one experimental design*
- Quality score settings  - *defines sequencing quality if simulation pipeline is used*
- Number of random seeds - *defines number of technical replicates*
- Reference files
- Tool settings

Experiments are automatically named using the format:  Ex<*experiment_number*>\_<*seed_index*>\_<*dataset_name*>

All pipeline parameters are defined in `config/worfklow_config.yaml`. 
This file controls datasets, simulation parameters, seeds, enabled tools and reporting behavior. 
No code changes are required to modify experimental setups.
The base version of Læmple comes with three configuration files: 

a. **workflow_config.yaml:** a minimal example running out of the box (see below)

b. **workflow_config_manuscript.yaml:** the config files for the analysis as presented in the manuscript. For this analysis to run, individual tools which are considered in the analysis need to be installed in `./bin`.

c. **workflow_config_commented.yaml:** a generic version of the config files with extensive comments about the scope and functionality of each parameter to set.

Additionally each subworkflow for any tool also contains tool-specific configuration if required in `rules/subworkflow_TOOLNAME/common.smk`. For more details for tool installation see [Developer Guide: Adding New Deconvolution Tools](#developer-guide-adding-new-deconvolution-tools)


## Run the Pipeline

From the project root directory:

```
python main.py
```

`main.py` is the top-level controller script for this project. It automates the execution of multiple Snakemake workflows by iterating over simulation parameters defined in `config/workflow_config.yaml` and launching complete pipeline for each experiment.

Instead of manually running each Snakemake workflow, `main.py`:
- Generates experiment names
- Randomizes seeds
- Updates `config/workflow_config.yaml` between each experiment
- Executes simulation, variant calling and lineage deconvolution workflows sequentially. 


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

- `--use-conda` ensures reproducible environments
- `--rerun-incomplete`restarts failed or interrupted jobs
- `--keep-going`continues independent jobs if one fails.
- `-c` describes number of cores to be used by each worfklow, default: 20

## Analyze Results

Knit the comparative report using `PostPredict_report.Rmd` for result visualization and interpretation. To this end, either use `Rstudio` or render the report directly from the command line using 

```
Rscript -e "rmarkdown::render('PostPredict_report.Rmd')"
```

# Advanced Customization

Læmple is designed to support custom benchmarking setups.

Common customization options include:

- Changing simulated variant decomposition
- Modifying simulated quality profiles & primer scheme
- Skipping simulation to use real world data (see)

All these customization are performed through configuration file `config/workflow_config.yaml`.

additional option are: 
- Changing tool specific settings or reference set
- Adding or removing tools

Please see more details in the specified subsections.

## Simulating lineage time courses

The simulation workflow supports two approaches for generating lineage abundance trajectories:

- Logistic growth simulation (default)
- Simulation based on a real time course using a user provided time course

### Configuring logistic growth simulation

By default, lineage abundance trajectories are generated using a logistic growth model. 
The lineages included in the simulation and their initial abundances are defined in `config/workflow_config.yaml`.

First, specify the reference sequences that should be available for simulation:

SIMULATION:
  REFERENCE_SEQUENCES:
    - references/reference_sequences.fa

A single FASTA file may contain one or multiple reference sequences. Each sequence must have a unique sequence identifier.

Next, define the simulated lineages within a dataset configuration:

SIMULATION:
  DATASET:
    SET_01:
      VARIANTS:
        variant_1:
          - 30
          - BA.2
        variant_2:
          - 230
          - XBB.1.5

The name of the variant must match an identifier present in the reference FASTA file. The number defines the initial timepoint when the lineage starts growing in the simulated time course.

**Multiple Datasets**

Several independent datasets can be simulated within a single workflow configuration by adding additional dataset sections:

SIMULATION:
  DATASET:
    SET_01:
      ...
    SET_02:
      ...
    SET_03:
      ...

Each dataset is simulated independently and can contain its own set of lineages, abundance distributions, and simulation parameters.

### Simulation based on a custom time course

To generate simulated sequencing data that follows an observed lineage abundance trajectory, enable the real time course option in `config/workflow_config.yaml`:

'''
USE_REAL_TIMECOURSE: true
TIMECOURSE_DATA: path/to/timecourse.tsv
'''

When enabled, lineage abundances are derived from the provided time-course file rather than being generated using the default logistic growth model.

Timecourse data needs to be provided using **tab-separated** file containing following columns:
- "sample_date" - sample date in YYYY-MM-DD format
- "lineage" - Name of lineage as string
- "relAbun" - realtive abundance data (0-1) as float

Reference sequences for all lineages in the custom timecourse need to be provided the same was as for the logistic groth model simulation. in `config/workflow_config.yaml`:

SIMULATION:
  REFERENCE_SEQUENCES:
    - references/reference_sequences.fa

## Using Real Sequencing Data

Instead of generating synthetic reads, users may start directly from existing sequencing datasets.
Following files are required to serve as ground truth for benchmarking:
    - sequencing data file
    - metadata file containing sample information
    - data file containing know relative abundance information

**Prepeare experimental directory**

Create a new experiment directory:

```
mkdir Ex01_01_realdata_example/data
```

**Sequencing data files requirements**

Place the sequencing files into the `data/` directory.

Supported input formats:
- Paired-end FASTQ files
- Paired BAM files

Example naming conventions:

- experimentName_sampleName-Timepoint_R1.fastq
- experimentName_sampleName-Timepoint_R2.fastq

or

- experimentName_sampleName-Timepoint_1.bam
- experimentName_sampleName-Timepoint_2.bam

**Metadata requirements**

Metadata needs to be provided as **tab-separated** tsv file named `experimentName_metadata.tsv`, containing following columns:

- "sample" - Complete sample name (experimentName_sampleName-Timepoint)
- "sample_date" - sample date in YYYY-MM-DD
- "timepoint" - integer value which timepoint in the timeseries the sample was taken

Place the corresponding metadata file into the same `data/` directory.

**Abundance data file requirements**

Abundance data needs to be provided as **comma-separated (,)** csv file named `experimentName_data.csv`, containing following columns:

- unnamed index column (optional)
- "*lineageName*" (for example "B.1.1.7") - abundance information 0-100 as float
- "*lineageName*" (for example "BA.2") - abundance information 0-100 as float
...
- "sample_date" - sample date in YYYY-MM-DD


**Configure the workflow**

In `config/workflow_config.yaml`: 

- Set `EXPERIMENT_NAME` to match the naming prefix used for all input files.
- Set `INCLUDE_SIMULATION_PIPELINE: false`

When this option is disabled, the pipeline skips the simulation stage and starts directly with variant calling and lineage deconvolution.

Note: The current variant-calling workflow is designed for amplicon-based sequencing data.

# Performing an Analysis

## Configure Workflow Parameters

Edit `config/workflow_config.yaml` to define:

- Datasets - *each dataset corresponce to one expermintal design*
- Quality score settings  - *defines sequencing quality if simulation pipeline is used*
- Number of random seeds - *defines number of technical replicates*
- Reference files
- Tool settings

### Experiment Naming Scheme

Experiments are automatically named using the format:  Ex<*experiment_number*>\_<*seed_index*>\_<*dataset_name*>

# Developer Guide: Adding New Deconvolution Tools

To compose a customized workflow, serving specific needs with respect to which tools/parameters should be compared, each lineage deconvolution tool needs to be installed independently and integrated into the pipeline in a standardized way.

To this end, the following components need to be set in place for each tool to be included:

1. Install/provide executables in `./bin`
2. Define `conda` environment
3. Compose Snakemake subworkflow
4. Compose output standardization script
5. Add configuration entry

### Install the tool itself

Follow the official installation instruction of the tools (e.g. GitHub README, documentation). Provide the scripts and executable of the tool encapsulated in a dedicated subdirectory under `./bin`.

#### Create a Conda environment

Create a conda yaml file named:  `envs/TOOLNAME.yaml`, which contains the tool (if applicable), all it dependencies and any required Python/R packages.

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
#### Add tool source code

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

This script should reformat the tool's output into a standardized post-prediction format used across the pipeline. 

The output table **must contain the following columns:**
- timepoint
- lineage 1
- lineage 2
- ...
- lineage N
- others
- sample_name
- tool_name

### Register the tool in the configuration

Finally, add the tool to `config/workflow_config.yaml` under the `TOOLS` section:

```
TOOLS:
  TOOLNAME_01:
    COLOUR_IN_REPORT: '#123456'
    INCLUDE_IN_ANALYSIS: true
    TOOL_NAME: toolname
    TOOL_LABEL: "ToolName version add.info"
```

**Configuration Fields**
- `COLOUR_IN_REPORT` Color used in plots and reports
- `INCLUDE_IN_ANALYSIS` Whether the tool is included in benchmarking
- `TOOL_NAME` Internal identifier (used by the pipeline)
- `TOOL_LABEL` Human-readable name for plots and tables

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
│   └── subworkflow_TOOL NAME           tool specific subworkflow
│       ├── Snakefile.smk               tool specific Snakemake workflow
│       ├── common.smk                  tool specific config file for subworkflow
│       └── prepareSummary.py           Python script to convert tool specfic output to standard output format
└──main.py              main script to start complete project
```

Additional folders that will be created after workflow completion:

```
├── logs
└── experiments
    ├── EXPERIMENT NAME 1               Format: <experiment name>_<replica number>_<timecourse name>
    │   ├── data                        contains simulated sequencing data in fastq format 
    │   ├── results                     
    │   │   ├── postPrediction          contains standardized output files from tool specific workflows
    │   │   ├── variantCall             contains all results from mutation calling pipeline
    │   │   └──TOOL NAME                tool specific output files 
    │   ├── simulation                  
    │   │   ├── abundances              abundances per samples in tsv format
    │   │   ├── QualityControl          QC report per sample
    │   │   ├── EXPERIMENT_NAME_data.csv        simulated abundances data over complete timecourse
    │   │   ├── EXPERIMENT_NAME_metadata.tsv    simulated metadata per sample 
    │   └── └── EXPERIMENT_NAME_plot.png        Plot of simulated timecourse      
    ├── EXPERIMENT NAME 2               
    │   ├── ...
    |   ┊
    ├── EXPERIMENT NAME 3               
    │   ├── ...
    |   ┊
    └── ...
```

# Notes & Best Practices

- Each tool is fully isolated via its Conda environment
- Tools can be enabled or disabled via configuration
- No changes to core workflows are required to add new tools
- Do not run multiple instances of main.py simultaneously
---


# Trouble shooting

* In case the individual conda environments can not be initiated 
   * make sure that the minimal version requirements for Snakemake (version ≥9.5.1) make and conda (version ≥24.7.1) are met.
   * disable the channel priority configuration for conda by specifying `conda config --set channel_priority disabled`
* In case you experience issues with the PostPredict_report.Rmd script, make sure that all required R packages are installed. Most conveniently, this can be done by open the script in Rstudio and allow Rstudio to make the required installations.

---

# Credits and third-party software

This repository includes third-party software developed by other authors, which is incorporated to enable the execution of the tool presented here and to support the reproduction of previously published results.

Detailed attribution, licensing information, commit identifiers, and references to the original publications are provided in the [credits](credits.md) file. Users are encouraged to cite the original works when making use of these components.
