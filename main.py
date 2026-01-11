# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 11:25:22 2023
This script starts overall snakemake workflows based on config file
@author: aschedl
"""
import random, yaml, os
from pathlib import Path

print("starting")

# get parameters from config file
with open("config/workflow_config.yaml", "r") as cfile:
    config = yaml.safe_load(cfile)

cores = config["CORES_AVAILABLE"]
quality_scores = config["SIMULATION"]["QUALITY_SCORES"]
seeds = config["SIMULATION"]["N_SEED"]
start_experiment = 1

for dataset in config["SIMULATION"]["DATASET"]:
    experimentNo = start_experiment
    start_experiment = start_experiment + len(config["SIMULATION"]["QUALITY_SCORES"])
    suffix = config["SIMULATION"]["DATASET"][dataset]["NAME"]

    # parameter modifications ares experiment name, quality score, seed
    experiment_param = {}

    for score in quality_scores:
        if len(str(experimentNo)) == 1:
            prefix = "Ex0"
        else: prefix = "Ex"

        # number of different seeds
        for i in range(1,seeds+1):
            name = f"{prefix}{experimentNo}_0{i}_{suffix}"
            experiment_param[name] = [score, random.randint(1,999)]
        experimentNo += 1

    for index, experiment in enumerate(experiment_param.keys()):
        print(experiment)
        print(f"Running Experiment {index+1} of", len(experiment_param.keys()))
        config["EXPERIMENT_NAME"] = experiment
        config["CURRENT_SET"] = dataset
        config["SIMULATION"]["AMPLICON_PSEUDOCOUNTS"] = experiment_param[experiment][0]
        config["APP"]["SEED"] = experiment_param[experiment][1]
        
        with open("config/workflow_config.yaml", "w") as cfile:
            yaml.dump(config, cfile)

        cli = ""
        cli += f"snakemake --snakefile rules/simulation.smk -c{cores} --use-conda --rerun-incomplete --rerun-triggers mtime"
        cli += " && "
        cli += f"snakemake --snakefile rules/variantCalling.smk -c{cores} --use-conda"
        cli += " && "
        cli += f"snakemake --snakefile rules/lineage_deconvolution.smk -c{cores} --use-conda --keep-going --rerun-incomplete --rerun-triggers mtime"
        cli += ""

        os.system(cli)
    
print("finished.")