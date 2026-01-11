# -*- coding: utf-8 -*-
"""
Created on Sun May 21 10:30:22 2023
@author: aschedl
Description: This script generates abundance files either based on logistic growth model
    or real timecourse data.
"""
# import libraries
import yaml, argparse, logging, ast, math, random, os, datetime
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from pathlib import Path
import PostPred_functions as func

## define functions
# logistic growth, exact shape based on parameters from config file
def logistic(t: int, init: int, b: int, c: int):
    """
    Function for logistic growth
    
    Parameters
    ----------
    t       :   time point
    init    :   inital value of cases
    b       :   growth rate, needs to be > 0
    c       :   maximum capacity
    
    formula: y(t) = c / (1 + a * e^-bt)
    initial value is: c / (1 + a)

    Returns
    -------
    float, rounded to next full number
    """
    return c / (1 + (c / init - 1) * math.exp(-b * t))

# generate growth line data set
def generate_growth_line(max_timepoints: int, start_timepoint: int, para: dict):
    """
    Function to generate logistic growth line
    
    Parameters
    ----------
    max_timepoints   :  maximal timepoints of simulation
    start_timepoint  :  start time of growth until max timepoint
    para    : dictionary for form parameter

    Returns
    -------
    list
    """

    x = np.linspace(start=0, stop=max_timepoints, num=max_timepoints+1)
    y = [logistic(i, para['init_cases'], para['growth_rate'], para['max_cases']) for i in x]
    growth_list = y[:((max_timepoints+1)-start_timepoint)]
    pad_list = []
    pad_list += [0] * (max_timepoints+1 - len(growth_list))

    return pad_list+growth_list

def draw_samples(df: pd.DataFrame, mode: str, s_number: int):
    """
    Function draws s_number of samples from dataframe

    Parameters
    ----------
    df          :   pandas dataframe to be drawn from
    mode        :   can be "all", "even" or "random" for interval of drawn samples
    s_number    :   number of samples to be drawn from dataframe

    Returns
    -------
    numpy ndarray that contains drawn sample timepoints
    """
    assert mode in ["random", "even", "all"], "Invalid Option for SAMPLE MODE, choose 'random'. 'even' or 'all."
    assert isinstance(s_number, int), "Please choose a interger as number of samples"

    t_s = np.linspace(start=0, stop=max_timepoints, num=max_timepoints+1, dtype=int)

    if mode == "random":
        random.shuffle(t_s)
        t_s = t_s[:s_number]
    elif mode == "even":
        t_s = np.linspace(start=0, stop=max_timepoints, num=s_number, dtype=int)

    # create tsv file for metadata
    sample_summary = ["sample\tsample_date\ttimepoint"]

    # create tsv file if timepoint is part of sampling list
    for index, row in df.iterrows():
        if index in t_s:
            sample_name = exp_name + "_simul-" + str(index+1) 
            data = []
            for column in df:
                if column != "sample_date":
                    if row[column] != 0.0:
                        data.append([column, row[column]])
            sample_summary.append(sample_name+"\t"+str(row["sample_date"])+"\t"+str(index+1)) 
            np.savetxt(fname=w_dir+"/abundances/" + sample_name +".tsv", X=data, delimiter="\t", fmt ='% s')

    with open(w_dir+"/" + exp_name + "_metadata.tsv", "w") as f:
        for line in sample_summary:
            f.write(line+"\n")

    # save dataframe
    df.to_csv(w_dir + "/" + exp_name + "_data.csv")

    print(df)
    return t_s

def plotAbundance(df: pd.DataFrame, mode: str, max_timepoints: int, timepoints_sampled: np.ndarray, real_timecourse: bool):
    """
    Function draws plot from timecourse datatable

    Parameters
    ----------
    df          :   pandas dataframe to be drawn from
    mode        :   can be "all", "even" or "random" for interval of drawn samples
    timepoints_sampled    :   timepoint where samples where taken

    Returns
    -------
    -
    """
    # plot abundances
    if real_timecourse:
        plt.figure().set_figwidth(15)
        plt.subplots_adjust(bottom=0.1, left=0.05, top=0.9)
    else:
        plt.figure()

    for col in df.columns:
        if col != "sample_date":
            plt.plot(df.index,df[col], label=col)

    plt.xticks(np.arange(start=0, stop=max_timepoints+1, step=10), rotation=90)
    if mode != "all":
        for index, s in enumerate(timepoints_sampled):
            plt.axvline(int(s), linestyle="dashed", color="gainsboro")
            plt.text(int(s), 110, "s"+str(index+1), rotation=90)
    if real_timecourse: plt.legend(frameon=False, bbox_to_anchor=(1, 1), loc='upper left')
    else: plt.legend()

    plt.savefig(w_dir + "/" + exp_name + "_plot.png")
    plt.close()

# add arguments
parser = argparse.ArgumentParser()
parser.add_argument("-o","--output_dir", type=str, default="")
parser.add_argument("-e","--experiment_name", type=str, default="")
parser.add_argument("-s","--seed", type=int, default="")
parser.add_argument("-m","--sample_mode", type=str, default="")
parser.add_argument("-n","--sample_number", type=int, default="")
parser.add_argument("-r","--real_timecourse", type=func.str2bool, default=False)
parser.add_argument("--real_timecourse_data", type=str, default="")
parser.add_argument("-t","--max_timepoints", type=int, default=365)
parser.add_argument("-v","--variants", type=ast.literal_eval, default={})
parser.add_argument("-f","--form", type=ast.literal_eval, default={})
parser.add_argument("-a","--constant_abundance", type=float, default=0.0)
parser.add_argument("-d","--start_date", type=str, default="")
args=parser.parse_args()

exp_name = args.experiment_name

if args.output_dir == "":
    output_dir = f"experiments/{exp_name}/simulation/abundances/"
else:
    output_dir = args.output_dir

# setup logging - send everything to stdout
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', 
                    level=os.environ.get("LOGLEVEL", "INFO"), 
                    datefmt='%Y-%m-%d %H:%M:%S')


w_dir = output_dir[:output_dir.rfind("/")]

logging.info("Starting simulation for {sample}.".format(sample=exp_name))
logging.info("Parameter settings are:")
logging.info("Seed: {s}".format(s=args.seed))
logging.info(yaml.dump(args))

# get experiment name
exp_name = args.experiment_name

# get sampling mode and number of samples
mode = args.sample_mode
s_number = args.sample_number

# create output directory and clean up from previous run
for folder in ["", "/abundances"]:
    if not os.path.exists(w_dir + folder):
        os.makedirs(w_dir + folder)
    else:
        for file_name in os.listdir(w_dir + folder):
            if os.path.isfile(w_dir+folder + "/" + file_name):
                os.remove(w_dir+folder + "/" + file_name)

if args.real_timecourse:
    f = Path.cwd() / Path(args.real_timecourse_data)
    df = pd.read_csv(f, sep="\t")

    df = df[["kw", "Pango_lineage", "r"]]
    df.rename(columns={"kw": "sample_date"}, inplace=True)

    df = df.pivot(index='sample_date', columns='Pango_lineage', values='r')
    df = df.fillna(0).round(decimals=2).reset_index()
    print(df.head())

    max_timepoints = len(df)

else:
    # define total length of simulated timecourse
    max_timepoints = args.max_timepoints

    # variant name and start timepoint from yaml file
    variants = args.variants
    constant = args.constant_abundance
    
    # handle constant variant abundance option
    try: constant = float(constant)
    except ValueError: raise AssertionError("constant variant abundance level needs to be a number")
    assert 0.0 <= constant <= 1.0, "constant variant abundance level needs to be between 0 and 1"

    # sort variants by ascending timepoints and generate dict
    variant_dict = {}

    for v_key in sorted(variants.items(), key=lambda x:x[1]):
        variant_dict[v_key[0]] = generate_growth_line(max_timepoints, v_key[1][0], args.form)

    # create dataframe from dict, round to next integer
    df = pd.DataFrame(variant_dict).round(0)

    # assuming later variants overtake previous variants
    for i, var in enumerate(variant_dict.keys()):
        df[var+"_rel"] = np.nan
        try: 
            next_var = list(variant_dict.keys())[i+1]
            df[var] = df.apply(lambda row: row[var] - row[next_var], axis=1)
        except IndexError:
            pass

    # total case numbers of all variants combined, minus constant
    df['sum'] = df.sum(axis=1) / (1 - constant)

    # create relative abundances from total case numbers
    for index, row in df.iterrows():
        percent = row["sum"] / 100.0 # is one percent

        for col in reversed(list(df.columns.values)):
            if '_rel' not in col and col != "sum":
                df.loc[index, [col+"_rel"]] = row[col] / percent

    # only create "others" if constant option was chosen
    if constant != 0.0:
        df['others_rel'] = [constant*100] * (max_timepoints+1)
    df['sum'] = df[[col for col in df.columns if "_rel" in col]].sum(axis=1)

    assert df['sum'].mean() == 100, "simulation does not sum to 100% - process aborted"

    cols_to_export= []
    for col in df.columns:
        if "_rel" in col: cols_to_export.append(col)

    df = df[cols_to_export].round(decimals=2)
    df.columns = df.columns.str.replace("_rel", "")

    rename = {}
    for var in variants.keys(): rename[var] = variants[var][1]
    df.rename(columns=rename, inplace=True)

    # add sample dates
    start_date = datetime.datetime.strptime(args.start_date, "%Y-%m-%d")
    dates=[]
    for i in range(len(df)):
        dates.append((start_date + datetime.timedelta(days=i+1)).date())
    
    df["sample_date"] = dates
    
timepoints_sampled = draw_samples(df, mode, s_number)
plotAbundance(df, mode, max_timepoints, timepoints_sampled, args.real_timecourse)

logging.info('\nSimulation done.')