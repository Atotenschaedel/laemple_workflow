# -*- coding: utf-8 -*-
"""
Created on 19-08-2023
Author: aschedl
Description: This script collect different function needed to create summary files for Post Prediction analysis.
"""

import yaml, argparse
import numpy as np
import pandas as pd


# surpress performance warning
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

# return empty dataframe if tool did not finish
def return_empty(sim_df: pd.DataFrame, tool_name: str):
    """
    Function retruns a standard format with bare minimum info, if no datafile is found
    
    Parameters
    ----------
    sim_df      :   dataframe containing data about simulation, mainly number of timepoint that should exist
    tool_name   :   name of the tool that failed to produce datafile

    Returns
    ----------
    pandas dataframe
    """

    df = pd.DataFrame({'timepoint': pd.Series(dtype=int)})

    print("No Results found for ", tool_name, ".")
    df["timepoint"] = sim_df[["timepoint"]].copy()
    df["tool_name"] = tool_name
    df["sample_name"] = df["tool_name"] + "-" + str(df["timepoint"])

    return df

# filter result dataframes from prediction pipeline per tool
def filter_dataframe(df: pd.DataFrame, threshold, tool_name):
    """
    Function filters dataframe in a standardized way
    
    Parameters
    ----------
    df          :   dataframe containing unfiltered result data
    threshold   :   threshold of relative abundance, everything below will be summarised under "ohters" category
    tool_name   :   name of the tool 

    Returns
    ----------
    pandas dataframe
    """

    assert "timepoint" in df.columns, f"no timepoint column in dataframe for {tool_name}"

    # change all columns to float except timepoint
    for col in df.columns:
        if col != "timepoint":
            try: 
                df[col] = df[col].astype(float)
            except: pass

    # add missing columns if necessary
    if not "others" in df:
        df["others"] = np.nan
    else: pass
    if not "sample_name" in df:
        df["sample_name"] = (tool_name + "-") + df["timepoint"].astype(str)
    else: pass
    if not "tool_name" in df:
        df["tool_name"] = tool_name
    else: pass

    # filter dataset
    # sum of all lineages under threshold
    for index, row in df.iterrows():
        other_sum = 0
        for col in df.select_dtypes(include=['float64']).columns:
            if isinstance(row[col], float) and row[col] <= threshold and not np.isnan(row[col]): 
                other_sum = other_sum + row[col]
                df.loc[index, col] = np.nan
        df.loc[index, "others"] = row["others"] + other_sum

    # add to 100% if below, add difference to "others"
    for index, row in df.iterrows():
        sum_total = 0
        for col in df.select_dtypes(include=['float64']).columns:
            if isinstance(row[col], float) and not np.isnan(row[col]):
                sum_total = sum_total + row[col]
        if sum_total < 1.0:
            diff = 1.0 - sum_total
            df.loc[index, "others"] = row["others"] + diff

    df = df.round(2)
    df.dropna(axis=1, how='all', inplace=True)
    return df

# translate WHO names to pangolin lineages they were first described on
def who_pango_translate(who_file: str, who_name: str):
    """
    Function returns closest associated pangloin lineage to any give WHO name, based on associated file.
    Parameters
    ----------
    who_file   :   config file that defines which pangolin lineage is associated to which WHO name
    who_name   :   WHO name to be translated

    Returns
    ----------
    string
    """
    with open(who_file, "r") as file:
            who_voc = yaml.safe_load(file)

    if who_name in who_voc.keys():
        return who_voc[who_name][0]
    else:
        return ""

def str2bool(flag: str):
    """
    Function return bool value based on string value
    Parameters
    ----------
    flag    :   string value to be translated to boolean value

    Returns
    ----------
    bool
    """
    if flag.lower() in ("true", "t"):
        return True
    elif flag.lower() in ("false", "f"):
        return False
    else:
        raise argparse.ArgumenttypeError("Boolean value expected.")