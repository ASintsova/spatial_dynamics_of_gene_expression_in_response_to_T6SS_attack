import configparser
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from matplotlib.patches import Ellipse
import os
import scipy as ss
import scipy.stats

import itertools
import seaborn as sns
from sklearn.decomposition import PCA

import visualization_functions as vis

# col1 = "#f93c3c"
# col2 = "#009a7f"
# col3 = "#9fa8ab"
# col5 = "#00490c"
# col4 = "#f4d895"


col1 = "#f34236"
col2 = "#d6c571"
col3 = "#88bc67"
col4 = "#2e8174"
col5 = "#143969"

colors = [col1, col2, col3, col4, col5]


wt_lfc = "WT Log2 Fold Change"
mut_lfc = "9C1 Log2 Fold Change"
log2FoldChange = "Log2 Fold Change"
padj = "Adjusted P value"

wt = "WT"
mut = "9C1"

mut_L = "Dienes Line"
wt_L = "Merge"

bh_mut = "Behind Dienes Line"
bh_wt = "Behind Merge"

thirty = "30 min"
zero = "0 min"
four = "4 hrs"

samples = {"Case1": "{} | {}".format(wt_L, zero),
           "Case2": "{} | {}".format(wt_L, four),
           "Case3": "{} | {}".format(mut_L, zero),
           "Case4": "{} | {}".format(mut_L, four),
           "Case5":"{} | {}".format(wt_L, thirty),
           "Case6":"{} | {}".format(wt_L, zero),
            "Case7":"{} | {} | {}".format(wt, bh_mut, thirty),
          "Case8":"{} | {} | {}".format(mut, bh_mut, thirty),
          "Case9": "{} | {}".format(mut_L, zero),
          "Case10": "{} | {}".format(mut_L, thirty),
          "Case11": "{} | {} | {}".format(mut, bh_mut, zero),
           "Case12":"{} | {} | {}".format(wt, bh_mut, zero),
          "Case13": "{} | {} | {}".format(wt, bh_wt,zero ),
           "Case14":"{} | {} | {}".format(wt, bh_wt, thirty)}


def process_config(config_file=''):
    """
    by default looks for config file in the same directory as the script
    :param config_file:
    :return:
    """
    if not config_file:
        config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config")
    print(config_file)
    config = configparser.ConfigParser()
    config.read(config_file)
    config_dict = {}
    for section in config.sections():
        config_dict[section] = {name: value for name, value in config.items(section)}
    return config_dict


def get_rid_of_genbank(df, column_name = "Function", genbank = "GenBank"):
    fx = df[column_name]
    fx = fx.str.replace("\({}\)".format(genbank), "")
    df[column_name] = fx
    return df
