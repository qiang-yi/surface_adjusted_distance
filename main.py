#import modules
import sys
print (sys.version)
sys.path.append("C:/Users/yi/git/TerrainMetrics_conda2/Update")
# next line loads packages installed for my user account
sys.path.append("C:/ProgramData/Anaconda3/lib/site-packages")
import time
import itertools
from osgeo import gdal
#import my modules
import surfaceAdjusted
import cProfile, pstats, StringIO
import os.path
import numpy as np
import pandas as pd
import csv


# for testing, just use the larger resolutions
resolution_L = [100, 1000] # [3,10,30,100,1000]

# NN interpolation does not work, leave commented out
methods = ['clos', 'wavg', 'biLin', 'biQua', 'biQub', 'TIN', 'p2p', 'NN']

# for testing, just use one state
study_areas = ["Colorado"] #["Colorado", "Nebraska", "Louisiana", "Washington", "NC", "Texas"]

data_dir = r"D:/SAD/Modified_2/"

def expandgrid(*itrs):
    """
    Generate all possible combinations of elements in multiple lists

    Args:
        - lists separated by commas

    Example:
        >>> expandgrid([0, 1], [2, 3, 4])
        >>> [[0, 0, 0, 1, 1, 1], [2, 3, 4, 2, 3, 4]]

    Returns:
        - a list of lists, with one element in each inner list for each
          combination of elements in the input lists
    """
    product = list(itertools.product(*itrs))
    return [[x[i] for x in product] for i in range(len(itrs))]

# Determine all of the cases to compute (area X resolution X transect X method combinations)
cases_list = []
for area in study_areas:
    area_start_time = time.time()
    area_path = data_dir + area + '/simulation/'
    area_transects = np.genfromtxt(area_path + 'tran_sim_pts.csv', delimiter=",")

    for resolution in resolution_L:

        n_transects = int(area_transects.shape[0] / 2)
        transect_indices = [i for i in range(n_transects)]

        if resolution == 3:
            # subset 3m resolution to "clos" method only
            cases = expandgrid(transect_indices, ["clos"], [resolution], [area_path], [area])
        else:
            # determine all possible combinations of transects and methods
            cases = expandgrid(transect_indices, methods, [resolution], [area_path], [area])

        n_cases = len(cases[0])

        df = pd.DataFrame(cases).transpose()
        df.columns = ["transect", "method", "resolution", "path", "area"]
        cases_list.append(df)

cases_df = pd.concat(cases_list)
print(cases_df.describe())

# subset transects
n_transects = 2
cases_df = cases_df.loc[lambda df: df.transect < n_transects, :]
cases_df

# try just one case
surfaceAdjusted.distance(
                cases_df['transect'].tolist()[0],
                cases_df['method'].tolist()[7],
                cases_df['resolution'].tolist()[0],
                cases_df['path'].tolist()[0])

# Calculate all cases
surfaceAdjusted.distance(
                cases_df['transect'].tolist(),
                cases_df['method'].tolist(),
                cases_df['resolution'].tolist(),
                cases_df['path'].tolist())
