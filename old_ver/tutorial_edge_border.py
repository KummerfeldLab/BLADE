from artifacts_detection import Artifact_detect
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy 
import os
import gzip
from statistics import mean, stdev

### STEP 1: Include your directory of tissue
#   Note: Must include "/filtered_feature_bc_matrix" and "spatial" folder of 10X format
#
#
###

dir = "/data_delivery/lniedern/umgc/2023-q2/230329_A00223_1029_AHTGGKDRX2/Niedernhofer_Project_039/Analysis/spaceranger_V12D05-320_D1"

test_class = Artifact_detect(dir)

df_border = test_class.get_border()
df_edge = test_class.get_edge()
df_concat = pd.concat([df_edge,df_border])
df_inner = test_class.get_inner(df_concat)

### STEP 2: Using ttest for comparing between reads on inner spots and target class of spots (edge, boarder, edge and boarder combined) 

test1 = scipy.stats.ttest_ind(df_border.gene_count, df_inner.gene_count,equal_var = False)
if len(df_boarder ) == 0:
    test1 = "Noboarder"

test2 = scipy.stats.ttest_ind(df_edge.gene_count, df_inner.gene_count,equal_var = False)

print(test1.pvalue)
print(test2.pvalue)



dir = "/Users/wan00232/Documents/UMNTMC-spatial/Batch5/V12F28_066_A1/spatial/tissue_positions.csv"
df = pd.read_csv(dir)