from malfunction import Malfunction
import pandas as pd
from pandas import Series
import numpy as np
import scipy 
from scipy import stats
import os
import gzip
from statistics import mean, stdev

### STEP1 : Read data
#
#   In this totorial we assume there are four tissues generated in one batch (A, B, C, and D).
#   Users need to adjust according to their batch size.
###

dir_A = "Your directory"
dir_B = "Your directory"
dir_C = "Your directory"
dir_D = "Your directory"


T_A = Malfunction(dir = dir_A)
df_sum_A = T_A.get_sum()
T_B = Malfunction(dir = dir_B)
df_sum_B = T_B.get_sum()
T_C = Malfunction(dir = dir_C)
df_sum_C = T_C.get_sum()
T_D = Malfunction(dir = dir_D)
df_sum_D = T_D.get_sum()

 

### STEP 2: outlier detection
#
###
outlier_A = Malfunction.outlier(df_sum_A)
outlier_B = Malfunction.outlier(df_sum_B)
outlier_C = Malfunction.outlier(df_sum_C)
outlier_D = Malfunction.outlier(df_sum_D)


### STEP 3: generate list of shared outlier spots
#
###

outlier_list = pd.concat([outlier_A.barcode, outlier_B.barcode, outlier_C.barcode, outlier_D.barcode ])
uq_list = pd.Series.unique(outlier_list)

occur_list = []
for i in uq_list:
    oc1 = i in list(outlier_A.barcode)
    oc2 = i in list(outlier_B.barcode)
    oc3 = i in list(outlier_C.barcode)
    oc4 = i in list(outlier_D.barcode)
    if oc1 + oc2 + oc3 + oc4 >=3:
        occur_list.append(i)
    


### STEP 4: inspect a certain tissue and find the max outlier spot cluster size
#   Here we use the second tissue as sample
###
clusters = T_B.GET_cluster(occur_list)

size = []
for i in clusters:
    size.append(len(i))
    
print("###########")
print(max(size))
