from Artifact_remove import Tissue_obj
from Artifact_remove import Artifact_detect
from Artifact_remove import Artifact_remove

import pandas as pd
import numpy as np
import scipy 
import os
from statistics import mean, stdev
import copy

import pickle
import sys

if not sys.argv: 
    print("need to pass the directory")

else: 
    dir = sys.argv[1]
type(dir)
test = Artifact_remove(dir) 

test.remove_border()
test.remove_edge(distance=3)
test.remove_malfunction()
#test.remove_edge(distance=6)


test.review_removing()
test.simple_cleanser()

#test.remove_border()
#test.remove_edge(distance=3)
#test.remove_malfunction()
#test.remove_edge(distance=6)
test.save(test.dir)



