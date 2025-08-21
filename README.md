test


## **BLADE**: a package for the purpose of tissue artifacts detection (Edge, Border, abnormal spots). 
It has both Python and R implementations. For detailed usage tutorial refer to usage tutorial section. 

### Main use:  Detecting edge artifacts, border artifacts, and location malfunction artifacts

### Context of use: Spatial transcriptomics data, especially Visium data

### Pitfalls: 
- Use of this method does not guarantee that there are no other types of artifacts.

- For different technologies or tissue samples, some parameters may need to be changed by hand to ensure optimal performance.

- Artifacts may not be detected if they: (1) occur only as gradual changes over larger pieces of tissue, (2) are low in strength, (3) occur in only small portions of the border or tissue edge.

- BLADE has not been thoroughly evaluated for data from technologies other than Visium.

### Principle of operation 

- Edge artifacts and border artifacts are detected by first identifying edge spots and border spots based on the image, and then t-tests are used to compare those spots to interior spots.

- Location malfunctions are identified at the batch level by identifying outlier spots in each image in the batch, identifying spot locations that are outliers in most or all of the images in the batch, and then testing whether there is a connected group of locations that is larger than should occur by chance.


### Theoretical properties and empirical evidence
- Edge and border artifacts detection inherits the theoretical properties from the t-test used.

- Location artifact detection inherits the theoretical properties of the outlier spot detection method. The test for detecting if the largest outlier cluster is implausibly large relies on simulating random spatial distributions tailored to the shape of the tissue samples in the images, and is asymptotically correct as the number of simulated tissue samples increases.

### Best practices

- Images tested with BLADE should be investigated in other ways (e.g. heatmaps, plotting sequencing depth by edge distance, etc.) to assess if any of the pitfall conditions might apply.

- BLADE should be used on Visium and other spatial transcriptomics data prior to other forms of analysis, to prevent bias
----------------------

## Current version 

## Installation 

### Python
```
pip install git+https://github.com/KummerfeldLab/artifactRemoval
```


### R
Due to many existing software are already developed in R. We provided a tutorial of applying package **'reticulate'** in R for run the same procedure as in python:  
```
install.packages("reticulate")
library(reticulate)
py_install("git+https://github.com/KummerfeldLab/artifactRemoval")
```



## Usage Tutorial


Here is an example of we read in a tissue sample and produce a cleaned up "tissue_position.csv" file. 

#### Explanation to major functions 
*remove_border()*: This function automatically remove any border spots. Return a message about how many border spots are removed.

*remove_edge(distance = 3)*: This function automatically remove any edge spots. It takes one parameter "distance" indicates the deepth of edge we are going to remove. Return a message about how many edge spots are removed.

*malfunction()*: This function automatically remove all malfunction spots. Return a message about how many malfunction spots are removed. 

NOTE: There might be overlap across the class of edge, border and malfunction spots. The final number of "in_tissue" spot returned might not necessarily equal to the sum of spots removal in three above functions. 
 

*review_removing()*: This function print the procedure of artifacts removal in the order of user run the removal functions. 

*simple_cleanser()*: This function restart the whole removal procedure. 

save(dir): This function specify a "tissue_position.csv" file in the dir specified by user.
WARNNING: for save function, user need to specify a "dir" that is different from the current location of "tissue_position.csv" file, since the save function might overwrite the "tissue_position.csv" file if one use the same address.  

#### Python Version:
```
dir = "your tissue directory"
test = Artifact_remove(dir = dir)

test.remove_border() # remove border spots
test.remove_edge(distance=3) # remove edge spots, distance means all spots within 3 spots from edge
test.remove_malfunction() # remove malfunction spots
test.remove_edge(distance=6) # for illustration for repeatedly remove more edge spots.


test.review_removing() # review what we have done
test.simple_cleanser() # back to original tissue sample
test.save(dir) # save the cleaned tissue "tissue_position.csv" file to the director you set

```

#### R Version:
```
use_python("~/.virtualenvs/r-reticulate/bin/python")
use_virtualenv("r-reticulate")
py_available()
py_config()
py_install("git+https://github.com/KummerfeldLab/artifactRemoval")
py_install("scipy")
py_install("numpy")
py_install("pickle")
#version <- "3.12.3"
#install_python(version)

artifactRemoval <- reticulate::import("artifactsRemoval")
dir = "your tissue directory"


test <- artifactRemoval$Artifact_remove(dir)
test$remove_border() # remove border spots
test$remove_edge(distance=3) # remove edge spots, distance means all spots within 3 spots from edge
test$remove_malfunction() # remove malfunction spots
test$remove_edge(distance=6) # for illustration for repeatedly remove more edge spots.


test$review_removing() # review what we have done
test$simple_cleanser() # back to original tissue sample
test$save(dir) # save the cleaned tissue "tissue_position.csv" file to the director you set
```
