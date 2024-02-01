# New implementation in GPcounts package
Note: some of the R-packages that DEtime relies on are no longer supported. 
The model has been re-implemented in the GPcounts python package and can be used with a Gaussian or counts (negative binomial) likelihood. 
See https://github.com/ManchesterBioinference/GPcounts/blob/master/demo_notebooks/Branching_GPcounts.ipynb for an example application there to single-cell pseudotime branching. 
If using the DEtime kernel in GPcounts then please cite the GPcounts paper https://doi.org/10.1093/bioinformatics/btab486 as well as the original DEtime paper

# DEtime
DEtime is an R package for two-sample time series analysis using Gaussian process methods. 
This package implements the Gaussian regression framework for perturbation time point inferrence in a two sample case. 

The paper describing this package is available at DOI: https://doi.org/10.1093/bioinformatics/btw329 and arXiv: http://arxiv.org/abs/1602.01743. Please refer to the Jupyter notebook [**DEtime_illustration.ipynb**](https://github.com/ManchesterBioinference/DEtime/blob/master/DEtime_illustration.ipynb) for R codes about how to run the package.

## Installation

There are two ways to install DEtime:
* use `devtools` package: 
  + Install and load the `devtools` package to be able to directly install R packages hosted on github :

   ```R
   install.packages("devtools")
   library(devtools)
   ```
  + To install `DEtime` type:

   ```R
   install_github("ManchesterBioinference/DEtime")
   ```

* download the tarball file [**DEtime_1.0.tar.gz**](https://github.com/ManchesterBioinference/DEtime/blob/master/DEtime_1.0.tar.gz). 
  + Install the dependent packages, **spam** and **gptk**, first:
  ```R
  install.packages("spam")
  install.packages("gptk")
  ```
  
  + Install the tarball
   ```R
   install.packages("DEtime_1.0.tar.gz", repos=NULL, type="source")
   ```

## Getting Started
The package contains two main functions: **DEtime_infer** and **DEtime_rank**. 

* **DEtime_infer** is the main function for perturbation point inference.

* **DEtime_rank** is used to filter out these silent genes before any focused perturbation point inference work. 

The user is required to provide _times_, _ControlData_, _PerturbedData_ etc to use
these two functions. For explanation of these arguments, please refer to the 
[vignettes](https://github.com/ManchesterBioinference/DEtime/tree/master/vignettes/vignettes.pdf) asscoiated with this package.

## Examples

+ For straightforward perturbation time point inference without ranking, 

```R
library(DEtime)

### inport simulated dataset
data(SimulatedData)

### go on with the perturbation time point inference
res <- DEtime_infer(ControlTimes = ControlTimes, ControlData = ControlData, PerturbedTimes = PerturbedTimes, PerturbedData = PerturbedData)

### Print a summary of the results
print_DEtime(res)
### plot results for all the genes
plot_DEtime(res)
}
```

+ If ranking is needed, 
```R
library(DEtime)

### inport simulated dataset
data(SimulatedData)

### calculating the loglikelihood ratio for these tested genes. the result is saved into DEtime_rank.txt

res_rank <- DEtime_rank(ControlTimes = ControlTimes, ControlData = ControlData, PerturbedTimes, PerturbedData=PerturbedData, savefile=TRUE)
 
### get the index of these data with loglikelihood ratio larger than 1
idx <- which(res_rank[,2]>1)

### go on with the perturbation time inference if some of the data has passed the threshould test 
if (length(idx)>0){
     res <- DEtime_infer(ControlTimes = ControlTimes, ControlData = ControlData[idx,], PerturbedTimes = PerturbedTimes, PerturbedData=PerturbedData[idx,])
     ### Print a summary of the results
     print_DEtime(res)
     ### plot results for all the genes
     plot_DEtime(res)
  }
```
