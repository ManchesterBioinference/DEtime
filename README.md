# DEtime
DEtime is an R package for two-sample time series analysis using Gaussian process methods. 
This package implements the Gaussian regression framework for perturbation time point inferrence in a two sample case. 

The paper describing this package is available at DOI: https://doi.org/10.1093/bioinformatics/btw329 and arXiv: http://arxiv.org/abs/1602.01743 

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
