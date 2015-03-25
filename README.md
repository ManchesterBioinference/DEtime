# DEtime
DEtime is a R package for two-sample time series analysis using Gaussian process methods. 
This package implements the Gaussian regression framework for perturbation time point inferrence in a two sample case. 

## Installation

There are different ways to install DEtime:
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

* download the tarball file **DEtime_1.0.tar.gz** from https://github.com/ManchesterBioinference/DEtime/ page and install it:

   ```R
   install.packages("DEtime_1.0.tar.gz", repos=NULL)
   ```

## Getting Started
The package contains two main functions: **DEtime_infer** and **DEtime_rank**. 

* **DEtime_infer** is the main function for perturbation point inference.

* **DEtime_rank** is used to filter out these silent genes before any focused perturbation point inference work. 

The user is required to provide _times_, _ControlData_, _PerturbedData_, _replicate_no_, _gene_no_, etc to use
these two functions. For explanation of these arguments, please refer to the 
[vignettes](https://github.com/ManchesterBioinference/DEtime/tree/master/vignettes/vignettes.pdf) or 
[package helpfile](https://github.com/ManchesterBioinference/DEtime/tree/master/man/DEtime.Rd) asscoiated with this package.

## Examples

```R
### inport simulated dataset
data(SimulatedData)

### calculating the loglikelihood ratio for these tested genes. the result is saved into DEtime_rank.txt
res_rank <- DEtime_rank(times = times, ControlData = ControlData, PerturbedData=PerturbedData,
                   replicate_no=replicate_no, gene_no=gene_no, gene_ID=gene_ID, savefile=TRUE)
 
### get the index of these data with loglikelihood ratio larger than 4
idx <- which(res_rank[,2]>4)

### go on with the perturbation time inference if some of the data has passed the threshould test 
if (length(idx)>0){
     res <- DEtime_infer(times = times, ControlData = ControlData[idx,], PerturbedData=
         PerturbedData[idx,], replicate_no=replicate_no, gene_no=length(idx),times_test=times, gene_ID=gene_ID[idx])

### Print a summary of the results
print_DEtime(res)
### plot results for all the genes
plot_DEtime(res)
}
```
