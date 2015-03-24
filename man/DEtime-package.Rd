\name{DEtime-package}
\alias{DEtime-package}
\docType{package}
\title{
Inferring the perturbation time from biological time course data
}
\description{This package implements the Gaussian regression framework for perturbation time point inferrence in a two sample case. The package contains two main functions: DEtime and DEtime_rank. DEtime is the main function for perturbation point inference and DEtime_rank is used to filter out these silent genes before any focused perturbation point inference work. The package works on the time course data from a wild-type and a perturbed system. Acting upon pre-defined testing perturbation time, the function goes over these perturbation time candidates and derives their likelihoods. From Bayes' theory, under a uniform prior assumption, the posterior distribution of the tested perturbation time is derived from their corresponding likeliooods. Maximum a posterior (MAP), mean or median of the posterior distribution can be taken as the solution to the estimated perturbation time point.
}
\details{
\tabular{ll}{
Package: \tab DEtime\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-11\cr
License: \tab GPL-3\cr
}
}
\author{
Jing Yang

Maintainer: Jing Yang <ynnjing@gmail.com>
}
\keyword{ package }
\seealso{
 \code{\link{DEtime}} \code{\link{DEtime_rank}} \code{\link{print_DEtime}} \code{\link{plot_DEtime}}
}
\examples{
## read simulated example data
data(SimulatedData)

## carrying out perturbation point inference without filtering
res <- DEtime(times = times, ControlData = ControlData, PerturbedData=PerturbedData,replicate_no=replicate_no, gene_no=gene_no, times_test=times, gene_ID=gene_ID)
## Print a summary of the results
print_DEtime(res)
## plot results for all the genes
plot_DEtime(res)

## carrying out perturbation point inference with filtering (use threshold of 4)
res_rank <- DEtime_rank(times = times, ControlData = ControlData, PerturbedData=PerturbedData,replicate_no=replicate_no, gene_no=gene_no, gene_ID=gene_ID, savefile=TRUE)
# get the index of these data with loglikelihood ratio larger than 4
idx <- which(res_rank[,2]>4)

if (length(idx)>0){
  res <- DEtime(times = times, ControlData = ControlData[idx,], PerturbedData=PerturbedData[idx,],replicate_no=replicate_no, gene_no=length(idx), times_test=times, gene_ID=gene_ID[idx])
## Print a summary of the results
  print_DEtime(res)
## plot results for all the genes
  plot_DEtime(res)}
}
