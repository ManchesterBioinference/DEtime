#' Inferring the perturbation time from biological time course data
#'
#' This package implements the Gaussian regression framework for perturbation time point inferrence in a two sample case. The package contains two main functions: DEtime and DEtime_rank. DEtime is the main function for perturbation point inference and DEtime_rank is used to filter out these silent genes before any focused perturbation point inference work. The package works on the time course data from a wild-type and a perturbed system. Acting upon pre-defined testing perturbation time, the function goes over these perturbation time candidates and derives their likelihoods. From Bayes' theory, under a uniform prior assumption, the posterior distribution of the tested perturbation time is derived from their corresponding likeliooods. Maximum a posterior (MAP), mean or median of the posterior distribution can be taken as the solution to the estimated perturbation time point.
#' @examples
#' ### Import simulated dataset
#' data(SimulatedData)
#'
#' ### Carrying out perturbation point inference for the first two genes in the
#' ### data with filtering by a threshold of 45 for the loglikelihood ratio.
#' ### This threshold is arbitrarily big and would not normally be used in practice.
#' ### We adopt it here in order to reduce the running time of this example.
#'
#' res_rank <- DEtime_rank(times = times, ControlData = ControlData[1:2,], 
#'   PerturbedData=PerturbedData[1:2,], replicate_no=replicate_no, gene_no=2, 
#'   gene_ID=gene_ID[1:2], savefile=TRUE)
#'
#' ### Get the index of these data with loglikelihood ratio larger than 4
#' idx <- which(res_rank[,2]>45)
#' 
#' if (length(idx)>0){
#'     res <- DEtime_infer(times = times, ControlData = ControlData[idx,], PerturbedData=
#'         PerturbedData[idx,], replicate_no=replicate_no, gene_no=length(idx), 
#'         times_test=times, gene_ID=gene_ID[idx])
#'      ### Print a summary of the results
#'      print_DEtime(res)
#'      ### Plot results of the first gene
#'      plot_DEtime(res,gene_ID[idx[1]])
#'  }
#' @docType package
#' @name DEtime
NULL

