#' @title Print the results from DEtime function
#' @param DEtimeOutput The returned value from \code{\link{DEtime_infer}} function
#' @return
#'    A table where the associated gene_ID, MAP, mean, median, 5 percentile and 95 percentile of the posterior distrbution of inferred perturbation time points are listed for each gene.
#' @description
#'    The function prints the results returned from \code{\link{DEtime_infer}} function, which will show the gene_ID associated with MAP, mean, median, ptl5 (lower 5 percentile) and ptl95 (upper 5 percentile) of the posterior distribution of inferred perturbation time points.
#' @seealso \code{\link{DEtime_infer}} \code{\link{plot_DEtime}}
#' @examples
#' data(SimulatedData)
#' res <- DEtime_infer(times = times, ControlData = ControlData, PerturbedData=PerturbedData,
#'        replicate_no=replicate_no, gene_no=gene_no, times_test=times_test, gene_ID=gene_ID)
#' print_DEtime(res)


print_DEtime <-
function(DEtimeOutput){
## This function is to print the statistical estimation results from the output of DEtime.R function. The MAP, mean, meadian, lower 5 percentile and upper 5 percentile of the estimation to perturbation time point will be printed, 
## ARG DEtimeOutput: the output from DEtime.R function

  cat('Perturbation point inference results from DEtime package: \n')
  cat('=======================================\n')
  
  print(DEtimeOutput$result, zero.print = ".")
  cat('=======================================\n')
  
  #cat('MAP of the perturbation point: ', DEtimeOutput$MAP, '\n')
  #cat('Mean of the perturbation point: ', DEtimeOutput$mean, '\n')
  #cat('Median of the perturbation point: ', DEtimeOutput$median, '\n')
  #cat('5ptl of the posterior distribution of inferred perturbation point: ', DEtimeOutput$fiveptl, '\n')
  #cat('95plt of the posterior distribution of inferred perturbation point: ', DEtimeOutput$nintyfiveptl, '\n')
}
