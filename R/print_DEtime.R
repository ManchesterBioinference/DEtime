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
