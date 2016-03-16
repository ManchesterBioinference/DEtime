#' @title Update a DEtime kernel structure with new parameters
#' @param kern The DEtime kernel structure whose parameters to be expanded
#' @param params Vector of parameters
#' @return
#'    kern Updated DEtime kernel structure
#' @description
#'    Update a DEtime kernel structure with new parameters.
#' @details
#'    \code{DEtimeKernExpandParam} returns a DEtime kernel structure filled with the parameters in the given vector. This is used as a helper function to enable parameters to be optimised in, for example, the optimisation functions.
#' @seealso
#'    \code{\link{DEtimeKernExtractParam}}
#' @examples
#' kern <- list()
#' kern <- DEtimeKernParamInit(kern)
#' K <- DEtimeKernCompute(kern, as.matrix(3:8))
#' ### the inverseWidth, variance and testing perturbation point for DEtime Kernel function
#' params <- c(2, 1, 20)  
#' Knew <- DEtimeKernExpandParam(K, params)
#' @export
DEtimeKernExpandParam <-
function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  kern$inverseWidth <- params[1]	## linear domain params, i.e. untransformed inverse-width and signal variance
  kern$variance <- params[2]
  kern$xp <- params[3]
  
  return (kern)
}
