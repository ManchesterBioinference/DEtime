#' @title Retrieve the diagonal of the DEtime kernel
#' @param kern DEtime kernel structure to be computed
#' @param x Depending on the number of inputs, x can be the input data matrix (rows are data points) to the kernel computation, or the first input matrix to the kernel computation (forms the rows of the kernel)
#' @return 
#'    Kd Vector containing computed diagonal elements of the kernel structure
#' @details
#'    \code{Kd <- DEtimekernDiagCompute(kern, x)} computes the diagonal of a DEtime kernel matrix for the given kernel.
#' @description
#'    Retrieve the diagonal of the DEtime kernel matrix.
#' @examples
#' kern <- list()
#' kern <- DEtimeKernParamInit(kern)
#' K <- DEtimeKernCompute(kern, as.matrix(3:8))
#' Kd <- DEtimeKernDiagCompute(kern, as.matrix(3:8))


DEtimeKernDiagCompute <-
function (kern, x) {
  #k <- diag(DEtimeKernCompute(kern,x))
  k <- matrix(kern$variance, dim(as.array(x))[1], 1)
  
  if ("isNormalised" %in% names(kern) && kern$isNormalised)
    k <- k * sqrt(kern$inverseWidth/(2*pi))

  return (k)
}
DEtimeKernExpandParam <-
function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  kern$inverseWidth <- params[1]	## linear domain params, i.e. untransformed inverse-width and signal variance
  kern$variance <- params[2]
  kern$xp <- params[3]
  
  return (kern)
}

