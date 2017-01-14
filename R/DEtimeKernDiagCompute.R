#' @title Retrieve the diagonal of the DEtime kernel
#' @param kern DEtime kernel structure to be computed
#' @param X A two-column matrix where the first column of this matrix is the time points for control and perturbed conditions and the second column uses '1' to represent time points from control condition and '2' to represent time points from perturbed condition
#' @return 
#'    Kd Vector containing computed diagonal elements of the kernel structure
#' @details
#'    \code{Kd <- DEtimekernDiagCompute(kern, X)} computes the diagonal of a DEtime kernel matrix for the given kernel.
#' @description
#'    Retrieve the diagonal of the DEtime kernel matrix.
#' @examples
#' kern <- list()
#' kern <- DEtimeKernParamInit(kern)
#' X <- matrix(c(seq(3:8),seq(4:8),rep(1,6),rep(2,5)),ncol=2)
#' Kd <- DEtimeKernDiagCompute(kern, X)
#' @export

DEtimeKernDiagCompute <-
function (kern, X) {
  #k <- diag(DEtimeKernCompute(kern,x))
  Kd <- matrix(kern$variance, dim(X)[1], 1)
  
  if ("isNormalised" %in% names(kern) && kern$isNormalised)
    Kd <- Kd * sqrt(kern$inverseWidth/(2*pi))

  return (Kd)
}

