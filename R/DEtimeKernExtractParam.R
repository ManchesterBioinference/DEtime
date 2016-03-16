#' @title Extract the parameters used in a DEtime kernel
#' @param kern The DEtime kernel structure whose parameters to be extracted from
#' @param only.values A BOOLEAN parameter with value TRUE or FALSE. When it is TRUE, only the values of the paramter vector are returned, when it is FALSE, both the names and the values of the paramter vector are returned
#' @param untransformed.values A BOOLEAN parameter with value TRUE or FALSE. When set to TRUE, the untransformed values of the parameter vector are returned
#' @return 
#'    params The parameters used in the DEtime kernel
#' @description
#'    Extract the parameters from a DEtime kernel.
#' @details
#'    \code{DEtimeKernExtractParam} function returns the hyperparameters used in a DEtime kernel structure.
#' @seealso
#'    \code{\link{DEtimeKernExpandParam}}
#' @examples
#' kern <- list()
#' kern <- DEtimeKernParamInit(kern)
#' params <- DEtimeKernExtractParam(kern)
#' @export

DEtimeKernExtractParam <-
function (kern, only.values=TRUE,
                                 untransformed.values=TRUE) {
  params <- c(kern$inverseWidth, kern$variance, kern$xp)

  if ( !only.values )
    names(params) <- c("inverseWidth", "variance","xp")

  return (params)
}

