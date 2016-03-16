#' @title DEtime kernel parameter initialisation
#' @param kern The DEtime kernel structure of which the parameters will be initialised
#' @return 
#'    kern The DEtime kernel structure with the parameters initialised
#' @description
#'    initialises the parameters of a DEtime kernel
#' @seealso
#'    \code{\link{DEtimeKernCompute}}
#' @examples
#' kern <- list()
#' K <- DEtimeKernParamInit(kern)
#' @export

DEtimeKernParamInit <-
function (kern) {
  kern$xp <- -1000
  kern$inverseWidth <- 1
  kern$variance <- 1
  kern$nParams <- 3
  kern$paramNames <- c("inverseWidth", "variance", "xp")
  
  kern$isStationary <- TRUE

  if ("options" %in% names(kern) && "isNormalised" %in% names(kern$options) && kern$options$isNormalised)
    kern$isNormalised <- TRUE
  else
    kern$isNormalised <- FALSE

  if ("options" %in% names(kern) && "inverseWidthBounds" %in% names(kern$options)) {
    kern$transforms <- list(list(index=1, type="bounded"),
                            list(index=2, type="positive"))
    kern$transformArgs <- list()
    kern$transformArgs[[1]] <- kern$options$inverseWidthBounds
    kern$inverseWidth <- mean(kern$options$inverseWidthBounds)
  } else if ("options" %in% names(kern) && "varianceBounds" %in% names(kern$options)) {
	  kern$transforms <- list(list(index=1, type="positive"),
							  list(index=2, type="bounded"))
	  kern$transformArgs <- list()
	  kern$transformArgs[[2]] <- kern$options$varianceBounds
	  kern$variance <- mean(kern$options$varianceBounds)
  } else if ("options" %in% names(kern) && "inverseWidthBounds" %in% names(kern$options) && "varianceBounds" %in% names(kern$options)) {
	  kern$transforms <- list(list(index=1, type="bounded"),
							  list(index=2, type="bounded"))
	  kern$transformArgs <- list()
	  kern$transformArgs[[1]] <- kern$options$inverseWidthBounds
	  kern$transformArgs[[2]] <- kern$options$varianceBounds
	  kern$inverseWidth <- mean(kern$options$inverseWidthBounds)
	  kern$variance <- mean(kern$options$varianceBounds)
  } else {
	  kern$transforms <- list(list(index=c(1,2), type="positive"))
  }

  return (kern)
}

