\name{DEtimeKernExpandParam}
\Rdversion{1.0}
\alias{DEtimeKernExpandParam}

\title{Update a DEtime kernel structure with new parameters.}
\description{
  Update a DEtime kernel structure with new parameters.
}
\usage{
DEtimeKernExpandParam(kern, params)
}
\arguments{
  \item{kern}{the DEtime kernel structure to be updated.}
  \item{params}{vector of parameters.}
}
\value{
  \item{kern}{updated DEtime kernel structure.}
}
\details{This function returns a DEtime kernel structure filled with the parameters in the given
  vector. This is used as a helper function to enable parameters to be
  optimised in, for example, the optimisation functions.
}
\seealso{
\code{
  \link{DEtimeKernExtractParam}}
}
\examples{
kern <- list()
kern <- DEtimeKernParamInit(kern)
K <- DEtimeKernCompute(kern, as.matrix(3:8))
params <- c(2, 1, 20)  ## the inverseWidth, variance and testing perturbation point for DEtime Kernel function
Knew <- DEtimeKernExpandParam(K, params) 
     
}
\keyword{kernel}
