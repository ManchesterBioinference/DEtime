\name{DEtimeKernDiagCompute}
\Rdversion{1.0}
\alias{DEtimeKernDiagCompute}
\title{Retrieve the diagonal of the DEtime kernel.}
\description{
  Retrieve the diagonal of the DEtime kernel matrix.
}
\usage{
DEtimeKernDiagCompute(kern, x)
}
\arguments{
  \item{kern}{DEtime kernel structure to be computed.}
  \item{x}{depending on the number of inputs, x can be the input data matrix
    (rows are data points) to the kernel computation, or
    the first input matrix to the kernel computation (forms the rows of the kernel).}
}
\details{
  \code{Kd <- DEtimekernDiagCompute(kern, x)}
  computes the diagonal of a DEtime kernel matrix for the given kernel.
}
\value{
  \item{Kd}{vector containing computed diagonal elements of the kernel structure.}
}
\examples{
kern <- list()
kern <- DEtimeKernParamInit(kern)
K <- DEtimeKernCompute(kern, as.matrix(3:8))
Kd <- DEtimeKernDiagCompute(kern, as.matrix(3:8))
}
\keyword{kernel}
