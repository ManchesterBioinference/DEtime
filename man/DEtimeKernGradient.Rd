\name{DEtimeKernGradient}
\Rdversion{1.0}
\alias{DEtimeKernGradient}
\title{Compute the gradient with respect to the kernel parameters.}
\description{
  Compute the gradient with respect to the kernel parameters.
}
\usage{
DEtimeKernGradient(kern, x, x2, covGrad)
## kernGradient(kern, x, partial)
## kernGradient(kern, x, x1, x2, partial_)
}
\arguments{
  \item{kern}{the DEtime kernel structure for which the gradients are being computed.}
  \item{x}{when x2 is provided, x is the input locations associated with the rows of the kernel matrix; when x2 is missing, x is the input locations associated with both the rows and the columns of the kernel matrix} 
  \item{x2}{the input locations associated with the columnss of the kernel matrix}
  \item{covGrad}{a matrix of partial derivatives of the function of interest with respect to the kernel matrix. the matrix should have the same number of rows as x1 and the same number of columns as x2 has rows.}
}
\details{
  \code{g <- kernGradient(kern, x, covGrad)}
computes the gradient of functions with respect to the kernel parameters. As well as the kernel structure and the input positions, the user provides a matrix covGrad which gives the partial derivatives of the function with respect to the relevant elements of the kernel matrix.

  \code{g <- kernGradient(kern, x1, x2, covGrad)}
  computes the derivatives as above, but input locations are now provided in two matrices associated with rows and columns of the kernel matrix.
}
\value{
  \item{g}{gradients of the function of interest with respect to the kernel parameters. The ordering of the vector should match that provided by the function DEtimeKernExtractParam.}
}
\seealso{
\code{\link{DEtimeKernCompute}}, \code{\link{DEtimeKernExtractParam}}.
}
\examples{
kern <- list()
kern <- DEtimeKernParamInit(kern)
g <- DEtimeKernGradient(kern, as.matrix(c(1, 4)), array(1, c(2, 2)))
}
\keyword{kernel}
