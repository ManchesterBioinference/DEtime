% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEtimeKernGradient.R
\name{DEtimeKernGradient}
\alias{DEtimeKernGradient}
\title{Compute the gradient with respect to the kernel parameters}
\usage{
DEtimeKernGradient(kern, X, X2, covGrad)
}
\arguments{
\item{kern}{The DEtime kernel structure for which the gradients are being computed}

\item{X}{When \code{X2} is provided, X is the input locations associated with the rows of the kernel matrix; when \code{X2} is missing, \code{X} is the input locations associated with both the rows and the columns of the kernel matrix}

\item{X2}{The input locations associated with the columnss of the kernel matrix}

\item{covGrad}{A matrix of partial derivatives of the function of interest with respect to thekernel matrix. the matrix should have the same number of rows as X1 and the same number of columns as \code{X2} has rows}
}
\value{
g Gradients of the function of interest with respect to the kernel parameters. The ordering of the vector should match that provided by \code{\link{DEtimeKernExtractParam}}.
}
\description{
Compute the gradient with respect to the kernel parameters.
}
\details{
\code{g <- kernGradient(kern, X, covGrad)} computes the gradient of functions with respect to the kernel parameters. As well as the kernel structure and the input positions, the user provides a matrix covGrad which gives the partial derivatives of the function with respect to the relevant elements of the kernel matrix.
    \code{g <- kernGradient(kern, X1, X2, covGrad)} computes the derivatives as above, but input locations are now provided in two matrices associated with rows and columns of the kernel matrix.
}
\examples{
kern <- list()
kern <- DEtimeKernParamInit(kern)
X <- matrix(c(seq(0,4),seq(0,4), rep(1,5),rep(2,5)),ncol=2)
g <- DEtimeKernGradient(kern, X, array(1, c(10, 10)))
}
\seealso{
\code{\link{DEtimeKernCompute}}, \code{\link{DEtimeKernExtractParam}}
}

