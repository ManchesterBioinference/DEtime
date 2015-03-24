\name{DEtimeKernParamInit}
\Rdversion{1.0}
\alias{DEtimeKernParamInit}
\title{Kernel parameter initialisation.}
\description{
  initialises the parameters of a DEtime kernel.
}
\usage{
DEtimeKernParamInit(kern)
}
\arguments{
  \item{kern}{the DEtime kernel structure for which the parameters will be initialised.}
}
\value{
  \item{kern}{the DEtime kernel structure with the parameters initialised.}
}
\seealso{
\code{\link{DEtimeKernCompute}}.
}
\examples{
kern <- list()
K <- DEtimeKernParamInit(kern)

}
\keyword{kernel}
