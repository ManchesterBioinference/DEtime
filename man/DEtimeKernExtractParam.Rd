\name{DEtimeKernExtractParam}
\Rdversion{1.0}
\alias{DEtimeKernExtractParam}

\title{Extract the parameters used in a DEtime kernel.}
\description{
  Extract the parameters from a DEtime kernel.
}
\usage{
DEtimeKernExtractParam(kern, only.values = TRUE, untransformed.values = TRUE)
}
\arguments{
  \item{kern}{the DEtime kernel structure whose parameters to be extracted from.}
  \item{only.values}{ a BOOLEAN parameter with value TRUE or FALSE. When it is TRUE, only the values of the paramter vector are returned, when it is FALSE, both the names and the values of the paramter vector are returned}
  \item{untransformed.values}{ a BOOLEAN parameter with value TRUE or FALSE. When set to TRUE, the untransformed values of the parameter vector are returned}
}
\value{
  \item{params}{The parameters used in the DEtime kernel.}
}
\details{This function returns the hyperparameters used in a DEtime kernel structure.
}
\seealso{
\code{
  \link{DEtimeKernExpandParam}}
}
\examples{
kern <- list()
kern <- DEtimeKernParamInit(kern)
params <- DEtimeKernExtractParam(kern) 
     
}
\keyword{kernel}
