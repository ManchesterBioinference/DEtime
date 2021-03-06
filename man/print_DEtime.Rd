% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_DEtime.R
\name{print_DEtime}
\alias{print_DEtime}
\title{Print the results from DEtime function}
\usage{
print_DEtime(DEtimeOutput)
}
\arguments{
\item{DEtimeOutput}{The returned value from \code{\link{DEtime_infer}} function}
}
\value{
A table where the associated gene_ID, MAP, mean, median, 5 percentile and 95 percentile of the posterior distrbution of inferred perturbation time points are listed for each gene.
}
\description{
The function prints the results returned from \code{\link{DEtime_infer}} function, which will show the gene_ID associated with MAP, mean, median, ptl5 (lower 5 percentile) and ptl95 (upper 5 percentile) of the posterior distribution of inferred perturbation time points.
}
\examples{
data(SimulatedData)
res <- DEtime_infer(ControlTimes = ControlTimes, ControlData = ControlData, 
PerturbedTimes = PerturbedTimes, PerturbedData = PerturbedData)
print_DEtime(res)
}
\seealso{
\code{\link{DEtime_infer}} \code{\link{plot_DEtime}}
}

