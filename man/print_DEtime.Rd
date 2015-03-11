\name{print_DEtime}
\alias{print_DEtime}

\title{
print the results from DEtime function
}
\description{
The function prints the results returned from DEtime() function, which will show the gene_ID associated with MAP, mean, median, ptl5 (lower 5 percentile) and ptl95 (upper 5 percentile) of the posterior distribution of inferred perturbation time points.
}
\usage{
print_DEtime(DEtimeOutput)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{DEtimeOutput}{the returned value from DEtime function}
}
\value{
A table where the associated gene_ID, MAP, mean, median, 5 percentile and 95 percentile of the posterior distrbution of inferred perturbation time points are listed for each gene.
}
\author{
Jing Yang <ynnjing@gmail.com>
}

\seealso{
\code{\link{DEtime}}
\code{\link{plot_DEtime}}
}

\examples{
## read simulated example data
## load simulated data set
data(SimulatedData)
res <- DEtime(times = times, ControlData = ControlData, PerturbedData=PerturbedData,replicate_no=replicate_no, gene_no=gene_no, times_test=times_test, gene_ID=gene_ID)
print_DEtime(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Gaussian Regression }
\keyword{print }% __ONLY ONE__ keyword per line
