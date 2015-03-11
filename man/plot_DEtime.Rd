\name{plot_DEtime}
\alias{plot_DEtime}

\title{
plot the results from DEtime function
}
\description{
This function plots the results returned from DEtime() function. The produced figures show the the posterior distribution of inferred perturbation time points on the upper panel and Gaussian Regression of the original data on the lower panel. 
}
\usage{
plot_DEtime(DEtimeOutput, plot_gene_ID = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DEtimeOutput}{
the returned value from DEtime() function
}
  \item{plot_gene_ID}{
the gene_IDs of those genes whose GP regression and posterior distribution of the perturbation time points will be plotted. If not supplied, all the genes will be plotted.
}
}
\value{
A figure for each gene whhich illustrates the GP regression (lower panel) as well as posterior distribution of its perturbation points (upper panel).
}
\author{
Jing Yang <ynnjing@gmail.com>
}
\note{
Please note that the MAP solution of the perturbation point is taken as the optimized estimate to the perturbation point and Gaussian Regression is derived based upon this estimated perturabtion point.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DEtime}}
\code{\link{print_DEtime}}
}

\examples{
## read simulated example data
## load simulated data set
data(SimulatedData)
res <- DEtime(times = times, ControlData = ControlData, PerturbedData=PerturbedData,replicate_no=replicate_no, gene_no=gene_no, times_test=times_test, gene_ID=gene_ID)
plot_DEtime(res,plot_gene_ID='gene1')
    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Gaussian Regression }
\keyword{plot }% __ONLY ONE__ keyword per line
