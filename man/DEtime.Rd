\name{DEtime}
\alias{DEtime}

\title{
Inferring perturbation time point from biological time course data
}
\description{ This is the main function in DEtime Package, which applies a mixedGP kernel to time course data under control and perturbed conditions. It returns the posterior distribution of these predefined perturbation time candidates and relevant statistical estimations of the inferred perturbation time point.  
}
\usage{
DEtime(times, ControlData, PerturbedData, replicate_no, gene_no, times_test = times, gene_ID = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{
Experimental time point at which time course biological data are measured}
  \item{ControlData}{
Time course data measured under control condition
}
  \item{PerturbedData}{
Time course data measured under perturbed condition
}
  \item{replicate_no}{
The number of replicates used in the studies experiment
}
  \item{gene_no}{
The number of genes addressed in this study
}
  \item{times_test}{
The predefined evenly spaced time points upon which perturbation will be evaluated
}
  \item{gene_ID}{
ID of these genes addressed in this study}
}
\details{
Both control and perturbed data have to be measured at the same time points with the same number of replicates. Replicates are required to be obtained across all time points. ControlData and PerturbedData are two matrices where each row represents the time course data for one particular gene under either control or perturbed condition. The columns for both ControlData and  PerturbedData are ordered by the time sequencing followed by replicates.
}
\value{ The function will return a DEtimeOutput object which includes:
\item{result}{
   the statistical estimation for the inferred perturbation time}
item{____$MAP}{
   maximum a posterior solution to the inferred perturbation time}  
\item{____$mean}{
   mean of the posterior distribution of the inferred perturbation time }
\item{____$median}{
   median of the posterior distribution of the inferred perturbation time}
\item{____$ptl5}{
   5 percentile of the posterior distribution of the inferred perturbation time}  
\item{_____$ptl95}{
    95 percentile of the posterior distribution of the inferred perturbation time}  
\item{posterior}{
    posterior distribution of the tested perturbation time points }
\item{model}{
    optimized GP model which will be used for later GP regression work}   
\item{best_param}{
    optimized hyperparameter for the optimized GP model } 
\item{originaltimes}{
    original experimental time points which will be used for future print or plot functions}  
\item{originaldata}{
   original measured time course data which will be used for future print or plot functions}  
\item{times_test}{
    tested perturbation time points}
\item{gene_ID}{the ID of genes for the data}
}
\author{
Jing Yang <ynnjing@gmail.com>
}
\seealso{
\code{\link{plot_DEtime}}
\code{\link{print_DEtime}}
}
\examples{
data(SimulatedData)
res <- DEtime(times = times, ControlData = ControlData, PerturbedData=PerturbedData, replicate_no=replicate_no, gene_no=gene_no, times_test=times, gene_ID=gene_ID)
  }
\keyword{Gaussian Process }
\keyword{Perturbation }% __ONLY ONE__ keyword per line
