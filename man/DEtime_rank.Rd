\name{DEtime_rank}
\alias{DEtime_rank}

\title{
Estimating the log-likelihood ratio of the biological time course data
}
\description{DEtime_rank function in DEtime Package intends to rank biological time course data measured under control and perturbed conditions. In the function, an independent GP model and a combined GP model are used to model the data separately, the difference between the log-likelihood is used as the rank factor for the dataset. A high rank factor normally indicates better differential expression. 
}
\usage{
DEtime_rank(times, ControlData, PerturbedData, replicate_no, gene_no, gene_ID = NULL, savefile=TRUE)
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
\item{gene_ID}{
ID of these genes addressed in this study}
\item{savefile}{
a BOOLEAN argument which is used to indicate if the ranking list will be saved or not}
}
\details{
Both control and perturbed data have to be measured at the same time points with the same number of replicates. Replicates are required to be obtained across all time points. ControlData and PerturbedData are two matrices where each row represents the time course data for one particular gene under either control or perturbed condition. The columns for both ControlData and  PerturbedData are ordered by the time sequencing followed by replicates.
}
\value{ The function will return a dataframe object which includes:
\item{res}{ a dataframe object, whose first column is the gene_ID and second column is the Loglikelihood_ratio of the named gene. }
}
\author{
Jing Yang <ynnjing@gmail.com>
}
\seealso{
\code{\link{DEtime}}
}
\examples{
data(SimulatedData)
res <- DEtime_rank(times = times, ControlData = ControlData, PerturbedData=PerturbedData, replicate_no=replicate_no, gene_no=gene_no, gene_ID=gene_ID)
  }
\keyword{Gaussian Process }
\keyword{Perturbation }% __ONLY ONE__ keyword per line
