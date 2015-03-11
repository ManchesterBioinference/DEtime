\name{RealData}
\alias{RealData}
\docType{data}
\title{experimental time-series data used for testing}
\description{In this experiment, the aim is to study the transcriptional change occuring in Arabidopsis following inoculation with P. syringae pv. tomato DC3000 (PtoDC3000) versus the disarmed strain Pto DC3000hrpA. Three time series data are included: a control time series comprising challenge with 10 mM MgCl2 (control condition), infection of Arabidopsis with
virulent Pseudomonas syringage pv. tomato DC3000, which leads to disease developmen (perturbed condition 1) and infection of Arabidopsis with the disarmed strain DC3000hrpA (perturbed condition 2)}

\usage{data(RealData)}
\format{a list contains 100 sets of gene_ID, ControledData (perturbed condition 1) and PerturbedData (perturbed condition 2)}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\examples{
## Original experimental time points
times <- matrix(c(0,2,3,4,6,7,8,10,11,12,14,16,17.5), ncol=1)
## number of replicates 
replicate_no  <- 4
## number of genes that will be studied in this example
gene_no <- 4
## read data from file
data(RealData)
## define the perturbation time point candidates
times_test <- seq(min(times),max(times),by=(max(times)-min(times))/(4*length(times)))
## gene_ID

## inferring the perturbation time point
res <- DEtime(times = times, ControlData = ControlData[1:gene_no,], PerturbedData=PerturbedData[1:gene_no,], replicate_no=replicate_no, gene_no=gene_no, times_test=times_test, gene_ID=gene_ID[1:gene_no])
}
\keyword{datasets}
