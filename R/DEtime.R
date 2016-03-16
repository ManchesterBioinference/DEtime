#' Inferring perturbation time point from biological time course data
#'
#' @param times Experimental time point at which time course biological data are measured
#' @param ControlData Time course data measured under control condition
#' @param PerturbedData Time course data measured under perturbed condition
#' @param times_test The predefined evenly spaced time points upon which perturbation will be evaluated
#' @param gene_ID ID of these genes addressed in this study
#' @param kernel covariance kernel used in this study, default is DEtime kernel, which is based on RBF kernel, and you can also choose DEtimeMatern32 or DEtimeMatern12 which are based on Matern32 and Matern12, respectively 
#' @return The function will return a DEtimeOutput object which includes:
#' \enumerate{
#'    \item result: the statistical estimation for the inferred perturbation time
#'    \item ____$MAP: maximum a posterior solution to the inferred perturbation time
#'    \item ____$mean: mean of the posterior distribution of the inferred perturbation time
#'    \item ____$median: median of the posterior distribution of the inferred perturbation time
#'    \item ____$ptl5: 5 percentile of the posterior distribution of the inferred perturbation time
#'    \item _____$ptl95: 95 percentile of the posterior distribution of the inferred perturbation time
#'    \item posterior: posterior distribution of the tested perturbation time points
#'    \item model: optimized GP model which will be used for later GP regression work
#'    \item best_param: optimized hyperparameter for the optimized GP model
#'    \item originaltimes: original experimental time points which will be used for future print or plot functions
#'    \item originaldata: original measured time course data which will be used for future print or plot functions
#'    \item times_test: tested perturbation time points
#'    \item gene_ID: the ID of genes for the data
#' }
#' @description 
#' This is the main function in DEtime Package, which applies a mixedGP kernel to time course data under control and perturbed conditions. It returns the posterior distribution of these predefined perturbation time candidates and relevant statistical estimations of the inferred perturbation time point.
#' @details
#' Both control and perturbed data have to be measured at the same time points with the same number of replicates. Replicates are required to be obtained across all time points. ControlData and PerturbedData are two matrices where each row represents the time course data for one particular gene under either control or perturbed condition. The columns for both ControlData and  PerturbedData are ordered by the time sequencing followed by replicates.
#' @examples
#' ### import simulated data
#' data(SimulatedData)
#' ### start perturbation time inference
#' res <- DEtime_infer(times = times, ControlData = ControlData, PerturbedData=PerturbedData)

DEtime_infer <- function(times,ControlData,PerturbedData,times_test=times,gene_ID=NULL) {

## DEtime_infer fits a mixed GP model on two dataset composed of both control data and perturbed data which diverges at a single perturbation point. The two dataset before the presumed perturbation point are fitted with a single GP and after the perturbation point, they are fitted with two independent GPs. The posterior distribution of the tested perturbation points is then obtained from which all thesestatistical info, for example, the MEDIAN, MODE and 5-95 percentile of the distribution can be derived.      
## ARG times: The time points of the measured data. At the moment, it is assumed the control and perturbed data are measured at exactly the same time points.
## ARG times_test: The time points which will be tested for perturbation effects. The original time points will be used if it is not provided.
## ARG ControlData : The control dataset used in the model, which is ordered by the first replicate for all time points then the second replicate for all time points and so on
## ARG PerturbedData : The perturbed dataset used in the model, which is ordered by the first replicate for all time points then the second replicate for all time points and so on
## COPYRIGHT: Jing Yang, 2015
##

if (is.null(times)) {
    stop("Time points for the measurements are not provided, pls provide times.")
}

if (is.null(times_test)) {
    times_test <- seq(min(times),max(times),length=50)
    dim(times_test) <- c(length(times_test),1)
}

if (is.null(ControlData)){
    stop("Control data are not provided, pls provide ControlData.")
}

if (is.null(PerturbedData)){
    stop("Perturbed data are not provided, pls provide PerturbedData.")
}

if (is.null(dim(ControlData))){
    ControlData <- array(ControlData)
}

if (is.null(dim(PerturbedData))){
    PerturbedData <- array(PerturbedData)
}

if (dim(PerturbedData)[2]!=length(times)) {
   if (dim(PerturbedData)[1]==length(times)) { PerturbedData = t(PerturbedData)}
   else {
   stop("Dimension of the perturbed data is not correct")
 }}

if (dim(ControlData)[2]!=length(times)) {
   if (dim(ControlData)[1]==length(times)) { ControlData = t(ControlData)}
   else {
   stop("Dimension of the control data is not correct")
 }}


if ((is.null(dim(ControlData)))||(is.null(dim(PerturbedData)))||(!all.equal((dim(ControlData)),dim(PerturbedData)))){
    stop("Dimension of the input data is not correct, please note ControlData and PerturbedData are both matrix and have to be in the same size.")
}

if (is.null(gene_ID)){
    gene_ID <- as.character(seq(1,dim(ControlData)[1]))
}

if ((max(times_test)>max(times))||(min(times_test)<min(times))){
    stop("input perturbation test points is our of the range of original times, pls change times_test")
}

### interpolation points
gene_no <- dim(PerturbedData)[1]
len_times <- length(times)
len_test <- length(times_test)
Data <- array(cbind(ControlData, PerturbedData),dim=c(gene_no,2*len_times))
times_rep2 <- c(times, times)

DEtimeOutput <- list()
MAP_DEtime <- matrix(0,nrow=gene_no,ncol=1)
mean_DEtime <- matrix(0,nrow=gene_no,ncol=1)
median_DEtime <- matrix(0,nrow=gene_no,ncol=1)
ptl5_DEtime <- matrix(0,nrow=gene_no,ncol=1)
ptl95_DEtime <- matrix(0,nrow=gene_no,ncol=1)

model_DEtime <- vector("list",gene_no)

posterior <- matrix(0,nrow=gene_no,ncol=len_test)
best_param <- matrix(0,nrow=gene_no,ncol=4)
model_DEtime <- list()
#### model and parameters initialization
for (idx in seq(1,gene_no)){
model <- list() ## Allocate space for model.
### Initialize hyperparameters used in the GP model
param <- rep(0, 5)  ## Hyperparameters
### Likelihood for each tested perturbation point of each gene
likelihood <- rep(0, len_test) ##  the likelihood

x <- matrix(times_rep2, ncol=1)
y <- matrix(scale((Data[idx,]),center=TRUE,scale=TRUE), ncol=1)

options=gpOptions(approx="ftc")
options$kern = list(type="cmpnd",comp=list(list(type=kernel,options=list(inverseWidthBounds=c(1/(2*max(times)),1/(min(times)+0.25*(max(times)-min(times)))),varianceBounds=c(max(min(y),0.5),max(max(y),5)))),list(type="white")))
#options$kern = list(type="cmpnd",comp=list(list(type="DEtime",list(type="white"))))
  if (sum(is.nan(y)) > (length(y)/2)) {
    cat('Majority of points in profile are NaN.\n')
    next
    }

options$isMissingData <- any(is.nan(y))
options$isSpherical <- !any(is.nan(y))
stdy = sd(c(y[!is.nan(y)])) ## Profile variance.

## Optimise GP log likelihoods.
model0 <- gpCreate(dim(x)[2], dim(y)[2], x, y, options)
model0 <- gpOptimise(model0,0)
param0 <- gpExtractParam(model0)  ## get optimized hyperparameters

param <- matrix(rep(param0,times=len_test), nrow = len_test, byrow=TRUE)
  
## estimate the likelihood for each tested time point

for (i in seq_along(times_test)){
    param0[3] <- times_test[i]
    model[[i]] = gpExpandParam(model0,param0)
    
    likelihood[i] <- exp(gpLogLikelihood(model[[i]]))
}
    

posterior[idx,] <- likelihood/sum(likelihood)
cum_posterior <- cumsum(likelihood/sum(likelihood))
MAP_DEtime[idx] <- times_test[which.max(likelihood)]
median_DEtime[idx] <- times_test[which(cum_posterior>0.50)[1]]
mean_DEtime[idx] <- sum(times_test * posterior[idx,])
ptl5_DEtime[idx] <- times_test[which(cum_posterior>0.05)[1]]
ptl95_DEtime[idx] <- times_test[which(cum_posterior>0.95)[1]]
param0[3] <- MAP_DEtime[idx,] ## use MAP as the best perturbation point
best_param[idx,] <- as.matrix(param0, nrow=1)
model_DEtime[[idx]] <- model0
}

#tstar <- matrix(seq(min(times)-(2*(max(times)-min(times))/10), max(times), length=200), ncol=1)
#  tstar2 <- matrix(rep(tstar,2, bycol=TRUE), ncol=1)

#if (TRUE) {
#    param0[3] <- times_test[which(likelihood==max(likelihood))]
#    model0 <- gpExpandParam(model0,param0)
#    Kx = kernCompute(model0$kern, x, tstar2)
#    Ktrain = kernCompute(model0$kern, x)
#    invKtrain = .jitCholInv(Ktrain, silent=TRUE)$invM
#    yPred = t(Kx) %*% invKtrain %*% y
#    yVar =diag(abs(kernCompute(model0$kern, tstar2) - t(Kx) %*% invKtrain %*% Kx))

    #gpPlot_DEtime(model0, tstar2, yPred, yVar, cbind(times_test,posterior), title =paste('GP regression plot of ',gene_ID,' with perturbation time at ',format(param0[3],digits=4),sep=""))
#    gpPlot_DEtime(model0, tstar2, yPred, yVar, cbind(times_test,posterior), title =paste('GPR result for ',gene_ID, ' by DEtime package', sep=""))

#}

DEtimeOutput$posterior <- posterior
DEtimeOutput$best_param <- best_param
DEtimeOutput$result <- data.frame(gene_ID = gene_ID, MAP= format(MAP_DEtime, digits=4, nsmall=2), mean = format(mean_DEtime, digits=4, nsmall=2), median = format(median_DEtime, digits=4, nsmall=2), ptl5 = format(ptl5_DEtime, digits=4, nsmall=2), ptl95 = format(ptl95_DEtime, digits=4, nsmall=2))

DEtimeOutput$originaltimes <- x
DEtimeOutput$originaldata <- Data
DEtimeOutput$times_test <- times_test
DEtimeOutput$model <- model_DEtime
DEtimeOutput$gene_ID <- gene_ID
DEtimeOutput$gene_no <- gene_no

return(DEtimeOutput)
}

