#' Inferring perturbation time point from biological time course data
#'
#' @param ControlTimes Experimental time point at which time course biological data for the control case are measured, they have to be repeated if there are replicated measurements
#' @param ControlData Time course data measured under control condition
#' @param PerturbedTimes Experimental time point at which time course biological data for the perturbed case are measured, they have to be repeated if there are replicated measurements
#' @param PerturbedData Time course data measured under perturbed condition
#' @param TestTimes The predefined evenly spaced time points upon which perturbation will be evaluated. If undefined, TestTimes <- seq(min(c(ControlTimes, PerturbedTimes)), max(c(ControlTimes, PerturbedTimes)),length=50)
#' @param gene_ID ID of these genes addressed in this study. If undefinied, numbers will be used instead
#' @param bound.lengthscale bounds for the lengthscale used in the DEtime RBF kernel. When not provided,bound.lengthscale <- c(min(ControlTimes,PerturbedTimes), 4*max(c(ControlTimes,PerturbedTimes)))
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
#'    \item ControlTimes: original experimental time points for the control case which will be used for future print or plot functions
#'    \item ControlData: original measured time course data for the control case which will be used for future print or plot functions
#'    \item PerturbedTimes: original experimental time points for the perturbed case which will be used for future print or plot functions
#'    \item PerturbedData: original measured time course data for the control case which will be used for future print or plot functions
#'    \item TestTimes: tested perturbation time points
#'    \item gene_ID: the ID of genes for the data
#' }
#' @description 
#' This is the main function in DEtime Package, which applies a mixedGP kernel to time course data under control and perturbed conditions. It returns the posterior distribution of these predefined perturbation time candidates and relevant statistical estimations of the inferred perturbation time point.
#' @details
#' ControlTimes and PerturbedTimes can be ordered by either time series, for instance time1, time1, time2, time2, time3, time3 ... or replicate sequences, for instance: time1, time2, time3, time1, time2, time3. ControlData and PerturbedData are two matrices where each row represents the time course data for one particular gene under either control or perturbed condition. The orders of the ControlData and PeruturbedData have to match those of the ControlTimes and PerturbedTimes, respectively. 

#' @examples
#' ### import simulated data
#' data(SimulatedData)
#' ### start perturbation time inference
#' res <- DEtime_infer(ControlTimes = ControlTimes, ControlData = ControlData, 
#' PerturbedTimes = PerturbedTimes, PerturbedData = PerturbedData)
#' @import gptk
#' @import stats
#' @export

DEtime_infer <- function(ControlTimes,ControlData,PerturbedTimes,PerturbedData,TestTimes=NULL,gene_ID=NULL,bound.lengthscale=NULL) {

if (is.null(ControlTimes)) {
    stop("Time points for the control measurements are not provided, pls provide times.")
}
ControlTimes <- as.numeric(ControlTimes)
 
if (is.null(PerturbedTimes)) {
    stop("Time points for the perturbed measurements are not provided, pls provide times.")
}
PerturbedTimes <- as.numeric(PerturbedTimes)

if (is.null(ControlData)){
    stop("ControlData are not provided, pls provide controlData.")
  }
else if ((length(ControlData)%%length(ControlTimes))>0) {
     stop("Dimension of the ControlData is not correct.")
  }
 else if (is.null(dim(ControlData))) {
    dim(ControlData) <- c(length(ControlData)%/%length(ControlTimes), length(ControlTimes))
  }
  else if (dim(ControlData)[2] != length(ControlTimes)){
    if(dim(ControlData)[1]==length(ControlTimes)) { 
     cat("ControlData is in the wrong dimension, it is transposed now.\n")
     ControlData = t(ControlData)}
  }
 else{
  cat("ControlData is accepted\n")}
  

  if (is.null(PerturbedData)){
    stop("PerturbedData are not provided, pls provide PerturbedData.")
  }
 else if ((length(PerturbedData)%%length(PerturbedTimes))>0) {
     stop("Dimension of the PerturbedData is not correct.")
  }
  else if (is.null(dim(PerturbedData))) {
    dim(PerturbedData) <- c(length(PerturbedData)%/%length(PerturbedTimes), length(PerturbedTimes))
  }
  else if (dim(PerturbedData)[2] != length(PerturbedTimes)){
    if(dim(PerturbedData)[1]==length(PerturbedTimes)) {   
    cat("PerturbedData is in the wrong dimension, it is transposed now.\n")
    PerturbedData = t(PerturbedData)}} 
 else{
 cat("PerturbedData is accepted\n")}

if (is.null(gene_ID)){
    cat("gene IDs are not provided. Numbers are used instead.\n")
    gene_ID <- as.character(seq(1,dim(ControlData)[1]))
}
  else if(length(gene_ID)<dim(ControlData)[1]){
    cat("Insufficient gene IDs are provided, numbers are used for remainning ones\n")
    gene_ID <- c(gene_ID, as.character(seq(length(gene_ID)+1,dim(ControlData)[1])))
}

times <- c(ControlTimes,PerturbedTimes)

if (is.null(TestTimes)){
    cat("Testing perturbation time points are not provided. Default one is used.\n")
    TestTimes <- seq(min(times),max(times),length=50)
    dim(TestTimes) <- c(length(TestTimes),1)
}

if ((max(TestTimes)>max(times))||(min(TestTimes)<min(times))){
    cat("input testing perturbation time points are our of the range of original times, its range is changed accordingly\n")
   TestTimes <-TestTimes[which(TestTimes<max(times))]
   TestTimes <-TestTimes[which(TestTimes>min(times))]
   }

len_test <- length(TestTimes)

if ((max(ControlTimes)<min(PerturbedTimes)) || (min(ControlTimes>max(PerturbedTimes)))){
  stop("There are no intersection between the ControlData and the PerturbedData. The program has to be stopped here") 
}
if (is.null(bound.lengthscale)){
    bound.lengthscale <- c(max(diff(times)),4*max(times))
}
else if ((length(bound.lengthscale)!=2) || (bound.lengthscale[1]>bound.lengthscale[2])){
   cat("bound.lengthscale is not right. The default one is used instead\n")
   bound.lengthscale <- c(max(diff(times)),4*max(times))
}


gene_no <- dim(ControlData)[1]

times <- matrix(c(ControlTimes,PerturbedTimes),ncol=1)
data <- matrix(c(ControlData,PerturbedData),nrow=gene_no,ncol=dim(times)[1])

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
param <- rep(0, 3)  ## Hyperparameters
### Likelihood for each tested perturbation point of each gene
likelihood <- rep(0, len_test) ##  the likelihood

x <- matrix(cbind(times,c(rep(1,length(ControlTimes)),rep(2,length(PerturbedTimes)))), ncol=2)
y <- matrix(scale((data[idx,]),center=TRUE,scale=TRUE), ncol=1)

options=gpOptions(approx="ftc")
#options$kern = list(type="cmpnd",comp=list(list(type="DEtime",options=list(inverseWidthBounds=c(1/(2*max(times)),1/(min(times)+0.25*(max(times)-min(times)))),varianceBounds=c(max(min(y),0.5),max(max(y),5)))),list(type="white")))
options$kern = list(type="cmpnd",comp=list(list(type="DEtime",options=list(inverseWidthBounds=1.0/bound.lengthscale)),list(type="white")))
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

for (i in seq_along(TestTimes)){
    param0[3] <- TestTimes[i]
    model[[i]] = gpExpandParam(model0,param0)
    likelihood[i] <- exp(gpLogLikelihood(model[[i]]))
}

posterior[idx,] <- likelihood/sum(likelihood)
cum_posterior <- cumsum(likelihood/sum(likelihood))
MAP_DEtime[idx] <- TestTimes[which.max(likelihood)]
median_DEtime[idx] <- TestTimes[which(cum_posterior>0.50)[1]]
mean_DEtime[idx] <- sum(TestTimes * posterior[idx,])
ptl5_DEtime[idx] <- TestTimes[which(cum_posterior>0.05)[1]]
ptl95_DEtime[idx] <- TestTimes[which(cum_posterior>0.95)[1]]
param0[3] <- MAP_DEtime[idx,] ## use MAP as the best perturbation point
best_param[idx,] <- as.matrix(param0, nrow=1)
model_DEtime[[idx]] <- model0

cat(paste('gene', gene_ID[idx], 'is done\n', sep=' '))
}

#tstar <- matrix(seq(min(times)-(2*(max(times)-min(times))/10), max(times), length=200), ncol=1)
#  tstar2 <- matrix(rep(tstar,2, bycol=TRUE), ncol=1)

#if (TRUE) {
#    param0[3] <- TestTimes[which(likelihood==max(likelihood))]
#    model0 <- gpExpandParam(model0,param0)
#    Kx = kernCompute(model0$kern, x, tstar2)
#    Ktrain = kernCompute(model0$kern, x)
#    invKtrain = .jitCholInv(Ktrain, silent=TRUE)$invM
#    yPred = t(Kx) %*% invKtrain %*% y
#    yVar =diag(abs(kernCompute(model0$kern, tstar2) - t(Kx) %*% invKtrain %*% Kx))

    #gpPlot_DEtime(model0, tstar2, yPred, yVar, cbind(TestTimes,posterior), title =paste('GP regression plot of ',gene_ID,' with perturbation time at ',format(param0[3],digits=4),sep=""))
#    gpPlot_DEtime(model0, tstar2, yPred, yVar, cbind(TestTimes,posterior), title =paste('GPR result for ',gene_ID, ' by DEtime package', sep=""))

#}

DEtimeOutput$posterior <- posterior
DEtimeOutput$best_param <- best_param
DEtimeOutput$result <- data.frame(gene_ID = gene_ID, MAP= format(MAP_DEtime, digits=4, nsmall=2), mean = format(mean_DEtime, digits=4, nsmall=2), median = format(median_DEtime, digits=4, nsmall=2), ptl5 = format(ptl5_DEtime, digits=4, nsmall=2), ptl95 = format(ptl95_DEtime, digits=4, nsmall=2))

DEtimeOutput$ControlTimes <- ControlTimes
DEtimeOutput$PerturbedTimes <- PerturbedTimes
DEtimeOutput$ControlData <- ControlData
DEtimeOutput$PerturbedData <- PerturbedData
DEtimeOutput$TestTimes <- TestTimes
DEtimeOutput$model <- model_DEtime
DEtimeOutput$gene_ID <- gene_ID
DEtimeOutput$gene_no <- gene_no

cat('DEtime inference is done.\nPlease use print_DEtime or plot_DEtime to view the results.\n')
return(DEtimeOutput)
}

