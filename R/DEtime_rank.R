#' @title Calculating the log-likelihood ratio of the biological time course data
#'
#' @param ControlTimes Experimental time point at which time course biological data for the control case are measured, they have to be repeated if there are replicated measurements
#' @param ControlData Time course data measured under control condition
#' @param PerturbedTimes Experimental time point at which time course biological data for the perturbed case are measured, they have to be repeated if there are replicated measurements
#' @param PerturbedData Time course data measured under perturbed condition
#' @param gene_ID ID of these genes addressed in this study
#' @param bound.lengthscale bounds for the lengthscale used in the DEtime RBF kernel. When not provided,bound.lengthscale <- c(max(diff(c(ControlTimes,PerturbedTimes))), 4*max(c(ControlTimes,PerturbedTimes)))
#' @param savefile A BOOLEAN argument which is used to indicate if the ranking list will be saved or not
#' @return 
#'   \code{DEtime_rank} returns a dataframe object whose first column is the gene_ID and second column is the Loglikelihood_ratio of the named gene.
#' @description 
#'   \code{DEtime_rank} intends to rank biological time course data measured under control and perturbed conditions. In the function, an independent GP model and a combined GP model are used to model the data separately, the difference between the log-likelihood is used as the rank factor for the dataset. A high rank factor normally indicates better differential expression.
#' @details
#' ControlTimes and PerturbedTimes can be ordered by either time series, for instance time1, time1, time2, time2, time3, time3 ... or replicate sequences, for instance: time1, time2, time3, time1, time2, time3. ControlData and PerturbedData are two matrices where each row represents the time course data for one particular gene under either control or perturbed condition. The orders of the ControlData and PeruturbedData have to match those of the ControlTimes and PerturbedTimes, respectively. 
#' @examples
#' ### import simulated data
#' data(SimulatedData)
#' ### calculating loglikelihood ratio
#' res <- DEtime_rank(ControlTimes = ControlTimes, ControlData = ControlData, 
#'        PerturbedTimes = PerturbedTimes,PerturbedData=PerturbedData)
#' @import gptk
#' @import stats
#' @import utils
#' @export

DEtime_rank <- function(ControlTimes,ControlData,PerturbedTimes,PerturbedData,gene_ID=NULL, bound.lengthscale=NULL, savefile=TRUE) {
  
  ## DEtime fits a mixed GP model on two dataset composed of both control data and perturbed data which diverges at a single perturbation point. The two dataset before the presumed perturbation point are fitted with a single GP and after the perturbation point, they are fitted with two independent GPs. The posterior distribution of the tested perturbation points is then obtained from which all thesestatistical info, for example, the MEDIAN, MODE and 5-95 percentile of the distribution can be derived.      
  ## ARG ControlTimes: The time points of the measured data under control condition. This does not have to be the same as the time points for perturbed condition
  ## ARG ControlData : The control dataset used in the model, which is the matched measurements at corresponding control time points
  ## ARG PerturbedTimes: The time points of the measured data under perturbed condition. This does not have to be the same as the time points for control condition
  ## ARG PerturbedData : The perturbed dataset used in the model, which is the matched measurements at corresponding perturbed time points
  ## COPYRIGHT: Jing Yang, 2015
  ##
  
if (is.null(ControlTimes)) {
    stop("Time points for the control measurements are not provided, pls provide times.")
}

if (is.null(PerturbedTimes)) {
    stop("Time points for the perturbed measurements are not provided, pls provide times.")
}

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
     warning("ControlData is in the wrong dimension, it is transposed now.")
     ControlData = t(ControlData)}
  }
 else{
  print("ControlData is accepted")}


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
    warning("PerturbedData is in the wrong dimension, it is transposed now.")
    PerturbedData = t(PerturbedData)}}
 else{
 print("PerturbedData is accepted")}

if (is.null(gene_ID)){
    print("gene IDs are not provided. Numbers are used instead")
    gene_ID <- as.character(seq(1,dim(ControlData)[1]))
}
  else if(length(gene_ID)<dim(ControlData)[1]){
    warning("Insufficient gene IDs are provided, numbers are used for remainning ones")
    gene_ID <- c(gene_ID, as.character(seq(length(gene_ID)+1,dim(ControlData)[1])))
}

times <- matrix(c(ControlTimes,PerturbedTimes),ncol=1)

if ((max(ControlTimes)<min(PerturbedTimes)) || (min(ControlTimes>max(PerturbedTimes)))){
  stop("There are no intersection between the ControlData and the PerturbedData. The program has to be stopped here")
}
if (is.null(bound.lengthscale)){
    bound.lengthscale <- c(max(diff(times)),4*max(times))
}

  gene_no <- dim(ControlData)[1]
  Data <- matrix(c(ControlData,PerturbedData),nrow=gene_no,ncol=dim(times)[1])

  ### interpolation points
  len_times <- length(times) 
  rank <- matrix(0,nrow=gene_no,ncol=1)
  
  #### model and parameters initialization
  for (idx in seq(1,gene_no)){
    model <- list() ## Allocate space for model.
    ### Initialize hyperparameters used in the GP model
    param <- rep(0, 5)  ## Hyperparameters
    ### Likelihood for each tested perturbation point of each gene
    loglikelihood_minusinfty <- numeric(0) ##  the likelihood
    loglikelihood_infty <- numeric(0) ##  the likelihood
    x <- matrix(cbind(times,c(rep(1,length(ControlTimes)),rep(2,length(PerturbedTimes)))), ncol=2)
    y <- matrix(scale((Data[idx,]),center=TRUE,scale=TRUE), ncol=1)
    
    options=gpOptions(approx="ftc")
    options$kern = list(type="cmpnd",comp=list(list(type="DEtime",options=list(inverseWidthBounds=1.0/bound.lengthscale)),list(type="white")))

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
    param0[3] <- -1000.00
    model0 <- gpExpandParam(model0,param0)
    model0 <- gpOptimise(model0,0)
    loglikelihood_minusinfty <- gpLogLikelihood(model0)
    
    param0[3] <- 1000.00
    model1 <- gpExpandParam(model0,param0)
    model1 <- gpOptimise(model1,0)
    loglikelihood_plusinfty <- gpLogLikelihood(model1)
    
    rank[idx] <- loglikelihood_minusinfty - loglikelihood_plusinfty
  }
  res <- data.frame(gene_ID=gene_ID,Loglikelihood_ratio=rank)
  #colnames(res) <- c('gene_ID', 'Loglikelihood_ratio')
  if (savefile){
    write.table(res, "DEtime_result.txt", sep="\t", row.names = FALSE, col.names = TRUE)
    cat('rank list saved in DEtime_rank.txt\n')
  }
  else{ print(res)}
  
  return(res)
  
}

