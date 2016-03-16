#' @title Calculating the log-likelihood ratio of the biological time course data
#'
#' @param times Experimental time point at which time course biological data are measured
#' @param ControlData Time course data measured under control condition
#' @param PerturbedData Time course data measured under perturbed condition
#' @param gene_ID ID of these genes addressed in this study
#' @param savefile A BOOLEAN argument which is used to indicate if the ranking list will be saved or not
#' @return 
#'   \code{DEtime_rank} returns a dataframe object whose first column is the gene_ID and second column is the Loglikelihood_ratio of the named gene.
#' @description 
#'   \code{DEtime_rank} intends to rank biological time course data measured under control and perturbed conditions. In the function, an independent GP model and a combined GP model are used to model the data separately, the difference between the log-likelihood is used as the rank factor for the dataset. A high rank factor normally indicates better differential expression.
#' @details
#'    Both control and perturbed data have to be measured at the same time points with the same number of replicates. Replicates are required to be obtained across all time points. ControlData and PerturbedData are two matrices where each row represents the time course data for one particular gene under either control or perturbed condition. The columns for both ControlData and  PerturbedData are ordered by the time sequencing followed by replicates.
#' @examples
#' ### import simulated data
#' data(SimulatedData)
#' ### calculating loglikelihood ratio
#' res <- DEtime_rank(times = times, ControlData = ControlData, 
#'   PerturbedData=PerturbedData, gene_ID=gene_ID)
#' @import gptk
#' @export

DEtime_rank <- function(times,ControlData,PerturbedData,gene_ID=NULL, savefile=TRUE) {
  
  ## DEtime fits a mixed GP model on two dataset composed of both control data and perturbed data which diverges at a single perturbation point. The two dataset before the presumed perturbation point are fitted with a single GP and after the perturbation point, they are fitted with two independent GPs. The posterior distribution of the tested perturbation points is then obtained from which all thesestatistical info, for example, the MEDIAN, MODE and 5-95 percentile of the distribution can be derived.      
  ## ARG times: The time points of the measured data. At the moment, it is assumed the control and perturbed data are measured at exactly the same time points.
  ## ARG replicate_no: The replicate nos of the meaured data. At the moment, same number of replicates are accepted for this model. If this is not provided, the value will be obtained by comparing the size of the measured data and the length of the time points
  ## ARG ControlData : The control dataset used in the model, which is ordered by the first replicate for all time points then the second replicate for all time points and so on
  ## ARG PerturbedData : The perturbed dataset used in the model, which is ordered by the first replicate for all time points then the second replicate for all time points and so on
  ## COPYRIGHT: Jing Yang, 2015
  ##
  
  
  
  if (is.null(times)) {
    stop("Time points for the measurements are not provided, pls provide times.")
  }
  
  if (is.null(ControlData)){
    stop("Control data are not provided, pls provide ControlData.")
  }
  else if ((length(ControlData)%%length(times))>0) {
     stop("Dimension of the Control data is not correct.")
  }
  else if (is.null(dim(ControlData))) {
    dim(ControlData) <- c(length(ControlData)%/%length(times), length(times))
  }
  else if (dim(ControlData)[2] != length(times)){
    if(dim(ControlData)[1]==length(times)) { ControlData = t(ControlData)}}
  
  if (is.null(PerturbedData)){
    stop("Perturbed data are not provided, pls provide ControlData.")
  }
  else if ((length(PerturbedData)%%length(times))>0) {
     stop("Dimension of the perturbed data is not correct.")
  }
  else if (is.null(dim(PerturbedData))) {
    dim(PerturbedData) <- c(length(PerturbedData)%/%length(times), length(times))
  }
  else if (dim(PerturbedData)[2] != length(times)){
    if(dim(PerturbedData)[1]==length(times)) { PerturbedData = t(PerturbedData)}}
 
  
  if (!all.equal((dim(ControlData)),dim(PerturbedData))){
    stop("Dimension of the input data is not correct, please note ControlData and PerturbedData are both matrix and have to be in the same size.")
  }
  
  
  if (is.null(gene_ID)){
    gene_ID <- as.character(seq(1,length(PerturbedData)%/%length(times)))
  }
  
  ### interpolation points
  Data <- cbind(PerturbedData,ControlData)
  times2 <- c(times, times)
  len_times <- length(times) 
  gene_no <- dim(ControlData)[1]
  dim(Data) <- c(gene_no, 2*len_times) 
  rank <- matrix(0,nrow=gene_no,ncol=1)
  
  #### model and parameters initialization
  for (idx in seq(1,gene_no)){
    model <- list() ## Allocate space for model.
    ### Initialize hyperparameters used in the GP model
    param <- rep(0, 5)  ## Hyperparameters
    ### Likelihood for each tested perturbation point of each gene
    loglikelihood_minusinfty <- numeric(0) ##  the likelihood
    loglikelihood_infty <- numeric(0) ##  the likelihood
    x <- matrix(times2, ncol=1)
    y <- matrix(scale((Data[idx,]),center=TRUE,scale=TRUE), ncol=1)
    
    options=gpOptions(approx="ftc")
    options$kern = list(type="cmpnd",comp=list(list(type="DEtime",options=list(inverseWidthBounds=c(1/(2*max(times)),1/(min(times)+0.25*(max(times)-min(times)))),varianceBounds=c(min(min(y),0.5),max(max(y),5)))),list(type="white")))
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

