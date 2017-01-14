#' @title Plot the results from DEtime function
#' @param DEtimeOutput The returned value from \code{\link{DEtime_infer}} function
#' @param BestPerturbPos The statistical factor used for the optimal estimate to the perturbation point. You can choose "mean", "median" or "MAP". The corresponding statistical inference of the posterior distribution of the pertubation points will be adopted as the best estimate to the perturbation. "MAP" will be used if this parameter is not provided.
#' @param plot_gene_ID The gene_IDs of those genes whose GP regression and posterior distribution of the perturbation time points will be plotted. If not supplied, all the genes will be plotted.
#' @return
#'    A figure for each gene whhich illustrates the GP regression (lower panel) as well as posterior distribution of its perturbation points (upper panel).   
#' @import gptk
#' @description
#'    This function plots the results returned from \code{\link{DEtime_infer}} function. The produced figures show the the posterior distribution of inferred perturbation time points on the upper panel and Gaussian Regression of the original data on the lower panel.
#' @seealso \code{\link{DEtime_infer}} \code{\link{print_DEtime}}
#' @examples
#' data(SimulatedData)
#' res <- DEtime_infer(ControlTimes = ControlTimes, ControlData = ControlData, 
#'        PerturbedTimes = PerturbedTimes, PerturbedData = PerturbedData)
#' plot_DEtime(res,BestPerturbPos="mean", plot_gene_ID=c('1','2'))
#' @import graphics
#' @import grDevices
#' @export

plot_DEtime <-
function(DEtimeOutput, BestPerturbPos=NULL, plot_gene_ID=NULL){
## plot the results from DEtime.R function.
## ARG DEtimeOutput : the return from DEtime.R function, which contains necessary information for plotting GP regression of the data and posterior distribution of the perturbation points
## ARG plot_gene_ID: a vector which represents the list of gene_IDs whose results will be plotted. If not provided, all genes will be plotted


  ControlTimes <- DEtimeOutput$ControlTimes
  PerturbedTimes <- DEtimeOutput$PerturbedTimes
  TestTimes <- DEtimeOutput$TestTimes
  posterior <- DEtimeOutput$posterior
  best_param <- DEtimeOutput$best_param
  ControlData <- DEtimeOutput$ControlData
  PerturbedData <- DEtimeOutput$PerturbedData
  gene_ID <- DEtimeOutput$gene_ID
  gene_no <- DEtimeOutput$gene_no

  
  param0 <-  as.character(DEtimeOutput$result$MAP)
  if (is.null(BestPerturbPos)){}
   else if (toupper(BestPerturbPos) == "MEAN"){
     cat('mean of the posterior distribution of tested time points is used as the best perturbation point in GP plotting\n')
     param0 <- as.character(DEtimeOutput$result$mean)
   }
   else if (toupper(BestPerturbPos) == "MEDIAN") {
     cat('median of the posterior distribution of tested time points is used as the best perturbation point in GP plotting\n')
     param0 <- as.character(DEtimeOutput$result$median)
  }
   else if (toupper(BestPerturbPos) == "MAP") {
     cat('Option for BestPerturbPos is wrong, MAP will be used instead\n')
  }
   else{
     cat('option for BestPerturbPos is wrong, MAP of the posterior distribution of tested time points is used as the best perturbation point in GP plotting\n')
}
  
  if (is.null(plot_gene_ID)){
    seq_genes <- seq(1,gene_no)
    cat('All genes will be plotted \n')
    }
  else {
    seq_genes <- match(plot_gene_ID, gene_ID)
    }
    
  if (!length(seq_genes)) {
   stop(" the gene_ID supplied for plotting is not in the range of original gene_IDs")
   }
   
  times <- c(ControlTimes,PerturbedTimes)
  X <- matrix(cbind(times,c(rep(1,length(ControlTimes)),rep(2,length(PerturbedTimes)))), ncol=2)

  gene_no <- dim(data)[1]
  tstar <- matrix(seq(min(times)-(2*(max(times)-min(times))/10), max(times), length=100), ncol=1)
  tstar2 <- matrix(rep(tstar,2, bycol=TRUE), ncol=1)
  X2 <- matrix(c(tstar2,c(rep(1,dim(tstar)[1]),rep(2,dim(tstar)[1]))), ncol=2)
  data <- cbind(ControlData,PerturbedData)
  for (idx in seq_genes){
    y <- matrix(scale((data[idx,]),center=TRUE,scale=TRUE), ncol=1)  
    
    best_param[idx,3] <- as.numeric(param0[idx])
    model <- gpExpandParam(DEtimeOutput$model[[idx]],best_param[idx,])
    Kx <- kernCompute(model$kern, X, X2)
    Ktrain <- kernCompute(model$kern, X)
    invKtrain <- .jitCholInv(Ktrain, silent=TRUE)$invM
    yPred <- t(Kx) %*% invKtrain %*% y
    yPred <- yPred*sd(data[idx,]) + mean(data[idx,]) 
    yVar <- diag(abs(kernCompute(model$kern, X2) - t(Kx) %*% invKtrain %*% Kx))
    yVar <- yVar * var(data[idx,])
    #yPred <- yPred * sd(data[idx,]) + mean(data[idx,])
    dev.new()
    #pdf(paste0(gene_ID[idx],".pdf"), width=7, height=5)
    cat(paste(gene_ID[idx], 'is plotted\n'))
    .gpPlot_DEtime(model, X2, ControlTimes, ControlData[idx,], PerturbedTimes, PerturbedData[idx,],yPred, yVar, cbind(TestTimes,posterior[idx,]), title =paste('GPR result for gene ',gene_ID[idx], ' by DEtime package', sep=""))
    #dev.off()
    }
}
