#' @title Plot the results from DEtime function
#' @param DEtimeOutput The returned value from \code{\link{DEtime_infer}} function
#' @param plot_gene_ID The gene_IDs of those genes whose GP regression and posterior distribution of the perturbation time points will be plotted. If not supplied, all the genes will be plotted.
#' @return
#'    A figure for each gene whhich illustrates the GP regression (lower panel) as well as posterior distribution of its perturbation points (upper panel).   
#' @import gptk
#' @description
#'    This function plots the results returned from \code{\link{DEtime_infer}} function. The produced figures show the the posterior distribution of inferred perturbation time points on the upper panel and Gaussian Regression of the original data on the lower panel.
#' @seealso \code{\link{DEtime_infer}} \code{\link{print_DEtime}}
#' @examples
#' data(SimulatedData)
#' res <- DEtime_infer(times = times, ControlData = ControlData, 
#' PerturbedData=PerturbedData, times_test=times_test, gene_ID=gene_ID)
#' plot_DEtime(res,plot_gene_ID='gene1')
#' @export
plot_DEtime <-
function(DEtimeOutput, plot_gene_ID=NULL){
## plot the results from DEtime.R function.
## ARG DEtimeOutput : the return from DEtime.R function, which contains necessary information for plotting GP regression of the data and posterior distribution of the perturbation points
## ARG plot_gene_ID: a vector which represents the list of gene_IDs whose results will be plotted. If not provided, all genes will be plotted


  times <- DEtimeOutput$originaltimes
  times_test <- DEtimeOutput$times_test
  posterior <- DEtimeOutput$posterior
  best_param <- DEtimeOutput$best_param
  data <- DEtimeOutput$originaldata
  gene_ID <- DEtimeOutput$gene_ID
  gene_no <- DEtimeOutput$gene_no
 
  if (is.null(plot_gene_ID)){
    seq_genes <- seq(1,gene_no)
    print(seq_genes)
    cat('All genes will be plotted \n')
    }
  else {
    seq_genes <- which(plot_gene_ID==gene_ID)
    }
    
  if (!length(seq_genes)) {
   stop(" the gene_ID supplied for plotting is not in the range of original gene_IDs")
   }
   
  gene_no <- dim(data)[1]
  tstar <- matrix(seq(min(times)-(2*(max(times)-min(times))/10), max(times), length=200), ncol=1)
  tstar2 <- matrix(rep(tstar,2, bycol=TRUE), ncol=1)
  for (idx in seq_genes){
    y <- matrix(scale((data[idx,]),center=TRUE,scale=TRUE), ncol=1)  
    model <- gpExpandParam(DEtimeOutput$model[[idx]],best_param[idx,])
    Kx <- kernCompute(model$kern, times, tstar2)
    Ktrain <- kernCompute(model$kern, times)
    invKtrain <- .jitCholInv(Ktrain, silent=TRUE)$invM
    yPred <- t(Kx) %*% invKtrain %*% y
    yVar <- diag(abs(kernCompute(model$kern, tstar2) - t(Kx) %*% invKtrain %*% Kx))
    #dev.new()
    #pdf(paste0(gene_ID[idx],".pdf"), width=7, height=5)
    cat(paste(gene_ID[idx], 'is plotted\n'))
    .gpPlot_DEtime(model, tstar2, yPred, yVar, cbind(times_test,posterior[idx,]), title =paste('GPR result for gene ',gene_ID[idx], ' by DEtime package', sep=""))
    #dev.off()
    }
}
