#' @title Compute the gradient with respect to the kernel parameters
#' @param kern The DEtime kernel structure for which the gradients are being computed
#' @param x When \code{x2} is provided, x is the input locations associated with the rows of the kernel matrix; when \code{x2} is missing, \code{x} is the input locations associated with both the rows and the columns of the kernel matrix
#' @param x2 The input locations associated with the columnss of the kernel matrix
#' @param covGrad A matrix of partial derivatives of the function of interest with respect to thekernel matrix. the matrix should have the same number of rows as x1 and the same number of columns as \code{x2} has rows   
#' @return
#'    g Gradients of the function of interest with respect to the kernel parameters. The ordering of the vector should match that provided by \code{\link{DEtimeKernExtractParam}}.
#' @description
#'    Compute the gradient with respect to the kernel parameters.
#' @details
#'    \code{g <- kernGradient(kern, x, covGrad)} computes the gradient of functions with respect to the kernel parameters. As well as the kernel structure and the input positions, the user provides a matrix covGrad which gives the partial derivatives of the function with respect to the relevant elements of the kernel matrix.
#'     \code{g <- kernGradient(kern, x1, x2, covGrad)} computes the derivatives as above, but input locations are now provided in two matrices associated with rows and columns of the kernel matrix.
#' @seealso
#'     \code{\link{DEtimeKernCompute}}, \code{\link{DEtimeKernExtractParam}}
#' @examples
#' kern <- list()
#' kern <- DEtimeKernParamInit(kern)
#' g <- DEtimeKernGradient(kern, as.matrix(c(1, 4)), array(1, c(2, 2)))


DEtimeKernGradient <-
function (kern, x, x2, covGrad) {
  
  if (dim(x)[1]==1){ x<-t(x)}
  replicate_no <- length(which(x==x[1]))/2
  if (replicate_no<1){
      replicate_no <- 1
  }
 
  l <- dim(x)[1]/(2*replicate_no)
  xp <- kern$xp
  x0 <- x[1:l]
  if ( nargs() == 3 ) {
    replicate_no2 <- replicate_no
    k <- DEtimeKernCompute(kern, x)
    dist2xx <- .dist2(x0,x0)
    dist2xxp1 <- .dist2(x0,xp)
    dist2xxp2 <- .dist2(xp,x0)
    l2 <- l
    covGrad <- x2
    if (xp >= max(x0)) { pos1 <- length(x0)}
    else { pos1 <- max(min(which((x0<xp)==FALSE))-1, 0)}
    pos2 <- pos1
    flag <- FALSE      
    
  } else if (nargs() == 4) {
    if (dim(x2)[1]==1){ x2<-t(x2)}
    replicate_no2 <- length(which(x2==x2[1]))/2
    if (replicate_no2<1){
      replicate_no2 <- 1
    }
    k <- DEtimeKernCompute(kern, x, x2)
    l2 <- dim(x2)[1]/(2*replicate_no2)
    x1 <- x2[1:l2]
    sorted_x <- sort(rbind(x0,x1),index.return=TRUE)
    xnew <- sorted_x$x
    ind <- sorted_x$ix
    dist2xx <- .dist2(xnew,xnew)
    dist2xxp1 <- .dist2(xnew,xp)
    dist2xxp2 <- .dist2(xp,xnew)
    if (xp >= max(xnew)) { pos1 <- length(xnew)}
    else {pos1 <- max(min(which((xnew<xp)==FALSE))-1,0)}
    pos2 <- pos1
    flag <- TRUE
  }
  
  k_tmp <- matrix(0,nrow =(2*l),ncol=(2*l2))
  k_tmp[1:l,1:l2] <- dist2xx
  k_tmp[(l+1):(2*l),(l2+1):(2*l2)] <- dist2xx
  k_tmp[1:l,(l2+1):(2*l2)] <- dist2xxp1 %*% dist2xxp2
  k_tmp[(l+1):(2*l),1:l2] <- t(t(dist2xxp2) %*% t(dist2xxp1))
  
  if ((pos1>0) & (xp>0)) {
     k_tmp[(l+1):(l+pos1),] <- k_tmp[1:pos1,]
     if (pos2 >0 ){
        k_tmp[,(l2+1):(l2+pos2)] <- k_tmp[,1:pos2]
        k_tmp[(l+1):(l+pos1),(l2+1):(l2+pos2)] <- k_tmp[1:pos1,1:pos2]
     }
  }
  
  kr <- k_tmp
  if (flag){
  for (i in seq(2*l)){
    for (j in seq(2*l2)){
      kr[ind[i],ind[j]] <- k_tmp[i,j]
    }
  }
  }
  
  k00 <- kr[rep(1:l,replicate_no), rep(1:l2,replicate_no2)]
  k01 <- kr[rep(1:l,replicate_no), rep((l2+1):(2*l2),replicate_no2)]
  k10 <- kr[rep((l+1):(2*l),replicate_no), rep(1:l2,replicate_no2)]
  k11 <- kr[rep((l+1):(2*l),replicate_no), rep((l2+1):(2*l2),replicate_no2)]  
  
  
  k_a <- rbind(cbind(k00,k01),cbind(k10,k11))
  
  
  g <- array()
  if ("isNormalised" %in% names(kern) && kern$isNormalised) {
    g[1] <- -0.5*sum(covGrad*k*k_a) +
      0.5 * sum(covGrad*k)/kern$inverseWidth
  }
  else {
    g[1] <- -0.5*sum(covGrad*k*k_a)
  }
  g[2] <- sum(covGrad*k)/kern$variance

  g[3] <- 0
  if ( any(is.nan(g)) )
    warning("g is NaN.")

  return (g)
}


