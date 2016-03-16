#' @title Compute the DEtime kernel given the parameters and X
#' @param kern DEtime kernel structure to be computed
#' @param x Depending on the number of inputs, x can be the input data matrix (rows are data points) to the kernel computation, or the first input matrix to the kernel computation (forms the rows of the kernel)
#' @param x2 Second input matrix to the kernel computation (forms the columns of the kernel)
#' @return 
#'    DEtime kernel structure.
#' @description
#'    Compute the DEtime kernel given the parameters and X.
#' @details
#'    \code{K <- DEtimeKernCompute(kern, x)} computes a DEtime kernel matrix given an input data matrix.
#'    \code{K <- DEtimeKernCompute(kern, x1, x2)} computes a DEtime kernel matrix for the given kernel type given two input data matrices, one for the rows and one for the columns.
#' @examples
#' kern <- list()
#' kern <- DEtimeKernParamInit(kern)
#' K <- DEtimeKernCompute(kern, as.matrix(3:8))
#' @export

DEtimeKernCompute <-
function (kern, x, x2) {  
  xp <- kern$xp
  if (dim(x)[1]==1){ x<-t(x)} 
  replicate_no <- length(which(x==x[1]))/2
  if (replicate_no<1){
      replicate_no <- 1
  }
  l <- dim(x)[1]/(2*replicate_no)
  x0 <- x[1:l]

  if ( nargs() < 3 ) {
    replicate_no2 <- replicate_no
    n2 <- .dist2(x0,x0)
    m21 <- .dist2(x0,xp)
    m22 <- .dist2(xp,x0)
    
    l2 <- l
    if (xp>=max(x0)) { 
      pos1 <- length(x0)}
    else {
      pos1 <- max(min(which((x0<=xp)==FALSE))-1, 1)}
    
    pos2 <- pos1

    flag <- TRUE
    
  } else {
    
    replicate_no2 <- length(which(x2==x2[1]))/2   
    l2 <- dim(x2)[1]/(2*replicate_no2)
    x1 <- x2[1:l2]
    sorted_x <- sort(c(x0,x1),index.return=TRUE)
    xnew <- sorted_x$x
    ind <- sorted_x$ix
    n2 <- .dist2(xnew,xnew)
    m21 <- .dist2(xnew,xp)
    m22 <- .dist2(xp,xnew)
    if (xp >= max(xnew)) { pos1 <- length(xnew)}
    else {
    pos1 <- max(min(which((xnew<=xp)==FALSE))-1)}
    pos2 <- pos1
    flag <- FALSE
  }
  
  wi2 <- 0.5*kern$inverseWidth
  k0 <- matrix(, nrow = 2*l, ncol= 2*l2)
  
  
  if (flag){
    
    k_0 <- matrix(, nrow = 2*l, ncol= 2*l2)
    k_0[1:l,1:l2] <- kern$variance*exp(-n2*wi2)
    k_0[(l+1):(2*l),(l2+1):(2*l2)] <- k_0[1:l,1:l2]
    k_0[1:l,(l2+1):(2*l2)] <- kern$variance*exp(-wi2*m21)%*%exp(-wi2*m22)
    k_0[(l+1):(2*l),1:l2] <- k_0[1:l,(l2+1):(2*l2)]
    
    #browser()  
    k0 <- k_0
    if ((pos1 > 0) & (xp>0)) {
    
      k0[(l+1):(l+pos1),] <- k_0[1:pos1,]
      if (pos2 > 0){
        k0[,(l+1):(l+pos2)] <- k_0[,1:pos2]
        k0[(l+1):(l+pos1),(l+1):(l+pos2)] <- k_0[1:pos1,1:pos2]
      }
    }
    
    }
  else {
    
    l0 <- l+l2
    k_0 <- matrix(, nrow = 2*l0, ncol= 2*l0)
    k_0[1:l0,1:l0] <- kern$variance*exp(-n2*wi2)
    k_0[(l0+1):(2*l0),(l0+1):(2*l0)] <- k_0[1:l0,1:l0]
    k_0[1:l0,(l0+1):(2*l0)] <- kern$variance*exp(-wi2*m21)%*%exp(-wi2*m22)
    k_0[(l0+1):(2*l0),1:l0] <- kern$variance*t(exp(-wi2*m22))%*%t(exp(-wi2*m21))
    
    k_1 <- k_0
    if ((pos1 > 0) & (xp>0)) {
      k_1[(l0+1):(l0+pos1),] <- k_0[1:pos1,]
      if (pos2 > 0){
        k_1[,(l0+1):(l0+pos2)] <- k_0[,1:pos2]
        k_1[(l0+1):(l0+pos1),(l0+1):(l0+pos2)] <- k_0[1:pos1,1:pos2]
      }
    }
    k0 <- matrix(0,nrow=(2*l),ncol=(2*l2))
    
    k_r <- matrix(0, nrow=2*l0,ncol=2*l0)

    
    for (i in seq(l0)){
      for (j in seq(l0)){
        k_r[ind[i],ind[j]] <- k_1[i,j]
      }
    }
    
    for (i in seq((l0+1),(2*l0))){
      for (j in seq(l0)){
        k_r[ind[i-l0]+l0,ind[j]] <- k_1[i,j]
      }
    }
    
    for (i in seq(l0)){
      for (j in seq((l0+1),(2*l0))){
        k_r[ind[i],ind[j-l0]+l0] <- k_1[i,j]
      }
    }
    
    for (i in seq((l0+1),(2*l0))){
      for (j in seq((l0+1),(2*l0))){
        k_r[ind[i-l0]+l0,ind[j-l0]+l0] <- k_1[i,j]
      }
    }
    
    
    k0[1:l,1:l2] <- k_r[1:l,(l+1):(l+l2)]
    k0[1:l,(l2+1):(2*l2)] <- k_r[1:l,(l0+l+1):(l0+l+l2)]
    k0[(l+1):(2*l),1:l2] <- k_r[(l0+1):(l0+l),(l+1):(l+l2)]
    k0[(l+1):(2*l),(l2+1):(2*l2)] <- k_r[(l0+1):(l0+l),(l0+l+1):(l0+l+l2)]

    
   
  }

 
  k00 <- k0[rep(1:l,replicate_no), rep(1:l2,replicate_no2)]
  k01 <- k0[rep(1:l,replicate_no), rep((l2+1):(2*l2),replicate_no2)]
  k10 <- k0[rep((l+1):(2*l),replicate_no), rep(1:l2,replicate_no2)]
  k11 <- k0[rep((l+1):(2*l),replicate_no), rep((l2+1):(2*l2),replicate_no2)]  
  
  k <- rbind(cbind(k00,k01),cbind(k10,k11))
  
  if ("isNormalised" %in% names(kern) && kern$isNormalised)
    k <- k * sqrt(kern$inverseWidth/(2*pi))

  return(k)

}

