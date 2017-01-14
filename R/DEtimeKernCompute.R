#' @title Compute the DEtime kernel given the parameters and X
#' @param kern DEtime kernel structure to be computed
#' @param X A two-column matrix where the first column of this matrix is the time points for control and perturbed conditions and the second column uses '1' to represent time points from control condition and '2' to represent time points from perturbed condition
#' @param X2 Second input matrix to the kernel computation (forms the columns of the kernel)
#' @return 
#'    DEtime kernel structure.
#' @description
#'    Compute the DEtime kernel given the parameters and X.
#' @details
#'    \code{K <- DEtimeKernCompute(kern, X)} computes a DEtime kernel matrix given an input data matrix.
#'    \code{K <- DEtimeKernCompute(kern, X, X2)} computes a DEtime kernel matrix for the given kernel type given two input data matrices, one for the rows and one for the columns.
#' @examples
#' kern <- list()
#' kern <- DEtimeKernParamInit(kern)
#' X <- matrix(c(seq(0,4),seq(0,4), rep(1,5),rep(2,5)),ncol=2)
#' K <- DEtimeKernCompute(kern, X)
#' @import AtmRay
#' @export

DEtimeKernCompute <-
function (kern,X,X2) {  
  xp <- kern$xp
  x <- matrix(X[,1],ncol=1)
  lfx <- sum(X[,2]==1)
  lx <- dim(x)[1]
  lgx = lx-lfx
  fx <- x[(1:lfx),]
  gx <- x[-(1:lfx),]
  if ( nargs() < 3 ) {
    x2 <- x
    ly <- lx
    lfy <- lfx
    lgy <- lgx
    fy <- x[(1:lfy),]
    gy <- x[-(1:lfy),]
  }
  else {
    x2 <- matrix(X2[,1],ncol=1)
    ly <- dim(x2)[1]
    lfy <- sum(X2[,2]==1)
    lgy <- ly - lfy
    fy <- x2[(1:lfy),]
    gy <- x2[-(1:lfy),]
  }

  wi2 <- 0.5*kern$inverseWidth
  k <- matrix(, nrow = lx, ncol= ly)
  k00 <- matrix(, nrow = lfx, ncol= lfy)
  k01 <- matrix(, nrow = lfx, ncol= lgy)
  k10 <- matrix(, nrow = lgx, ncol= lfy)
  k11 <- matrix(, nrow = lgx, ncol= lgy)
  
  n2ff <- .dist2(fx,fy)
  n2gg <- .dist2(gx,gy)
  n2fg1 <- .dist2(fx,gy)
  n2fg2 <- .dist2(gx,fy)
  m2fxp <- .dist2(fx,xp)
  m2fyp <- .dist2(fy,xp)
  m2gxp <- .dist2(gx,xp)
  m2gyp <- .dist2(gy,xp)
    
  k00 <- exp(-n2ff*wi2)
  k11 <- exp(-wi2*m2gxp)%*%t(exp(-wi2*m2gyp))
  k01 <- exp(-wi2*m2fxp)%*%t(exp(-wi2*m2gyp))
  k10 <- exp(-wi2*m2gxp)%*%t(exp(-wi2*m2fyp)) 
  
  idx_gx <- which(gx<xp)
  idx_gy <- which(gy<xp)
 
  idx_gx1 <- which(gx>=xp)
  idx_gy1 <- which(gy>=xp)
 
  idx1 <- meshgrid(1:lfx,idx_gy)
  idx2 <- meshgrid(idx_gx,1:lfy)
  idx3 <- meshgrid(idx_gx,idx_gy)
  idx4 <- meshgrid(idx_gx1,idx_gy1)
  

  m_idx1 <- matrix(c(idx1$x,idx1$y),ncol=2)
  m_idx2 <- matrix(c(idx2$x,idx2$y),ncol=2)
  m_idx3 <- matrix(c(idx3$x,idx3$y),ncol=2)
  m_idx4 <- matrix(c(idx4$x,idx4$y),ncol=2)
  
  k01[m_idx1] <- exp(-n2fg1[m_idx1]*wi2)
  k10[m_idx2] <- exp(-n2fg2[m_idx2]*wi2)

  k11[m_idx3] <- exp(-n2gg[m_idx3]*wi2)
  k11[m_idx4] <- exp(-n2gg[m_idx4]*wi2)
   
  k <- rbind(cbind(k00,k01),cbind(k10,k11))
  
  k <- kern$variance*k
  if ("isNormalised" %in% names(kern) && kern$isNormalised)
    k <- k * sqrt(kern$inverseWidth/(2*pi))

  return(k)

}

