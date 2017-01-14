#' @title Compute the gradient with respect to the kernel parameters
#' @param kern The DEtime kernel structure for which the gradients are being computed
#' @param X When \code{X2} is provided, X is the input locations associated with the rows of the kernel matrix; when \code{X2} is missing, \code{X} is the input locations associated with both the rows and the columns of the kernel matrix
#' @param X2 The input locations associated with the columnss of the kernel matrix
#' @param covGrad A matrix of partial derivatives of the function of interest with respect to thekernel matrix. the matrix should have the same number of rows as X1 and the same number of columns as \code{X2} has rows   
#' @return
#'    g Gradients of the function of interest with respect to the kernel parameters. The ordering of the vector should match that provided by \code{\link{DEtimeKernExtractParam}}.
#' @description
#'    Compute the gradient with respect to the kernel parameters.
#' @details
#'    \code{g <- kernGradient(kern, X, covGrad)} computes the gradient of functions with respect to the kernel parameters. As well as the kernel structure and the input positions, the user provides a matrix covGrad which gives the partial derivatives of the function with respect to the relevant elements of the kernel matrix.
#'     \code{g <- kernGradient(kern, X1, X2, covGrad)} computes the derivatives as above, but input locations are now provided in two matrices associated with rows and columns of the kernel matrix.
#' @seealso
#'     \code{\link{DEtimeKernCompute}}, \code{\link{DEtimeKernExtractParam}}
#' @examples
#' kern <- list()
#' kern <- DEtimeKernParamInit(kern)
#' X <- matrix(c(seq(0,4),seq(0,4), rep(1,5),rep(2,5)),ncol=2)
#' g <- DEtimeKernGradient(kern, X, array(1, c(10, 10)))
#' @import AtmRay
#' @export

DEtimeKernGradient <-
function (kern, X, X2, covGrad) {
  
  xp <- kern$xp
  x <- matrix(X[,1],ncol=1)
  lfx <- sum(X[,2]==1)
  lx <- dim(x)[1]
  lgx = lx-lfx
  fx <- x[(1:lfx),]
  gx <- x[-(1:lfx),]

  if ( nargs()==3 ) {
    ly <- lx
    lfy <- lfx
    lgy <- lgx
    fy <- x[(1:lfy),]
    gy <- x[-(1:lfy),]
    covGrad <- X2
    x2 <- x
    X2 <- X
  }
  if ( nargs()==4 ) {
    x2 <- matrix(X2[,1],ncol=1)
    ly <- dim(x2)[1]
    lfy <- sum(X2[,2]==1)
    lgy <- ly-lfy
    fy <- x2[(1:lfy),]
    gy <- x2[-(1:lfy),]
  }

  k <- matrix(, nrow = lx, ncol= ly)
  k <- DEtimeKernCompute(kern,X,X2)

  k_a <- matrix(, nrow = lx, ncol= ly)
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

  k00 <- n2ff
  k11 <- m2gxp%*%t(m2gyp)
  k01 <- m2fxp%*%t(m2gyp)
  k10 <- m2gxp%*%t(m2fyp)
 
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

  k01[m_idx1] <- n2fg1[m_idx1]
  k10[m_idx2] <- n2fg2[m_idx2]

  k11[m_idx3] <- n2gg[m_idx3]
  k11[m_idx4] <- n2gg[m_idx4]


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


