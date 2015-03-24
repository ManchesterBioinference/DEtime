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
    k_0[1:l,(l2+1):(2*l2)] <- exp(-wi2*m21)%*%exp(-wi2*m22)/kern$variance
    k_0[(l+1):(2*l),1:l2] <- k_0[1:l,(l2+1):(2*l2)]
      
    k0 <- k_0
    if (pos1 > 0) {
    
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
    k_0[1:l0,(l0+1):(2*l0)] <- exp(-wi2*m21)%*%exp(-wi2*m22)/kern$variance
    k_0[(l0+1):(2*l0),1:l0] <- t(exp(-wi2*m22))%*%t(exp(-wi2*m21))/kern$variance
    
    k_1 <- k_0
    if (pos1 > 0) {
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
DEtimeKernDiagCompute <-
function (kern, x) {
  #k <- diag(DEtimeKernCompute(kern,x))
  k <- matrix(kern$variance, dim(as.array(x))[1], 1)
  
  if ("isNormalised" %in% names(kern) && kern$isNormalised)
    k <- k * sqrt(kern$inverseWidth/(2*pi))

  return (k)
}
DEtimeKernExpandParam <-
function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  kern$inverseWidth <- params[1]	## linear domain params, i.e. untransformed inverse-width and signal variance
  kern$variance <- params[2]
  kern$xp <- params[3]
  
  return (kern)
}

DEtimeKernExtractParam <-
function (kern, only.values=TRUE,
                                 untransformed.values=TRUE) {
  params <- c(kern$inverseWidth, kern$variance, kern$xp)

  if ( !only.values )
    names(params) <- c("inverseWidth", "variance","xp")

  return (params)
}

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
  
  if (pos1>0 ) {
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

DEtimeKernParamInit <-
function (kern) {
  kern$xp <- -1000
  kern$inverseWidth <- 1
  kern$variance <- 1
  kern$nParams <- 3
  kern$paramNames <- c("inverseWidth", "variance", "xp")
  
  kern$isStationary <- TRUE

  if ("options" %in% names(kern) && "isNormalised" %in% names(kern$options) && kern$options$isNormalised)
    kern$isNormalised <- TRUE
  else
    kern$isNormalised <- FALSE

  if ("options" %in% names(kern) && "inverseWidthBounds" %in% names(kern$options)) {
    kern$transforms <- list(list(index=1, type="bounded"),
                            list(index=2, type="positive"))
    kern$transformArgs <- list()
    kern$transformArgs[[1]] <- kern$options$inverseWidthBounds
    kern$inverseWidth <- mean(kern$options$inverseWidthBounds)
  } else if ("options" %in% names(kern) && "varianceBounds" %in% names(kern$options)) {
	  kern$transforms <- list(list(index=1, type="positive"),
							  list(index=2, type="bounded"))
	  kern$transformArgs <- list()
	  kern$transformArgs[[2]] <- kern$options$varianceBounds
	  kern$variance <- mean(kern$options$varianceBounds)
  } else if ("options" %in% names(kern) && "inverseWidthBounds" %in% names(kern$options) && "varianceBounds" %in% names(kern$options)) {
	  kern$transforms <- list(list(index=1, type="bounded"),
							  list(index=2, type="bounded"))
	  kern$transformArgs <- list()
	  kern$transformArgs[[1]] <- kern$options$inverseWidthBounds
	  kern$transformArgs[[2]] <- kern$options$varianceBounds
	  kern$inverseWidth <- mean(kern$options$inverseWidthBounds)
	  kern$variance <- mean(kern$options$varianceBounds)
  } else {
	  kern$transforms <- list(list(index=c(1,2), type="positive"))
  }

  return (kern)
}


.gpPlot_DEtime <-
function(model,Xstar,mu,S,likelihood,simpose=NULL,xlim=NULL,ylim=NULL,xlab='',ylab='',col='blue',title='') {
    ## GPPLOT Plots the GP mean and variance.
    
    if (missing(model) || missing(Xstar)) {
        stop('Missing GP model or points of prediction Xstar.')
    } else {
        if (missing(mu) || missing(S)) {
            meanVar = gpPosteriorMeanVar(model, Xstar, varsigma.return=TRUE)
            mu = meanVar$mu; S = meanVar$varsigma
        }
    }
    
    lstar <- dim(Xstar)[1]/2
    l <- dim(model$X)[1]/2
    #f = c(mu+2*sqrt(abs(S)), rev(mu-2*sqrt(abs(S))))
    f1  = c(mu[1:lstar]+2*sqrt(abs(S[1:lstar])), rev(mu[1:lstar]-2*sqrt(abs(S[1:lstar]))))
    f2 = c(mu[(lstar+1):(2*lstar)]+2*sqrt(abs(S[(lstar+1):(2*lstar)])), rev(mu[(lstar+1):(2*lstar)]-2*sqrt(abs(S[(lstar+1):(2*lstar)]))))
    
    if (is.null(xlim))
    xlim = range((model$X))
    if (is.null(ylim))
    ylim = range(c(f1,f2))
    
    #   par(pty="s")
    old_par <- par(no.readonly=TRUE)
    #attach(mtcars)
    #layout(2,1, widths=4,heights=c(1,3))
    #par(mar=c(0,0,0,0))
    par(fig=c(0.1,1.0,0.8,1),mar=c(0.0,2.0,2.0,2.0))
    plot(likelihood, type="n", xaxt="n", yaxt="n", xlim=xlim, bty="l")
    xx <- c(likelihood[,1],rev(likelihood[,1]))
    yy <- c(rep(0,nrow(likelihood)), rev(likelihood[,2]))
    polygon(xx,yy, col=rgb(0,0,1,0.4), border=NA, ylim=range(xx))
    par(fig=c(0.1,1.0,0.1,0.8), mar=c(2.0,2.0,0.0,2.0),new=TRUE)
    plot(1, type="n", xlim=xlim, ylim=ylim, cex.axis=.5,cex.lab=.5, cex.main=.5, cex.sub = .5, bty="l") ## Empty plot basis.
    
    if (col=='blue') shade = rgb(0,0,1,alpha=.1)
    else if (col=='red') shade = rgb(255,0,0,alpha=.1)
    else shade = 'gray'
    
    lstar <- dim(Xstar)[1]/2
    l <- dim(model$X)[1]/2
    
    replicate_no <- length(which(model$X==model$X[1]))/2
    size_x <- length(model$X)/(2*replicate_no)
    mu_f <- matrix(0,nrow = size_x, ncol =1)
    mu_g <- matrix(0,nrow = size_x, ncol =1)
    
    for (i in seq(replicate_no)){
        mu_f <- mu_f + model$y[((i-1)*size_x+1):(i*size_x)]
        mu_g <- mu_g + model$y[((i-1)*size_x+1+replicate_no*size_x):(i*size_x+replicate_no*size_x)]
    }
    mu_f <- mu_f/replicate_no
    mu_g <- mu_g/replicate_no
    
    polygon(c(Xstar[1:lstar,], rev(Xstar[1:lstar,])), f1, col = rgb(0.9,1.0,1.0,1.0), border = shade)	## Confidence intervals.
    polygon(c(Xstar[(lstar+1):(2*lstar),], rev(Xstar[(lstar+1):(2*lstar),])), f2, col = rgb(1,0.95,0.7,.8), border = shade)	## Confidence intervals.
    points(model$X[1:l,], model$y[1:l,], pch = 0, cex = .5, lwd=.5, col = 'blue')	## Training points.
    points(model$X[(l+1):(2*l),], model$y[(l+1):(2*l),], pch = 4, cex = .5, lwd=.5, col = 'red')	## Training points.
    lines(Xstar[1:lstar,], mu[1:lstar,], col='blue', lwd=.5, lty = 1)	## Mean function.
    lines(Xstar[(lstar+1):(2*lstar),], mu[(lstar+1):(2*lstar),], col='red', lwd=.5, lty = 1)	## Mean function.
    #lines(model$X[1:size_x,], mu_f[1:size_x,], col='black', lwd=.5, lty = 2)  ## Mean function
    #lines(model$X[1:size_x,], mu_g[1:size_x,], col='blue', lwd=.5, lty = 2)  ## Mean function .
    #par(mar=c(0, 0, 0, 0))
    #legend("bottom", c("Original control", "Original perturbed", "Mean of control", "Mean of perturbed", "Estimated control", "Estimated perturbed"), col = c('black', 'blue', 'black','darkblue','black','blue'),text.col = "black", lty = c(NA,NA,2, 2, 1,1), pch = c(0,1, NA,NA,NA,NA), bg = "gray90", ncol=3,bty="n", cex=0.8)
    
    legend("bottom", c("Original control", "Original perturbed", "Estimated control", "Estimated perturbed"), col = c('blue', 'red','blue','red'),text.col = "black", lty = c(NA,NA,1,1), pch = c(0,4, NA,NA), bg = "gray90", ncol=3,bty="n", cex=0.8)
    
    mtext(title, side=3, line=3, cex=1, col="black")
    mtext("Gene expression", side=2, line=2, padj=0.5, cex=1, col="black")
    mtext("Time", side=1, line=2, cex=1, col="black")
    if (!is.null(simpose)) {
        y = mu[simpose] + rnorm(6, 0, exp(model$params$xmin[3]/2))
        points(simpose, y, pch = 4, cex = 0.5, lwd=3, col = col)
    }
    par(old_par)
    #zeroAxes()
}

.jitCholInv <-
function ( M, Num=10, silent=FALSE ) {
  jitter <- 0
  jitter1 <- abs(mean(diag(M)))*1e-6
  eyeM <- diag( 1, nrow=length(M[,1]), ncol=length(M[1,]) )

  for ( i in 1:Num ) {

    ## clear the last error message
    try(stop(""),TRUE)

    Ch <- try( chol( M + jitter*eyeM ), silent=TRUE )

    nPos <- grep("not positive definite",  geterrmessage())

    if ( length(nPos) != 0 ) {
      jitter1 <- jitter1*10
      jitter <- jitter1

      if (! silent) {
        warnmsg <- paste("Matrix is not positive definite, adding",
                         signif(jitter,digits=4), "jitter!")
        warning(warnmsg)
      }
    }
    else break
  }

  invCh <- try (solve( Ch, eyeM ), silent=TRUE)

  if ( class(invCh) == "try-error" ) {
    return (NaN)
  }
  else {
    invM <- invCh %*% t(invCh)

    if ( jitter == 0 ) {
      ans <- list(invM=invM, jitter=jitter, chol=Ch)
    }
    else ans <- list(invM=invM, jitM=M+jitter*eyeM , jitter=jitter, chol=Ch)

    return (ans)
  }
}

.dist2 <-
function (x, x2) {
  xdim <- dim(as.matrix(x))
  x2dim <- dim(as.matrix(x2))

  xMat <- array(apply(as.matrix(x*x),1,sum), c(xdim[1], x2dim[1]))
  x2Mat <- t(array(apply(as.matrix(x2*x2),1,sum), c(x2dim[1], xdim[1])))

  if ( xdim[2] != x2dim[2] )
    stop("Data dimensions are not matched.")

  n2 <-   xMat+x2Mat-2*tcrossprod(x, x2)

  return (n2)
}


