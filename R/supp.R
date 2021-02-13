## For simultaneous confidence intervals
Indicator <- function(x, ...)
{
  return(as.numeric(x > ... ))
}

#since our goal is to optimize z, this function returns the Confidence interval
#for a particular z

CIz <- function(z, p1 , p2, theta.hat, phi, ci.sigma.mat, n, mean = TRUE)
{
  se_ests <- sqrt(diag(ci.sigma.mat)/n)
  
  if (mean == FALSE)
  {
    lower.ci.p2 <- phi - z * se_ests
    upper.ci.p2 <- phi + z * se_ests 
    return (list("lower.ci" = lower.ci.p2, "upper.ci" = upper.ci.p2))  
  }
  p <- p1 + p2
  
  lower.ci.p1 <- theta.hat - z * se_ests[1:p1]
  upper.ci.p1 <- theta.hat + z * se_ests[1:p1] 
  lower.ci.p2 <- phi - z * se_ests[(p1+1):p]
  upper.ci.p2 <- phi + z * se_ests[(p1+1):p]
  
  return (list("lower.ci" = c(lower.ci.p1, lower.ci.p2), "upper.ci" = c(upper.ci.p1, upper.ci.p2)))
}

"set.mfrow" <-function (Nparms = 1) 
{
  ## Set up dimensions of graphics window: 
  ##	1 x 1	if Nparms = 1 
  ##	1 X 2 	if Nparms = 2 
  ##	2 X 2 	if Nparms = 3 or 4 
  ##	3 X 2 	if Nparmss = 5 or 6 or 10 - 12 
  ##	3 X 3 	if Nparms = 7 - 9 or >= 13 
    ## One plot per variable
  mfrow <- switch(min(Nparms,13),
                    c(1,1),
                    c(1,2),
                    c(2,2),
                    c(2,2),
                    c(3,2),
                    c(3,2),
                    c(3,3),
                    c(3,3),
                    c(3,3),
                    c(3,2),
                    c(3,2),
                    c(3,2),
                    c(3,3))
  return(mfrow)
}


addCI <- function(data, CIs, component = 1, bord = NA, 
                  mean = TRUE,mean.color = 'plum4', quan.color = 'lightsteelblue3'
                  , opaq = 0.7, ...)
{
  x = ts(data[, component])
  mean.color = adjustcolor(mean.color, alpha.f = opaq)
  quan.color = adjustcolor(quan.color, alpha.f = opaq)
  mn <- CIs$mean.est[component]
  quans <- CIs$xi.q[ , component]
  mcil = CIs$lower.ci.mean[component]
  mciu = CIs$upper.ci.mean[component]
  qcil = CIs$lower.ci.mat[, component]
  qciu = CIs$upper.ci.mat[, component]
  if(mean){
    dum1 <- density(x, from = mcil, to = mciu)
    polygon(c(mcil, dum1$x, mciu), c(0, dum1$y, 0), col = mean.color, border = bord )
  }
  for(j in 1:length(quans))
  {
    dum1 <- density(x, from = qcil[j], to = qciu[j])
    polygon(c(qcil[j], dum1$x, qciu[j]), c(0, dum1$y, 0), col = quan.color, border = bord)
  }
  if(mean){
    segments(mn,0,mn, density(x, from = mn, to = mn, n = 1 )$y,...)
  }
  for(j in 1:length(quans))
  {
    segments(quans[j],0,quans[j], density(x, from = quans[j], to = quans[j], n = 1 )$y,...)
  }
}

plot.CIs <- function(x,dimn, CIs, bord = NULL, mean.color, quan.color, rug,
  mean = TRUE, auto.layout, ask, ...)
{
  pars <- NULL
  if(is.null(attributes(x)$varnames)) 
  {
    varnames <- as.character(1:dim(x)[2])
  }else{
    varnames <- attributes(x)$varnames
  }
  
  if(dimn < 4) {
    par(mfrow = c(dimn,1))
    for(i in 1:dimn)
    {
      beta = ts(x[, i])
      main1 = paste("Density of ",varnames[i])
      plot(density(beta), main = main1, ...)
      if(rug == TRUE) rug(beta, ticksize=0.03, side=1, lwd=0.5)
      addCI(x, CIs, component = i, bord = bord, 
        mean.color = mean.color, quan.color = quan.color, 
        mean = mean, ...)
    }
  }
  else if(dimn == 4) {
    par(mfrow = c(2,2))
    for(i in 1:4)
    {
      beta = ts(x[, i])
      main1 = paste("Density of ",varnames[i])
      plot(density(beta), main = main1, ...)
      if(rug == TRUE) rug(beta, ticksize=0.03, side=1, lwd=0.5)
      addCI(x, CIs, component = i, bord = bord, 
            mean.color = mean.color, quan.color = quan.color, 
            mean = mean, ...)
    }
  }else if(dimn == 5||dimn == 6){
    par(mfrow = c(3,2))
    for(i in 1:dimn)
    {
      beta = ts(x[, i])
      main1 = paste("Density of ",varnames[i])
      plot(density(beta), main = main1, ...)
      if(rug == TRUE) rug(beta, ticksize=0.03, side=1, lwd=0.5)
      addCI(x, CIs,component = i, bord = bord, 
            mean.color = mean.color, quan.color = quan.color, 
            mean = mean, ...)
    }
  }else {
    on.exit(par(pars))
    
    
    if (auto.layout) {
      mfrow <- set.mfrow(Nparms = 6)
      pars <- par(mfrow = mfrow)
    }
    
    for(i in 1:dimn)
    {
      beta = ts(x[, i])
      plot(density(beta), main = main1, ...)
      if(rug == TRUE) rug(beta, ticksize=0.03, side=1, lwd=0.5)
      main1 = paste("Density of ",varnames[i])
      addCI(x, CIs, component = i, bord = bord, 
            mean.color = mean.color, quan.color = quan.color, 
            mean = mean, ...)
      if (i == 1)
        pars <- c(pars, par(ask=ask))
    }
  }
  on.exit(par(pars, ask=FALSE,mfrow=c(1,1)))
  
}

#Q is list of quantiles to be estimated. It will be estimated for each component, m is the length of Q
makeCI <- function(x, Q = c(0.1, 0.9), alpha = 0.05, thresh = 0.001,
                   iid = FALSE, mean = TRUE, b.size = NULL) 
{
  
  if(is.null(b.size)) 
    b.size <- batchSize(x)
  
  mq <- length(Q)
  n <- dim(x)[1]
  
  # p1 is the dimension of g(x)
  # p2 is defined this because p1+p2 will ncols in lambda and sigma
  # v is the vector of all means and quantiles to be estimated
  p1 <- length(x[1,])
  p2 <- mq*length(x[1,])
  theta.hat <- colMeans(x) #g bar
  xi.q <- apply(x, 2, quantile, Q)
  xi.q <- as.matrix(xi.q)
  if(mq==1) xi.q <- t(xi.q)
  #phi is the vector of all quantiles
  phi <- rep(0, p2)
  for(i in 1:mq)
  {
    phi[((i-1)*(p2/mq) + 1):(i*(p2/mq))] = xi.q[i,]
  }
  
  fs <- rep(0, p2)
  for(j in 1:mq)
  {
    for(i in 1:(p2/mq))
    {
      fs[(j-1)*(p2/mq) + i] <- density(x[, i], from = xi.q[j, i], to = xi.q[j, i], n = 1 )$y 
    }
  }
  
  I.flat <- rep(1, p1)
  
  #since p2 was m*dim(h(x))
  lower.ci.mat <- matrix(0, nrow = mq, ncol = p2/mq)
  upper.ci.mat <- matrix(0, nrow = mq, ncol = p2/mq) 
  indis <- matrix(0,nrow = n,ncol = p2)
  for(i in 1:mq)
  {
    
    indi <- (apply(x, 1, Indicator, xi.q[i,]))
    if(p2 > 1)
    {
      indi <- t(indi)
    }
    indis[,((i-1)*(p2/mq) + 1):(i*(p2/mq))] <- indi 
  }
  if(mean == FALSE) Y <- indis else Y <- cbind(x, indis)
  
  if(iid == FALSE)sigma.mat <- mcse.multi(Y, size = b.size)$cov else sigma.mat <- cov(Y)
  
  if(mean == FALSE) lambda <- 1/fs else (lambda <- 1/c(I.flat, fs))
  
  ci.sigma.mat <- (t(t(sigma.mat)*lambda))*lambda
  p <- p1 + p2
  if(mean == FALSE) p = p2
  
  z1 <- qnorm(1 - alpha/2)
  z2 <- qnorm(1 - alpha/(2*p))
  foo1 <- CIz(z1, p1, p2, theta.hat, phi,ci.sigma.mat, n, mean)
  foo2 <- CIz(z2, p1, p2, theta.hat, phi,ci.sigma.mat, n, mean)
  if(mean == FALSE) v <- phi else v <- c(theta.hat, phi)
  
  count <- 0
  prob1 <- pmvnorm(lower = foo1$lower.ci, upper = foo1$upper.ci, mean = v, sigma = (ci.sigma.mat/n))[1]
  prob2 <- pmvnorm(lower = foo2$lower.ci, upper = foo2$upper.ci, mean = v, sigma = (ci.sigma.mat/n))[1]
  
  while(prob2 - prob1 > thresh)
  {
    count <- count + 1
    z.star <- (z1 + z2)/2
    foo.star <- CIz(z.star, p1, p2, theta.hat, phi, ci.sigma.mat, n, mean)
    prob.star <- pmvnorm(lower = foo.star$lower.ci, upper = foo.star$upper.ci, mean = v, sigma = (ci.sigma.mat/n))[1]
    if(prob.star > 1- alpha) 
    {
      z2 <- z.star
      prob2 <- prob.star
    }else
    {
      z1 <- z.star
      prob1 <- prob.star
    }
    if(abs(prob1 - (1 - alpha)) < thresh)
    {
      temp <- CIz(z1, p1, p2, theta.hat, phi,ci.sigma.mat, n, mean)
      break
    }
  }
  temp <- CIz(z1, p1, p2, theta.hat, phi,ci.sigma.mat, n, mean)
  
  for(i in 1:mq)
  {
    if(mean == FALSE)
    {
      lower.ci.mat[i, ] <- temp$lower.ci[((i-1)*(p2/mq)+1):(i*(p2/mq) )]
      upper.ci.mat[i, ] <- temp$upper.ci[((i-1)*(p2/mq) +1):(i*(p2/mq) )]
    }
    else{
      lower.ci.mat[i, ] <- temp$lower.ci[((i-1)*(p2/mq) + p1 + 1):(i*(p2/mq) + p1)]
      upper.ci.mat[i, ] <- temp$upper.ci[((i-1)*(p2/mq) + p1 + 1):(i*(p2/mq) + p1)]
    }
  }
  
  row.names(lower.ci.mat) <- Q
  row.names(upper.ci.mat) <- Q
  if(mean)
  {
    lower.mean <- temp$lower.ci[1:p1]
    upper.mean <- temp$upper.ci[1:p1]
    
    foo3 <- list("lower.ci.mean" = lower.mean, "upper.ci.mean" = upper.mean, "lower.ci.mat" = lower.ci.mat, "upper.ci.mat" = upper.ci.mat, "mean.est" = theta.hat, "xi.q" = xi.q)
  } else foo3 <- list("lower.ci.mat" = lower.ci.mat, "upper.ci.mat" = upper.ci.mat, "mean.est" = theta.hat, "xi.q" = xi.q)
  
  return(foo3)
  
}
#####




## For boxplots
plot.boxx <- function(x, dimn, CIs, mean.color, quan.color, range, width, varwidth, notch, outline, plot, border,
                      col, ann, horizontal, add,...) 
{
  if(is.null(attributes(x)$varnames)) 
  {
    varnames <- as.character(1:dim(x)[2])
  }else{
    varnames <- attributes(x)$varnames
  }
  boxplot.matrix(x, ..., range = range, width = width, varwidth = varwidth, notch = notch, outline = outline,
                 names=varnames, plot = plot, border = border, col = col, ann = ann, horizontal = horizontal, add = add)
  for(i in 1:dimn) {
    boxCI(x, CI = CIs, component = i,dimn = dimn, mean.color=mean.color, 
          quan.color = quan.color, horizontal = horizontal)
  }
}

boxCI <- function(x,CI,components = c(1),dimn = 1,mean.color = 'plum4', quan.color = 'lightsteelblue3',horizontal = FALSE) 
{
  quans = CI$xi.q
  mn = CI$mean.est
  quansi <- quans[, component]
  qcil = CI$lower.ci.mat[, component]
  qciu = CI$upper.ci.mat[, component]
  i <- component
  for(j in 1:length(quansi))
  {
    if(horizontal==TRUE) {
      polygon(c(qcil[j],qcil[j],qciu[j],qciu[j]), c(i-(0.2*min(dimn,2)),i+(0.2*min(dimn,2)),i+(0.2*min(dimn,2)), i-(0.2*min(dimn,2))),  col = quan.color, border = FALSE)
    }else {
      polygon(c(i-(0.2*min(dimn,2)),i+(0.2*min(dimn,2)),i+(0.2*min(dimn,2)), i-(0.2*min(dimn,2))), c(qcil[j],qcil[j],qciu[j],qciu[j]), col = quan.color, border = FALSE)
    } 
  } 
  for(j in 1:length(quansi))
  {
    if(horizontal==TRUE) {
      segments(quansi[j],i-(0.2*min(dimn,2)), quansi[j], i+(0.2*min(dimn,2)))
    }else {
      segments(i-(0.2*min(dimn,2)), quansi[j], i+(0.2*min(dimn,2)), quansi[j])
    }
    
  }
}

chain_stacker <- function(x) {
  m <- length(x)
  
  if(class(x) != "list")
    stop("must be list of chains")
  
  if(is.null(x))
    stop("Chains are null")
  
  n <- as.integer(nrow(x[[1]]))
  p <- as.integer(ncol(x[[1]]))
  
  b <- 0
  for(i in 1:m) {
    b <- b + attributes(x[[i]])$size
  }
  b.final <- floor(b/m)
  a <- floor(n/b.final)
  ab <- a*b.final
  trash <- n-ab
  big.chain <- matrix(0,ncol = p,nrow = ab*m)
  for (i in 1:m) {
    big.chain[((i-1)*ab+1):(i*ab),] <- x[[i]][-(1:trash),]
  }
  
  return(list("b.size" = b.final, "new.data" = big.chain)) 
}
