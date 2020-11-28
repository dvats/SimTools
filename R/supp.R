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

"set.mfrow" <-function (Nparms = 1, nplots = 1) 
{
  ## Set up dimensions of graphics window: 
  ## If only density plots OR trace plots are requested, dimensions are: 
  ##	1 x 1	if Nparms = 1 
  ##	1 X 2 	if Nparms = 2 
  ##	2 X 2 	if Nparms = 3 or 4 
  ##	3 X 2 	if Nparmss = 5 or 6 or 10 - 12 
  ##	3 X 3 	if Nparms = 7 - 9 or >= 13 
  ## If both density plots AND trace plots are requested, dimensions are: 
  ##	1 x 2	if Nparms = 1 
  ##	2 X 2 	if Nparms = 2 
  ##	3 X 2 	if Nparms = 3, 4, 5, 6, 9 
  ##	4 x 2	if Nparms otherwise 
  if (nplots==1) {
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
  }
  else {
    ## Two plot per variable
    ##
    mfrow <- switch(min(Nparms, 13),
                    c(1,2),
                    c(2,2),
                    c(3,2),
                    c(3,2),
                    c(3,2),
                    c(3,2),
                    c(4,2),
                    c(4,2),
                    c(3,2),
                    c(4,2),
                    c(4,2),
                    c(4,2),
                    c(4,2))
  }
  return(mfrow)
}

den.plot <- function(x, mn, quans, mcil,mciu, qcil, qciu, mean = TRUE,bord = NA,mean.color, quan.color,main = "Component-wise Densities",...)
{
  plot(density(x, ...),main = main,...)
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
    segments(mn,0,mn, density(x, from = mn, to = mn, n = 1 )$y)
  }
  for(j in 1:length(quans))
  {
    segments(quans[j],0,quans[j], density(x, from = quans[j], to = quans[j], n = 1 )$y)
  }
}

plot.CIs <- function(x,dimn, CIs, bord = NULL, mean.color, quan.color, mean = TRUE, mn, quans,auto.layout,ask,...)
{
  if(dimn < 4) {
    par(mfrow = c(dimn,1))
    for(i in 1:dimn)
    {
      beta = ts(x[, i])
      den.plot(x = beta, mn = mn[i], quans = quans[, i], mcil = CIs$lower.ci.mean[i], mciu = CIs$upper.ci.mean[i], qcil = CIs$lower.ci.mat[, i],qciu = CIs$upper.ci.mat[, i],bord = bord,mean.color = mean.color, quan.color = quan.color,mean = mean,...)
    }
  }
  else if(dimn == 4) {
    par(mfrow = c(2,2))
    for(i in 1:4)
    {
      beta = ts(x[, i])
      den.plot(x = beta, mn = mn[i], quans = quans[, i], mcil = CIs$lower.ci.mean[i], mciu = CIs$upper.ci.mean[i], qcil = CIs$lower.ci.mat[, i],qciu = CIs$upper.ci.mat[, i],bord = bord,mean.color = mean.color, quan.color = quan.color,mean = mean,...)
    }
  }else if(dimn == 5||dimn == 6){
    par(mfrow = c(3,2))
    for(i in 1:dimn)
    {
      beta = ts(x[, i])
      den.plot(x = beta, mn = mn[i], quans = quans[, i], mcil = CIs$lower.ci.mean[i], mciu = CIs$upper.ci.mean[i], qcil = CIs$lower.ci.mat[, i],qciu = CIs$upper.ci.mat[, i],bord = bord,mean.color = mean.color, quan.color = quan.color,mean = mean,...)
    }
  }else {
    pars <- NULL
    on.exit(par(pars))
    
    
    if (auto.layout) {
      mfrow <- set.mfrow(Nparms = 6, nplots = 1)
      pars <- par(mfrow = mfrow)
    }
    
    for(i in 1:dimn)
    {
      beta = ts(x[, i])
      den.plot(x = beta, mn = mn[i], quans = quans[, i], mcil = CIs$lower.ci.mean[i], mciu = CIs$upper.ci.mean[i], qcil = CIs$lower.ci.mat[, i],qciu = CIs$upper.ci.mat[, i],bord = bord,mean.color = mean.color, quan.color = quan.color,mean = mean,...)
      if (i == 1)
        pars <- c(pars, par(ask=ask))
    }
  }
  
}

#Q is list of quantiles to be estimated. It will be estimated for each component, m is the length of Q
error.est <- function(x, Q = c(0.1, 0.9), alpha = 0.05, thresh, mean = TRUE,iid, ...) 
{
  (m = length(Q))
  x <- as.matrix(x)
  (n <- length(x[, 1]))
  if(!is.matrix(x) && !is.data.frame(x))
    stop("'x' must be a matrix or data frame.")
  #p1 is the dimension of g(x)
  #p2 is defined this because p1+p2 will ncols in lambda and sigma
  #v is the vector of all means and quantiles to be estimated
  (p1 <- length(x[1,]))
  (p2 <- m*length(x[1,]))
  (theta.hat <- colMeans(x)) #g bar
  (xi.q <- apply(x, 2, quantile, Q))
  
  #phi is the vector of all quantiles
  phi <- rep(0, p2)
  for(i in 1:m)
  {
    phi[((i-1)*(p2/m) + 1):(i*(p2/m))] = xi.q[i,]
  }
  
  fs <- rep(0, p2)
  for(j in 1:m)
  {
    for(i in 1:(p2/m))
    {
      (fs[(j-1)*(p2/m) + i] <- density(x[, i], from = xi.q[j, i], to = xi.q[j, i], n = 1 )$y )
    }
  }
  
  (I.flat <- rep(1, p1))
  
  #since p2 was m*dim(h(x))
  lower.ci.mat <- matrix(0, nrow = m, ncol = p2/m)
  upper.ci.mat <- matrix(0, nrow = m, ncol = p2/m)
  indis <- matrix(0,nrow = n,ncol = p2)
  for(i in 1:m)
  {
    
    (indi <- (apply(x, 1, Indicator, xi.q[i,])))
    if(p2 > 1)
    {
      indi <- t(indi)
    }
    indis[,((i-1)*(p2/m) + 1):(i*(p2/m))] <- indi 
  }
  if(mean == FALSE) Y <- indis else Y <- cbind(x, indis)
  
  if(iid == FALSE)sigma.mat <- mcse.multi(Y)$cov else sigma.mat <- cov(Y)
  
  if(mean == FALSE) lambda <- 1/fs else (lambda <- 1/c(I.flat, fs))
  
  (ci.sigma.mat <- (t(t(sigma.mat)*lambda))*lambda)
  (p <- p1 + p2)
  if(mean == FALSE) p = p2
  (z1 <- qnorm(1 - alpha/2))
  (z2 <- qnorm(1 - alpha/(2*p)))
  (foo1 <- CIz(z1, p1, p2, theta.hat, phi,ci.sigma.mat, n, mean))
  (foo2 <- CIz(z2,p1, p2, theta.hat, phi,ci.sigma.mat, n, mean))
  if(mean == FALSE) v <- phi else (v <- c(theta.hat, phi))
  
  count <- 0
  (prob1 <- pmvnorm(lower = foo1$lower.ci, upper = foo1$upper.ci, mean = v, sigma = (ci.sigma.mat/n))[1])
  (prob2 <- pmvnorm(lower = foo2$lower.ci, upper = foo2$upper.ci, mean = v, sigma = (ci.sigma.mat/n))[1])
  
  while(prob2 - prob1 > thresh)
  {
    count <- count + 1
    (z.star <- (z1 + z2)/2)
    (foo.star <- CIz(z.star, p1, p2, theta.hat, phi, ci.sigma.mat, n, mean))
    (prob.star <- pmvnorm(lower = foo.star$lower.ci, upper = foo.star$upper.ci, mean = v, sigma = (ci.sigma.mat/n))[1])
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
  (temp <- CIz(z1, p1, p2, theta.hat, phi,ci.sigma.mat, n, mean))
  
  for(i in 1:m)
  {
    if(mean == FALSE)
    {
      lower.ci.mat[i, ] <- temp$lower.ci[((i-1)*(p2/m)+1):(i*(p2/m) )]
      upper.ci.mat[i, ] <- temp$upper.ci[((i-1)*(p2/m) +1):(i*(p2/m) )]
    }
    else{
      lower.ci.mat[i, ] <- temp$lower.ci[((i-1)*(p2/m) + p1 + 1):(i*(p2/m) + p1)]
      upper.ci.mat[i, ] <- temp$upper.ci[((i-1)*(p2/m) + p1 + 1):(i*(p2/m) + p1)]
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

plot.Smcmc <- function(x, Q = c(0.1, 0.9), alpha = 0.05,thresh = 0.001, iid = FALSE, plot = TRUE, auto.layout = TRUE, ask = dev.interactive(),mean = TRUE, border = NA, mean.col = 'plum4', quan.col = 'lightsteelblue3', opaq = 0.7, ...)
{
  foo3 <- error.est(x, Q, alpha,thresh = thresh,iid = iid,mean = mean, ...)
  if(plot == TRUE)
  {
    plot.CIs(x, dimn = length(x[1,]), CIs = foo3, bord = border, mean.color = adjustcolor(mean.col, alpha.f = opaq), quan.color = adjustcolor(quan.col, alpha.f = opaq), mn = foo3$mean.est,mean = mean, quans = foo3$xi.q,auto.layout = auto.layout,ask = ask,...)
  }
  return(foo3)
}
