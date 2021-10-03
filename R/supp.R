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



plot.CIs <- function(x, 
                    dimn, 
                    CIs, 
                    bord = NULL, 
                    mean.color, 
                    quan.color, 
                    rug, 
                    mean = TRUE, 
                    auto.layout, 
                    ask, ...)
{
  
  pars <- NULL
  if(class(x) == "Smcmc")
  {
    if(is.null(x$varnames)) 
    {
      varnames <- as.character(1:dim(x$stacked)[2])
    }else{
      varnames <- x$varnames
    }
    data <- x$stacked
  } else{
    if(class(x) == "Siid")
    {
      if(is.null(x$varnames)) 
      {
        varnames <- as.character(1:dim(x$data)[2])
      }else{
        varnames <- x$varnames
      }
      data <- x$data
    }
  }
  
  setLayout(dimn)
  for(i in 1:dimn)
  {
    beta = ts(data[, i])
    
    main1 = paste("Density of ",varnames[i])
    plot(density(beta), main = main1, ...)
    if(rug == TRUE) rug(beta, ticksize=0.03, side=1, lwd=0.5)
    main1 = paste("Density of ",varnames[i])
    addCI(x, CIs, component = i, bord = bord, 
          mean.color = mean.color, quan.color = quan.color, 
          mean = mean, ...)
    if(dimn>6)
      if (i == 1)
        pars <- c(pars, par(ask=ask))
  }
  on.exit(par(pars, ask=FALSE,mfrow=c(1,1)))
  
}



## For boxplots
plot.boxx <- function(x, dimn, CIs, quan.color, range, width, varwidth, notch, outline, plot, border,
                      col, ann, horizontal, add,...) 
{
  if(is.null(x$varnames)) 
  {
    varnames <- as.character(1:dim(x$data)[2])
  }else{
    varnames <- x$varnames
  }
  boxplot.matrix(x$data, ..., range = range, width = width, varwidth = varwidth, notch = notch, outline = outline,
                 names=varnames, plot = plot, border = border, col = col, ann = ann, horizontal = horizontal, add = add)
  for(i in 1:dimn) {
    boxCI(x, CI = CIs, component = i,dimn = dimn,  
          quan.color = quan.color, horizontal = horizontal)
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
    b <- b + batchSize(x[[i]])
  }
  b.final <- floor(b/m)
  a <- floor(n/b.final)
  ab <- a*b.final
  trash <- n-ab
  big.chain <- matrix(0,ncol = p, nrow = ab*m)
  if(ab != n)
  {
    for (i in 1:m) {
      big.chain[((i-1)*ab+1):(i*ab),] <- x[[i]][-(1:trash),]
    }
  }else{
    big.chain <- Reduce("rbind", x)
  }
  return(list("b.size" = b.final, "stacked.data" = big.chain))
}

setLayout <- function(p=1, auto.layout = TRUE,pars = NULL)
{
  if(p < 4) {
    par(mfrow = c(dimn,1))
  }
  else if(p == 4) {
    par(mfrow = c(2,2))
  }
  else if(p == 5||p == 6){
    par(mfrow = c(3,2))
  }
  else {
    on.exit(par(pars))
    
    if (auto.layout) {
      mfrow <- set.mfrow(Nparms = 6)
      pars <- par(mfrow = mfrow)
    }
  }  
  
}

acf.Smcmc <- function(x,
                      lag.max     = NULL,
                      type        = c("correlation", "covariance", "partial"),
                      plot        = TRUE,
                      na.action   = na.fail,
                      auto.layout = TRUE,
                      ask         = dev.interactive(),
                      avgcol      = "blue",
                      chain.col   = "red",
                      ...)
{
  
  if(is.null(x))
    stop("data can not be null")
  
  if(!is.list(x))
    stop("x must be a list of chains")
  
  if (is.null(lag.max))
    lag.max <- floor(10 * (log10(n)))
  
  lag.max <- as.integer(min(lag.max, n - 1L))
  m <- length(x)
  dimn <- dim(x[[1]])
  n <- dimn[1]
  p <- dimn[2]
  
  #x <- as.Smcmc(x)
  
  if(plot) 
  {
    pars <- NULL
    setLayout(p)
  }
  
    for(i in 1:p)
    {
      xi <- matrix(data = 0, nrow = n, ncol = m)
      for (j in 1:m) {
        xi[,j] <- x[[j]][,i]
      }
      xi <- xi - mean(xi)
      
      acf <- list()
      
      for(j in 1:m) {
        acf[[j]] <- acf(xi[,j], type = type, plot = FALSE, demean = FALSE, lag.max = lag.max)
      }
      avgf <- acf[[1]]
      for(j in 2:m) {
        avgf$acf <- avgf$acf + acf[[j]]$acf
      }
      avgf$acf <- avgf$acf/m
      
      if(plot)
      {
        plot(avgf, main = i)
        
        lines(0:lag.max, as.matrix(avgf$acf), type = "l", col = adjustcolor(avgcol, alpha.f = .5), lwd = 2, lty = 2, yaxt = 'n', xaxt = 'n')
        
        for(j in 1:m)
        {
          lines(0:lag.max, as.matrix(acf[[j]]$acf), type = "l", col = adjustcolor(chain.col, alpha.f = .5), lwd = 2, lty = 2, yaxt = 'n', xaxt = 'n')
        }
        if(p > 6)
          if (i == 1)
            pars <- c(pars, par(ask=ask))
        
      }
      
    }
    on.exit(par(pars, ask=FALSE,mfrow=c(1,1)))
    
}
