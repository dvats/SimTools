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

## For plotting density plots and Histograms
plot.CIs <- function(x, 
                     dimn, 
                     CIs, 
                     bord = NULL, 
                     mean.color = 'plum4' , 
                     quan.color = 'lightsteelblue3', 
                     rug, 
                     main,
                     mean = TRUE, 
                     auto.layout, 
                     ask, ...)
{
  if(class(x) == "Smcmc")
  {
    if(is.null(x$varnames)) 
    {
      varnames <- as.character(1:dim(x$stacked)[2])
    }else{
      varnames <- x$varnames
    }
    data <- x$stacked
  }else{
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
  if(is.null(main)){par(oma = c(2,0,2,0))}else{par(oma = c(2,0,4,0))}
  setLayout_trace(dimn)
  for(i in 1:dimn)
  {
    ylab = paste("Density of ",i)
    par(mar = c(1.2, 4.1,1.2, 2.1))
    beta = data[, i]
    if (max(abs(beta - floor(beta))) == 0 || bndw(beta) == 0 || length(unique(beta)) == 1)
    {
      beta = as.vector(beta)
      h = hist(beta, plot = F,..) # or hist(x,plot=FALSE) to avoid the plot of the histogram
      h$density = h$counts/sum(h$counts)
      plot(h,freq=FALSE,ylab=ylab, xlab = NA, main = NA, col = "lightgreen")
    }
    
    else
    {
      beta = ts(beta)
      plot(density(beta,...),main = NA, xlab = NA, ylab = ylab)
      if(rug == TRUE) rug(beta, ticksize=0.03, side=1, lwd=0.5)
      addCI(x, CIs, component = i, bord = bord, 
            mean.color = mean.color, quan.color = quan.color, 
            mean = mean)
    }
  }
  if(!is.null(main)){mtext(main, side = 3, line = 1,outer = TRUE, cex = 1.5)}
  on.exit(par(ask = FALSE,mfrow=c(1,1)))
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  par(fig = c(0, 1, 0 , 1))
  par(oma = c(0, 0, 0, 0))
  
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
  
  b.final <- floor(mean(sapply(x, batchSize))) # mean batch size

  a <- floor(n/b.final)
  ab <- a*b.final
  trash <- n-ab
  big.chain <- matrix(0,ncol = p, nrow = ab*m)
  colnames(big.chain) <- colnames(x[[1]])
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

## setLayout for all plots
setLayout_trace <- function(p, ask = FALSE)
{
  if(p%%2 == 0)
  {
    mfrow = c(p/2,2)
  }
  
  else
  {
    mfrow = c((p+1)/2,2)
  }
  par(ask = FALSE, mfrow = mfrow)
  k = list(ask,mfrow)
  names(k) = c("ask","mfrow")
  return(k)
}

## for calculating Bandwidth
bndw <- function(x) {
  x <- x[!is.na(as.vector(x))]
  return(1.06 * min(sd(x), IQR(x)/1.34) * length(x)^-0.2)
}

