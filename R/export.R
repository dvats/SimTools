#' @title Add simultaneous confidence interval to existing plot.
#'
#' @description Adds simultaneous confidence intervals for quantiles and means to an existing plot.
#'
#' @name addCI
#' @aliases addCI
#' @usage addCI(x, CIs, component = 1, bord = NA, mean = TRUE, mean.color = 'plum4',
#'                quan.color = 'lightsteelblue3', opaq = 0.7, ...)
#' @param x : a `Smcmc' class object
#' @param CIs : the output from the `getCI` function
#' @param component : numeric indicating which component to draw the confidence intervals for
#' @param bord : logical for whether a border is desired around the confidence intervals
#' @param mean : logical argument whether the mean is to be plotted
#' @param mean.color : color for the mean confidence interval
#' @param quan.color : color for the quantile confidence intervals
#' @param opaq : opacity of \code{mean.col} and \code{quan.col}. A value of 0 is transparent and 1 is completely opaque.
#' @param ... : arguments passed on to the boundaries of the confidence intervals in `segments`
#' @return adds segments for confidence intervals into an already existing plot environment
#'
#' @examples
#' chain <- matrix(0, ncol = 1, nrow = 1e3)
#' chain[1,] <- 0
#' err <- rnorm(1e3)
#' for(i in 2:1e3)
#' {
#'   chain[i,] <- .3*chain[i-1,] + err[i]
#' }
#' chain <- Smcmc(list(chain))
#' plot(density(chain$stacked[,1]))
#' CIs <- getCI(chain)
#' addCI(chain, CIs, component = 1)
#' @export
addCI <- function(x, 
                  CIs, 
                  component   = 1, 
                  bord        = NA, 
                  mean        = TRUE,
                  mean.color  = 'plum4', 
                  quan.color  = 'lightsteelblue3',
                  opaq        = 0.7, ...)
{
  if(class(x) == "Smcmc")  obj <- ts(x$stacked[, component])
  if(class(x) == "Siid")  obj <- ts(x$data[, component])
  mean.color = adjustcolor(mean.color, alpha.f = opaq)
  quan.color = adjustcolor(quan.color, alpha.f = opaq)
  mn <- CIs$mean.est[component]
  quans <- CIs$xi.q[ , component]
  mcil = CIs$lower.ci.mean[component]
  mciu = CIs$upper.ci.mean[component]
  qcil = CIs$lower.ci.mat[, component]
  qciu = CIs$upper.ci.mat[, component]
  if(mean){
    dum1 <- density(obj, from = mcil, to = mciu)
    polygon(c(mcil, dum1$x, mciu), c(0, dum1$y, 0), col = mean.color, border = bord )
  }
  for(j in 1:length(quans))
  {
    dum1 <- density(obj, from = qcil[j], to = qciu[j])
    polygon(c(qcil[j], dum1$x, qciu[j]), c(0, dum1$y, 0), col = quan.color, border = bord)
  }
  if(mean){
    segments(mn,0,mn, density(obj, from = mn, to = mn, n = 1 )$y,...)
  }
  for(j in 1:length(quans))
  {
    segments(quans[j],0,quans[j], density(obj, from = quans[j], to = quans[j], n = 1 )$y,...)
  }
}


#' @title Calculates simultaneous confidence intervals.
#'
#' @description Calculates simultaneous confidence intervals for means and
#'              quantiles as indicated for the desired MCMC output
#'
#' @name getCI
#' @aliases getCI
#' @usage getCI(x, Q = c(0.1, 0.9), alpha = 0.05, thresh = 0.001, iid = FALSE,
#'             mean = TRUE)
#' @param x : a `Smcmc' class object
#' @param Q : vector of quantiles 
#' @param alpha : confidence levels of the simulatenous intervals
#' @param thresh : threshold for the optimization methodology that calculates the simultaneous CIs
#' @param iid : logical argument for constructing density plot for iid samples. Defaults to \code{FALSE}
#' @param mean : logical  indicating whether mean is to be plotted
#' @return adds segments for confidence intervals into an already existing plot environment
#'
#' @examples
#' chain <- matrix(0, ncol = 1, nrow = 1e3)
#' chain[1,] <- 0
#' err <- rnorm(1e3)
#' for(i in 2:1e3)
#' {
#'   chain[i,] <- .3*chain[i-1,] + err[i]
#' }
#' chain <- Smcmc(list(chain))
#' plot(density(chain$stacked[,1]))
#' CIs <- getCI(chain)
#' addCI(chain, CIs, component = 1)
#'
#' @references
#' Robertson, N., Flegal, J. M., Vats, D., and Jones, G. L., 
#' “Assessing and Visualizing Simultaneous Simulation Error”, 
#' Journal of Computational and Graphical Statistics,  2020. 
#'
#' @export
getCI <- function(x,
                  Q      = c(0.1, 0.9), 
                  alpha  = 0.05, 
                  thresh = 0.001,
                  iid    = FALSE, 
                  mean   = TRUE) 
{
  
  if(class(x) == "Smcmc")
  {
    if(is.null(x$size)) {
      b.size <- 0
      for (i in length(x$chains)) {
        b.size <- b.size + batchSize(x$chains[[i]])
      }
      b.final <- floor(b.size/length(x$chains))
    } 
    else b.size <- x$size
    
    x <- x$stacked
  }else{
    if(class(x) == "Siid")
    {
      b.size <- 1
      x <- x$data
    }
  }
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
  
  if(iid == FALSE) suppressWarnings(sigma.mat <- mcse.multi(Y, size = b.size)$cov) else sigma.mat <- cov(Y)
  
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



#' @title Add simultaneous confidence interval to existing boxplot
#'
#' @description Adds simultaneous confidence intervals for quantiles to an existing boxplot.
#'
#' @name boxCI
#' @aliases boxCI
#' @usage boxCI(x, CI, component = c(1), dimn = 1, 
#'                quan.color = 'lightsteelblue3', horizontal = FALSE)
#' @param x : a `Smcmc' class object
#' @param CI : the output from the `getCI` function with `iid = TRUE`
#' @param component : vector indicating which components to draw the confidence intervals for
#' @param dimn : numeric for how many components are being plotted
#' @param quan.color : color for the quantile confidence intervals
#' @param horizontal : logical for whether boxplots are horizontal
#' @return adds segments for confidence intervals into an already existing plot environment
#'
#' @examples
#' output <- matrix(rnorm(3*1e3), nrow = 1e3, ncol = 3)
#'
#' @export
boxCI <- function(x,
                  CI,
                  component = 1,
                  dimn       = 1,
                  quan.color = 'lightsteelblue3',
                  horizontal = FALSE) 
{ 
  quans = CI$xi.q
  mn = CI$mean.est
  quansi <- quans[, component]
  qcil = CI$lower.ci.mat[, component]
  qciu = CI$upper.ci.mat[, component]
  i <- component
  if(dimn == 1) i <- 1
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







#' @title ACF Plot for Markov chain Monte Carlo
#'
#' @description Autocorrelation function plots for MCMC data (including multiple chains)
#'
#'
#' @name ACF
#' @usage ACF(x,component = NULL, type = c("correlation", "covariance"),
#'              plot= TRUE, lag.max = NULL, avg.col = "blue", chain.col   = "red",
#'              na.action   = na.fail, auto.layout = TRUE, ask = dev.interactive()) 
#'          
#' @param x : an `Smcmc' class object or a list of Markov chains or a Markov chain matrix
#' @param component : a vector of integers indicating which components' ACF plots are needed. By default all components are drawn.
#' @param type : the kind of ACF plot: "correlation" or "covariance"
#' @param plot : TRUE if plots are required. If FALSE, raw values are returned
#' @param lag.max : Maximum lag for the ACF plot
#' @param avg.col  : color for the overall ACF of each component
#' @param chain.col : color for the ACF of the individual chains.
#' @param na.action :  function to be called to handle missing values. ‘na.pass’ can be used.
#' @param auto.layout : logical argument for an automatic layout of plots
#' @param ask : activating interactive plots


#' @return returns the autocorrelation function plots of the Markov chains. Uses the
#'         more accurate globally-centered ACFs.
#' @examples
#' # Producing Markov chain
#' chain <- matrix(0, ncol = 1, nrow = 1e3)
#' chain[1,] <- 0
#' err <- rnorm(1e3)
#' for(i in 2:1e3)
#' {
#'   chain[i,] <- .3*chain[i-1,] + err[i]
#' }
#' chain <- Smcmc(list(chain))
#' ACF(chain)
#'
#' @references
#' Agarwal, M., and Vats, D., “Globally-centered autocovariances in MCMC”, 
#' arxiv - 2009.01799,  2020. 
#'
#' @export

ACF <- function(x,
                component   = NULL,
                type        = c("correlation", "covariance", "partial"),
                plot        = TRUE,
                lag.max     = NULL,
                avgcol      = "blue",
                chain.col   = "red",
                na.action   = na.fail,
                auto.layout = TRUE,
                ask         = dev.interactive())
{
  
  if(is.matrix(x))
  {
    x <- list(x)
  }
  if(class(x) == "Smcmc")
  {
    x <- x$chains
  }else if(!is.list(x))
  {
    stop("x must be a matrix, list or an Smcmc  object")
  }
  dimn <- dim(x[[1]])
  n <- dimn[1]
  p <- dimn[2]

  which <- as.numeric(component)
  if(is.null(component)) which <- 1:p

  varnames <- colnames(x[[1]])
  if (is.null(lag.max)) lag.max <- floor(10 * (log10(n)) )
  
  lag.max <- as.integer(min(lag.max, n - 1L))
  m <- length(x)

  
  if(plot) 
  {
     setLayout(length(which))
  }
  
  for(i in which)
  {
    xi <- matrix(data = 0, nrow = n, ncol = m)

    for (j in 1:m) {
      xi[,j] <- x[[j]][,i]
    }
    xi <- xi - mean(xi) # global mean
    
    chain.acf <- list(length = m)
    
    for(j in 1:m) {
      chain.acf[[j]] <- acf(xi[,j], type = type, plot = FALSE, demean = FALSE, lag.max = lag.max)
    }

    avgf <-  chain.acf[[1]]
    k = 100
    avgf$acf <-  0
    for(j in 1:m) 
    {
      avgf$acf <- avgf$acf + chain.acf[[j]]$acf
      k = min(k,min(chain.acf[[j]]$acf))
    }
    avgf$acf <- avgf$acf/m
    
    if(plot)
    {
      plot(avgf, ci = 0, main = varnames[i], lwd = .2, ylim = c( min(- 1.96/sqrt(n),min(avgf$acf),k), max(avgf$acf)))
      
      for(j in 1:m)
      {
        lines(0:lag.max, as.matrix(chain.acf[[j]]$acf), type = "l", col = adjustcolor(chain.col, alpha.f = .5), lwd = 1, lty = 1, yaxt = 'n', xaxt = 'n')
      }
      lines(0:lag.max, as.matrix(avgf$acf), type = "l", col = adjustcolor(avgcol, alpha.f = .6), lwd = 1, lty = 1, yaxt = 'n', xaxt = 'n')
      
    }
    
  }
on.exit( par (ask = FALSE,mfrow=c(1,1)))
invisible(list("combined" = avgf, "individual" = chain.acf))    
}





#' @title Trace Plot for Markov chain Monte Carlo
#' @description traceplot is a graphical tool commonly used in Bayesian statistics and Markov Chain Monte Carlo(MCMC) methods to diagnose the convergence and mixing properties of a chain.
#' @name traceplot
#' @usage function(x, fast = NULL,which = NULL, xlim = NULL, ylim = NULL,
#'                   xlab = "Iteration",ylab = NULL, auto.layout = TRUE,
#'                   alpha.f = 0.9, ask = dev.interactive(),
#'                   col = c("palevioletred3","steelblue3","tan3","lightsteelblue3",
#'                           "springgreen2","skyblue3","khaki3",
#'                           "lightpink1","palegreen3","thistle3","aquamarine3","dimgrey",
#'                           "tomato3"))
#'
#'@param x : an `Smcmc' class object or a list of Markov chains or a Markov chain matrix or a vector.
#'@param fast : a Boolean argument that will be auto set to TRUE when chain size is larger than 1e6
#'@param which : if we want full size trace plots of specific dimensions of chain, we can pass a vector of respective dimension/components. 
#'@param xlim : range of x-axis
#'@param ylim : range of y-axis
#'@param main : usual heading for plot
#'@param xlab : labels of x-axis
#'@param ylab : labels of y-axis(it should be a vector of length equal to dimension of chain)
#'@param auto.layout : logical argument for an automatic layout of plots
#'@param ask : activating interactive plots
#'@param col : color vector for multiple chains
#'@param alpha.f : To fix the opacity of lines as per user convenience, by default it is 0.9.
#'
#' @return Returns the Trace Plots of Markov Chain(s)
#' @examples
#' # example code
#' # Defining a function to produce Markov chain with dimension p and size n
#' MakeChain <- function(p, n , h = .5)
#' {
#'   chain <- matrix(0, nrow = n,ncol = p)
#'   for (i in 2:n) 
#'   {
#'     prop <- chain[i-1, ] + rnorm(p, mean = 0, sd = h)
#'     log.ratio <- sum(dnorm(prop, log = TRUE) - dnorm(chain[i-1, ], log = TRUE))
#'     if(log(runif(1)) < log.ratio) 
#'     {chain[i, ] <- prop}
#'     else{chain[i, ] <- chain[i - 1, ]}
#'   }
#'   v = vector(length = p)
#'   for(i in 1:p){v[i] = paste("Comp ",i)}
#'   colnames(chain) = v
#'   return(chain)
#' }
#'
#'
#' chain1 <- MakeChain(p=4,n=1000)
#' chain2 <- MakeChain(p=4,n=1000)
#' chain3 <- MakeChain(p=4,n=1000)
#' out <- Smcmc(list(chain1,chain2,chain3))
#' traceplot(out)
#'
#'
#' chain1 <- MakeChain(p=6,n=1000)
#' chain2 <- MakeChain(p=6,n=1000)
#' chain3 <- MakeChain(p=6,n=1000)
#' out <- Smcmc(list(chain1,chain2,chain3))
#' traceplot(out)

#'
#'
#'
#'
#' @export
traceplot <- function(x, fast = NULL, which = NULL, xlim = NULL, ylim = NULL,main = NULL,
                      xlab = "Iteration",ylab = NULL,alpha.f = 0.9, auto.layout = TRUE,
                      ask = dev.interactive(),
                      col = c("palevioletred3","steelblue3","tan3","lightsteelblue3",
                              "springgreen2","skyblue3","khaki3",
                              "lightpink1","palegreen3","thistle3","aquamarine3","dimgrey",
                              "tomato3"))
{
  ## Check for the Data Type
  if(TRUE %in% (class(x) == "Smcmc")){x <- x$chains}else if(is.list(x))
  {
    for(i in 1:length(x)) { x[[i]] = as.matrix(x[[i]]) }
  }else if(is.vector(x))
  {
    x <- as.matrix(x)
    x <- list(x)
  }else if(is.matrix(x))
  {x <- list(x)}else
      {stop("x must be a matrix, list or an Smcmc object")}
  
  dimn <- dim(x[[1]])
  n <- dimn[1]
  p <- dimn[2]
  m <- length(x)
  
  if(m>5){stop("Maximum 5 chains are allowed")}
  
  if(is.null(ylab) | !is.vector(ylab) | length(ylab)!=p)
  {
    if(!is.null(ylab) && (!is.vector(ylab) | length(ylab) != p ))
    {message("ylab should be NULL, or vector of length no. of dimension. Default value of ylab applied.")}
    ylab <- vector(length = p)
    for(i in 1:p) { ylab[i] <- paste("Comp ",i) }
  }
  
  if(dim(x[[1]])[1] > 100000 && is.na(fast)){ fast = TRUE }
  if(!is.null(fast) && isTRUE(fast))
  {
    if(n < 1e3){index = 1:n}
    else
    {
      index <- sample(1:n,1000,replace = FALSE)
      index <- sort(index)
    }
  }else{index = 1:n}
  
  xlim <- if (is.null(xlim)) c(0,n) else xlim
  
  maxi <- rep(-1e9,p)
  mini <- rep(1e9,p)
  vec <- vector(length = m)
  for(j in 1:m)
  {
    if(is.null(ylim)){
      mat <- as.matrix(x[[j]])
      tempx <- rbind(apply(mat,2,max), maxi)
      tempn <- rbind(apply(mat,2,min),mini)
      maxi <- apply(tempx,2,max)
      mini <- apply(tempn,2, min)
    }
    else{mini = rep(ylim[1],p)
         maxi = rep(ylim[2],p)}
      vec[j] <- paste("Chain",j)
  }
  
  ## this if condition is for plotting traceplots of chains, 
  ## when which not equal to NULL and Dimension of chain exceeds 12.
  if(!is.null(which)| p > 10)
  {
    if(p > 10 && is.null(which) )
    {
      which <- 1:p
      message("Number of dimension of chain(s) is more than 10. Series of plot returned!")
      
    }
    lay <- par()
    leg <- lay$mfrow[1]*lay$mfrow[2]
    for(i in which)
    {
       par(mar = c(4, 4, 3, 2))
       j <- 1  
       ylim <- c(mini[i],maxi[i])
       plot(x = index, y = x[[1]][index,i], 
            xlab = xlab, ylab = ylab[i],
            type = "n", lwd = 1.3, lty = 1, ylim = ylim, xlim = xlim,
            col = adjustcolor(col[1],alpha.f = alpha.f))
       while(j <= m)
       {
        lines( x= index, y=x[[j]][index,i], type = "l", 
              col = adjustcolor(col[j],alpha.f = alpha.f),
              ylim = ylim, xlim = xlim, lwd = 1.3, lty = 1, 
              yaxt = 'n', xaxt = 'n')
         j  <- j + 1
        }
      
       ##Setting of legend
       if(i%%leg == 0|i == which[length(which)])
       {
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
        legend("top",legend = vec, col = col[1:m],lty = 1,lwd =2, xpd = TRUE, 
                horiz = TRUE,cex = 1.25, seg.len= 1, bty = 'n')
        par(mfrow = lay$mfrow)
        }
      }
       on.exit(par(ask = FALSE, mfrow = c(1,1)))
       par(mar = c(5.1, 4.1, 4.1, 2.1))
       par(fig = c(0, 1, 0, 1))
       par(oma = c(0, 0, 0, 0))
      
  }
  else
  {  
    par(oma = c(4,0,3,0))
    axs = p-1
    if(p <= 4){axs = p-1}
    else {axs = p-2}
    setLayout_trace(p)
    for(i in 1:p)
    {
      j = 2
      par(mar = c(0, 4.1,0, 2.1))
      ylim <-  c(mini[i],maxi[i])
      plot(x= index, y = x[[1]][index,i], xlab = if(i>axs) xlab else "" , ylab =ylab[i],
           type = "l", lwd = 1.3, lty = 1, ylim = ylim, xlim = xlim,
           col = adjustcolor(col[1],alpha.f = alpha.f), xaxt = if(i>axs) 's' else 'n')
      while(j<=m)
      {
        lines(x= index, y=x[[j]][index,i], type = "l", 
              col = adjustcolor(col[j],alpha.f = alpha.f),
              ylim = ylim, xlim = xlim, lwd = 1.3, lty = 1, 
              yaxt = 'n', xaxt = 'n')
        j = j + 1
      }
      if(p%%2 != 0 && i == p-1 && p>4){mtext("Iteration", side = 1, line = 2.3,cex = 1)}
    }
    if(p<=4){mtext("   Iteration", side = 1, line = 2.3, outer = TRUE, cex= 1)}
    else{mtext("Iteration", side = 1, line = 2.3,at = 0.3,outer = TRUE,cex = 1)
      if(p%%2 == 0){mtext("Iteration", side = 1, line = 2.3,at = 0.8,outer = TRUE,cex = 1)}}
    ##Setting of legend
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
      legend("top",legend = vec[1:m], col = col[1:m],lty = 1,lwd =2, xpd = TRUE, 
             horiz = TRUE,cex = 1.25, seg.len= 1, bty = 'n')
      
      
  }
  if(isTRUE(fast) && n > 100000){message("Chain size is very large, fast changed to TRUE.")}
  on.exit(par(ask = FALSE, mfrow = c(1,1)))
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  par(fig = c(0, 1, 0, 1))
  par(oma = c(0, 0, 0, 0))
}
