## usethis namespace: start
#' @importFrom grDevices adjustcolor dev.interactive
#' @importFrom graphics boxplot par polygon segments boxplot.matrix
#' @importFrom stats cov density qnorm quantile ts
#' @importFrom mcmcse mcse.multi batchSize
#' @importFrom mvtnorm pmvnorm
## usethis namespace: end

#' @title Smcmc class
#'
#' @description Smcmc class for simulated data using Markov chain Monte Carlo
#'
#' @name Smcmc
#' @aliases Smcmc as.Smcmc as.Smcmc.default is.mcmc 
#' @usage Smcmc(data, batch.size = TRUE, stacked = TRUE, varnames = colnames(data))
#' @param data : a list of MCMC output matrices each with `nsim` rows and `p` columns
#' @param batch.size : logical argument, if true, calculates the batch size appropriate for this Markov chain. Setting to TRUE saves time in future steps.
#' @param stacked : recommended to be `TRUE`. logical argument, if true, stores a carefully stacked version of the MCMC output for use later.
#' @param varnames : a character string equal to the number of columns in \code{data}
#'
#' @return an Smcmc class object
#' @examples
#' # Producing Markov chain
#' chain <- matrix(0, nrow = 1e3, ncol = 1)
#' chain[1,] <- 0
#' err <- rnorm(1e3)
#' for(i in 2:1e3)
#' {
#'   chain[i,] <- .3*chain[i-1,] + err[i]
#' }
#' smcmc.obj <- Smcmc(chain)
#' @export
Smcmc <- function(data,
                    batch.size = TRUE, 
                    stacked = TRUE,
                    varnames = NULL) # make Smcmc object
{
  if(missing(data))
    stop("Data must be provided.")
  
  
  if(!is.list(data))
  {
    data = as.matrix(data)
    data <- list(data)
  }else
  {
    for(i in 1:length(data))
    {
      data[[i]] = as.matrix(data[[i]])
    }
  }
  
  nsim <- dim(data[[1]])[1]
  if(is.null(varnames))
  {
    for(i in 1:length(data)){if(!is.null(colnames(data[[i]]))){varnames = colnames(data[[i]]);break }}
  }
  
  if(stacked == TRUE)
  {
    foo <- chain_stacker(data)
    stacked.chain <- foo$stacked.data
    
    if(batch.size == TRUE)
    {
      size <- foo$b.size
    }
    else{
      size <- NULL
    }
  }
  
  out <- list( chains = data,
               stacked  = stacked.chain,
               b.size   = size,
               nsim     = nsim,
               varnames = varnames)
  
  class(out) <- "Smcmc"
  return(out)
}

"is.Smcmc" <- function (x) 
{
  if (inherits(x, "Smcmc")) 
    return(TRUE)
  return(FALSE)
}

#' @export
"as.Smcmc" <- function (x, ...) 
  UseMethod("as.Smcmc")

#' @export
"as.Smcmc.default" <- function (x, ...) 
  if (is.Smcmc(x)) x else Smcmc(x)




#' @title density plot form Smcmc class
#'
#' @description Density plots with simultaenous error bars around means and quantiles
#'  for MCMC data. The error bars account for the correlated nature of the process.
#'
#'
#' @name densityplot
#' @usage densityplot(x, Q = c(0.1, 0.9), alpha = 0.05, thresh = 0.001, iid = FALSE,
#'                             plot = TRUE, mean = TRUE, border = NA, mean.col = 'plum4', 
#'                             quan.col = 'lightsteelblue3',rug = TRUE, opaq = 0.7, 
#'                             auto.layout = TRUE, ask = dev.interactive(),...)    
#' @param x : a `Smcmc' class object
#' @param Q : vector of quantiles
#' @param alpha : confidence level of simultaneous confidence intervals 
#' @param thresh : numeric typically less than .005 for the accuracy of the simulteaneous procedure
#' @param main : To add main heading
#' @param iid : logical argument for constructing density plot for iid samples. Defaults to \code{FALSE}
#' @param plot :  logical argument for is plots are to be returned 
#' @param mean : logical argument whether the mean is to be plotted
#' @param which : A vector of components, if you want plots of specific components.
#' @param border : whether a border is required for the simultaneous confidence intervals
#' @param mean.col : color for the mean confidence interval
#' @param quan.col : color for the quantile confidence intervals
#' @param rug : logical indicating whether a rug plot is desired
#' @param opaq : opacity of \code{mean.col} and \code{quan.col}. A value of 0 is transparent and 1 is completely opaque.
#' @param auto.layout : logical argument for an automatic layout of plots
#' @param ask : activating inter active plots
#' @param ... : arguments passed on to the \code{density} plot in base R
#' @return returns a plot of the univariate density estimates with simultaneous
#'			confidence intervals wherever asked. If \code{plot == FALSE} a list of
#'			estimates and simultaneous confidence intervals.
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
#' densityplot(chain)
#'
#' @references
#' Robertson, N., Flegal, J. M., Vats, D., and Jones, G. L., 
#' “Assessing and Visualizing Simultaneous Simulation Error”, 
#' Journal of Computational and Graphical Statistics,  2020. 
#'
#' @export


"densityplot" <- function(x, 
                          Q        = c(0.1, 0.9), 
                          alpha    = 0.05, 
                          thresh   = 0.001, 
                          main     = NA,
                          iid      = FALSE, 
                          plot     = TRUE,  
                          mean     = TRUE,
                          which    = NULL,
                          border   = NA, 
                          mean.col = 'plum4', 
                          quan.col = 'lightsteelblue3',
                          rug      = FALSE, 
                          opaq     = 0.7, 
                          auto.layout = TRUE, 
                          ask      = dev.interactive(), ...)
{
  
  x <- as.Smcmc(x)
  out <- getCI(x, Q, alpha, thresh = thresh, iid = iid, mean = mean)
  if(plot == TRUE)
  {
    plot.CIs(x, CIs = out, bord = border, 
             mean.color = adjustcolor(mean.col, alpha.f = opaq), 
             quan.color = adjustcolor(quan.col, alpha.f = opaq), 
             mean = mean, auto.layout = auto.layout, rug = rug,
             ask = ask, main = main, which= which, ...)
  }
  invisible(out)
}


#' @title Summary plot function for Smcmc objects
#'
#' @description Plots traceplot, acfplot and densityplot of all the dimensions of 
#' chains in  Smcmc object
#'
#' @name plot.Smcmc
#' @usage plot.Smcmc(x)    
#' @param x : a `Smcmc' class object
#' @return return plot(s) of the all the dimensions of Smcmc object
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
#' plot(chain)
#'
#'@export

"plot.Smcmc" <- function(x,which = NULL,...)
{
  if(class(x)!="Smcmc"){stop("Argument must be Smcmc object")}
  y = x$chains
  dimn <- dim(y[[1]])
  n <- dimn[1]
  p <- dimn[2]
  m <- length(x)
  p2 <- 0
  p3 <- 0
  if(!is.null(which)){p = length(which)}else{which = 1:p}
  if(p>12){stop("Maximum allowed dimension is 12")}
  if(p>4){p2 = p-4;p = 4}
  if(p2 >4)
  {
    p3 <- p2-4
    p2 <- 4
  }
  par(mfrow = c(3,p),oma = c(0,0,0,0),mar = c(2.2,4,1,1))
  traceplot(x,which = which[1:p],legend = F,xlab = NA,...)
  acfplot(x,which = which[1:p],xlab = NA,...)
  densityplot(x,which = which[1:p],...)
  
  if(p2>0)
  {
    par(mfrow = c(3,p2),mar = c(2,4,1,2))
    traceplot(x,which = which[5:(4+p2)],legend = F,xlab = NA)
    acfplot(x,which = which[5:(4+p2)],xlab = NA)
    densityplot(x,which = which[5:(4+p2)])
  }
  
  if(p3>0)
  {
    par(mfrow = c(3,p3),mar = c(2,4,1,2))
    traceplot(x,which = which[9:(8+p3)],legend = F,xlab = NA)
    acfplot(x,which = which[9:(8+p3)],xlab = NA)
    densityplot(x,which = which[9:(8+p3)])
  }
  on.exit(par(ask = FALSE,mfrow=c(1,1)))
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  par(fig = c(0, 1, 0 , 1))
  par(oma = c(0, 0, 0, 0))
}


#' @title Summary function for Smcmc objects
#' @name summary.Smcmc
#' @description To show different statistics of the Smcmc object
#' @usage \method{summary}{Smcmc}(x)    
#' @param x : a `Smcmc' class object
#' @return return statistics of the all the dimensions(& chains) in Smcmc object
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
#' summary(chain)
#'
#'@export

"summary.Smcmc" <- function (object ,...)
{
  object <- as.Smcmc(object)
  object.class <- class(object)
  Batch_Size <- object$b.size
  Smcmc_output <- object$chains[[1]]
  stacked <- object$stacked
  
  chains <- object$chains
  n <- dim(chains[[1]])[1]
  p <- dim(chains[[1]])[2]
  columns <- ncol(Smcmc_output)
  m <- length(chains)
  dimen <- vector(length = columns)
  if(!is.null(object$varnames)){dimen <- object$varnames }else{for(i in 1:columns){dimen[i] <- paste("Component",i)}}
  colname <- c("Mean","SD","MCSE","ESS","Gelman-Rubin")
  
  statistics <- matrix(0,ncol = 5, nrow = p)
  statistics[ ,1] <- apply(stacked,2,mean)
  statistics[ ,2] <- apply(stacked,2,sd)
  statistics[ ,3] <- (round(mcmcse::mcse.mat(stacked,size = Batch_Size), 4))[ ,2]
  mp_ess <- mcmcse::multiESS(stacked, size = Batch_Size)
  foo <- p*log(mp_ess/n)
  ms <- cov(chains[[1]])
  if(m>1)
  {
    for(i in 2:m)
    {
      ms <- ms + cov(chains[[i]])
    }
  }
  ms <- ms/m
  correction <- sum(log(eigen(ms, symmetric = TRUE, 
                              only.values = TRUE)$values)) - sum(log(eigen(cov(stacked), symmetric = TRUE, 
                                                                           only.values = TRUE)$values))
  foo <- (foo + correction)/p
  n_ess <- n*exp(foo)
  multiess <- floor(n_ess)
  multigelmann <- round(sqrt((n-1)/n + m/n_ess),5)
  
  s <- vector(length = p)
  for(i in 1:p)
  {
    ms <- var(chains[[1]][,i])
    if(m>1)
    {
      for(j in 2:m)
      {
        ms <- ms + var(chains[[j]][,i])
      }
    }
    s[i] <- ms
  }
  
  s <- s/m
  statistics[ ,4] <- s/(statistics[,3]^2)
  statistics[ ,5] <- round(sqrt(1 + m/statistics[,4]),5)
  statistics[,4] <- floor(statistics[,4])
  statistics <- as.data.frame(statistics)
  stacked_rows <- dim(stacked)[1]
  colnames(statistics) <- colname
  rownames(statistics) <- dimen
  statistics <- drop(statistics)
  summary_list <- list(nsim = n,
                       Dimensions = p,
                       no.chains = m,
                       stacked_rows = stacked_rows,
                       Statistics = statistics,
                       MultiESS = multiess,
                       MultiGelmann = multigelmann)
  class(summary_list) <- "summary.Smcmc"
  return(summary_list)
}


#'@rdname summary.Smcmc
#'@export

"print.summary.Smcmc" <-function (x, ...) 
{
  cat("\n", "No. of Iterations = ",x$nsim," samples\n", sep = "")
  cat("No. of Dimensions = ",x$Dimensions,"\n", sep = "")
  cat("No. of chains = ",x$no.chains,"\n", sep = "")
  cat("Number of Iterations Considered=",x$stacked_rows, "\n" )
  cat("\n1. Statistics for each variable,\n")
  print(x$Statistics, ...)
  cat("\n")
  cat("MultiESS value =", x$MultiESS, "\n")
  cat("Multi Gelman-Rubin =", x$MultiGelmann)
  cat("\n")
  invisible(x)
}



#' @title Covert to Smcmc Object
#'
#' @description To covert different MCMC objects to Smcmc object
#'
#' @name convert2Smcmc
#' @usage convert2Smcmc(x)    
#' @param x : a object belongs from any of "mcmc.list", "stanfit", "rstan", "array", "matrix" classes.
#' @return return Smcmc object having same chain(s)
#' @examples
#' 
#' library(StanHeaders)
#' library(rstan)
#'
#' stanmodelcode <- "
#' data { 
#'  int<lower=0> n;
#'  vector[2] x[n];
#'  }
#'  parameters {
#'  vector[2] mu;
#'  vector<lower=0>[2] lambda;
#'  real<lower=-1,upper=1> r;
#'  } 
#'  transformed parameters {
#'  vector<lower=0>[2] sigma;
#'  cov_matrix[2] T;
#'  sigma[1] <- inv_sqrt(lambda[1]);
#'  sigma[2] <- inv_sqrt(lambda[2]);
#'  T[1,1] <- square(sigma[1]);
#'  T[1,2] <- r * sigma[1] * sigma[2];
#'  T[2,1] <- r * sigma[1] * sigma[2];
#'  T[2,2] <- square(sigma[2]);
#'  }
#'  model {
#'  mu ~ normal(0, inv_sqrt(.001));
#'  lambda ~ gamma(.001, .001);
#'  x ~ multi_normal(mu, T);
#'  }"
#'
#'  x <- matrix(c(0.8,102,1.0,98,0.5,100,0.9,105,0.7,103,0.4,110, 1.2,99, 1.4,87,0.6,113,1.1,89,
#'                 1.3,93), nrow=11, ncol=2, byrow=TRUE) 
#'                 
#'  n <- nrow(x) # number of people/units measured
#'
#' data <- list(x=x, n=n) # to be passed on to Stan
#' myinits <- list(list(r=0))
#'
#' # parameters to be monitored: 
#' parameters <- c("r")
#'
#' # The following command calls Stan with specific options.
#' # For a detailed description type "?rstan".
#' fit <- stan(model_code=stanmodelcode,data=data,init=myinits,iter=1000,chains=1,thin=1)
#'
#' convert2Smcmc(fit)
#'
#'@export


convert2Smcmc <- function(x)
{
  if(class(x) =="mcmc.list")
  {
    temp = list(as.matrix(x[[1]]))
    
    for(i in 2:length(x))
    {
      append(temp,as.matrix(x[[i]]))
    }
    return(as.Smcmc(temp))
  }
  
  else if("stanfit"%in%class(x) || "rstan"%in%class(x))
  {
    x = as.array(x)
    dim_samples <- dim(x)
    # Reshape the multidimensional array into separate chains
    num_chains <- dim_samples[1]  # Number of chains
    num_iterations <- dim_samples[2]  # Number of iterations
    num_variables <- dim_samples[3]  # Number of variables
    
    chains = list(as.matrix(x[,1 , ]))
    if(num_iterations > 1){
      for (i in 2:num_iterations) {
        append(chains, as.matrix(x[ ,i, ]))
      }
    }
    return(as.Smcmc(chains))
  }
  
  else if("array"%in%class(x) || "matrix"%in%class(x) ||"list"%in%class(x))
  {
    return(as.Smcmc(x))
  }
}
