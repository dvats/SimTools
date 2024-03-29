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
    data <- list(data)

  nsim <- dim(data[[1]])[1]
  if(is.null(varnames)) varnames <- colnames(data[[1]])  
  
  if(stacked == TRUE)
  {
    foo <- chain_stacker(data)
    stacked.chain <- foo$stacked.data
    
    if(batch.size == TRUE)
    {
      size <- foo$b.size
    }
    else{
      size <- batch.size
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




#' @title Plot Smcmc
#'
#' @description Density plots with simultaenous error bars around means and quantiles
#'  for MCMC data. The error bars account for the correlated nature of the process.
#'
#'
#' @name plot.Smcmc
#' @usage \method{plot}{Smcmc}(x, Q = c(0.1, 0.9), alpha = 0.05, thresh = 0.001, iid = FALSE,
#'                             plot = TRUE, mean = TRUE, border = NA, mean.col = 'plum4', 
#'                             quan.col = 'lightsteelblue3',rug = TRUE, opaq = 0.7, 
#'                             auto.layout = TRUE, ask = dev.interactive(),...)    
#' @param x : a `Smcmc' class object
#' @param Q : vector of quantiles
#' @param alpha : confidence level of simultaneous confidence intervals 
#' @param thresh : numeric typically less than .005 for the accuracy of the simulteaneous procedure
#' @param iid : logical argument for constructing density plot for iid samples. Defaults to \code{FALSE}
#' @param plot :  logical argument for is plots are to be returned 
#' @param mean : logical argument whether the mean is to be plotted
#' @param border : whether a border is required for the simultaneous confidence intervals
#' @param mean.col : color for the mean confidence interval
#' @param quan.col : color for the quantile confidence intervals
#' @param rug : logical indicating whether a rug plot is desired
#' @param opaq : opacity of \code{mean.col} and \code{quan.col}. A value of 0 is transparent and 1 is completely opaque.
#' @param auto.layout : logical argument for an automatic layout of plots
#' @param ask : activating interactive plots
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
#' plot(chain)
#'
#' @references
#' Robertson, N., Flegal, J. M., Vats, D., and Jones, G. L., 
#' “Assessing and Visualizing Simultaneous Simulation Error”, 
#' Journal of Computational and Graphical Statistics,  2020. 
#'
#' @export
"plot.Smcmc" <- function(x, 
                         Q        = c(0.1, 0.9), 
                         alpha    = 0.05, 
                         thresh   = 0.001, 
                         iid      = FALSE, 
                         plot     = TRUE,  
                         mean     = TRUE, 
                         border   = NA, 
                         mean.col = 'plum4', 
                         quan.col = 'lightsteelblue3',
                         rug      = TRUE, 
                         opaq     = 0.7, 
                         auto.layout = TRUE, 
                         ask      = dev.interactive(), ...)
{
  
  x <- as.Smcmc(x)
  out <- getCI(x, Q, alpha, thresh = thresh, iid = iid, mean = mean)
  if(plot == TRUE)
  {
    plot.CIs(x, dimn = length(x$stacked[1,]), CIs = out, bord = border, 
    	mean.color = adjustcolor(mean.col, alpha.f = opaq), 
    	quan.color = adjustcolor(quan.col, alpha.f = opaq), 
    	mean = mean, auto.layout = auto.layout, rug = rug,
    	ask = ask, ...)
  }
  invisible(out)
}



