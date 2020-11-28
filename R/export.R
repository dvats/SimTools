#' @title Plot for Smcmc class
#'
#' @description Density plots with simultaenous error bars around means and quantiles
#'  for MCMC data. The error bars account for the correlated nature of the process
#'
#' @name plot.Smcmc
#' @usage plot.Smcmc(x, Q = c(0.1, 0.9), alpha = 0.05, thresh = 0.001, iid = FALSE, plot = TRUE, auto.layout = TRUE, ask = dev.interactive(),mean = TRUE, border = NA, mean.col = 'plum4', quan.col = 'lightsteelblue3', opaq = 0.7, ...)
#' @param x : multivariate or univariate Markov chain  \eqn{n \times p} matrix where \eqn{n} is the number of Monte Carlo samples and \eqn{p} is the dimension of the Markov chain.
#' @param Q : vector of quantiles
#' @param alpha : confidence level of simultaneous confidence intervals 
#' @param thresh : numeric typically less than .005 for the accuracy of the simulteaneous procedure
#' @param iid : logical argument for constructing density plot for iid samples. Defaults to \code{FALSE}
#' @param plot :  logical argument for is plots are to be returned 
plot.Smcmc <- function(x, Q = c(0.1, 0.9), alpha = 0.05, thresh = 0.001, iid = FALSE, plot = TRUE, auto.layout = TRUE, ask = dev.interactive(),mean = TRUE, border = NA, mean.col = 'plum4', quan.col = 'lightsteelblue3', opaq = 0.7, ...)
{
  out <- error.est(x, Q, alpha,thresh = thresh,iid = iid,mean = mean, ...)
  if(plot == TRUE)
  {
    plot.CIs(x, dimn = length(x[1,]), CIs = out, bord = border, mean.color = adjustcolor(mean.col, alpha.f = opaq), quan.color = adjustcolor(quan.col, alpha.f = opaq), mn = out$mean.est,mean = mean, quans = out$xi.q,auto.layout = auto.layout,ask = ask,...)
  }
  return(out)
}
