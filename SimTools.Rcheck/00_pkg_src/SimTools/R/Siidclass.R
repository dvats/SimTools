#' @title Siid class
#'
#' @description Class for independent and identically distributed (iid) samples
#'
#' @name Siid
#' @aliases Siid as.Siid as.Siid.default is.iid 
#' @usage Siid(data, varnames = colnames(data))
#' @param data : an iid output matrix with nsim rows and p columns
#' @param varnames : a character string equal to the number of columns in \code{data}
#'
#' @return an Siid class object
#' @examples
#' # Generating iid data
#' chain <- matrix(rnorm(3*1e3), nrow = 1e3, ncol = 3)
#' siid.obj <- Siid(chain)
#'
#' @export
"Siid" <- function(data, varnames = colnames(data)) # make Siid object
{
  if(missing(data))
    stop("Data must be provided.")

  if(is.vector(data))
    data <- as.matrix(data, ncol = 1)

  nsim <- dim(data)[1]

  attr(data,"nsim") <- nsim
  attr(data,"varnames") <- varnames
  attr(data,"class") <- "Siid"
  return(data)
}

"is.Siid" <- function (x) 
{
  if (inherits(x, "Siid")) 
    return(TRUE)
  return(FALSE)
}

#' @export
"as.Siid" <- function (x, ...) 
  UseMethod("as.Siid")

#' @export
"as.Siid.default" <- function (x, ...) 
  if (is.Siid(x)) x else Siid(x)


#' @title Boxplot for Siid
#'
#' @description Boxplots with simultaenous error bars around all quantiles for iid data.
#'
#'
#' @name boxplot.Siid
#' @usage \method{boxplot}{Siid}(x, ...,  alpha = 0.05, thresh = 0.001,  mean.col = 'plum4',
#'                      quan.col = 'lightsteelblue3', opaq = .6, range = 1.5, width = NULL, varwidth = FALSE,
#'                      outline = TRUE, plot = TRUE, border = par("fg"), col = 'white',
#'                      ann = !add, horizontal = FALSE, add = FALSE)
#' @param x : a `Siid' class object
#' @param ... : arguments sent to boxplot  
#' @param alpha : confidence level of simultaneous confidence intervals 
#' @param thresh : numeric typically less than .005 for the accuracy of the simulteaneous procedure
#' @param mean.col : color for the mean confidence interval
#' @param quan.col : color for the quantile confidence intervals
#' @param opaq : opacity of \code{mean.col} and \code{quan.col}. A value of 0 is transparent and 1 is completely opaque.
#' @param range :  as defined for base \code{boxplot}
#' @param width : as defined for base \code{boxplot}
#' @param varwidth : as defined for base \code{boxplot}
#' @param outline : as defined for base \code{boxplot}
#' @param plot : logical indicating whether the plot is to be constructed
#' @param border : as defined for base \code{boxplot}
#' @param col : as defined for base \code{boxplot}
#' @param ann : as defined for base \code{boxplot}  
#' @param horizontal : as defined for base \code{boxplot}  
#' @param add : as defined for base \code{boxplot}    
#' @return returns the base \code{boxplot} with simultaneous confidence intervals around all quantiles
#' @examples
#' # Generating iid data
#' chain <- matrix(rnorm(3*1e3), nrow = 1e3, ncol = 3)
#' siid.obj <- Siid(chain)
#' boxplot(chain)
#'
#' @references
#' Robertson, N., Flegal, J. M., Vats, D., and Jones, G. L., 
#' “Assessing and Visualizing Simultaneous Simulation Error”, 
#' Journal of Computational and Graphical Statistics,  2020. 
#'
#' @export

boxplot.Siid <- function(x, ...,  alpha = 0.05, thresh = 0.001, mean.col = 'plum4',
                      quan.col = 'lightsteelblue3', opaq = .6, range = 1.5, width = NULL, varwidth = FALSE,
                      outline = TRUE, plot = TRUE, border = par("fg"), col = 'white',
                      ann = !add, horizontal = FALSE, add = FALSE)
{
  x <- as.Siid(x)
  Q <- c(0.25, 0.50,0.75)
  notch <- FALSE
  foo3 <- error.est(x, Q, alpha, thresh = thresh, iid = TRUE)
  plot.boxx(x, dimn = length(x[1,]), CIs = foo3, mean.color = mean.col, quan.color = adjustcolor(quan.col, alpha.f = opaq), mn = foo3$mean.est, 
            quans = foo3$xi.q, range = range, width = width, varwidth = varwidth, notch = notch, outline = outline,
            plot = plot, border = border, col = col, ann = ann, horizontal = horizontal, add = add)
}
