pkgname <- "SimTools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SimTools')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ACF")
### * ACF

flush(stderr()); flush(stdout())

### Name: ACF
### Title: ACF Plot for Markov chain Monte Carlo
### Aliases: ACF

### ** Examples

# Producing Markov chain
chain <- matrix(0, ncol = 1, nrow = 1e3)
chain[1,] <- 0
err <- rnorm(1e3)
for(i in 2:1e3)
{
  chain[i,] <- .3*chain[i-1,] + err[i]
}
chain <- Smcmc(list(chain))
ACF(chain)




cleanEx()
nameEx("Siid")
### * Siid

flush(stderr()); flush(stdout())

### Name: Siid
### Title: Siid class
### Aliases: Siid as.Siid as.Siid.default is.iid

### ** Examples

# Generating iid data
chain <- matrix(rnorm(3*1e3), nrow = 1e3, ncol = 3)
siid.obj <- Siid(chain)




cleanEx()
nameEx("Smcmc")
### * Smcmc

flush(stderr()); flush(stdout())

### Name: Smcmc
### Title: Smcmc class
### Aliases: Smcmc as.Smcmc as.Smcmc.default is.mcmc

### ** Examples

# Producing Markov chain
chain <- matrix(0, nrow = 1e3, ncol = 1)
chain[1,] <- 0
err <- rnorm(1e3)
for(i in 2:1e3)
{
  chain[i,] <- .3*chain[i-1,] + err[i]
}
smcmc.obj <- Smcmc(chain)



cleanEx()
nameEx("addCI")
### * addCI

flush(stderr()); flush(stdout())

### Name: addCI
### Title: Add simultaneous confidence interval to existing plot.
### Aliases: addCI

### ** Examples

chain <- matrix(0, ncol = 1, nrow = 1e3)
chain[1,] <- 0
err <- rnorm(1e3)
for(i in 2:1e3)
{
  chain[i,] <- .3*chain[i-1,] + err[i]
}
chain <- Smcmc(list(chain))
plot(density(chain$stacked[,1]))
CIs <- getCI(chain)
addCI(chain, CIs, component = 1)



cleanEx()
nameEx("boxCI")
### * boxCI

flush(stderr()); flush(stdout())

### Name: boxCI
### Title: Add simultaneous confidence interval to existing boxplot
### Aliases: boxCI

### ** Examples

output <- matrix(rnorm(3*1e3), nrow = 1e3, ncol = 3)




cleanEx()
nameEx("boxplot.Siid")
### * boxplot.Siid

flush(stderr()); flush(stdout())

### Name: boxplot.Siid
### Title: Boxplot for Siid
### Aliases: boxplot.Siid

### ** Examples

# Generating iid data
chain <- matrix(rnorm(3*1e3), nrow = 1e3, ncol = 3)
siid.obj <- Siid(chain)
boxplot(siid.obj)




cleanEx()
nameEx("getCI")
### * getCI

flush(stderr()); flush(stdout())

### Name: getCI
### Title: Calculates simultaneous confidence intervals.
### Aliases: getCI

### ** Examples

chain <- matrix(0, ncol = 1, nrow = 1e3)
chain[1,] <- 0
err <- rnorm(1e3)
for(i in 2:1e3)
{
  chain[i,] <- .3*chain[i-1,] + err[i]
}
chain <- Smcmc(list(chain))
plot(density(chain$stacked[,1]))
CIs <- getCI(chain)
addCI(chain, CIs, component = 1)




cleanEx()
nameEx("plot.Siid")
### * plot.Siid

flush(stderr()); flush(stdout())

### Name: plot.Siid
### Title: Plot Siid
### Aliases: plot.Siid

### ** Examples

# Generating iid data
chain <- matrix(rnorm(3*1e3), nrow = 1e3, ncol = 3)
siid.obj <- Siid(chain)
plot(siid.obj)




cleanEx()
nameEx("plot.Smcmc")
### * plot.Smcmc

flush(stderr()); flush(stdout())

### Name: plot.Smcmc
### Title: Plot Smcmc
### Aliases: plot.Smcmc

### ** Examples

# Producing Markov chain
chain <- matrix(0, ncol = 1, nrow = 1e3)
chain[1,] <- 0
err <- rnorm(1e3)
for(i in 2:1e3)
{
  chain[i,] <- .3*chain[i-1,] + err[i]
}
chain <- Smcmc(list(chain))
plot(chain)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
