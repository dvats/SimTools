pkgname <- "SimTools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SimTools')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Smcmc")
### * Smcmc

flush(stderr()); flush(stdout())

### Name: Smcmc
### Title: Plot Smcmc
### Aliases: Smcmc as.Smcmc as.Smcmc.default is.mcmc

### ** Examples

# Producing Markov chain
chain <- numeric(length = 1e3)
chain[1] <- 0
err <- rnorm(1e3)
for(i in 2:1e3)
{
  chain[i] <- .3*chain[i-1] + err[i]
}
smcmc.obj <- Smcmc(chain)



cleanEx()
nameEx("plot.Smcmc")
### * plot.Smcmc

flush(stderr()); flush(stdout())

### Name: plot.Smcmc
### Title: Plot Smcmc
### Aliases: plot.Smcmc

### ** Examples

# Producing Markov chain
chain <- numeric(length = 1e3)
chain[1] <- 0
err <- rnorm(1e3)
for(i in 2:1e3)
{
  chain[i] <- .3*chain[i-1] + err[i]
}
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
