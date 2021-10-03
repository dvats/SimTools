##########################################
## This document is a running file
## in place of a vignette for now.
##
## We run an MCMC example and an IID
## Monte Carlo example and demonstrate the 
## plots that are available.
##########################################
set.seed(1)


############################
# Smcmc class
############################

# Making a Markov chain
# RW MH for a standard normal multivariate target
MakeChain <- function(p = 4, n = 1e3, h = .5)
{
  chain <- matrix(0, nrow = n,ncol = p)

  for (i in 2:n) 
  {
    prop <- chain[i-1, ] + rnorm(p, mean = 0, sd = h)
    log.ratio <- sum(dnorm(prop, log = TRUE) - dnorm(chain[i-1, ], log = TRUE))
    
    if(log(runif(1)) < log.ratio) 
    {
      chain[i, ] <- prop
    }
    else {
      chain[i, ] <- chain[i - 1, ]
    }
  }
  colnames(chain) <- c("Comp 1", "Comp 2", "Comp 3", "Comp 4")
  return(chain)
}

# Make chain
chain <- MakeChain()


# install SimTools package
# download the tarball from GitHub and then
# on the terminal run: R CMD INSTALL SimTools_1.0-0.tar.gz
library(SimTools)


# Make ACF plot
ACF(chain)

# First, we tell the package that chain is an MCMC process
# using the class `Smcmc`. A few things, setting:
# batch.size = TRUE calculates the batch.size for this
#         chain and uses it in the future, and not wasting
#         time in the future
# varnames = colnames(chain) by default

out.one <- Smcmc(chain)

# the chain is accessed via
head(out.one$chains[[1]])

# the batch size chosen
# for variance estimation
out.one$b.size


# the stacked chain stores a cleaned 
# and "stacked" version of all chains
# it stacks multiple chains so that 
# n = ab for all chain
head(out.one$stacked)


# Density plots are called by `plot`
# This is because plot(Smcmc object) calls an underlying
# plot.Smcmc function
plot(out.one)

# rug plots can be removed
plot(out.one, rug = FALSE)

# flexibility in this function
plot(out.one, Q = c(.05, .20))  # change quantiles
plot(out.one, mean = FALSE)  # no mean



### Multiple chains
# The package now supports multiple chain implementations

chain1 <- MakeChain()
chain2 <- MakeChain()
chain3 <- MakeChain()

out.three <- Smcmc(list(chain1, chain2, chain3))
str(out.three)
# access list of chains by out.three$chains 


# Make ACF plot
ACF(out.three)


# calculated batch size is the average of 
# individual chain batch sizes
out.three$b.size

# dimension is typically not #3k \times 4
dim(out.three$stacked)

# plots the density from the stacked estimates
plot(out.three)


# plotting individual chain densities versus combined
# using the getCI and the addCI functions
# we can add appropriate simulateneous intervals
par(mfrow = c(2,2))

foo1 <- Smcmc(chain1)
plot(density(foo1$stacked[,1]), main = "Markov chain 1")
CIs <- getCI(foo1)
addCI(foo1, CIs, component = 1)

foo2 <- Smcmc(chain2)
plot(density(foo2$stacked[,1]), main = "Markov chain 2")
CIs <- getCI(foo2)
addCI(foo2, CIs, component = 1)

foo3 <- Smcmc(chain3)
plot(density(foo3$stacked[,1]), main = "Markov chain 3")
CIs <- getCI(foo3)
addCI(foo3, CIs, component = 1)

plot(density(out.three$stacked[,1]), main = "Combined")
CIs <- getCI(out.three)
addCI(out.three, CIs, component = 1)


# Typically, we would expect people to work
# with the stacked chains. So separate plots
# for each component can be made

# get CIs for each component
CIs <- getCI(out.three)

for(i in 1:4)
{
  plot(density(out.three$stacked[,i]), main = paste("Component", i), lty = 2)
  addCI(out.three, CIs, component = i)
}



############################
# Siid class
############################
n <- 5e2
p <- 5
sim <- matrix(rnorm(n*p), ncol = p, nrow = n) 

out <- Siid(sim)  # make into Siid object)
plot(out)  # since Siid object, mcmcse is not used in variance estimation
boxplot(out)


# Picks up the column names if present
colnames(sim) <- paste("Comp", 1:5)
out <- Siid(sim)  # make into Siid object)
plot(out)  # since Siid object, mcmcse is not used in variance estimation
boxplot(out)

# only first component
boxplot(sim[,1], col = "white")
boxCI(out, getCI(out, Q = c(.25, .5, .75), iid = TRUE), component = 1)
# need to change getCI so that there is an automatic Q option in it


