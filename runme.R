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
    ratio <- prod((dnorm(prop)/dnorm(chain[i-1, ])))
    alpha <- min(1, ratio)
    if(runif(1) <= alpha) {
      chain[i, ] <- prop
    }
    else {
      chain[i, ] <- chain[i - 1, ]
    }
  }
  return(chain)
}

# Make chain
chain <- MakeChain()


# install SimTools package
# On the terminal run: R CMD INSTALL SimTools_1.0-0.tar.gz
library(SimTools)


# First, we tell the package that chain is an MCMC process
# using the class `Smcmc`. A few things, setting:
# batch.size = TRUE calculates the batch.size for this
#         chain and uses it in the future, and not wasting
#         time in the future
# varnames = colnames(chain) by default

out1 <- Smcmc(chain)
out2 <- Smcmc(chain, batch.size = TRUE)

# View the attributes of these classes
attributes(out1)
attributes(out2)   # attribute $size is stored in this


# Density plots are called by `plot`
# This is because plot(Smcmc object) calls an underlying
# plot.Smcmc function
plot(out1)
plot(out2)

# flexibility in this function
plot(out2, Q = c(.05, .20))  # change quantiles
plot(out2, mean = FALSE)  # no mean
plot(out2, lty = 2, lwd = 2)  # change line type
plot(out2, main = "Density")  # change title


# Repeating this again for a larger size of chain
chain <- MakeChain(p = 13, n = 1e5)


out3 <- Smcmc(chain)
out4 <- Smcmc(chain, batch.size = TRUE) # takes a second longer

# Interactive plot windows
plot(out3)
plot(out4) # slightly faster everytime I call this, since batchsize is not being calculated



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




