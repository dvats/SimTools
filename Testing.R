##########################################
## This document is a running file
## in place of a vignette for now.
##
## We run an MCMC example and an IID
## Monte Carlo example and demonstrate the 
## plots that are available.
##########################################
# install SimTools package
# download the tarball from GitHub and then
# on the terminal run: R CMD INSTALL SimTools_1.0-0.tar.gz
library(SimTools)
set.seed(1)
############################
# Smcmc class
############################

# Making a Markov chain
# RW MH for a standard normal multivariate target
MakeChain <- function(p, phi = diag(rep(.99, p)), nsim = 1e5,  omega = diag(rep(1, p)), last = NULL)
{
  # sig <- diag(rep(rho, p))
  chain <- matrix(0, nrow = nsim, ncol = p)	
  decomp <- svd(omega)
  omega.root <- as.matrix(decomp$v) %*% as.matrix(diag((decomp$d)^(1/2), nrow = p, ncol = p)) %*% as.matrix(t(decomp$u))
  if(is.null(last)){
    chain[1, ] <- rnorm(p)
  } else {
    chain[1, ] <- last
  }
  # sigma.root%*%rnorm(p)
  for(i in 2:nsim)
  {
    chain[i, ] <- as.matrix(phi) %*% as.matrix(chain[i-1, ]) + omega.root %*% as.matrix(rnorm(p))
  }
  return(chain)
}

chain <- MakeChain(p = 4)

# Make ACF plot
acfplot(chain)

# First, we tell the package that chain is an MCMC process
# using the class `Smcmc`. A few things, setting:
# batch.size = TRUE calculates the batch.size for this
#         chain and uses it in the future, and not wasting
#         time in the future
# varnames = colnames(chain) by default


out.one <- Smcmc(chain)
acfplot(out.one)
l = acfplot(out.one,lag.max = 600, plot=F)
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

## For making traceplot of chains, 
## it has various features

traceplot(out.one)
traceplot(out.one, legend = F)
traceplot(out.one, main = "Traceplots")
traceplot(out.one, main = "Traceplots", legend = F)


# Density plots are called by `densityplot`
# This is because densityplot(Smcmc object) calls an underlying
# densityplot function
densityplot(out.one)

# rug plots can also be seen
densityplot(out.one, rug = TRUE)

# flexibility in this function
densityplot(out.one, Q = c(.05, .20))  # change quantiles
densityplot(out.one, mean = FALSE)  # no mean

### Summary plot
plot(out.one)
summary(out.one)
###********************************************************************############
###********************************************************************###########

### Multiple chains
# The package now supports multiple chain implementations

chain1 <- MakeChain(p = 4,n=1e4)
chain2 <- MakeChain(p=4,n=1e4)
chain3 <- MakeChain(p =4,n=1e4)
chain4 <- MakeChain(p = 4,n=1e4)
colnames(chain3) <- c("alpha","beta","gamma","delta")
out.four <- Smcmc(list(chain1, chain2, chain3,chain4))

str(out.four)
acfplot(out.four)
traceplot(out.four, fast = F)
traceplot(out.four, main = "Hui Hui Hui")
densityplot(out.four)
plot(out.four)

# access list of chains by out.three$chains 
 chain1 <- MakeChain(p = 5,n=1e3)
 chain2 <- MakeChain(p=5,n=1e3)
 chain3 <- MakeChain(p =5,n=1e3)
 chain4 <- MakeChain(p = 5,n=1e3)
 out.five <- Smcmc(list(chain1, chain2, chain3,chain4))

# Make ACF plot, with maximum lag = 300
acfplot(out.five,lag.max=300)

# calculated batch size is the average of 
# individual chain batch sizes
out.five$b.size

# dimension is typically not #3k \times 4
dim(out.five$stacked)

# plots the density from the stacked estimates
densityplot(out.five, main = "Density Plots of All variate")

plot(out.five, which = 1)
summary(out.five)


par(mfrow = c(3,2))
## Testing of Layouts
traceplot(out.five)
traceplot(out.five,which = 1:5,main = "Traceplot")

acfplot(out.five,which = 1:5)
acfplot(out.five,which = 1:5,main = "ACF plot")

densityplot(out.five,which = 1:5)
densityplot(out.five,which = 1:5,main = "Density Plot")

par(mfrow = c(1,1))
traceplot(out.five,which = 1:2,main = "Traceplot")
acfplot(out.five,which = 1:2,main = "ACF plot")
densityplot(out.five,which = 1:2,main = "Density Plot")
par(mfrow = c(1,1))
plot(out.five)


chain1 <- MakeChain(p = 7,n=1e3)
chain2 <- MakeChain(p=7,n=1e3)
chain3 <- MakeChain(p =7,n=1e3)
chain4 <- MakeChain(p = 7,n=1e3)
out.sev <- Smcmc(list(chain1, chain2, chain3,chain4),varnames = c("GJ", "DV","JF","KG","GJ","SP","Unkonwn"))
acfplot(out.sev,main = "ACF Plot")
traceplot(out.sev)
densityplot(out.sev)

### Testing of large chains
chain7 <- MakeChain(p =10,nsim=1e5)
chain8 <- MakeChain(p =10,nsim=1e5)
chain9 <- MakeChain(p =10,nsim=1e5)
chain10 <- MakeChain(p=10,nsim=1e5)
chain11 <- MakeChain(p =10,nsim=1e5)
outlarge.five <- Smcmc(list(chain7, chain8, chain9,chain10,chain11))
densityplot(outlarge.five,main = "Density Plot")
traceplot(outlarge.five,main = "Trace Plot", which = 1:4)
acfplot(outlarge.five,main = "ACF Plot",which = 1:4)
plot(outlarge.five)


### Testing of large chains with dimesions greater than 10
chain11 <- MakeChain(p =12,n=1e5)
chain12 <- MakeChain(p =12,n=1e5)
chain13 <- MakeChain(p =12,n=1e5)
chain14 <- MakeChain(p=12,n=1e5)
chain15 <- MakeChain(p =12,n=1e5)
out12.five <- Smcmc(list(chain12, chain13, chain14,chain15,chain11))
densityplot(out12.five,main = "Density Plot")
densityplot(out12.five,main = "Density Plot", fast = F)
traceplot(out12.five,main = "Trace Plot")
acfplot(out12.five,main = "ACF Plot")
plot(out12.five)
summary(out12.five)



chain11 <- MakeChain(p =15,n=7*1e4)
chain12 <- MakeChain(p =15,n=7*1e4)
chain13 <- MakeChain(p =15,n=7*1e4)
chain14 <- MakeChain(p=15,n=7*1e4)
chain15 <- MakeChain(p =15,n=7*1e4)
out15.five <- Smcmc(list(chain12, chain13, chain14,chain15,chain11))
densityplot(out15.five,main = "Density Plot")
traceplot(out15.five,main = "Trace Plot")
acfplot(out15.five,main = "ACF Plot")
summary(out15.five)



chain11 <- MakeChain(p =14,n=1e5)
chain12 <- MakeChain(p =14,n=1e5)
chain13 <- MakeChain(p =14,n=1e5)
chain14 <- MakeChain(p=14,n=1e5)
chain15 <- MakeChain(p =14,n=1e5)
out14.five <- Smcmc(list(chain12, chain13, chain14,chain15,chain11))
densityplot(out14.five,main = "Density Plot")
traceplot(out14.five,main = "Trace Plot")
acfplot(out14.five,main = "ACF Plot")
summary(out14.five)

##Testing of Discrete State Space compatibility
c1 <- sample(1:100, 10000 ,replace = T)
c2 <- sample(1:100, 10000 ,replace = T)
out <- as.Smcmc(list(cbind(c1,c2)))
densityplot(out)


c1 <- as.matrix(rbinom(50,n= 10000, p = 0.7))
c2 <- as.matrix(rbinom(50,n= 10000,p=0.8))
c3 <- as.matrix(rpois(10000,lambda=10))
c4 <- as.matrix(rpois(10000,lambda=10))
c5 <- as.matrix(rbinom(50,n= 100000, p = 0.7))
c6 <- as.matrix(rbinom(50,n= 100000,p=0.8))
c7 <- as.matrix(rpois(10000,lambda=10))
c8 <- as.matrix(rpois(10000,lambda=10))
out <- as.Smcmc(cbind(c1,c2,c3,c4),cbind(c5,c6,c7,c8))
densityplot(out)
plot(out)
summary(out)


## Mixed chain(2 dimensions Continuous and 1 Dimensions Discrete)
chain1 = MakeChain(p=2, nsim=10000)
chain2 = MakeChain(p=2, nsim=10000)
chain3 = MakeChain(p=2, nsim=10000)
chain4 = MakeChain(p=2, nsim=10000)

out <- as.Smcmc(list(cbind(chain1,c1),cbind(chain2,c2),cbind(chain3,c3),cbind(chain4,c4)))
densityplot(out)
plot(out)
summary(out)
