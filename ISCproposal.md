# SimTools: Variability assessment for simulation methods in R
[Dootika Vats](http://home.iitk.ac.in/~dootika/), [James Flegal](https://faculty.ucr.edu/~jflegal/), [Galin Jones](http://users.stat.umn.edu/~galin/)

An R package for implementing the Portkey two-coin Bernoulli factory of Vats et. al. (2020). The main function is `portkey` which is meant to be used within an MCMC algorithm to determine whether the move should be accepted or rejected. The regular two-coin algorithm can also be implemented using `twocoin` which just calls `portkey` for $\beta = 1$. 


# Installation
This R package is not on CRAN and is hosted on GitHub only

To download this development repo,  through the the `devtools` package:

```{r}
# install.packages("devtools")
library(devtools)
devtools::install_github("dvats/portkey")
```

