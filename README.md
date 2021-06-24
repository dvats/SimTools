# SimTools: Variability assessment for simulation methods in R

Leadership:
[Dootika Vats](http://dvats.github.io/), [James Flegal](https://faculty.ucr.edu/~jflegal/), [Galin Jones](http://users.stat.umn.edu/~galin/)

Contributors:
Gunjan Jalori

The development of this package has begun. Updates will be described in this README. For communication regarding issues in the code and other suggestions,  we can try and use the GitHub Issues interface.


**Run**: `runme.R`

---------------------
 February 24, 2021

Added functions `getCI`, `addCI`, and `boxCI` to get and add simultaneous CIs to existing plots.

Additionally, now the `Smcmc` class supports multiple chains. The
input is a list of chains, and calling `Smcmc` does some processing
to make a `stacked` version of the chains.

---------------------
 December 7, 2020 

`Siid` class created. This corresponds to "IID simulation objects". The `boxplot` and `plot` commands work directly on this and produces boxplots and density plots. 

---------------------
 December 4, 2020 

`Smcmc` class created. This corresponds to "MCMC simulation objects". The `plot` command works directly on this and produces density plots with the desired quantiles and means indicated.

---------------------

## Funding

Grateful for support by SERB, DST, India (MSC/2020/000165).
