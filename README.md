# SimTools: Variability assessment for simulation methods in R

Leadership:
[Dootika Vats](http://dvats.github.io/), [James Flegal](https://faculty.ucr.edu/~jflegal/), [Galin Jones](http://users.stat.umn.edu/~galin/)

Contributors:
Gunjan Jalori, Siddharth Pathak

Updates will be described in this README. For communication regarding issues in the code and other suggestions,  please employ the issues feature on GitHub.


**Run**: `runme.R`

---------------------
Install from this branch using:
`devtools::install_github("dvats/SimTools", ref = "Siddharth-Pathak")`

---------------------

---------------------
Next update goal:

---------------------
 August 25, 2023

The following were accomplished
- multiple chain compatibility
- summary function for output
- compatibility for discrete state space
- make efficient trace plots
- New functions to use are: `densityplot()`, `acfplot()`, `traceplot()`, `plot()`, `summary()`.

---------------------
 October 3, 2021

Added function `ACF` that plots ACF plots for Markov chains for each
component. The function plots globally-centered ACFs, and can thus be
used for multiple chain setups

Other clean-up to the code has also been done.

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
