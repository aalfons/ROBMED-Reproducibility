# A Robust Bootstrap Test for Mediation Analysis


## About ROBMED

The robust bootstrap test ROBMED for mediation analysis is less sensitive to deviations from model assumptions (such as outliers or heavily tailed distributions) than the standard bootstrap test of Preacher & Hayes ([2004](http://dx.doi.org/10.3758/BF03206553), [2008](http://dx.doi.org/10.3758/BRM.40.3.879)).  ROBMED utilizes the robust MM-regression estimator ([Yohai, 1987](https://projecteuclid.org/euclid.aos/1176350366)) instead of the OLS estimator for regression, and runs bootstrap tests with the fast and robust bootstrap methodology ([Salibián-Barrera & Zamar, 2002](https://projecteuclid.org/euclid.aos/1021379865); [Salibián-Barrera & Van Aelst, 2008](https://doi.org/10.1016/j.csda.2008.05.007)).

More information can be found in our manuscript:

Alfons, A., Ates, N.Y., & Groenen, P.J.F. (2021). A Robust Bootstrap Test for
Mediation Analysis. Organizational Research Methods, accepted for publication.

An earlier version of our manuscript is available as an 
[internal research report](https://repub.eur.nl/pub/109594/) of the 
Erasmus Research Institute of Management (ERIM).


## Reproduce results

This repository provides a collection of [R](https://CRAN.R-project.org/) 
scripts to reproduce all examples, simulations and figures in our manuscript 
and the supplementary report.  

The easiest way to reproduce the results is to clone this repository with 
[RStudio](https://rstudio.com/products/rstudio/download/).  Running the 
scripts within the resulting RStudio project ensures that there are no issues 
with file paths for storing or reading results, or for producing files 
containing plots.  In addition, the RStudio project uses 
[packrat](https://rstudio.github.io/packrat/) to make sure that the correct 
versions of all required R packages are used.  After opening the RStudio 
project for the fist time, please type `packrat::restore()` on the R command 
line to retrieve the correct versions of all required packages.

Please note that this repository is rather large because it also contains R 
data files with all simulation results.  This way, if you only want to quickly
reproduce the figures with simulation results, you do not actually need to run 
the simulations first.
