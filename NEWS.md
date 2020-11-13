# pedprobr 0.4.0

## Breaking changes

* `likelihood()` has been refactored, and the functionality for calculating the likelihood of two linked markers is moved to a separate function, `likelihood2()`. 

## New features

* Theta-correction is implemented in likelihood calculations, through the argument `theta` of `likelihood()`.

* `merlin()` and `likelihoodMerlin()` have been overhauled.

* A new function `checkMerlin()` checks if MERLIN is installed and available.

* `merlin()` now has an argument `linkageMap` facilitating analysis of linked markers.

* `likelihoodMerlin()` gains arguments `rho` and `logbase`.

* New function `lumpAlleles()`.


# pedprobr 0.3.0

## New features

* The `likelihood()` function now handles multiple marker inputs. For example, 
the call `likelihood(x, 1:2)` results in a vector of length 2 with the likelihoods 
of the first two markers attached to `x`. 

* Recombination parameter `theta` is renamed to `rho` everywhere, to align 
with other ped suite packages. `theta` still works, though.

## Other changes

* Many input checks have been added, giving more useful error messages.

* The data structure used in likelihood calculations is simplified and more 
efficient.

* Likelihood computation of two linked markers (and therefore also 
`twoMarkerDistribution()`) is much faster now.


# pedprobr 0.2.0

## New features

* New general MERLIN wrapper, `merlin()`.

* Both `likelihoodMerlin()` and the new `merlin()` now accepts ped lists.

## Bug fixes

* Fixed a bug affecting lumped mutation models.

# pedprobr 0.1.0

* Initial CRAN release.
