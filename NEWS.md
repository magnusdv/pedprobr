# pedprobr 0.9.6

* Methods for allele lumping have been rewritten and expanded. In particular, markers with unlumpable mutation models (in the Kemeny-Snell sense) can now, in certain cases, be lumped using the *special lumping* recently implemented in **pedmut**. To activate this, set `special = TRUE` in `likelihood()` and `likelihood2()`. By default, `special` is FALSE in the former and TRUE in the latter.

* `likelihood()` and `likelihood2()` also gain a new argument `alleleLimit`, limiting the number of alleles in unlumpable markers. If the limit is exceeded, the mutation model is replaced with the simpler `equal` model, which always permits lumping. By default, this behaviour is disabled (i.e. `alleleLimit = Inf`).

* Improve structure of `likelihood()` and `likelihood2()`, to avoid redundant calculations. 

* Skip unneeded checks in `...MarkerDistribution()`

# pedprobr 0.9.5

* `oneMarkerDistribution()` now fully supports pedigrees with multiple components.

* `oneMarkerDistribution()` has a new argument `output` taking values `"array"` (default, as before), `"table"`, or `"sparse"`. Both `"table"` and `"sparse"` produce data frames where each row corresponds to a genotype combination. The `"sparse"` format only includes combinations with non-zero probability.

* In `oneMarkerDistribution()`, the argument `partialmarker` has been renamed to `marker`. The old name still works as an alias, but will be removed in a future version. While the old argument had no default value, the new defaults to the first attached marker. This simplifies the call in many cases, for example `singleton("A") |> addMarker() |> oneMarkerDistribution("A")`.

* In `twoMarkerDistribution()`, the arguments `partialmarker1` and `partialmarker2` have been renamed to `marker1` and `marker2`, respectively.

* Fixed bug affecting likelihood calculations in pedigrees with partial genotypes (e.g. `"1/-"`) in founders.

* The function `allGenotypes()` is ~2-4 times faster due to a better implementation.

* Updated dependencies: **pedtools** v2.6.0, **pedmut** v0.7.0.


# pedprobr 0.9.4

* Fixed rare bug in the peeling algorithm manifesting with reversed peeling order.

* **pedtools** v2.2.0 is now required.

# pedprobr 0.9.3

* Fix bug in `likelihood()` affecting singleton with partial genotype (e.g. "1/-").

* Minor optimisations in `reduceAlleles()` and `startdata_M()`.

* More efficient `oneMarkerDistribution()` and `twoMarkerDistribution()` in cases without mutation models.

* README updates.


# pedprobr 0.9.2

The version implements several improvements in the peeling algorithm, mostly invisible to end users (except that likelihoods calculations will be noticeably faster, and use much less memory, in many cases).

* Less memory footprint for untyped founders with only 1 child: Sample alleles directly, avoid complete set of genotypes.

* Avoid costly 3-dimensional arrays in peeling calculations.

* More efficient genotype elimination in likelihood calculations involving markers without mutation modelling. As a result, this is now enforced (rather than optional) everywhere, and the `eliminate` argument is deprecated.

* Added citation info.

# pedprobr 0.8.0

## New features

* Mutation models are now supported in `likelihood2()`, providing likelihood calculations for pairwise linked markers.

* New functions `haldane()` and `kosambi()`.

## Other

* Update README.

* Update package doc.

* Include license (GPL >= 2).


# pedprobr 0.7.1

* `reduceAlleles()` is faster in some cases, due to the new `pedmut::lumpedModel()`.

* A few minor speedups and code improvements.


# pedprobr 0.7.0

## New features

* `merlin()` gains argument `checkpath`. Set this to FALSE to avoid redundant checks of MERLIN availability.

* `likelihoodMerlin()` gains argument `perChrom`, used to parse the chromosome-wise likelihoods from the MERLIN output.

* `lumpAlleles()` gains argument `always`. By default lumping is skipped for markers where all individuals are genotyped, but this can cause problems e.g. with MERLIN which operates with an upper limit of alleles.

* Some efforts are done to check for (and warn about) underflow/overflow in MERLIN results.


# pedprobr 0.6.1

## Bug fixes

* Fixed a regression error involving the external program MINX (part of MERLIN). The new version of `checkMerlin()` by default checks that both `merlin` and `minx` are available on the system, and that they both come from the latest MERLIN version.


# pedprobr 0.6.0

## Breaking changes

* `pedprobr` now depends on R 4.1 and `pedtools` 1.1.0.

* The deprecated argument `loop_breakers` has been completely removed (renamed to `loopBreakers`)

## New features

* `likelihood.ped()` gains a new argument `lump`, which activates allele lumping. This is TRUE by default, and should remain so in most cases.

* Using theta correction in pedigrees with inbred founders now gives an error.

* The internal `peelingProcess()` now checks for unbroken loops in the pedigree. (Of interest for developers mainly.)

## Bug fixes

* Fixed bug affecting likelihood computations with theta correction.

## Other changes

* When calling MERLIN, the window-centric and unnecessary `.exe` extensions has been removed. 

* Many examples have been rewritten in simpler code, taking advantage of `pedtools::addMarker` and the pipe `|>`.


# pedprobr 0.5.0

## New features

* Removing/disabling mutation models is now easier, with `setMutationModel(..., model = NULL)`.

## Bug fixes

* Fixed bug in `setMutationModel()` affecting lists of multiple pedigrees.


# pedprobr 0.4.0

## Breaking changes

* `likelihood()` has been refactored, moving the treatment of (two) linked markers to a separate function, `likelihood2()`. 

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
with other pedsuite packages.

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
