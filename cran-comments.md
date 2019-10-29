## Resubmission
This is a resubmission. In this version I have

As requested:
* Removed unconditional print()/cat() statements [in likelihood.ped() and merlin()]

In addition:
* Added a URL in Description file
* Replaced T/F with TRUE/FALSE everywhere
* Fixed a bug in oneMarkerDistribution()/twoMarkerDistribution (and added unit test)

I have re-run devtools::check_win_devel() and rhub::check_for_cran(), 
which gave 1 NOTE (New submission).

## Resubmission
This is a resubmission. In this version I have

* In DESCRIPTION, single-quoted the package name 'pedprobr', as requested.

* Added \value fields to allGenotypes.Rd and genoCombinations.Rd, as requested.
I also did some modifications of code and documentation of these two functions.
All checks were re-run.


## Test environments
* local Windows 10 install, R 3.6.1
* devtools::check_win_devel()
* rhub::check_for_cran()

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
