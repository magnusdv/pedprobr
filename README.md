
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pedprobr <img src="man/figures/logo.png" align="right" height=140 />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/pedprobr)](https://CRAN.R-project.org/package=pedprobr)
[![](https://cranlogs.r-pkg.org/badges/grand-total/pedprobr?color=yellow)](https://cran.r-project.org/package=pedprobr)
[![](https://cranlogs.r-pkg.org/badges/last-month/pedprobr?color=yellow)](https://cran.r-project.org/package=pedprobr)
<!-- badges: end -->

## Introduction

The main content of **pedprobr** is an implementation of the
Elston-Stewart algorithm for pedigree likelihoods given marker
genotypes. It is part of the
[pedsuite](https://magnusdv.github.io/pedsuite/), a collection of
packages for pedigree analysis in R.

The **pedprobr** package does much of the hard work in several other
pedsuite packages:

- [**forrel**](https://github.com/magnusdv/forrel) : relatedness
  analysis and forensic pedigree analysis
- [**dvir**](https://github.com/magnusdv/dvir) : disaster victim
  identification
- [**paramlink2**](https://github.com/magnusdv/paramlink2) : parametric
  linkage analysis
- [**pedbuildr**](https://github.com/magnusdv/pedbuildr) : pedigree
  reconstruction

The workhorse of **pedprobr** is the `likelihood()` function, supporting
a variety of situations:

- autosomal and X-linked markers
- a single marker or two linked markers
- complex inbred pedigrees
- pedigrees with inbred founders
- mutation models

## Installation

To get the current official version of **pedprobr**, install from CRAN
as follows:

``` r
install.packages("pedprobr")
```

Alternatively, you can obtain the latest development version from
GitHub:

``` r
# install.packages("devtools") # install devtools if needed
devtools::install_github("magnusdv/pedprobr")
```

## Getting started

``` r
library(pedprobr)
#> Loading required package: pedtools
```

To set up a simple example, we first use **pedtools** utilities to
create a pedigree where two brothers are genotyped with a single SNP
marker. The marker has alleles `a` and `b`, with frequencies 0.2 and 0.8
respectively, and both brothers are heterozygous `a/b`.

``` r
# Pedigree with SNP marker
x = nuclearPed(nch = 2) |> 
  addMarker(geno = c(NA, NA, "a/b", "a/b"), afreq = c(a = 0.2, b = 0.8))

# Plot with genotypes
plot(x, marker = 1)
```

<img src="man/figures/README-pedplot-1.png" style="display: block; margin: auto;" />

The pedigree likelihood, i.e., the probability of the genotypes given
the pedigree, is obtained as follows:

``` r
likelihood(x, marker = 1)
#> [1] 0.1856
```

## Genotype probability distributions

Besides `likelihood()`, other important functions in **pedprobr** are:

- `oneMarkerDistribution()`: the joint genotype distribution at a single
  marker, for any subset of pedigree members
- `twoMarkerDistribution()`: the joint genotype distribution at two
  linked markers, for a single person

In both cases, the distributions are computed conditionally on any known
genotypes at the markers in question.

To illustrate `oneMarkerDistribution()` we continue our example from
above, and consider the following question: **What is the joint genotype
distribution of the parents, conditional on the genotypes of the
children?**

The answer is found as follows:

``` r
oneMarkerDistribution(x, ids = 1:2, verbose = F)
#>            a/a        a/b       b/b
#> a/a 0.00000000 0.01724138 0.1379310
#> a/b 0.01724138 0.13793103 0.2758621
#> b/b 0.13793103 0.27586207 0.0000000
```

The output confirms the intuitive result that the parents cannot both be
homozygous for the same allele. The most likely combination is that one
parent is heterozygous `a/b`, while the other is homozygous `b/b`.

The argument `output` controls how the output of
`oneMarkerDistribution()` is formatted. Instead of the default matrix
(or multidimensional array, if more than 2 individuals), we can also get
the distribution in table format:

``` r
oneMarkerDistribution(x, ids = 1:2, verbose = F, output = "table")
#>     1   2       prob
#> 1 a/a a/a 0.00000000
#> 2 a/b a/a 0.01724138
#> 3 b/b a/a 0.13793103
#> 4 a/a a/b 0.01724138
#> 5 a/b a/b 0.13793103
#> 6 b/b a/b 0.27586207
#> 7 a/a b/b 0.13793103
#> 8 a/b b/b 0.27586207
#> 9 b/b b/b 0.00000000
```

A third possibility is `output = "sparse"`, which gives a table similar
to the above, but with only the rows with non-zero probability.
