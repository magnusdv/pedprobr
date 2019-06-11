#' Probability distribution for a single marker
#'
#' Computes the genotype probability distribution of one or several pedigree
#' members, conditional on known genotypes for the marker.
#'
#' @param x A `ped` object.
#' @param ids A numeric with ID labels of one or more pedigree members.
#' @param partialmarker Either a `marker` object or the name (or index) of a
#'   marker attached to `x`.
#' @param grid.subset (Optional; not relevant for most users.) A numeric matrix
#'   describing a subset of all marker genotype combinations for the `ids`
#'   individuals. The matrix should have one column for each of the `ids`
#'   individuals, and one row for each combination: The genotypes are described
#'   in terms of the matrix `M = allGenotypes(n)`, where `n` is the number of
#'   alleles for the marker. If the entry in column `j` is the integer `k`, this
#'   means that the genotype of individual `ids[j]` is row `k` of `M`.
#' @param loop_breakers (Only relevant if the pedigree has loops). A vector with
#'   ID labels of individuals to be used as loop breakers. If NULL (default)
#'   loop breakers are selected automatically. See [breakLoops()].
#' @param eliminate A non-negative integer, indicating the number of iterations
#'   in the internal genotype-compatibility algorithm. Positive values can save
#'   time if `partialmarker` has many alleles.
#' @param verbose A logical.
#' @return A named `k`-dimensional array, where `k = length(ids)`, with the
#'   joint genotype distribution for the `ids` individuals. The probabilities
#'   are conditional on the known genotypes and the allele frequencies of
#'   `partialmarker`.
#' @author Magnus Dehli Vigeland
#' @seealso [twoMarkerDistribution()]
#'
#' @examples
#' library(pedtools)
#' x = nuclearPed(2)
#'
#' # Empty SNP marker
#' snp = marker(x, alleles=1:2)
#' oneMarkerDistribution(x, ids=3:4, partialmarker=snp)
#'
#' #### Different example for the same pedigree. A marker with 4 alleles:
#' m2 = marker(x, `3`=c('C','D'), `4`=c('C','D'), alleles=LETTERS[1:4])
#' oneMarkerDistribution(x, ids=1, partialmarker=m2)
#'
#' # Same as above, but computing only the cases where individual 1 is heterozygous.
#' # (The numbers 5:10 refer to the 6 last rows of allGenotypes(4),
#' # which contain the heterozygous genotypes.)
#' oneMarkerDistribution(x, ids=1, partialmarker=m2, grid.subset=matrix(5:10, ncol=1))
#'
#' #### Expanding on the previous example:
#' # Joint genotype probabilities of the parents, but including only the combinations
#' # where the father is heterozygous and the mother is homozygous:
#' grid = expand.grid(5:10, 1:4)
#' oneMarkerDistribution(x, ids=1:2, partialmarker=m2, grid.subset=grid)
#'
#' #### Something else:
#' # The genotype distribution of an individual (id=5) whose half cousin (id=9) is
#' # homozygous for a rare allele.
#' y = halfCousinPed(1)
#' snp = marker(y, `9`="a", alleles=c("a", "b"), afreq=c(0.01, 0.99))
#' oneMarkerDistribution(y, ids=5, partialmarker=snp)
#'
#' #### X-linked example: TODO
#'
#' @export
oneMarkerDistribution <- function(x, ids, partialmarker, grid.subset = NULL,
                                  loop_breakers = NULL, eliminate = 0, verbose = TRUE) {
  if(!is.ped(x))
    stop2("Input is not a `ped` object")
  if(!isCount(eliminate, minimum = 0))
    stop2("`eliminate` must be a nonnegative integer")

  m = partialmarker

  if (!is.marker(m)) {
    if(length(m) != 1)
      stop2("`partialmarker` must have length 1")
    m = getMarkers(x, markers = m)[[1]]
  }

  if (!is.null(x$LOOP_BREAKERS))
    stop2("`ped` objects with pre-broken loops are not allowed as input to `oneMarkerDistribution()`")

  alleles = alleles(m)
  onX = is_Xmarker(m)

  if (verbose) {
    cat(sprintf("Partial marker (%s):\n", ifelse(onX, "X-linked", "autosomal")))
    print(m)
  }

  starttime = Sys.time()

  allgenos = allGenotypes(nAlleles(m))

  # Compute grid before loop breaking (works better with eliminate2)
  if (is.null(grid.subset))
    grid.subset = genoCombinations(x, m, ids, make.grid=T)
  else
    grid.subset = as.matrix(grid.subset)

  if (x$UNBROKEN_LOOPS) {
    x = breakLoops(setMarkers(x, m), loop_breakers = loop_breakers, verbose = verbose)
    m = x$markerdata[[1]]
  }

  # Ensure peeling order is set (to avoid redundant computation)
  if(is.null(attr(x, "PEELING_ORDER")))
    attr(x, "PEELING_ORDER") = peelingOrder(x)

  int.ids = internalID(x, ids)
  gt.strings = paste(alleles[allgenos[, 1]], alleles[allgenos[, 2]], sep = "/")

  geno.names = if(onX) list(alleles, gt.strings)[getSex(x, ids)]
               else rep(list(gt.strings), length(ids))

  marginal = likelihood(x, marker1 = m, eliminate = eliminate)
  if (marginal == 0)
      stop2("Partial marker is impossible")
  probs = array(0, dim = lengths(geno.names, use.names = F), dimnames = geno.names)
  probs[grid.subset] = apply(grid.subset, 1, function(allg_rows) {
      m[int.ids, ] = allgenos[allg_rows, ]
      likelihood(x, marker1 = m, eliminate = eliminate)
  })

  res = probs/marginal
  if (verbose) {
    cat("==============================\n\n")
    cat("Analysis finished in ", round(Sys.time() - starttime,2), " seconds\n")
    if(length(ids)==1)
      cat("\nGenotype probability distribution for individual ", ids, ":\n", sep="")
    else
      cat("\nJoint genotype probability distribution for individuals ",
          toString(ids), ":\n", sep = "")

    print(round(res, 4))
    return(invisible(res))
  }
  else
    res
}
