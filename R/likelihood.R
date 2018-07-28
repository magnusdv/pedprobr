#' Pedigree likelihood
#'
#' Calculates the likelihood of a pedigree (or a list of pedigrees) given
#' genotypes for a marker or a pair of linked markers.
#'
#' The `likelihood` function is the heart of `pedprobr`. It implements the
#' Elston-Stewart algorithm, and works in a variety of situations:
#'
#' * complex pedigrees with multiple layers inbreeding
#' * autosomal and X-linked markers
#' * a single marker or two linked markers
#' * markers with mutation models
#'
#'
#' @param x a `ped` object, a `singleton` object, or a list of such objects.
#' @param marker1 a [marker()] object compatible with `x`. If `x` is a list,
#'   then `marker1` must be a list of corresponding `marker` objects.
#' @param marker2 either NULL, or a [marker()] object compatible with `x`. See
#'   Details.
#' @param theta the recombination rate between `marker1` and `marker2`. To make
#'   biological sense `theta` should be between 0 and 0.5.
#' @param eliminate mostly for internal use: a non-negative integer indicating
#'   the number of iterations in the internal genotype-compatibility algorithm.
#'   Positive values can save time if the number of alleles is large.
#' @param logbase a numeric, or NULL. If numeric the log-likelihood is returned,
#'   with `logbase` as basis for the logarithm.
#' @param loop_breakers a vector of ID labels indicating loop breakers. If NULL
#'   (default), automatic selection of loop breakers will be performed. See
#'   [breakLoops()].
#' @param startdata for internal use.
#' @param verbose a logical
#' @param total a logical; if TRUE, the product of the likelihoods is returned,
#'   otherwise a vector with the likelihoods for each pedigree in the list.
#' @param \dots further arguments.

#' @return The likelihood of the data. If the parameter `logbase` is a
#' positive number, the output is `log(likelihood, logbase)`.

#' @author Magnus Dehli Vigeland
#'
#' @examples
#'
#' library(pedtools)
#'
#' # likelihood of inbred pedigree (grandfather/granddaughter incest)
#' x = nuclearPed(father="grandfather", children="a")
#' x = addDaughter(x, "a", id="granddaughter")
#' x = addChildren(x, father="grandfather", mother="granddaughter", nch=1, id="child")
#' m = marker(x, grandfather=1, child=1:2)
#'
#' plot(x, m)
#' lik = likelihood(x, m)
#'
#' stopifnot(lik==0.09375)
#'
#' @export
likelihood = function(x, ...) UseMethod("likelihood", x)


#' @export
#' @rdname likelihood
likelihood.ped = function(x, marker1, marker2 = NULL, theta = NULL, startdata = NULL,
                          eliminate = 0, logbase = NULL, loop_breakers = NULL, verbose=TRUE, ...) {

  twolocus = !is.null(marker2)
  if(twolocus && is.null(theta))
    stop2("Argument `theta` is missing")

  if(missing(marker1) || is.null(marker1))
    stop2("Argument `marker1` is missing")


  if (!is.marker(marker1)) {
    if(length(marker1) != 1)
      stop2("Length of `marker1` must be 1")
    marker1 = getMarkers(x, markers = marker1)[[1]]
  }

  if (twolocus && !is.marker(marker2)) {
    if(length(marker2) != 1)
      stop2("Length of `marker2` must be 1")
    marker2 = getMarkers(x, markers = marker2)[[1]]
  }

  marker1 = .reduce_alleles(marker1)
  marker2 = .reduce_alleles(marker2)  # unchanged if NULL

  if (x$UNBROKEN_LOOPS) {
    if(verbose) message("Tip: To optimize speed, consider breaking loops before calling 'likelihood'. See ?breakLoops.")
    m = list(marker1)
    m[[2]] = marker2  # no effect if NULL
    x = breakLoops(setMarkers(x, m), loop_breakers = loop_breakers, verbose = verbose)
    marker1 = x$markerdata[[1]]
    if (twolocus)
      marker2 = x$markerdata[[2]]
  }

  if (is.null(startdata)) {
    inform = .informative(x, marker1, marker2)
    inform_subnucs = inform$subnucs
    attr(x, "treat_as_founder") = inform$newfounders
    if(!twolocus)
      startdata = startdata_M(x, marker = marker1, eliminate = eliminate)
    else
      startdata = startdata_MM(x, marker1 = marker1, marker2 = marker2, eliminate = eliminate)
  }

  if (attr(startdata, "impossible"))
    return(ifelse(is.numeric(logbase), -Inf, 0))

  Xchrom = is_Xmarker(marker1)
  mutmat = attr(marker1, "mutmat")
  PEEL = choosePeeler(twolocus, theta, Xchrom, x$SEX, mutmat)

  dat = startdata

  if (is.null(dups <- x$LOOP_BREAKERS)) {
    for (sub in inform_subnucs) {
      dat = PEEL(dat, sub)
      if (sub$link > 0 && attr(dat, "impossible"))
        return(ifelse(is.numeric(logbase), -Inf, 0))
    }
    likelihood = dat
  }
  else {
    origs = dups[, 1]
    copies = dups[, 2]

    # Utility for comparing matrices, needed below:
    two2one = function(matr) {
      # If input is vector (i.e. X-linked male genotypes), return it unchanged
      if (is.matrix(matr)) 1000 * matr[1, ] + matr[2, ] else matr
    }

    # For each orig, find the indices of its haplos (in orig$hap) that also occur in its copy.
    # Then take cross product of these vectors.
    loopgrid = fast.grid(lapply(seq_along(origs), function(i) {
      ori = two2one(dat[[c(origs[i], 1)]])
      seq_along(ori)[ori %in% two2one(dat[[c(copies[i], 1)]])]
    }), as.list = TRUE)

    # Initialise likelihood
    likelihood = 0

    for (r in loopgrid) {
      # r is an index vector: r[i] gives a column number of the hap matrix of orig[i].
      dat1 = dat
      attr(dat1, "impossible") = FALSE

      for (i in seq_along(origs)) {
        orig.int = origs[i]
        copy.int = copies[i]
        orighap = dat[[orig.int]]$hap
        origprob = dat[[orig.int]]$prob
        hap = if (is.matrix(orighap))  orighap[, r[i], drop = F] else orighap[r[i]]
        prob = origprob[r[i]]
        if (sum(prob) == 0)
          print("Loop-loekke: Alle sannsynligheter er null. Magnus lurer paa om dette kan gi feilmelding.")
        dat1[[orig.int]] = list(hap = hap, prob = prob)
        dat1[[copy.int]] = list(hap = hap, prob = 1)
      }

      for (sub in inform_subnucs) {
        dat1 = PEEL(dat1, sub)

        # If impossible data - break out of ES-algorithm and go to next r in loopgrid.
        if (sub$link > 0 && attr(dat1, "impossible")) break

        # If pedigree traversed, add to totalt and go to next r
        if (sub$link == 0) likelihood = likelihood + dat1
      }
    }
  }
  if (is.numeric(logbase)) log(likelihood, logbase)
  else likelihood
}

#' @export
#' @rdname likelihood
likelihood.singleton = function(x, marker1, logbase = NULL, ...) {
  if (!is.marker(marker1)) {
    if(is.count(marker1) && marker1 < nMarkers(x))
      marker1 = x$markerdata[[marker1]]
    else
      stop2("Argument `marker1` must be either a `marker` object or the index of an attached marker: ", marker1)
  }
  if (is.null(marker1) || all(marker1 == 0))
    return(if (is.numeric(logbase)) 0 else 1)

  m = marker1
  afr = afreq(m)

  onX = is_Xmarker(m)
  if (onX && x$SEX == 1) {
    if (all(m > 0) && m[1] != m[2])
      stop2("Heterozygous genotype at X-linked marker in male singleton")
    res = afr[m[1]]
  }

  else if (0 %in% m) {
    p = afr[m[m != 0]]
    res = p^2 + 2 * p * (1 - p)
  }
  else {
    res = prod(afr[m]) * ifelse(m[1] != m[2], 2, 1)  # assumes HWE
  }
  if (is.numeric(logbase)) log(res, logbase) else res
}

#' @export
#' @rdname likelihood
likelihood.list = function(x, marker1, marker2 = NULL, logbase = NULL, total = TRUE, ...) {
  assert_that(is.pedList(x), is.count(marker1), is.null(marker2) || is.count(marker2))

  liks = vapply(x, function(xx) likelihood(xx, marker1, marker2, logbase=logbase, ...), numeric(1))

  if (total)
    if(!is.null(logbase)) sum(liks) else prod(liks)
  else liks
}



