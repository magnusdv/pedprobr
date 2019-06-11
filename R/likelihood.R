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
#' @param setup for internal use.
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
#' # likelihood of inbred pedigree (grandfather/granddaughter incest)
#' x = nuclearPed(father = "grandfather", child = "a")
#' x = addDaughter(x, "a", id = "granddaughter")
#' x = addChildren(x, father = "grandfather", mother = "granddaughter",
#'                 nch = 1, id = "child")
#' m = marker(x, grandfather = 1, child = 1:2)
#'
#' plot(x, m)
#' lik = likelihood(x, m)
#'
#' stopifnot(lik == 0.09375)
#'
#' @export
likelihood = function(x, ...) UseMethod("likelihood", x)


#' @export
#' @rdname likelihood
likelihood.ped = function(x, marker1, marker2 = NULL, theta = NULL, setup = NULL,
                          eliminate = 0, logbase = NULL, loop_breakers = NULL,
                          verbose = FALSE, ...) {

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

  marker1 = reduceAlleles(marker1)
  marker2 = reduceAlleles(marker2)  # unchanged if NULL

  if (x$UNBROKEN_LOOPS) {
    if(verbose) message("Tip: To optimize speed, consider breaking loops before calling 'likelihood'. See ?breakLoops.")
    m = list(marker1)
    m[[2]] = marker2  # no effect if NULL
    x = breakLoops(setMarkers(x, m), loop_breakers = loop_breakers, verbose = verbose)
    marker1 = x$markerdata[[1]]
    if (twolocus)
      marker2 = x$markerdata[[2]]
  }

  setup = setupData(x, marker1, marker2, eliminate, setup)

  dat = setup$startdata
  if (attr(dat, "impossible"))
    return(ifelse(is.numeric(logbase), -Inf, 0))

  Xchrom = is_Xmarker(marker1)
  mut = mutmod(marker1)
  PEEL = choosePeeler(twolocus, theta, Xchrom, x$SEX, mut)

  if (is.null(dups <- x$LOOP_BREAKERS)) {
    for (sub in setup$informativeNucs) {
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
    loopgrid = fastGrid(lapply(seq_along(origs), function(i) {
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

      for (sub in setup$informativeNucs) {
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
likelihood.singleton = function(x, marker1, marker2 = NULL, logbase = NULL, ...) {
  twolocus = !is.null(marker2)
  if(missing(marker1) || is.null(marker1))
    stop2("Argument `marker1` is missing")

  # If two markers: Linkage is irrelevant for singletons
  if (!is.null(marker2)) {
    lik1 = likelihood.singleton(x, marker1)
    lik2 = likelihood.singleton(x, marker2)
    total = lik1 * lik2
    if (is.numeric(logbase)) total = log(total, logbase)
    return(total)
  }

  if (!is.marker(marker1)) {
    if(length(marker1) != 1)
      stop2("Length of `marker1` must be 1")
    marker1 = getMarkers(x, markers = marker1)[[1]]
  }

  # Quick return if empty
  if(all(marker1 == 0))
    return(if (is.numeric(logbase)) 0 else 1)

  m = marker1
  afr = afreq(m)
  chromX = is_Xmarker(m)
  finb = founderInbreeding(x, chromType = if(chromX) "x" else "autosomal")

  if (chromX && x$SEX == 1) {
    if (all(m > 0) && m[1] != m[2])
      stop2("Heterozygous genotype at X-linked marker in male singleton")
    res = afr[m[1]]
  }
  else if (m[1] == 0 || m[2] == 0) {
    p = afr[m[m != 0]]
    res = p^2 + 2 * p * (1 - p)
  }
  else {
    res = HWprob(m[1], m[2], afr, finb)
  }
  if (is.numeric(logbase)) log(res, logbase) else res
}

#' @export
#' @rdname likelihood
likelihood.list = function(x, marker1, marker2 = NULL, logbase = NULL, total = TRUE, ...) {
  if(!is.pedList(x))
    stop2("Input is a list, but not a list of `ped` objects")

  nx = length(x)
  if (is.atomic(marker1))
    marker1 = rep(list(marker1), length = nx)
  if (is.atomic(marker2))
    marker2 = rep(list(marker2), length = nx)  # Note: NULL is atomic

  liks = vapply(1:nx,
                function(i) likelihood(x[[i]], marker1[[i]], marker2[[i]], logbase=logbase, ...),
                numeric(1))

  if (total)
    if(!is.null(logbase)) sum(liks) else prod(liks)
  else
    liks
}


setupData = function(x, marker1, marker2, eliminate, setup) {
  if (is.null(setup))
    setup = list()

  informativeNucs = setup$informativeNucs
  treatAsFounder  = setup$treatAsFounder
  startdata       = setup$startdata

  if(is.null(informativeNucs)) {
    inform = informativeSubnucs(x, marker1, marker2)
    informativeNucs = inform$subnucs
    treatAsFounder = inform$newfounders
  }

  if(is.null(startdata)) {
    if(is.null(marker2))
      startdata = startdata_M(x, marker = marker1, eliminate = eliminate,
                              treatAsFounder = treatAsFounder)
    else
      startdata = startdata_MM(x, marker1 = marker1, marker2 = marker2,
                               eliminate = eliminate, treatAsFounder = treatAsFounder)
  }

  list(informativeNucs = informativeNucs, treatAsFounder = treatAsFounder, startdata = startdata)
}
