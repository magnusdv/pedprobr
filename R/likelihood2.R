#' @rdname likelihood
#' @export
likelihood2 = function(x, ...) UseMethod("likelihood2", x)

#' @export
#' @rdname likelihood
likelihood2.ped = function(x, marker1, marker2, rho = NULL, peelOrder = NULL,
                          eliminate = 0, logbase = NULL, loopBreakers = NULL,
                          verbose = FALSE, ...) {

  if(hasInbredFounders(x))
    stop2("Likelihood of linked markers is not implemented in pedigrees with founder inbreeding.\n",
          "(Note that this is usually not well-defined)")

  if(hasSelfing(x))
    stop2("Likelihood of pedigrees with selfing is not implemented.\n",
          "Contact the maintainer if this is important to you.")

  if(is.vector(marker1) && !is.list(marker1) && length(marker1) == 1)
    marker1 = getMarkers(x, markers = marker1)[[1]]
  else if(!is.marker(marker1))
      stop2("Argument `marker1` must be a single marker. Received: ", class(marker1))

  if(is.atomic(marker2) && !is.list(marker2) && length(marker2) == 1)
    marker2 = getMarkers(x, markers = marker2)[[1]]
  else if(!is.marker(marker2))
    stop2("Argument `marker2` must be a single marker. Received: ", class(marker2))

  checkRho(rho)

  ### Quick return if singleton (linkage is then irrelevant)
  if(is.singleton(x)) {
    lik1 = likelihoodSingleton(x, marker1)
    lik2 = likelihoodSingleton(x, marker2)
    if(is.numeric(logbase))
      res = log(lik1, logbase) + log(lik2, logbase)
    else
      res = lik1 * lik2
    return(res)
  }

  # Allele lumping
  marker1 = reduceAlleles(marker1, verbose = verbose)
  marker2 = reduceAlleles(marker2, verbose = verbose)

  # Break unbroken loops TODO: move up (avoid re-attaching)
  if (x$UNBROKEN_LOOPS) {
    if(verbose)
      message("Tip: To optimize speed, consider breaking loops before calling 'likelihood'. See ?breakLoops.")
    x = breakLoops(setMarkers(x, list(marker1, marker2)), loopBreakers = loopBreakers, verbose = verbose)
    marker1 = x$MARKERS[[1]]
    marker2 = x$MARKERS[[2]]
  }

  ### Quick return if unlinked
  if(rho == 0.5) {
    if(verbose)
      message("Unlinked markers; computing likelihoods separately")
    lik1 = likelihood.ped(x, marker1, logbase = logbase, verbose = FALSE)
    lik2 = likelihood.ped(x, marker2, logbase = logbase, verbose = FALSE)
    if(is.numeric(logbase))
      res = log(lik1, logbase) + log(lik2, logbase)
    else
      res = lik1 * lik2
    return(res)
  }

  # Peeling order
  if(is.null(peelOrder))
    peelOrder = informativeSubnucs(x, mlist = list(marker1, marker2), peelOrder = peelingOrder(x))

  treatAsFou = attr(peelOrder, "treatAsFounder")

  # Autosomal or X?
  x1 = isXmarker(marker1)
  x2 = isXmarker(marker2)
  if(x1 != x2)
    stop2("Both markers must be either autosomal or X-linked")
  Xchrom = x1 && x2
  if(verbose)
    message("Chromosome type: ", if(Xchrom) "X" else "autosomal")

  # Peeler function
  if(!Xchrom)
    peeler = function(dat, sub) .peel_MM_AUT(dat, sub, rho, mut1 = attr(marker1, "mutmod"), mut2 = attr(marker2, "mutmod"))
  else
    peeler = function(dat, sub) .peel_MM_X(dat, sub, rho, x$SEX, mut1 = attr(marker1, "mutmod"), mut2 = attr(marker2, "mutmod"))

  # Startdata
  if(Xchrom)
    startdata = startdata_MM_X(x, marker1, marker2, eliminate, treatAsFou)
  else
    startdata = startdata_MM_AUT(x, marker1, marker2, eliminate, treatAsFou)

  res = peelingProcess(x, m = NULL, startdata = startdata, peeler = peeler, peelOrder = peelOrder)

  if(is.numeric(logbase)) log(res, logbase) else res
}


#' @export
#' @rdname likelihood
likelihood2.list = function(x, marker1, marker2, logbase = NULL, ...) {
  if(!is.pedList(x))
    stop2("Input is a list, but not a list of `ped` objects")

  if (!(is.vector(marker1) && !is.list(marker1)))
    stop2("`likelihood.list()` requires `marker1` to be a vector of marker names or indices. Received: ", class(marker1))
  if (!(is.vector(marker2) && !is.list(marker2)))
    stop2("`likelihood.list()` requires `marker2` to be a vector of marker names or indices. Received: ", class(marker2))

  likel = vapply(x, function(comp)
    likelihood2.ped(comp, marker1, marker2, logbase = logbase, ...),
    FUN.VALUE = 1)

  if(is.numeric(logbase)) sum(likel) else prod(likel)
}
