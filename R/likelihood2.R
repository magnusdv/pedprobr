#' @rdname likelihood
#' @export
likelihood2 = function(x, ...) UseMethod("likelihood2", x)

#' @export
#' @rdname likelihood
likelihood2.ped = function(x, marker1, marker2, rho = NULL, peelOrder = NULL,
                           lump = TRUE, special = TRUE, alleleLimit = Inf,
                           logbase = NULL, loopBreakers = NULL, verbose = FALSE, ...) {

  if(hasInbredFounders(x))
    stop2("Likelihood of linked markers is not implemented in pedigrees with founder inbreeding.\n",
          "(Note that this is usually not well-defined)")

  if(hasSelfing(x))
    stop2("Likelihood of pedigrees with selfing is not implemented.\n",
          "Contact the maintainer if this is important to you.")

  # Attach markers if not already
  if(is.marker(marker1) && is.marker(marker2))
    x = setMarkers(x, list(marker1, marker2))
  else if(length(marker1) == 1 && length(marker2) == 1)
    x = selectMarkers(x, c(marker1, marker2))
  else
    stop2("Unrecognised markers; received objects of class: ", class(marker1), ", ", class(marker2))

  # Check rho input
  checkRho(rho)

  # If both markers empty, return 1
  if(sum(x$MARKERS[[1]]) + sum(x$MARKERS[[2]]) == 0) {
    if(verbose) message("Both markers are empty; returning likelihood 1")
    return(1)
  }

  # Quick return if singleton (linkage is then irrelevant)
  if(is.singleton(x)) {
    lik1 = likelihoodSingleton(x, x$MARKERS[[1]])
    lik2 = likelihoodSingleton(x, x$MARKERS[[2]])
    res = if(is.numeric(logbase)) log(lik1, logbase) + log(lik2, logbase) else lik1 * lik2
    return(res)
  }

  # Allele lumping
  if(lump)
    x = lumpAlleles(x, always = FALSE, special = special, alleleLimit = alleleLimit,
                    verbose = verbose)

  # Break unbroken loops
  if(x$UNBROKEN_LOOPS)
    x = breakLoops(x, loopBreakers = loopBreakers, verbose = verbose)

  # Peeling order
  if(is.null(peelOrder))
    peelOrder = informativeSubnucs(x, peelOrder = peelingOrder(x))

  # Quick return if unlinked
  if(rho == 0.5) {
    if(verbose)
      message("Unlinked markers; computing likelihoods separately")
    lik1 = likelihood.ped(x, 1, peelOrder = peelOrder, lump = FALSE, logbase = logbase, verbose = FALSE)
    lik2 = likelihood.ped(x, 2, peelOrder = peelOrder, lump = FALSE, logbase = logbase, verbose = FALSE)
    res = if(is.numeric(logbase)) log(lik1, logbase) + log(lik2, logbase) else lik1 * lik2
    return(res)
  }

  treatAsFou = attr(peelOrder, "treatAsFounder")

  m1 = x$MARKERS[[1]]
  m2 = x$MARKERS[[2]]
  mut1 = mutmod(m1)
  mut2 = mutmod(m2)

  # Autosomal or X?
  x1 = isXmarker(m1)
  x2 = isXmarker(m2)
  if(x1 != x2)
    stop2("Both markers must be either autosomal or X-linked")
  Xchrom = x1 && x2
  if(verbose)
    message("Chromosome type: ", if(Xchrom) "X" else "autosomal")

  # Peeler function
  if(!Xchrom)
    peeler = function(dat, sub) .peel_MM_AUT(dat, sub, rho, mut1 = mut1, mut2 = mut2)
  else
    peeler = function(dat, sub) .peel_MM_X(dat, sub, rho, x$SEX, mut1 = mut1, mut2 = mut2)

  # Startdata
  pedInfo = .pedInfo(x, treatAsFounder = treatAsFou, Xchrom = Xchrom)
  startdata = startdata_MM(x, m1, m2, pedInfo = pedInfo)

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
