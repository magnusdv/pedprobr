#' @rdname likelihood
#' @export
likelihood2 = function(x, ...) UseMethod("likelihood2", x)

#' @export
#' @rdname likelihood
likelihood2.ped = function(x, marker1, marker2, rho, peelOrder = NULL,
                          eliminate = 0, logbase = NULL, loop_breakers = NULL,
                          verbose = FALSE, theta = NULL, ...) {

  if(!is.null(theta)) {
    message("Argument `theta` has been renamed to `rho`")
    rho = theta
  }

  if(is.null(rho))
    stop2("Argument `rho` is missing")

  if(hasInbredFounders(x))
    stop2("Likelihood of linked markers is not implemented in pedigrees with founder inbreeding.\n",
          "(Note that this is usually not well-defined)")

  if(hasSelfing(x))
    stop2("Likelihood of pedigrees with selfing is not implemented.\n",
          "Contact the maintainer if this is important to you.")

  if(!is.marker(marker1)) {
    if(is.atomic(marker1) && length(marker1) == 1)
      marker1 = getMarkers(x, markers = marker1)[[1]]
    else
      stop2("Argument `marker1` must be a single marker")
  }

  if(!is.marker(marker2)) {
    if(is.atomic(marker2) && length(marker2) == 1)
      marker2 = getMarkers(x, markers = marker2)[[1]]
    else
      stop2("Argument `marker2` must be a single marker")
  }

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
    x = breakLoops(setMarkers(x, list(marker1, marker2)), loop_breakers = loop_breakers, verbose = verbose)
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
    peeler = function(dat, sub) .peel_MM_AUT(dat, sub, rho)
  else
    peeler = function(dat, sub) .peel_MM_X(dat, sub, rho, x$SEX)

  # Startdata
  if(Xchrom)
    startdata = startdata_MM_X(x, marker1, marker2, eliminate, treatAsFou)
  else
    startdata = startdata_MM_AUT(x, marker1, marker2, eliminate, treatAsFou)

  res = peelingProcess(x, startdata = startdata, peeler = peeler, peelOrder = peelOrder)

  if(is.numeric(logbase)) log(res, logbase) else res
}


#' @export
#' @rdname likelihood
likelihood2.list = function(x, marker1, marker2, logbase = NULL, ...) {
  if(!is.pedList(x))
    stop2("Input is a list, but not a list of `ped` objects")

  if (!is.atomic(marker1))
    stop2("`likelihood.list()`requires `marker1` to be a vector or marker names or indices. Received: ", class(marker1))
  if (!is.atomic(marker2))
    stop2("`likelihood.list()`requires `marker2` to be a vector or marker names or indices. Received: ", class(marker2))

  likel = vapply(x, function(comp)
    likelihood2.ped(comp, marker1, marker2, logbase = logbase, ...),
    FUN.VALUE = 1)

  if(is.numeric(logbase)) sum(likel) else prod(likel)
}


# Internal function: likelihood of a single marker
peelingProcess2 = function(x, marker1, marker2, startdata, peeler, peelOrder = NULL) {

  if(is.null(peelOrder))
    peelOrder = informativeSubnucs(x, list(marker1, marker2))

  if(is.function(startdata))
    startdata = startdata(x, marker1, marker2, peelOrder)

  if(attr(startdata, "impossible"))
    return(0)

  dat = startdata

  if (is.null(x$LOOP_BREAKERS)) { # If no loops
    for (nuc in peelOrder) {
      dat = peeler(dat, nuc)
      if (nuc$link > 0 && attr(dat, "impossible"))
        return(0)
    }

    # Output of last peeling is the likelihood
    return(dat)
  }

  ### If broken loops
  LB = x$LOOP_BREAKERS

  # For each orig, find the indices of its genotypes that also occur in its copy.
  genoMatching = lapply(1:nrow(LB), function(i)
    matchDat(dat[[LB[[i,"orig"]]]], dat[[LB[[i,"copy"]]]]))

  # Then take cross product of these vectors.
  loopgrid = fastGrid(genoMatching, as.list = TRUE)

  # Initialise likelihood
  likelihood = 0

  for (r in loopgrid) {
    dat1 = dat
    attr(dat1, "impossible") = FALSE

    for (i in 1:nrow(LB)) { # Note: r[i] is a valid index of orig[i]$pat/mat
      origi = LB[[i, "orig"]]
      copyi = LB[[i, "copy"]]
      origDat = dat[[origi]]
      dat1[[origi]] = dat1[[copyi]] = lapply(origDat, function(vec) vec[r[i]])
      dat1[[copyi]]$prob = 1
      if (sum(dat1[[origi]]$prob) == 0)
        message("The likelihood algorithm reached a strange place. The maintainer would be grateful to see this example.")
    }

    for (nuc in peelOrder) {
      dat1 = peeler(dat1, nuc)

      # If impossible data - break out of ES-algorithm and go to next r in loopgrid.
      if (nuc$link > 0 && attr(dat1, "impossible")) break

      # If pedigree traversed, add to total and go to next r
      if (nuc$link == 0) likelihood = likelihood + dat1
    }
  }

  likelihood
}

