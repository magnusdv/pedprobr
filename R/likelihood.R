#' Pedigree likelihood
#'
#' The `likelihood()` and `likelihood2()` functions constitute the heart of
#' **pedprobr**. The former computes the pedigree likelihood for each indicated
#' marker. The latter computes the likelihood for a pair of linked markers
#' separated by a given recombination rate.
#'
#' The implementation is based on the peeling algorithm of Elston and Stewart
#' (1971). A variety of situations are covered; see the Examples section for
#' some demonstrations.
#'
#' * autosomal and X-linked markers
#' * complex inbred pedigrees
#' * markers with mutation models
#' * pedigrees with inbred founders
#' * single markers or two linked markers
#'
#' For more than two linked markers, see [likelihoodMerlin()].
#'
#' Allele lumping can significantly reduce computation time with highly
#' polymorphic STR markers and many untyped pedigree members. This is
#' particularly important in `likelihood2()` which is prone to run out of memory
#' without lumping. If a non-lumpable mutation model is used, specialised
#' lumping may still be possible in some situations. This is attempted if
#' `special = TRUE`, which is the default in `likelihood2()` but not in
#' `likelihood()`.
#'
#' @param x A `ped` object, a `singleton` object, or a list of such objects.
#' @param markers One or several markers compatible with `x`. Several input
#'   forms are possible:
#'
#'   * A `marker` object compatible with `x`.
#'
#'   * A list of `marker` objects.
#'
#'   * A vector of names or indices of markers attached to `x`. If `x` is a
#'   list, this is the only valid input.
#'
#' @param marker1,marker2 Single markers compatible with `x`.
#' @param rho The recombination rate between `marker1` and `marker2`. To make
#'   biological sense `rho` should be between 0 and 0.5.
#' @param lump Activate allele lumping, i.e., merging unobserved alleles. This
#'   is an important time saver, and should be applied in nearly all cases. (The
#'   parameter exists mainly for debugging purposes.) If any markers use a
#'   non-lumpable mutation model, the `special` argument may be used to apply
#'   more advanced methods.
#' @param special A logical indicating if special lumping procedures should be
#'   attempted if the mutation model is not generally lumpable. By default FALSE
#'   in `likelihood()` and TRUE in `likelihood2()`.
#' @param alleleLimit A positive number or `Inf` (default). If the mutation
#'   model is not generally lumpable, and the allele count exceeds this limit,
#'   switch to an `equal` model with the same rate and reapply lumping.
#' @param logbase Either NULL (default) or a positive number indicating the
#'   basis for logarithmic output. Typical values are `exp(1)` and 10.
#' @param loopBreakers A vector of ID labels indicating loop breakers. If NULL
#'   (default), automatic selection of loop breakers will be performed. See
#'   [pedtools::breakLoops()].
#' @param peelOrder For internal use.
#' @param allX For internal use; set to TRUE if all markers are X-chromosomal.
#' @param verbose A logical.
#' @param theta Theta correction.
#' @param \dots Further arguments.

#' @return A numeric with the same length as the number of markers indicated by
#'   `markers`. If `logbase` is a positive number, the output is
#'   `log(likelihood, logbase)`.
#'
#' @seealso [likelihoodMerlin()], for likelihoods involving more than 2 linked markers.
#'
#' @author Magnus Dehli Vigeland
#' @references Elston and Stewart (1971). _A General Model for the Genetic
#'   Analysis of Pedigree Data_. \doi{https://doi.org/10.1159/000152448}
#'
#' @examples
#'
#' ### Simple likelihoods ###
#' p = 0.1
#' q = 1 - p
#' afr = c("1" = p, "2" = q)
#'
#' # Singleton
#' s = singleton() |> addMarker(geno = "1/2", afreq = afr)
#' stopifnot(all.equal(likelihood(s), 2*p*q))
#'
#' # Trio
#' x = nuclearPed() |> addMarker(geno = c("1/1", "1/2", "1/1"), afreq = afr)
#' lik = likelihood(x, verbose = TRUE)
#' stopifnot(all.equal(lik, p^2 * 2*p*q * 0.5))
#'
#'
#' ### Example with inbred founder ###
#'
#' # Set 100% inbreeding for the father in the previous example
#' y = setFounderInbreeding(x, ids = 1, value = 1)
#'
#' # Plot (notice the inbreeding coefficient)
#' plot(y, marker = 1)
#'
#' stopifnot(all.equal(likelihood(y), p * 2*p*q * 0.5))
#'
#'
#' ### Example with two linked markers
#'
#' # Add a second marker, highly polymorphic
#' x = addMarker(x, geno = c(NA, NA, "1/1"), alleles = 1:10)
#'
#' # Likelihood assuming complete linkage
#' likelihood2(x, 1, 2, rho = 0, verbose = TRUE)
#'
#' @export
likelihood = function(x, ...) UseMethod("likelihood", x)

#' @export
#' @rdname likelihood
likelihood.ped = function(x, markers = NULL, peelOrder = NULL, lump = TRUE,
                          special = FALSE, alleleLimit = Inf,
                          logbase = NULL, loopBreakers = NULL, allX = NULL,
                          verbose = FALSE, theta = 0, ...) {

  if(theta > 0 && hasInbredFounders(x))
    stop2("Theta correction cannot be used in pedigrees with inbred founders")

  if(hasSelfing(x))
    stop2("Likelihood of pedigrees with selfing is not implemented.\n",
          "Contact the maintainer if this is important to you.")

  # Catch erroneous input
  if(is.ped(peelOrder))
    stop2("Invalid input for argument `peelOrder`. Received object type: ", class(peelOrder))

  if(!is.null(markers)) {
    if(is.vector(markers) && !is.list(markers))
      x = selectMarkers(x, markers = markers)
    else if(is.marker(markers) || is.markerList(markers))
      x = setMarkers(x, markers)
    else
      stop2("Invalid input for argument `markers`. Received object type: ", class(markers))
  }

  # Just pedigree
  ped = setMarkers(x, NULL)

  nM = length(x$MARKERS)
  if(verbose)
    message("Number of markers: ", nM)

  if(nM == 0)
    return(numeric(0))

  # Autosomal or X?
  if(length(allX) == 1 && is.logical(allX))
    Xchrom = allX
  else {
    isX = vapply(x$MARKERS, isXmarker, logical(1))
    if(!all(isX == isX[1]))
      stop2("Cannot mix autosomal and X-linked markers in the same likelihood calculation")
    Xchrom = isX[1]
  }

  if(verbose)
    message("Chromosome type: ", if(Xchrom) "X" else "autosomal")

  ### Quick calculations if singleton
  if(is.singleton(x)) {
    if(verbose) message("Passing to singleton method")
    liks = vapply(x$MARKERS, function(m)
      likelihoodSingleton(ped, m, theta = theta, Xchrom = Xchrom), FUN.VALUE = 1)
    return(if(is.numeric(logbase)) log(liks, logbase) else liks)
  }

  # Allele lumping
  if(lump)
    x = lumpAlleles(x, always = FALSE, special = special,
                    alleleLimit = alleleLimit, verbose = verbose)

  # Break unbroken loops
  if (x$UNBROKEN_LOOPS)
    x = breakLoops(x, loopBreakers = loopBreakers, verbose = verbose)

  # Peeling order: Same for all markers
  if(is.null(peelOrder))
    peelOrder = peelingOrder(x)
  if(theta == 0)
    peelOrder = informativeSubnucs(x, peelOrder = peelOrder)

  if(verbose)
    message("Informative nuclei: ", length(peelOrder))

  treatAsFou = attr(peelOrder, "treatAsFounder")

  # Select tools for peeling
  # TODO: Organise better, e.g., skip startdata if theta > 0

  pedInfo = .pedInfo(x, treatAsFounder = treatAsFou, Xchrom = Xchrom)
  starter = function(x, m) startdata_M(x, m, pedInfo = pedInfo)
  if(Xchrom)
    peeler = function(x, m) function(dat, sub) .peel_M_X(dat, sub, SEX = x$SEX, mutmat = mutmod(m))
  else
    peeler = function(x, m) function(dat, sub) .peel_M_AUT(dat, sub, mutmat = mutmod(m))

  # Loop over markers
  resList = lapply(x$MARKERS, function(m) {
    # If no data, return 1
    if (sum(m) == 0) return(1)

    # If theta correction, go to different function
    if(theta > 0)
      likTheta(x, m, theta, peeler(x,m), peelOrder)
    else
      peelingProcess(x, m, starter(x,m), peeler(x,m), peelOrder)
  })

  res = unlist(resList, recursive = FALSE)

  if(is.numeric(logbase)) log(res, logbase) else res
}



# Internal function: likelihood of a single marker
peelingProcess = function(x, m = x$MARKERS[[1]], startdata = NULL, peeler = NULL, peelOrder = NULL) {

  if(hasUnbrokenLoops(x))
    stop2("Peeling process cannot handle unbroken pedigree loops")

    # Default start
  if(is.null(startdata))
    startdata = function(x, m) startdata_M(x, m)

  # Default peeler
  if(is.null(peeler))
    if(isXmarker(m))
      peeler = function(dat, sub) .peel_M_X(dat, sub, SEX = x$SEX, mutmat = mutmod(m))
    else
      peeler = function(dat, sub) .peel_M_AUT(dat, sub, mutmat = mutmod(m))

  if(is.null(peelOrder))
    peelOrder = informativeSubnucs(x, m)

  if(is.function(startdata))
    startdata = startdata(x, m)

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
  nr = nrow(LB)

  # For each orig, find the indices of its genotypes that also occur in its copy.
  genoMatching = lapply(seq_len(nr), function(i)
    matchDat(dat[[LB[[i,"orig"]]]], dat[[LB[[i,"copy"]]]]))

  # Then take cross product of these vectors.
  loopgrid = fastGrid(genoMatching, as.list = TRUE)

  # Initialise likelihood
  likelihood = 0

  for (r in loopgrid) {
    dat1 = dat
    attr(dat1, "impossible") = FALSE

    for (i in seq_len(nr)) { # Note: r[i] is a valid index of orig[i]$pat/mat
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

      # If pedigree traversed, dat1 is a number. Add to total and goto next
      if (nuc$link == 0) likelihood = likelihood + dat1
    }
  }

  likelihood
}



#' @export
#' @rdname likelihood
likelihood.list = function(x, markers = NULL, logbase = NULL, ...) {
  if(!is.pedList(x))
    stop2("Input is a list, but not a list of `ped` objects")

  # Theta correction not implemented for lists
  if("theta" %in% names(list(...)))
    stop2("Theta correction is not implemented for lists")

  if(is.null(markers))
    markers = seq_len(nMarkers(x))
  else if (!(is.vector(markers) && !is.list(markers)))
    stop2("`likelihood.list()` requires `markers` to be a vector of marker names or indices. Received: ", class(markers))

  if(length(markers) == 0)
    return(numeric(0))

  liks = vapply(x, function(comp)
    likelihood.ped(comp, markers, logbase = logbase, ...),
    FUN.VALUE = numeric(length(markers)))

  if(length(markers) == 1)
    dim(liks) = c(1, length(x))

  if(is.numeric(logbase)) rowSums(liks) else apply(liks,1,prod)
}



# Utility for finding which genotypes in dat1 are also in dat2
matchDat = function(dat1, dat2) {
  twolocus = !is.null(dat1$allele2) || !is.null(dat1$mat2)
  if(!twolocus) {
    nseq = seq_along(dat1$mat)
    Xchrom = is.null(dat1$pat)
    if(Xchrom)
      nseq[dat1$mat %in% dat2$mat]
    else
      nseq[(dat1$pat + 1000*dat1$mat) %in% (dat2$pat + 1000*dat2$mat)]
  }
  else {
    nseq = seq_along(dat1$mat1)
    Xchrom = is.null(dat1$pat1)
    if(Xchrom)
      nseq[(dat1$mat1 %in% dat2$mat1) & (dat1$mat2 %in% dat2$mat2)]
    else {
      loc1 = (dat1$pat1 + 1000*dat1$mat1) %in% (dat2$pat1 + 1000*dat2$mat1)
      loc2 = (dat1$pat2 + 1000*dat1$mat2) %in% (dat2$pat2 + 1000*dat2$mat2)
      nseq[loc1 & loc2]
    }
  }
}

likelihoodSingleton = function(x, m, theta = 0, Xchrom = isXmarker(m)) {
  m1 = m[1]
  m2 = m[2]

  # Quick return if empty
  if(m1 + m2 == 0)
    return(1)

  afr = attr(m, "afreq")

  # Theta correction or founder inbreeding
  if(theta > 0)
    f = theta
  else
    f = founderInbreeding(x, chromType = if(Xchrom) "x" else "autosomal")

  # Male on X
  if (Xchrom && x$SEX == 1) {
    if (m1 != m2)
      stop2("Heterozygous genotype at X-linked marker in male singleton")
    return(afr[m1])
  }

  # One missing allele
  if (m1 == 0 || m2 == 0) {
    p = afr[m[m != 0]]
    res = p^2 + 2 * p * (1 - p)
    return(res)
  }

  # Otherwise: The usual HW formula, possibly with inbreeding correction
  HWprob(m1, m2, afr, f = f)
}
