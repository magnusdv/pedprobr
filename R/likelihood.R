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
#'
#' * 1 marker or 2 linked markers
#'
#' * complex inbred pedigrees
#'
#' * markers with mutation models
#'
#' * pedigrees with inbred founders
#'
#' For more than two linked markers, see [likelihoodMerlin()].
#'
#' @param x A `ped` object, a `singleton` object, or a list of such objects.
#' @param markers One or several markers compatible with `x`. Several input
#'   forms are possible:
#'
#'   * A [marker()] object compatible with `x`.
#'
#'   * A list of marker objects.
#'
#'   * A vector of names or indices of markers attached to `x`. If `x` is a
#'   list, this is the only valid input.
#'
#' @param marker1,marker2 Single markers compatible with `x`.
#' @param rho The recombination rate between `marker1` and `marker2`. To make
#'   biological sense `rho` should be between 0 and 0.5.
#' @param lump Activate allele lumping, i.e., merging unobserved alleles. This
#'   is an important time saver, and should be applied in nearly all cases. (The
#'   parameter exists mainly for debugging purposes.) The lumping algorithm will
#'   detect (and complain) if any markers use a non-lumpable mutation model.
#'   Default: TRUE.
#' @param eliminate Mostly for internal use: a non-negative integer indicating
#'   the number of iterations in the internal genotype-compatibility algorithm.
#'   Positive values can save time if the number of alleles is large.
#' @param logbase Either NULL (default) or a positive number indicating the
#'   basis for logarithmic output. Typical values are `exp(1)` and 10.
#' @param loopBreakers A vector of ID labels indicating loop breakers. If NULL
#'   (default), automatic selection of loop breakers will be performed. See
#'   [breakLoops()].
#' @param peelOrder For internal use.
#' @param dropout A single number (applied to all individuals) or a named
#'   vector. The probability of allelic dropout, implemented as in Method 3 of
#'   Dørum et al (2014).
#' @param verbose A logical.
#' @param theta Theta correction.
#' @param \dots Further arguments.

#' @return A numeric with the same length as the number of markers indicated by
#'   `markers`. If `logbase` is a positive number, the output is
#'   `log(likelihood, logbase)`.
#'
#' @seealso [likelihoodMerlin()], for likelihoods involving more than 2 linked
#'   markers.
#'
#' @author Magnus Dehli Vigeland
#'
#' @references
#'
#' * For a general overview of the implementation:
#' Vigeland (2021). _Pedigree analysis in R_. ISBN: 9780128244302.
#'
#' * Elston-Stewart algorithm:
#' Elston and Stewart (1971). _A General Model for the Genetic Analysis of
#' Pedigree Data_. \doi{https://doi.org/10.1159/000152448}
#'
#' * Dropout model:
#' Dørum et al (2014). _Models and implementation for relationship problems with
#' dropout_. \doi{https://doi.org/10.1007/s00414-014-1046-5}
#'
#' @examples
#'
#' ### Simple likelihood ###
#' p = 0.1
#' q = 1 - p
#'
#' # Singleton
#' s = singleton() |> addMarker(geno = "1/2", afreq = c("1" = p, "2" = q))
#'
#' stopifnot(all.equal(likelihood(s), 2*p*q))
#'
#' # Trio
#' t = nuclearPed() |>
#'   addMarker(geno = c("1/1", "1/2", "1/1"), afreq = c("1" = p, "2" = q))
#'
#' stopifnot(all.equal(likelihood(t), p^2 * 2*p*q * 0.5))
#'
#'
#' ### Example of calculation with inbred founders ###
#'
#' ### Case 1: Trio with inbred father
#' x = cousinPed(0, child = TRUE)
#' x = addSon(x, 5)
#' x = relabel(x, old = 5:7, new = c("father", "mother", "child"))
#'
#' # Add equifrequent SNP; father homozygous, child heterozygous
#' x = addMarker(x, father = "1/1", child = "1/2")
#'
#' # Plot with genotypes
#' plot(x, marker = 1)
#'
#' # Compute the likelihood
#' lik1 = likelihood(x, markers = 1)
#'
#'
#' ### Case 2: Using founder inbreeding
#' # Remove ancestry of father
#' y = subset(x, c("father", "mother", "child"))
#'
#' # Indicate that the father has inbreeding coefficient 1/4
#' founderInbreeding(y, "father") = 1/4
#'
#' # Plot (notice the inbreeding coefficient)
#' plot(y, marker = 1)
#'
#' # Likelihood should be the same as above
#' lik2 = likelihood(y, markers = 1)
#'
#' stopifnot(all.equal(lik1, lik2))
#'
#'
#'
#' @export
likelihood = function(x, ...) UseMethod("likelihood", x)

#' @export
#' @rdname likelihood
likelihood.ped = function(x, markers = NULL, peelOrder = NULL, lump = TRUE,
                          eliminate = 0, logbase = NULL, loopBreakers = NULL,
                          verbose = FALSE, theta = 0, dropout = NULL, ...) {

  if(theta > 0 && hasInbredFounders(x))
    stop2("Theta correction cannot be used in pedigrees with inbred founders")

  if(hasSelfing(x))
    stop2("Likelihood of pedigrees with selfing is not implemented.\n",
          "Contact the maintainer if this is important to you.")

  # Catch erroneous input
  if(is.ped(peelOrder))
    stop2("Invalid input for argument `peelOrder`. Received object type: ", class(peelOrder))

  if(is.null(markers))
    markers = x$MARKERS
  else if(is.vector(markers) && !is.list(markers))
    markers = getMarkers(x, markers = markers)
  else if(is.marker(markers))
    markers = list(markers)
  else if(!is.markerList(markers))
    stop2("Invalid input for argument `markers`. Received object type: ", class(markers))

  if(verbose)
    message("Number of markers: ", length(markers))

  if(length(markers) == 0)
    return(numeric(0))

  ### Quick calculations if singleton
  if(is.singleton(x)) {
    if(verbose)
      message("Passing to singleton method")
    dropout = if(x$ID %in% names(dropout)) dropout[x$ID] else dropout %||% 0
    liks = vapply(markers, function(m)
      likelihoodSingleton(x, m, theta = theta, dropout = dropout), FUN.VALUE = 1)
    return(if(is.numeric(logbase)) log(liks, logbase) else liks)
  }

  # Allele lumping
  if(lump)
    markers = lapply(markers, function(m) reduceAlleles(m, verbose = verbose))

  # Break unbroken loops TODO: move up (avoid re-attaching)
  if (x$UNBROKEN_LOOPS) {
    if(verbose)
      message("Tip: To optimize speed, consider breaking loops before calling 'likelihood'. See ?breakLoops.")
    x = breakLoops(setMarkers(x, markers), loopBreakers = loopBreakers, verbose = verbose)
    markers = x$MARKERS

  }

  # Dropout
  if(is.null(dropout)) {
    dropout = rep_len(0, length(x$ID))
  } else if(!is.null(dnms <- names(dropout))) {
    dropout = ifelse(x$ID %in% dnms, dropout[x$ID], 0)
  } else
    dropout = rep_len(as.numeric(dropout), length(x$ID))

  # Peeling order: Same for all markers
  if(is.null(peelOrder))
    peelOrder = peelingOrder(x)
  if(theta == 0)
    peelOrder = informativeSubnucs(x, mlist = markers, peelOrder = peelOrder)

  if(verbose)
    message(sprintf("%d informative %s", length(peelOrder), if(length(peelOrder) == 1) "nucleus" else "nuclei"))

  treatAsFou = attr(peelOrder, "treatAsFounder")

  # Autosomal or X?
  isX = vapply(markers, isXmarker, logical(1))
  Xchrom = all(isX)
  if(!Xchrom && any(isX))
    stop2("Cannot mix autosomal and X-linked markers in the same likelihood calculation")
  if(verbose)
    message("Chromosome type: ", if(Xchrom) "X" else "autosomal")

  # Select tools for peeling
  # TODO: Organise better, e.g., skip startdata if theta > 0
  if(Xchrom) {
    starter = function(x, m) startdata_M_X(x, m, eliminate = eliminate, treatAsFounder = treatAsFou)
    peeler = function(x, m) function(dat, sub) .peel_M_X(dat, sub, SEX = x$SEX, mutmat = mutmod(m))
  }
  else {
    starter = function(x, m) startdata_M_AUT(x, m, eliminate = eliminate, treatAsFounder = treatAsFou, dropout = dropout)
    peeler = function(x, m) function(dat, sub) .peel_M_AUT(dat, sub, mutmat = mutmod(m))
  }

  # Loop over markers
  resList = lapply(markers, function(m) {

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
peelingProcess = function(x, m, startdata, peeler, peelOrder = NULL) {

  if(hasUnbrokenLoops(x))
    stop2("Peeling process cannot handle unbroken pedigree loops")

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

      # If pedigree traversed, add to total and go to next r
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
  twolocus = length(dat1) > 4
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

likelihoodSingleton = function(x, m, theta = 0, dropout = 0) {
  m1 = m[1]
  m2 = m[2]

  # Quick return if empty
  if(m1 == 0 && m2 == 0)
    return(1)

  afr = afreq(m)
  chromX = isXmarker(m)

  # Theta correction or founder inbreeding
  if(theta > 0)
    f = theta
  else
    f = founderInbreeding(x, chromType = if(chromX) "x" else "autosomal")

  # Male on X
  if (chromX && x$SEX == 1) {
    if (all(m > 0) && m1 != m2)
      stop2("Heterozygous genotype at X-linked marker in male singleton")
    return(afr[m1])
  }

  # One missing allele
  if (m1 == 0 || m2 == 0) {
    p = afr[m[m != 0]]
    res = p^2 + 2 * p * (1 - p)
    return(res)
  }

  # Otherwise: The usual HW formula, possibly corrected for inbreeding and/or dropout
  HWprob(m1, m2, afr, f = f, dropout = dropout)
}
