#' Allele lumping
#'
#' Perform allele lumping (i.e., merging unobserved alleles) for all markers
#' attached to the input pedigree.
#'
#' @param x A `ped` object or a list of such.
#' @param markers A vector of names or indices referring to markers attached to
#'   `x`. (Default: All markers.)
#' @param always A logical. If TRUE, lumping is always attempted. By default
#'   (FALSE) lumping is skipped for markers where no individuals, or all
#'   individuals, are genotyped.
#' @param special A logical. If TRUE, special lumping procedures (depending on
#'   the pedigree) will be attempted if the marker has a mutation model that is
#'   not generally lumpable (in the Kemeny-Snell sense).
#' @param alleleLimit A positive number or `Inf` (default). If the mutation
#'   model is not generally lumpable, and the allele count exceeds this limit,
#'   switch to an `equal` model with the same rate and reapply lumping.
#' @param verbose A logical.
#'
#' @return An object similar to `x`, but whose attached markers have reduced
#'   allele set.
#'
#' @examples
#' x = nuclearPed() |> addMarker(geno = c(NA, NA, "1/1"), alleles = 1:5)
#'
#' # Before lumping
#' afreq(x, 1)
#'
#' # Lump
#' y = lumpAlleles(x, verbose = TRUE)
#' afreq(y, 1)
#'
#' # With lumpable mutation model
#' x2 = setMutmod(x, model = "equal", rate = 0.1)
#' mutmod(x2, 1)
#'
#' y2 = lumpAlleles(x2, verbose = TRUE)
#' mutmod(y2, 1)
#'
#' # Mutation model requiring special lumping
#' x3 = setMutmod(x, model = "random", rate = 0.1, seed = 1)
#' mutmod(x3, 1)
#'
#' # Lump
#' y3 = lumpAlleles(x3, verbose = TRUE)
#'
#' mutmod(y3, 1)
#'
#' stopifnot(all.equal(likelihood(x), likelihood(y)),
#'           all.equal(likelihood(x2), likelihood(y2)),
#'           all.equal(likelihood(x3), likelihood(y3)))
#'
#' @importFrom pedmut alwaysLumpable getParams isLumpable lumpedModel
#'   lumpMutSpecial
#' @export
lumpAlleles = function(x, markers = NULL, always = FALSE, special = TRUE,
                       alleleLimit = Inf, verbose = FALSE) {

  if(is.pedList(x))
    return(lapply(x, function(comp) lumpAlleles(comp, markers, verbose = verbose)))

  # For special lumping, pass along pedigree for info
  ped = setMarkers(x, NULL)

  # Index of markers to be reduced
  midx = whichMarkers(x, markers %||% seq_along(x$MARKERS))

  # Loop through
  for(i in midx) {
    m = x$MARKERS[[i]]

    if(sum(m) == 0 && !always) {
      if(verbose) message(sprintf("Marker %s: empty -- no lumping needed", name(m)))
      next
    }

    newm = .lumpMarker(m, always, special, ped, verbose = verbose)

    # If not regular lumpable, switch to fallback model if too many alleles
    if(nAlleles(newm) > alleleLimit && !alwaysLumpable(mutmod(m))) {
      if(verbose)
        message("  --- allele limit exceeded; switching to `equal` model")
      y = setMutmod(x, i, model = "equal", update = TRUE)
      newm = .lumpMarker(y$MARKERS[[i]], always, verbose = verbose)
    }

    x$MARKERS[[i]] = newm
  }

  x
}

.lumpMarker = function(marker, always = FALSE, special = TRUE, ped = NULL, verbose = FALSE) {

  attrs = attributes(marker)
  origAlleles = attrs$alleles
  nall = length(origAlleles)
  mut = attrs$mutmod

  # Internal allele indices as vector (speeds up things below)
  mInt = as.integer(marker)

  # Observed alleles
  isObs = seq_len(nall) %in% mInt
  nObs = sum(isObs)

  if(verbose) {
    modname = if(is.null(mut)) "no" else modelAbbr[getParams(mut, params = "model", format = 3)[[1]]]
    message(sprintf("Marker %s: %d of %d alleles, %s model. ",
                name(marker), nObs, nall, modname), appendLF = FALSE)
  }

  if (!always && all(mInt != 0)) {
    if(verbose) message("Lumping not needed (all members genotyped)")
    return(marker)
  }

  # No lumping if all, or all but one, are observed
  if(nObs >= nall - 1) {
    if(verbose) message("Lumping not needed")
    return(marker)
  }

  lump = if(!nObs) origAlleles else origAlleles[!isObs]
  mut = attrs$mutmod

  # No mutation model - or lumpable model
  if(is.null(mut) || isLumpable(mut, lump)) {
    marker[] = match(mInt, which(isObs), nomatch = 0L)
    keepAls = origAlleles[isObs]
    keepFreqs = attrs$afreq[isObs]
    attr(marker, "alleles") = c(keepAls, "lump")
    attr(marker, "afreq") = c(keepFreqs, 1 - sum(keepFreqs))

    if(!is.null(mut))
      attr(marker, "mutmod") = pedmut::lumpedModel(mut, lump = lump, check = FALSE)

    if(verbose)
      message(sprintf("Regular lumping: %d -> %d", nall, nObs + 1))
    return(marker)
  }

  # Attempt special lumping?
  if(!special || is.null(ped)) {
    if(verbose) message("Special lumping disabled")
    return(marker)
  }

  ### Special lumping

  # Signature of untyped individuals (F-depth, N-depth, F-width, N-width)
  usign = uSignature(ped, marker = marker)

  # Check if signature is implemented (better to do it in lumpMutSpecial?)
  if(usign[1] > 1 || usign[3] + usign[4] > 0) {
    if(verbose) message("Special lumping not available in this case")
    return(marker)
  }

  # Main calculation done in pedmut
  lumpedMut = lumpMutSpecial(mut, lump = lump, uSign = usign, verbose = FALSE)

  # Check if lumping was successful
  newAls = colnames(lumpedMut$female) # attr(lumpedMut, "afreq") may be NULL here
  if(length(newAls) == nall) {
    if(verbose) message("Special lumping ineffective")
    return(marker)
  }

  # Update marker
  marker[] = match(mInt, which(isObs), nomatch = 0L)

  attr(marker, "alleles") = newAls
  attr(marker, "afreq") = attr(lumpedMut, "afreq")
  attr(marker, "mutmod") = lumpedMut

  if(verbose)
    message(sprintf("Special lumping: %d -> %d", nall, length(newAls)))

  marker
}

# Signature of untyped individuals in pedigree x:
# (F-depth, F-width, N-depth, N-width)
uSignature = function(x, untyped = NULL, marker = NULL) {

  if(is.pedList(x)) {
    usList = lapply(x, function(comp) uSignature(comp, untyped, marker))
    return(do.call(pmax, usList))
  }

  # If marker is NULL, check all markers; return max scores
  if(is.null(untyped) && is.null(marker)) {

    if(!length(x$MARKERS))
      return(c(Fdep = 0, Ndep = 0, Fwid = 0, Nwid = 0))

    untyp = lapply(x$MARKERS, function(m) x$ID[m[,1] + m[,2] == 0])
    usigs = lapply(unique.default(untyp), function(u) uSignature(x, untyped = u))
    return(do.call(pmax, usigs))
  }

  if(!is.null(marker) && !is.marker(marker))
    marker = getMarkers(x, marker)[[1]]

  u = untyped %||% x$ID[marker[,1] + marker[,2] == 0]

  x = trim(x, u, verbose = FALSE)
  u = .myintersect(u, x$ID)

  if(is.null(x) || !length(u))
    return(c(Fdep = 0, Ndep = 0, Fwid = 0, Nwid = 0))

  isUntyped = x$ID %in% u
  names(isUntyped) = x$ID

  uF = .myintersect(u, founders(x))
  uN = .mysetdiff(u, uF)

  # F-depth: Longest untyped chain starting with founder
  d1 = descentPaths(x, uF) |> unlist(recursive = FALSE, use.names = FALSE)

  fchains = lapply(d1, function(p)
    p[seq_len(match(FALSE, isUntyped[p], nomatch = length(p) + 1) - 1)])

  Fdep = if(length(fchains)) max(lengths(fchains)) else 0

  # F-width: Max number of children of untyped founder
  Fwid = if(length(uF)) max(nChildren(x, uF)) else 0

  # N-depth: Longest untyped chain starting with nonfounder
  urem = .mysetdiff(uN, unlist(fchains, use.names = FALSE))

  d2 = descentPaths(x, urem) |>
    unlist(recursive = FALSE, use.names = FALSE)

  nchains = lapply(d2, function(p)
    p[seq_len(match(FALSE, isUntyped[p], nomatch = length(p) + 1) - 1)])

  Ndep = if(length(nchains)) max(lengths(nchains)) else 0

  # N-width: Max number of children of untyped nonfounder
  Nwid = if(length(uN)) max(nChildren(x, uN)) else 0

  # Return
  c(Fdep = Fdep, Fwid = Fwid, Ndep = Ndep, Nwid = Nwid)
}

modelAbbr = c(custom = "custom", dawid = "dawid", equal ="equal",
              proportional = "prop", random = "random", onestep = "onestep",
              stepwise = "stepw", trivial = "trivial")
