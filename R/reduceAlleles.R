#' Allele lumping
#'
#' Perform allele lumping (i.e., merging unobserved alleles) for all markers
#' attached to the input pedigree.
#'
#' @param x A `ped` object or a list of such.
#' @param markers A vector of names or indices referring to markers attached to
#'   `x`. (Default: All markers.)
#' @param always A logical. If TRUE, lumping is always attempted. By default
#'   (FALSE) lumping is skipped for markers where all individuals are genotyped.
#' @param verbose A logical.
#'
#' @return An object similar to `x`, but whose attached markers have reduced
#'   allele set.
#'
#' @examples
#' x = nuclearPed() |> addMarker(geno = c("1/1", NA, NA), alleles = 1:4)
#'
#' # Before lumping
#' afreq(x, 1)
#'
#' # Lump
#' y = lumpAlleles(x, verbose = TRUE)
#' afreq(y, 1)
#'
#' @export
lumpAlleles = function(x, markers = NULL, always = FALSE, verbose = FALSE) {
  markers = markers %||% seq_len(nMarkers(x))

  if(is.pedList(x))
    return(lapply(x, function(comp) lumpAlleles(comp, markers, verbose = verbose)))

  if(!is.ped(x))
    stop2("Input must be a `ped` object or a list of such")

  mlist = getMarkers(x, markers)
  mlistLumped = lapply(seq_along(mlist), function(i) {
    m = mlist[[i]]
    label = name(m)
    if(is.na(label)) label = i
    if(verbose) message("Marker ",label, ". ", appendLF = FALSE)
    reduceAlleles(m, always = always, verbose = verbose)
  })

  setMarkers(x, mlistLumped)
}

#' @importFrom pedmut isLumpable lumpedModel
reduceAlleles = function(marker, always = FALSE, verbose = FALSE) {

  if (is.null(marker)) {
    if(verbose) message("Lumping not needed - NULL marker")
    return(NULL)
  }

  # Internal allele indices as vector (speeds up things below)
  mInt = as.integer(marker)

  if (!always && all(mInt != 0)) {
    if(verbose) message("Lumping not needed - all members genotyped")
    return(marker)
  }

  attrs = attributes(marker)
  origAlleles = attrs$alleles
  nall = length(origAlleles)

  # Observed alleles
  isObs = seq_len(nall) %in% mInt
  nObs = sum(isObs)

  # No lumping if all, or all but one, are observed
  if(nObs >= nall - 1) {
    if(verbose)
      message(sprintf("Lumping not needed: %d of %d alleles observed", nObs, nall))
    return(marker)
  }

  lump = if(!nObs) origAlleles else origAlleles[!isObs]

  # No lumping if mutation model is present and not lumpable for this lump
  mut = attrs$mutmod
  if(!is.null(mut) && !isLumpable(mut, lump)) {
    if(verbose) message("Mutation model is not lumpable")
    return(marker)
  }

  marker[] = match(mInt, which(isObs), nomatch = 0L)
  attr(marker, "alleles") = c(origAlleles[isObs], "lump")

  presentFreq = attrs$afreq[isObs]
  attr(marker, "afreq") = c(presentFreq, 1 - sum(presentFreq))

  if (!is.null(mut))
    mutmod(marker) = lumpedModel(mut, lump = lump, check = FALSE)

  if(verbose)
    message(sprintf("Lumping: %d -> %d alleles", nall, nObs + 1))

  marker
}

