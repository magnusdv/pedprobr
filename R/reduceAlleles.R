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

#' @importFrom pedmut isLumpable lumpedMatrix mutationModel
reduceAlleles = function(marker, always = FALSE, verbose = FALSE) {

  if (is.null(marker)) {
    if(verbose) message("Lumping not needed - NULL marker")
    return(NULL)
  }

  if (!always && all(marker != 0)) {
    if(verbose) message("Lumping not needed - all members genotyped")
    return(marker)
  }

  attrs = attributes(marker)
  origAlleles = attrs$alleles

  # Index of observed alleles
  presentIdx = unique.default(marker[marker > 0])

  # No lumping if all, or all but one, are observed
  if(length(presentIdx) >= length(origAlleles) - 1) {
    if(verbose) message(sprintf("Lumping not needed: %d of %d alleles observed",
                                length(presentIdx), length(origAlleles)))
    return(marker)
  }

  presentIdx = sort.int(presentIdx, method = "shell")
  lump = if(!length(presentIdx)) origAlleles else origAlleles[-presentIdx]

  # No lumping if mutation model is present and not lumpable for this lump
  mut = attrs$mutmod
  if(!is.null(mut) && !isLumpable(mut, lump)) {
    if(verbose) message("Mutation model is not lumpable")
    return(marker)
  }

  marker[] = match(marker, presentIdx, nomatch = 0)
  attr(marker, "alleles") = c(origAlleles[presentIdx], "lump")

  presentFreq = attrs$afreq[presentIdx]
  attr(marker, "afreq") = c(presentFreq, 1 - sum(presentFreq))

  if (!is.null(mut)) {
    lumpedMale = lumpedMatrix(mut$male, lump)
    lumpedFemale = lumpedMatrix(mut$female, lump)
    mutmod(marker) = mutationModel(list(female = lumpedFemale, male = lumpedMale))
  }

  if(verbose) message(sprintf("Lumping: %d -> %d alleles",
              length(origAlleles), length(presentIdx) + 1))
  marker
}

