
#' @importFrom utils modifyList
#' @importFrom pedmut isLumpable lumpedMatrix mutationModel
reduceAlleles = function(marker) {
  if (all(marker != 0))
    return(marker)  # no reduction needed (OK!)

  attrs = attributes(marker)
  origAlleles = attrs$alleles

  # Index of observed alleles
  presentIdx = sort.int(unique.default(marker[marker > 0]))

  # No lumping if all, or all but one, are observed
  if (length(presentIdx) >= length(origAlleles) - 1)
    return(marker)

  lump = if(!length(presentIdx)) origAlleles else origAlleles[-presentIdx]

  # No lumping if mutation model is present and not lumpable for this lump
  allowsMut = allowsMutations(marker)
  if(allowsMut && !isLumpable(mutmod(marker), lump))
    return(marker)

  newMarker = match(marker, presentIdx, nomatch = 0)
  newAlleles = c(origAlleles[presentIdx], "lump")
  presentFreq = attrs$afreq[presentIdx]
  newFreq = c(presentFreq, 1 - sum(presentFreq))
  attributes(newMarker) = modifyList(attrs, list(alleles = newAlleles, afreq = newFreq))

  if (allowsMut) {
    mut = mutmod(marker)
    lumpedMale = lumpedMatrix(mut$male, lump)
    lumpedFemale = lumpedMatrix(mut$female, lump)
    mutmod(newMarker) = mutationModel(list(female = lumpedFemale, male = lumpedMale))
  }
  else {
    mutmod(newMarker) = NULL
  }

  newMarker
}

