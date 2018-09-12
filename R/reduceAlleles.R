
#' @importFrom utils modifyList
#' @importFrom pedmut isLumpable lumpedMatrix mutationModel
reduceAlleles = function(marker) {
  if (all(marker != 0))
    return(marker)  # no reduction needed (OK!)

  attrs = attributes(marker)
  orig_alleles = attrs$alleles

  # Index of observed alleles
  present_idx = unique.default(marker[marker > 0])

  # No lumping if all, or all but one, are observed
  if (length(present_idx) >= length(orig_alleles) - 1)
    return(marker)

  lump = if(!length(present_idx)) orig_alleles else orig_alleles[-present_idx]

  # No lumping if mutation model is present and not lumpable for this lump
  allowsMut = allowsMutations(marker)
  if(allowsMut && !isLumpable(mutmod(marker), lump))
    return(marker)

  new_marker = match(marker, present_idx, nomatch = 0)
  new_alleles = c(orig_alleles[present_idx], "lump")
  present_freq = attrs$afreq[present_idx]
  new_freq = c(present_freq, 1 - sum(present_freq))
  attributes(new_marker) = modifyList(attrs, list(alleles = new_alleles, afreq = new_freq))

  if (allowsMut) {
    mut = mutmod(marker)
    lumpedMale = lumpedMatrix(mut$male, lump)
    lumpedFemale = lumpedMatrix(mut$female, lump)
    mutmod(new_marker) = mutationModel(list(female = lumpedFemale, male = lumpedMale))
  }

  new_marker
}

