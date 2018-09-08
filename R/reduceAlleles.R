
#' @importFrom utils modifyList
#' @importFrom pedmut isLumpable lumpedMatrix
.reduce_alleles = function(marker) {
  if (all(marker != 0))
    return(marker)  # no reduction needed (OK!)

  attrs = attributes(marker)
  orig_alleles = attrs$alleles

  # Index of observed alleles
  present = unique.default(as.numeric(marker))
  present = sort.int(present[present > 0])

  if (length(present) >= length(orig_alleles) - 1)
    return(marker)  # return unchanged if all, or all but one, are observed

  redund = .mysetdiff(seq_along(orig_alleles), present) # TODO: avoid setdiff?

  allowsMut = allowsMutations(marker)
  mut = attrs$mutmod$male # TODO for now assuming equal. Include both sexes in mutmod!

  if(allowsMut && !isLumpable(mut, redund))
    return(marker)

  new_marker = match(marker, present, nomatch = 0)
  new_alleles = c(orig_alleles[present], "lump")
  present_freq = attrs$afreq[present]
  new_freq = c(present_freq, 1 - sum(present_freq))
  attributes(new_marker) = modifyList(attrs, list(alleles = new_alleles, afreq = new_freq))

  if (allowsMut) {
    newMut = lumpedMatrix(mut, redund)
    print(newMut)
    mutmod(new_marker) = list(male = newMut, female = newMut)
  }

  new_marker
}

