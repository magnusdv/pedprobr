
#' @importFrom utils modifyList
#' @importFrom pedmut isLumpable lumpedMatrix mutationModel
reduceAlleles = function(marker, verbose = FALSE) {

  if (is.null(marker)) {
    if(verbose) message("Allele lumping: Not needed - NULL marker")
    return(NULL)
  }

  if (all(marker != 0)) {
    if(verbose) message("Allele lumping: Not needed - all members genotyped")
    return(marker)
  }

  attrs = attributes(marker)
  origAlleles = attrs$alleles

  # Index of observed alleles
  presentIdx = unique.default(marker[marker > 0])

  # No lumping if all, or all but one, are observed
  if (length(presentIdx) >= length(origAlleles) - 1) {
    if(verbose) message("Allele lumping: Not needed - all (or all but one) alleles present")
    return(marker)
  }

  presentIdx = sort.int(presentIdx, method = "shell")
  lump = if(!length(presentIdx)) origAlleles else origAlleles[-presentIdx]

  # No lumping if mutation model is present and not lumpable for this lump
  mut = attrs$mutmod
  if(!is.null(mut) && !isLumpable(mut, lump)) {
    if(verbose) message("Allele lumping: mutation model is not lumpable")
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

  if(verbose) message(sprintf("Allele lumping: %d alleles -> %d alleles",
              length(origAlleles), length(presentIdx) + 1))
  marker
}

