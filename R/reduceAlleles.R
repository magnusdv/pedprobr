
#' @importFrom utils modifyList
.reduce_alleles = function(marker) {
  if (all(marker != 0))
    return(marker)  # no reduction needed (OK!)
  attrs = attributes(marker)

  allowsMut = allowsMutations(marker) 
  
  if (allowsMut) {
    malem = attrs$mutmat$male
    femalem = attrs$mutmat$female
    male_lump = identical(attr(malem, "lumpability"), "always")
    female_lump = identical(attr(femalem, "lumpability"), "always")
    if (!male_lump || !female_lump)
      return(marker)
  }
  orig_alleles = attrs$alleles

  # Index of observed alleles
  present = unique.default(as.numeric(marker))
  present = sort.int(present[present > 0])
  if (length(present) >= length(orig_alleles) - 1)
    return(marker)  # return unchanged if all, or all but one, are observed

  redund = .mysetdiff(seq_along(orig_alleles), present) # TODO: avoid setdiff?
  dummylab = paste(orig_alleles[redund], collapse = "_")

  if (length(present) == 0) {
    new_marker = rep.int(0, length(marker))
    attributes(new_marker) = modifyList(attrs, list(alleles = dummylab, afreq = 1))
    if (!is.null(attrs$mutmat)) {
      mm = matrix(1, dimnames = list(dummylab, dummylab))
      attr(new_marker, "mutmat") = list(male = mm, female = mm)
    }
    return(new_marker)
  }

  new_marker = match(marker, present, nomatch = 0)
  new_alleles = c(orig_alleles[present], dummylab)
  present_freq = attrs$afreq[present]
  new_freq = c(present_freq, 1 - sum(present_freq))
  n = length(present) + 1

  attributes(new_marker) = modifyList(attrs, list(alleles = new_alleles, afreq = new_freq))

  if (!is.null(attrs$mutmat)) { # TODO allowsMutation()
    if (male_lump) {
      mm = malem[c(present, redund[1]), c(present, redund[1])]
      mm[, n] = 1 - rowSums(mm[, -n, drop = F])
    }
    if (female_lump) {
      mf = femalem[c(present, redund[1]), c(present, redund[1])]
      mf[, n] = 1 - rowSums(mf[, -n, drop = F])
    }
    # for(i in 1:(n-1)) { present_i = present[i] #m_weight = f_weight = attrs$afreq[redund]
    # m_weight = malem[present_i, redund] f_weight = femalem[present_i, redund] mm[n, i] =
    # (m_weight/sum(m_weight)) %*% malem[redund, present_i] mf[n, i] = (f_weight/sum(f_weight))
    # %*% femalem[redund, present_i] }
    attr(new_marker, "mutmat") = list(male = mm, female = mf)
  }
  new_marker
}

