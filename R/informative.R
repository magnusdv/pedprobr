###### OTHER AUXILIARY FUNCTIONS

informativeSubnucs = function(x, mlist, peelOrder = peelingOrder(x)) {
  # Trim pedigree by removing leaves without genotypes.
  # Also remove nucs that are completely uninformative

  if(is.marker(mlist)) mlist = list(mlist)
  stationary = all(unlist(lapply(mlist, hasStationaryModel)))

  # Find untyped individuals
  nInd = pedsize(x)
  markermat = unlist(mlist)
  dim(markermat) = c(nInd, length(markermat)/nInd)
  untyped = rowSums(markermat) == 0
  typed = !untyped

  nucs = peelOrder

  # Quick return if x is nuclear
  if(length(nucs) == 1) {
    if (stationary) {
      sub = nucs[[1]]
      sub$children = sub$children[typed[sub$children]]
      goodNucs = list(sub)
    }
    else {
      goodNucs = nucs
    }
    return(list(subnucs = goodNucs, treatAsFounder = numeric(0)))
  }

  # Return unchanged if all are genotyped
  if (all(typed))
    return(list(subnucs = nucs, treatAsFounder = numeric(0)))

  treatAsFounder = numeric(0)
  goodNucs = list()
  NONFOU = nonfounders(x, internal = TRUE)
  LEAVES = leaves(x, internal = TRUE)

  untyped[x$LOOP_BREAKERS] = FALSE  # works (and quick) also if no loops.
  isUninfLeaf = isUninfFou = untyped

  isUninfLeaf[-LEAVES] = FALSE  # logical with T only if uninformative leaf
  isUninfFou[NONFOU] = FALSE  # logical with T only if uninformative founder

  for (sub in nucs) {
    fa = sub$father
    mo = sub$mother
    offs = sub$children
    link = sub$link
    isUninfLeaf[link] = FALSE # just in case link is an untyped child. TODO: necessary??

    # Keep only genotyped leaves
    keepOffs = offs[!isUninfLeaf[offs]]
    nkeep = length(keepOffs)

    # Skip nuc if: all kids are uninf leaves AND the non-link parent is uninf
    if (nkeep == 0 && link == fa && isUninfFou[mo])
      next
    if (nkeep == 0 && link == mo && isUninfFou[fa])
      next

    # Skip nuc if: Single child is link, both parents uninf founders
    # NB: Only skip if stationary mutation model!
    if (stationary && nkeep == 1 && link == keepOffs &&
        isUninfFou[fa] && isUninfFou[mo]) {
      treatAsFounder = c(treatAsFounder, link)
      next
    }
    sub$children = keepOffs

    goodNucs = c(goodNucs, list(sub))
    isUninfFou[link] = FALSE  #added in v0.8-1 to correct a bug marking certain 'middle' subnucs uninformative
  }
  list(subnucs = goodNucs, treatAsFounder = treatAsFounder)
}

