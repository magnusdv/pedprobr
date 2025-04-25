###### OTHER AUXILIARY FUNCTIONS

informativeSubnucs = function(x, mlist = x$MARKERS, peelOrder = peelingOrder(x)) {
  # Trim pedigree by removing leaves without genotypes.
  # Also remove nucs that are completely uninformative

  if(is.marker(mlist)) {
    stationary = hasStationaryModel(mlist)
    untyped = mlist[,1] + mlist[,2] == 0
  }
  else {
    stationary = all(unlist(lapply(mlist, hasStationaryModel),
                            recursive = FALSE, use.names = FALSE))

    # Find untyped individuals
    nInd = pedsize(x)
    markermat = unlist(mlist, recursive = FALSE, use.names = FALSE)
    dim(markermat) = c(nInd, length(markermat)/nInd)
    untyped = rowSums(markermat) == 0
  }

  typed = !untyped

  # Return unchanged if all are genotyped
  if (all(typed))
    return(peelOrder)

  # Quick return if x is nuclear
  if(length(peelOrder) == 1) {
    # TODO! Remove this?
    #if (!stationary)
    #  return(peelOrder)

    ch = peelOrder[[1]]$children
    peelOrder[[1]]$children = ch[typed[ch]]
    return(peelOrder)
  }


  treatAsFounder = numeric(0)
  goodNucs = list()
  NONFOU = nonfounders(x, internal = TRUE)
  LEAVES = leaves(x, internal = TRUE)

  untyped[x$LOOP_BREAKERS] = FALSE  # works (and quick) also if no loops.
  isUninfLeaf = isUninfFou = untyped

  isUninfLeaf[-LEAVES] = FALSE  # logical with T only if uninformative leaf
  isUninfFou[NONFOU] = FALSE  # logical with T only if uninformative founder

  for (nuc in peelOrder) {
    fa = nuc$father
    mo = nuc$mother
    offs = nuc$children
    link = nuc$link
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
    nuc$children = keepOffs

    goodNucs = c(goodNucs, list(nuc))
    isUninfFou[link] = FALSE  #added in v0.8-1 to correct a bug marking certain 'middle' nucs uninformative
  }

  attr(goodNucs, "treatAsFounder") = treatAsFounder

  goodNucs
}

