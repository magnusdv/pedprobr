###### OTHER AUXILIARY FUNCTIONS

informativeSubnucs = function(x, marker1, marker2 = NULL) {
  # Trim pedigree by removing leaves without genotypes.
  # Also remove nucs that are completely uninformative

  nucs = attr(x, "PEELING_ORDER")
  if(is.null(nucs))
    nucs = attr(x, "PEELING_ORDER") = peelingOrder(x)

  stationary = hasStationaryModel(marker1) &&
    (is.null(marker2) || hasStationaryModel(marker2))

  isMiss = marker1[, 1] == 0 & marker1[, 2] == 0
  if(!is.null(marker2))
    isMiss = isMiss & marker2[, 1] == 0 & marker2[, 2] == 0

  # Quick return if x is nuclear
  if(length(nucs) == 1) {
    if (stationary) {
      sub = nucs[[1]]
      sub$children = sub$children[!isMiss[sub$children]]
      goodNucs = list(sub)
    }
    else {
      goodNucs = nucs
    }
    return(list(subnucs = goodNucs, newfounders = numeric(0)))
  }

  # Return unchanged if all are genotyped
  if (!any(isMiss))
    return(list(subnucs = nucs, newfounders = numeric(0)))


  newfounders = numeric(0)
  goodNucs = list()
  NONFOU = nonfounders(x, internal = T)
  LEAVES = leaves(x, internal = T)

  isMiss[x$LOOP_BREAKERS] = F  # works (and quick) also if no loops.
  isUninfLeaf = isUninfFou = isMiss

  isUninfLeaf[-LEAVES] = F   # logical with T only if uninformative leaf
  isUninfFou[NONFOU] = F # logical with T only if uninformative founder

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

      newfounders = c(newfounders, link)
      next
    }
    sub$children = keepOffs

    goodNucs = c(goodNucs, list(sub))
    isUninfFou[link] = FALSE  #added in v0.8-1 to correct a bug marking certain 'middle' subnucs uninformative
  }
  list(subnucs = goodNucs, newfounders = newfounders)
}

