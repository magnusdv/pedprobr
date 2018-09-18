###### OTHER AUXILIARY FUNCTIONS

informativeSubnucs = function(x, marker1, marker2 = NULL) {
  # Trim pedigree by removing leaves without genotypes.
  # Also remove nucs that are completely uninformative

  nucs = attr(x, "PEELING_ORDER")
  if(is.null(nucs))
    nucs = attr(x, "PEELING_ORDER") = peelingOrder(x)

  if (!is.null(marker2))
    marker1 = marker1 + marker2
  if (all(marker1[, 1] > 0))
    return(list(subnucs = nucs, newfounders = numeric(0)))

  stationary = hasStationaryModel(marker1) &&
               (is.null(marker2) || hasStationaryModel(marker2))

  newfounders = numeric(0)
  good_nucs = list()
  NONFOU = nonfounders(x, internal=T)
  LEAVES = leaves(x, internal=T)

  is_miss = marker1[, 1] == 0 & marker1[, 2] == 0
  is_miss[x$LOOP_BREAKERS] = F  # works (and quick) also if no loops.
  is_uninf_leaf = is_uninf_fou = is_miss

  is_uninf_leaf[-LEAVES] = F   # logical with T only if uninformative leaf
  is_uninf_fou[NONFOU] = F # logical with T only if uninformative founder

  for (sub in nucs) {
    fa = sub$father
    mo = sub$mother
    offs = sub$children
    link = sub$link
    is_uninf_leaf[link] = FALSE # just in case link is an untyped child. TODO: necessary??

    # Keep only genotyped leaves
    keep_offs = offs[!is_uninf_leaf[offs]]
    nkeep = length(keep_offs)

    # Skip nuc if: all kids are uninf leaves AND the non-link parent is uninf
    if (nkeep == 0 && link == fa && is_uninf_fou[mo])
      next
    if (nkeep == 0 && link == mo && is_uninf_fou[fa])
      next

    # Skip nuc if: Single child is link, both parents uninf founders
    # NB: Only skip if stationary mutation model!
    if (stationary && nkeep == 1 && link == keep_offs &&
        is_uninf_fou[fa] && is_uninf_fou[mo]) {

      newfounders = c(newfounders, link)
      next
    }
    sub$children = keep_offs

    good_nucs = c(good_nucs, list(sub))
    is_uninf_fou[link] = FALSE  #added in v0.8-1 to correct a bug marking certain 'middle' subnucs uninformative
  }
  list(subnucs = good_nucs, newfounders = newfounders)
}

