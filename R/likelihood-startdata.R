#### FUNCTIONS FOR CREATING THE INTITIAL HAPLOTYPE COMBINATIONS W/PROBABILITIES.

startData = function(x, marker1, marker2, Xchrom, eliminate = 0, treatAsFounder = NULL) {

  if(is.null(marker2)) {
    if(Xchrom)
      startdata_M_X(x, marker1, eliminate, treatAsFounder)
    else
      startdata_M_AUT(x, marker1, eliminate, treatAsFounder)
  }
  else {
    if(Xchrom)
      startdata_MM_X(x, marker1, marker2, eliminate, treatAsFounder)
    else
      startdata_MM_AUT(x, marker1, marker2, eliminate, treatAsFounder)
  }
}

# TODO: Remove?
startdata_M = function(x, marker, eliminate = 0, treatAsFounder = NULL) {
  if(isXmarker(marker))
    startdata_M_X(x, marker, eliminate, treatAsFounder)
  else
    startdata_M_AUT(x, marker, eliminate, treatAsFounder)
}

startdata_M_AUT = function(x, marker, eliminate = 0, treatAsFounder = NULL) {

  glist = .buildGenolist(x, marker, eliminate, treatAsFounder)

  if (attr(glist, "impossible"))
    return(structure(list(), impossible = TRUE))

  FOU = founders(x, internal = TRUE)

  # Founder inbreeding: A vector of length pedsize(x), with NA's at nonfounders
  # Enables quick look-up e.g. FOU_INB[i].
  FOU_INB = rep(NA_real_, pedsize(x))
  FOU_INB[FOU] = founderInbreeding(x)

  # Add any members which should be treated as founders
  FOU = c(FOU, treatAsFounder)

  afr = afreq(marker)
  impossible = FALSE

  # Add probabilities to each genotype
  dat = lapply(1:pedsize(x), function(i) {

    # If impossible, speed through
    if(impossible)
      return(NULL)

    g = glist[[i]]
    if (i %in% FOU) {
      prob = HWprob(g$pat, g$mat, afr, f = FOU_INB[i])
      if (sum(prob) == 0){
        impossible = TRUE
        return(NULL)
      }
    }
    else
      prob = rep.int(1, length(g$mat))

    g$prob = as.numeric(prob)
    g
  })

  attr(dat, "impossible") = impossible
  dat
}

startdata_M_X = function(x, marker, eliminate = 0, treatAsFounder = NULL) {

  glist = .buildGenolistX(x, marker, eliminate, treatAsFounder = treatAsFounder)

  if (attr(glist, "impossible"))
    return(structure(list(), impossible = TRUE))

  FOU = founders(x, internal = TRUE)

  # Add any members which should be treated as founders
  FOU = c(FOU, treatAsFounder)

  sex = x$SEX
  afr = afreq(marker)
  impossible = FALSE

  dat = lapply(1:pedsize(x), function(i) {

    # If impossible, speed through
    if(impossible)
      return(NULL)

    g = glist[[i]]
    if (i %in% FOU) {
      prob = switch(sex[i], afr[g$mat], HWprob(g$pat, g$mat, afr))
      if (sum(prob) == 0) {
        impossible = TRUE
        return(NULL)
      }
    }
    else
      prob = rep.int(1, length(g$mat))

    g$prob = as.numeric(prob)
    g
  })

  attr(dat, "impossible") = impossible
  dat
}


##################################
# STARTDATA FOR TWO LINKED MARKERS
##################################

startprob_MM_AUT = function(g, afreq1, afreq2, founder) {
  if (founder) {
    p1 = g$pat1
    m1 = g$mat1
    p2 = g$pat2
    m2 = g$mat2

    # Multiply with two if heteroz for at least 1 marker.
    # If heteroz for both, then both phases are included in g, hence the factor 2 (not 4) in this case as well.
    hetfact = ((p1 != m1) | (p2 != m2)) + 1
    afreq1[p1] * afreq1[m1] * afreq2[p2] * afreq2[m2] * hetfact
  }
  else {
    rep(1, length(g$mat1))
  }
}

startprob_MM_X = function(g, afreq1, afreq2, sex, founder) {
  if (founder && sex == 1)
    afreq1[g$mat1] * afreq2[g$mat2]
  else
    startprob_MM_AUT(g, afreq1, afreq2, founder)
}


startdata_MM = function(x, marker1, marker2, eliminate = 0, treatAsFounder = NULL) {
  if(isXmarker(marker1))
    startdata_MM_X(x, marker1, marker2, eliminate = eliminate, treatAsFounder = treatAsFounder)
  else
    startdata_MM_AUT(x, marker1, marker2, eliminate = eliminate, treatAsFounder = treatAsFounder)
}

startdata_MM_X = function(x, marker1, marker2, eliminate = 0, treatAsFounder = NULL) {
  glist1 = .buildGenolistX(x, marker1, eliminate, treatAsFounder)
  glist2 = .buildGenolistX(x, marker2, eliminate, treatAsFounder)
  if (attr(glist1, "impossible") || attr(glist2, "impossible")) {
    dat = list()
    attr(dat, "impossible") = TRUE
    return(dat)
  }

  isFounder = logical(pedsize(x))
  isFounder[founders(x, internal = TRUE)] = TRUE
  isFounder[treatAsFounder] = TRUE

  sex = x$SEX
  impossible = FALSE

  dat = lapply(1:pedsize(x), function(i) {

    # If impossible, speed through
    if(impossible)
      return(NULL)

    sexi = sex[i]
    g1 = glist1[[i]]
    g2 = glist2[[i]]
    len1 = length(g1$mat)
    len2 = length(g2$mat)

    idx1 = rep(seq_len(len1), each = len2)
    idx2 = rep(seq_len(len2), times = len1)

    if (sexi == 1)
      g = list(mat1 = g1$mat[idx1], mat2 = g2$mat[idx2])
    else {
      pat1 = g1$pat[idx1]
      mat1 = g1$mat[idx1]
      pat2 = g2$pat[idx2]
      mat2 = g2$mat[idx2]

      if (isFounder[i]) {
        # Doubly heterozygous founders: Include the other phase as well.
        # (Since .buildGenolist() returns unordered genotypes for founders.)
        doublyhet = pat1 != mat1 & pat2 != mat2
        if (any(doublyhet)) {
          pat1 = c(pat1, pat1[doublyhet])
          mat1 = c(mat1, mat1[doublyhet])
          p2 = pat2; m2 = mat2
          pat2 = c(p2, m2[doublyhet])  # note switch
          mat2 = c(m2, p2[doublyhet])
        }
      }

      g = list(pat1 = pat1, mat1 = mat1, pat2 = pat2, mat2 = mat2)
    }

    prob = startprob_MM_X(g, afreq1 = afreq(marker1),
      afreq2 = afreq(marker2), sex = sexi, founder = isFounder[i])

    g$prob = prob

    keep = prob > 0
    if (!any(keep)){
      impossible = TRUE
      return(NULL)
    }

    if(!all(keep))
      g[] = lapply(g, function(vec) vec[keep])

    g
  })
  attr(dat, "impossible") = impossible
  dat
}

startdata_MM_AUT = function(x, marker1, marker2, eliminate = 0, treatAsFounder = NULL) {
  glist1 = .buildGenolist(x, marker1, eliminate, treatAsFounder)
  glist2 = .buildGenolist(x, marker2, eliminate, treatAsFounder)
  if (attr(glist1, "impossible") || attr(glist2, "impossible")) {
    dat = list()
    attr(dat, "impossible") = TRUE
    return(dat)
  }

  isFounder = logical(pedsize(x))
  isFounder[founders(x, internal = TRUE)] = TRUE
  isFounder[treatAsFounder] = TRUE
  impossible = FALSE

  dat = lapply(1:pedsize(x), function(i) {

    # If impossible, speed through
    if(impossible)
      return(NULL)

    g1 = glist1[[i]]
    g2 = glist2[[i]]
    len1 = length(g1$mat)
    len2 = length(g2$mat)
    idx1 = rep(seq_len(len1), each = len2)
    idx2 = rep(seq_len(len2), times = len1)

    pat1 = g1$pat[idx1]
    mat1 = g1$mat[idx1]
    pat2 = g2$pat[idx2]
    mat2 = g2$mat[idx2]

    if (isFounder[i]) {
      # Doubly heterozygous founders: Include the other phase as well.
      # (Since .buildGenolist() returns unordered genotypes for founders.)
      doublyhet = pat1 != mat1 & pat2 != mat2
      if (any(doublyhet)) {
        pat1 = c(pat1, pat1[doublyhet])
        mat1 = c(mat1, mat1[doublyhet])
        p2 = pat2; m2 = mat2
        pat2 = c(p2, m2[doublyhet])  # note switch
        mat2 = c(m2, p2[doublyhet])
      }
    }

    g = list(pat1 = pat1, mat1 = mat1, pat2 = pat2, mat2 = mat2)

    # Add probabilities
    prob = startprob_MM_AUT(g, afreq1 = afreq(marker1),
                            afreq2 = afreq(marker2), founder = isFounder[i])
    g$prob = prob

    keep = prob > 0
    if (!any(keep)){
      impossible = TRUE
      return(NULL)
    }

    if(!all(keep))
      g[] = lapply(g, function(vec) vec[keep])

    g
  })

  attr(dat, "impossible") = impossible
  dat
}


#### .buildGenolist and ELIMINATE

.genotypeHaploList = function(gt, n, unordered, complete = NULL) {
  nseq = seq_len(n)

  # The complete matrix can (should!) be supplied to avoid making it each time
  complete = complete %||% list(pat = rep(nseq, each = n), mat = rep.int(nseq, times = n))

  a = gt[1]
  b = gt[2]

  if(a == 0 && b == 0)
    g = complete
  else if (a == 0 && b > 0)
    g = list(pat = c(nseq, rep(b, n - 1)), mat = c(rep(b, n), nseq[-b]))
  else if (a > 0 && b == 0)
    g = list(pat = c(nseq, rep(a, n - 1)), mat = c(rep(a, n), nseq[-a]))
  else if (a == b)
    g = list(pat = a, mat = b)
  else
    g = list(pat = c(a, b), mat = c(b, a))

  if (unordered) { # TODO: condition on this earlier?
    keep = g$pat <= g$mat
    g$pat = g$pat[keep]
    g$mat = g$mat[keep]
  }

  g
}

.buildGenolist = function(x, marker, eliminate = 0, treatAsFounder = NULL) {
  n = nAlleles(marker)

  # Founders (except loop breaker copies) need only *unordered* genotypes
  founderNotLB = logical(pedsize(x))
  founderNotLB[founders(x, internal = TRUE)] = TRUE
  founderNotLB[treatAsFounder] = TRUE
  founderNotLB[x$LOOP_BREAKERS[, 2]] = FALSE

  # A matrix containing a complete set of ordered genotypes
  nseq = seq_len(n)
  COMPLETE = list(pat = rep(nseq, each = n), mat = rep.int(nseq, times = n))

  # Building a list of genotypes for each indiv.
  genolist = lapply(1:pedsize(x), function(i) {
    gt = marker[i, ]
    unordered = founderNotLB[i]
    .genotypeHaploList(gt, n, unordered, COMPLETE)
  })

  #als = alleles(marker)
  #genolist = lapply(genolist, function(g) {colnames(g) = paste(als[g[1,]], als[g[2,]], sep = "/"); g})

  attr(genolist, "impossible") = FALSE

  # If mutations, don't eliminate any genotypes
  if (allowsMutations(marker))
    return(genolist)

  .eliminate(x, genolist, n, repeats = eliminate, treatAsFounder = treatAsFounder)
}



.eliminate = function(x, genolist, nall, repeats = 0, treatAsFounder = NULL) {
  if (repeats == 0 || attr(genolist, "impossible"))
    return(genolist)
  N = pedsize(x)

  FOU = founders(x, internal = TRUE)
  NONFOU = nonfounders(x, internal = TRUE)

  # Adjust for "extra" founders
  if(length(treatAsFounder) > 0) {
    FOU = c(FOU, treatAsFounder)
    NONFOU = .mysetdiff(NONFOU, treatAsFounder)
  }

  offs = lapply(1:N, function(i) children(x, i, internal = TRUE))

  # Current number of genotypes of each
  Ngt = unlist(lapply(genolist, function(g) length(g$pat)))

  informative = logical(N)

  for (k in seq_len(repeats)) {
    ng = Ngt

    # Who are now informative (i.e. reduced genotype list)?
    informative[FOU] = (ng[FOU] < nall * (nall + 1)/2)
    informative[NONFOU] = (ng[NONFOU] < nall^2)

    # Loop through entire pedigree
    for (i in 1:N) {
      if (ng[i] == 1)
        next

      g = genolist[[i]]
      sx = x$SEX[i]

      if (i %in% NONFOU) {
        if(informative[fa <- x$FIDX[i]]) {
          kp = g$pat %in% genolist[[fa]]$pat | g$pat %in% genolist[[fa]]$mat
          if(!all(kp)) {
            g$pat = g$pat[kp]
            g$mat = g$mat[kp]
          }
        }
        if(informative[mo <- x$MIDX[i]]) {
          kp = g$mat %in% genolist[[mo]]$pat | g$mat %in% genolist[[mo]]$mat
          if(!all(kp)) {
            g$pat = g$pat[kp]
            g$mat = g$mat[kp]
          }
        }
      }

      barn = offs[[i]]
      for (b in barn[informative[barn]]) {
        hap = genolist[[b]][[sx]]
        kp = g$pat %in% hap | g$mat %in% hap
        if(!all(kp)) {
          g$pat = g$pat[kp]
          g$mat = g$mat[kp]
        }
      }

      genolist[[i]] = g
    }

    Ngt = unlist(lapply(genolist, function(g) length(g$pat)))
    if (any(Ngt == 0)) {
      attr(genolist, "impossible") = TRUE
      break
    }

    # If no change, break
    if (sum(Ngt) == sum(ng))
      break
  }

  genolist
}


#------------X-linked-------------------

.buildGenolistX <- function(x, marker, eliminate, treatAsFounder = NULL) {

  n = nAlleles(marker)
  nseq = seq_len(n)

  # Founders (except loop breaker copies) need only *unordered* genotypes
  founderNotLB = logical(pedsize(x))
  founderNotLB[founders(x, internal = TRUE)] = TRUE
  founderNotLB[treatAsFounder] = TRUE
  founderNotLB[x$LOOP_BREAKERS[, 2]] = FALSE

  genolist = vector(pedsize(x), mode = "list")

  # Males
  SEX = x$SEX
  a1 = marker[, 1]
  genolist[SEX == 1 & a1 == 0] = list(list(mat = nseq))
  genolist[SEX == 1 & a1 != 0] = lapply(a1[SEX == 1 & a1 != 0], function(a) list(mat = a))

  # Females
  COMPLETE = list(pat = rep(nseq, each = n), mat = rep.int(nseq, times = n))
  women = females(x, internal = TRUE)
  genolist[women] = lapply(women, function(i) {
    gt = marker[i, ]
    unordered = founderNotLB[i]
    .genotypeHaploList(gt, n, unordered, COMPLETE)
  })

  attr(genolist, "impossible") = FALSE

  # If mutations, don't eliminate any genotypes
  if (allowsMutations(marker))
    return(genolist)

  .eliminateX(x, genolist, n, eliminate, treatAsFounder)
}


.eliminateX = function(x, genolist, nall, repeats = 0, treatAsFounder = NULL) {
  if (repeats == 0 || attr(genolist, "impossible"))
    return(genolist)

  SEX = x$SEX
  FIDX = x$FIDX
  MIDX = x$MIDX
  FOU = founders(x, internal = TRUE)
  NONFOU = nonfounders(x, internal = TRUE)

  # Adjust for "extra" founders
  if(length(treatAsFounder) > 0) {
    FOU = c(FOU, treatAsFounder)
    NONFOU = .mysetdiff(NONFOU, treatAsFounder)
  }

  xsize = pedsize(x)
  males = males(x, internal = TRUE)
  females = females(x, internal = TRUE)
  femFou = .myintersect(females, FOU)
  femNonfou = .myintersect(females, NONFOU)

  isNonfou = (1:xsize) %in% NONFOU

  offs = lapply(1:xsize, function(i) children(x, i, internal = TRUE))

  informative = logical(xsize)

  Ngt = unlist(lapply(genolist, function(g) length(g$mat)))

  for (k in seq_len(repeats)) {
    ng = Ngt

    # Who are now informative (i.e. reduced genotype list)?
    informative[males] = (ng[males] < nall)
    informative[femFou] = (ng[femFou] < nall * (nall + 1)/2)
    informative[femNonfou] = (ng[femNonfou] < nall^2)

    for (i in males) {
      if (ng[i] == 1)
        next

      g = genolist[[i]]
      if (isNonfou[i] && informative[mo <- MIDX[i]]) {
        kp = g$mat %in% genolist[[mo]]$pat | g$mat %in% genolist[[mo]]$mat
        if(!all(kp))
          g$mat = g$mat[kp]
      }

      barn = offs[[i]]
      for (b in barn[informative[barn] & SEX[barn] == 2]) {
        kp = g$mat %in% genolist[[b]]$pat
        if(!all(kp))
          g$mat = g$mat[kp]
      }

      genolist[[i]] = g
    }

    for (i in females) {
      if (ng[i] == 1)
        next
      g = genolist[[i]]

      if(isNonfou[i]) {
        if(informative[fa <- x$FIDX[i]]) {
          kp = g$pat %in% genolist[[fa]]$mat
          if(!all(kp)) {
            g$pat = g$pat[kp]
            g$mat = g$mat[kp]
          }
        }
        if(informative[mo <- x$MIDX[i]]) {
          kp = g$mat %in% genolist[[mo]]$pat | g$mat %in% genolist[[mo]]$mat
          if(!all(kp)) {
            g$pat = g$pat[kp]
            g$mat = g$mat[kp]
          }
        }
      }

      barn = offs[[i]]
      for (b in barn[informative[barn]]) {
        hap = genolist[[b]]$mat
        kp = g$pat %in% hap | g$mat %in% hap
        if(!all(kp)) {
          g$pat = g$pat[kp]
          g$mat = g$mat[kp]
        }
      }

      genolist[[i]] = g
    }

    Ngt = unlist(lapply(genolist, function(g) length(g$mat)))

    if (any(Ngt == 0)) {
      attr(genolist, "impossible") = TRUE
      break
    }

    # If no change, break
    if (sum(Ngt) == sum(ng))
      break
  }

  genolist
}
