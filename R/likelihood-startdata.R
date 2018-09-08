#### FUNCTIONS FOR CREATING THE INTITIAL HAPLOTYPE COMBINATIONS W/PROBABILITIES.

startdata_M = function(x, marker, eliminate = 0) {
  if(is_Xmarker(marker))
    startdata_M_X(x, marker, eliminate)
  else
    startdata_M_AUT(x, marker, eliminate)
}

startdata_M_AUT = function(x, marker, eliminate = 0) {
  afr = afreq(marker)
  impossible = FALSE

  glist = .build_genolist(x, marker, eliminate)
  if (attr(glist, "impossible")) { #TODO why not simply return dat here?
    #dat = list()
    #attr(dat, "impossible") = TRUE
    return(dat)
  }

  FOU = founders(x, internal=T)

  # Founder inbreeding: A vector of length pedsize(x), with NA's at nonfounders
  # Enables quick look-up e.g. FOU_INB[i].
  FOU_INB = rep(NA_real_, pedsize(x))
  FOU_INB[FOU] = founder_inbreeding(x, ids=founders(x))

  # Add any members which should be treated as founders
  FOU = c(FOU, attr(x, "treat_as_founders"))

  dat = lapply(1:pedsize(x), function(i) {
    h = glist[[i]]
    if (i %in% FOU) {
      #prob = afr[h[1, ]] * afr[h[2, ]] * ((h[1, ] != h[2, ]) + 1)
      prob = HW_prob(h[1, ], h[2, ], afr, f = FOU_INB[i])
      if (sum(prob) == 0)
        impossible = TRUE
    }
    else prob = rep.int(1, ncol(h))
    list(hap = h, prob = as.numeric(prob))
  })

  attr(dat, "impossible") = impossible
  dat
}

startdata_M_X = function(x, marker, eliminate = 0) {
  sex = x$SEX
  afr = afreq(marker)

  impossible = FALSE

  glist = .build_genolist_X(x, marker, eliminate)
  if (attr(glist, "impossible")) {
    dat = list()
    attr(dat, "impossible") = TRUE
    return(dat)
  }

  FOU = founders(x, internal=T)

  dat = lapply(1:pedsize(x), function(i) {
    h = glist[[i]]
    if (i %in% FOU) {
      prob = switch(sex[i], afr[h], HW_prob(h[1, ], h[2, ], afr))
                  #  afr[h[1, ]] * afr[h[2, ]] * ((h[1, ] != h[2, ]) + 1))
      if (sum(prob) == 0)
        impossible = TRUE
    }
    else
      prob = rep.int(1, length(h)/sex[i])

    list(hap = h, prob = as.numeric(prob))
  })
  attr(dat, "impossible") = impossible
  dat
}


##################################
# STARTDATA FOR TWO LINKED MARKERS
##################################

startprob_MM_AUT = function(h, afreq1, afreq2, founder) {
  if (founder) {
    m1_1 = h[1, ]
    m1_2 = h[2, ]
    m2_1 = h[3, ]
    m2_2 = h[4, ]
    # multiply with two if heteroz for at least 1 marker. If heteroz for both, then both phases are included in h, hence the factor 2 (not 4) in this case as well.
    hetfact = ((m1_1 != m1_2 | m2_1 != m2_2) + 1)
    afreq1[m1_1] * afreq1[m1_2] * afreq2[m2_1] * afreq2[m2_2] * hetfact
  }
  else {
    rep(1, ncol(h))
  }
}

startprob_MM_X = function(h, afreq1, afreq2, sex, founder) {
  if (founder && sex == 1)
    afreq1[h[1, ]] * afreq2[h[2, ]]
  else
    startprob_MM_AUT(h, afreq1, afreq2, founder)
}


startdata_MM = function(x, marker1, marker2, eliminate = 0) {
  if(is_Xmarker(marker1))
    startdata_MM_X(x, marker1, marker2, eliminate=eliminate)
  else
    startdata_MM_AUT(x, marker1, marker2, eliminate=eliminate)
}

startdata_MM_X = function(x, marker1, marker2, eliminate = 0) {
  m1_list = .build_genolist_X(x, marker1, eliminate)
  m2_list = .build_genolist_X(x, marker2, eliminate)
  if (attr(m1_list, "impossible") || attr(m2_list, "impossible")) {
    dat = list()
    attr(dat, "impossible") = TRUE
    return(dat)
  }

  afreq1 = afreq(marker1)
  afreq2 = afreq(marker2)
  is_founder = logical(pedsize(x))
  is_founder[founders(x, internal=T)] = TRUE
  is_founder[attr(x, "treat_as_founder")] = TRUE

  sex = x$SEX
  impossible = FALSE

  dat = lapply(1:pedsize(x), function(i) {
    sexi = sex[i]
    h1 = m1_list[[i]]
    h2 = m2_list[[i]]
    if (sexi == 1)
      hap = rbind(rep(h1, each = length(h2)), rep(h2, times = length(h1)))  #matrix with two rows: m1, m2
    else {
      hl1 = dim(h1)[2]
      hl2 = dim(h2)[2]
      hap = rbind(h1[, rep(seq_len(hl1), each = hl2), drop = F],
                  h2[, rep(seq_len(hl2), times = hl1), drop = F])  #matrix with four rows: m1_1, m1_2, m2_1, m2_2
      if (is_founder[i]) {
        # Doubly heterozygous founders: Include the other phase as well.
        # (Since .build_genolist() returns unordered genotypes for founders.)
        doublyhet = hap[1, ] != hap[2, ] & hap[3, ] != hap[4, ]
        if (any(doublyhet))
          hap = cbind(hap, hap[c(1, 2, 4, 3), doublyhet, drop = F])
      }
    }
    prob = startprob_MM_X(hap, afreq1 = afreq1, afreq2 = afreq2, sex = sexi, founder = is_founder[i])
    keep = prob > 0
    if (!any(keep))
      impossible = TRUE
    list(hap = hap[, keep, drop = F], prob = as.numeric(prob[keep]))
  })
  attr(dat, "impossible") = impossible
  dat
}

startdata_MM_AUT = function(x, marker1, marker2, eliminate = 0) {
  m1_list = .build_genolist(x, marker1, eliminate)
  m2_list = .build_genolist(x, marker2, eliminate)
  if (attr(m1_list, "impossible") || attr(m2_list, "impossible")) {
    dat = list()
    attr(dat, "impossible") = TRUE
    return(dat)
  }

  afreq1 = afreq(marker1)
  afreq2 = afreq(marker2)
  is_founder = logical(pedsize(x))
  is_founder[founders(x, internal=T)] = TRUE
  is_founder[attr(x, "treat_as_founder")] = TRUE
  impossible = FALSE

  dat = lapply(1:pedsize(x), function(i) {
    h1 = m1_list[[i]]
    hl1 = dim(h1)[2]
    h2 = m2_list[[i]]
    hl2 = dim(h2)[2]
    hap = rbind(h1[, rep(seq_len(hl1), each = hl2), drop = F],
                h2[, rep(seq_len(hl2), times = hl1), drop = F])  #matrix with four rows: m1_1, m1_2, m2_1, m2_2
    if (is_founder[i]) {
      # Doubly heterozygous founders: Include the other phase as well.
      # (Since .build_genolist() returns unordered genotypes for founders.)
      doublyhet = hap[1, ] != hap[2, ] & hap[3, ] != hap[4, ]
      if (any(doublyhet))
        hap = cbind(hap, hap[c(1, 2, 4, 3), doublyhet, drop = FALSE])
    }
    prob = startprob_MM_AUT(hap, afreq1 = afreq1, afreq2 = afreq2, founder = is_founder[i])
    keep = prob > 0
    if (!any(keep))
      impossible = TRUE
    list(hap = hap[, keep, drop = F], prob = as.numeric(prob[keep]))
  })
  attr(dat, "impossible") = impossible
  dat
}


#### .BUILD_GENOLIST and ELIMINATE

.genotypeMatrix = function(gt, n, unordered, complete=NULL) {
  nseq = seq_len(n)

  # The complete matrix can (should!) be supplied to avoid making it each time
  if(is.null(complete))
    complete = rbind(rep(nseq, each = n), rep.int(nseq, times = n))

  a = gt[1]
  b = gt[2]

  if(a==0 && b==0)
    m = complete
  else if (a==0 && b>0)
    m = rbind(c(nseq, rep(b, n - 1)), c(rep(b, n), nseq[-b]))
  else if (a>0 && b==0)
    m = rbind(c(nseq, rep(a, n - 1)), c(rep(a, n), nseq[-a]))
  else if (a==b)
    m = cbind(c(a,b))
  else
    m = cbind(c(a,b), c(b, a))

  if (unordered) # TODO: condition on this earlier?
    m = m[, m[1, ] <= m[2, ], drop = FALSE]
  m
}

.build_genolist = function(x, marker, eliminate = 0) {
  n = nAlleles(marker)

  # Founders (except loop breaker copies) need only *unordered* genotypes
  founder_not_loopbreaker = logical(pedsize(x))
  founder_not_loopbreaker[founders(x, internal=T)] = TRUE
  founder_not_loopbreaker[attr(x, "treat_as_founder")] = TRUE
  founder_not_loopbreaker[x$LOOP_BREAKERS[, 2]] = FALSE

  # A matrix containing a complete set of ordered genotypes
  nseq = seq_len(n)
  COMPLETE = rbind(rep(nseq, each = n), rep.int(nseq, times = n))

  # Building a list of genotypes for each indiv.
  genolist = lapply(1:pedsize(x), function(i) {
    gt = marker[i, ]
    unordered = founder_not_loopbreaker[i]
    .genotypeMatrix(gt, n, unordered, COMPLETE)
  })

  attr(genolist, "impossible") = FALSE

  # If mutations, don't eliminate any genotypes
  if (allowsMutations(marker))
    return(genolist)

  .eliminate(x, genolist, n, repeats = eliminate)
}



.eliminate = function(x, genolist, nall, repeats = 0) {
  if (repeats == 0 || attr(genolist, "impossible"))
    return(genolist)
  N = pedsize(x)

  FOU = founders(x, internal=T)
  NONFOU = nonfounders(x, internal=T)

  # Adjust for "extra" founders
  treat_as_fou = attr(x, "treat_as_founder")
  if(length(treat_as_fou) > 0) {
    FOU = c(FOU, treat_as_fou)
    NONFOU = setdiff(NONFOU, treat_as_fou)
  }

  offs = lapply(1:N, function(i) children(x, i, internal=TRUE))
  ncols_ny = unlist(lapply(genolist, ncol))

  informative = logical(N)
  for (k in seq_len(repeats)) {
    ncols = ncols_ny
    informative[FOU] = (ncols[FOU] < nall * (nall + 1)/2)
    informative[NONFOU] = (ncols[NONFOU] < nall^2)
    for (i in 1:N) {
      if (ncols[i] == 1)
        next
      g = genolist[[i]]
      kjonn = x$SEX[i]
      if (i %in% NONFOU && informative[far <- x$FIDX[i]])
        g = g[, g[1, ] %in% genolist[[far]][1, ] | g[1, ] %in% genolist[[far]][2, ],
              drop = F]
      if (i %in% NONFOU && informative[mor <- x$MIDX[i]])
        g = g[, g[2, ] %in% genolist[[mor]][1, ] | g[2, ] %in% genolist[[mor]][2, ],
              drop = F]
      barn = offs[[i]]
      for (b in barn[informative[barn]]) {
        g = g[, g[1, ] %in% genolist[[b]][kjonn, ] | g[2, ] %in% genolist[[b]][kjonn, ],
              drop = F]
      }
      genolist[[i]] = g
    }
    ncols_ny = unlist(lapply(genolist, ncol))
    if (any(ncols_ny == 0)) {
      attr(genolist, "impossible") = TRUE
      return(genolist)
    }
    if (sum(ncols_ny) == sum(ncols))
      return(genolist)
  }
  genolist
}


#------------X-linked-------------------

.build_genolist_X <- function(x, marker, eliminate) {

  n = nAlleles(marker)
  nseq = seq_len(n)

  # Founders (except loop breaker copies) need only *unordered* genotypes
  founder_not_loopbreaker = logical(pedsize(x))
  founder_not_loopbreaker[founders(x, internal=T)] = TRUE
  founder_not_loopbreaker[attr(x, "treat_as_founder")] = TRUE
  founder_not_loopbreaker[x$LOOP_BREAKERS[, 2]] = FALSE

  genolist = vector(pedsize(x), mode="list")

  # Males
  SEX = x$SEX
  genolist[SEX == 1 & marker[, 1] == 0] = list(nseq)
  genolist[SEX == 1 & marker[, 1] != 0] = marker[SEX == 1 & marker[, 1] != 0, 1]

  # Females
  COMPLETE = rbind(rep(nseq, each = n), rep.int(nseq, times = n))
  women = females(x, internal=T)
  genolist[women] = lapply(women, function(i) {
    gt = marker[i, ]
    unordered = founder_not_loopbreaker[i]
    .genotypeMatrix(gt, n, unordered, COMPLETE)
  })

  attr(genolist, "impossible") = FALSE

  # If mutations, don't eliminate any genotypes
  if (allowsMutations(marker))
    return(genolist)

  .eliminate_X(x, genolist, n, eliminate)
}


.eliminate_X = function(x, genolist, nall, repeats = 0) {
  if (repeats == 0 || attr(genolist, "impossible"))
    return(genolist)

  SEX = x$SEX
  FIDX = x$FIDX
  MIDX = x$MIDX
  FOU = founders(x, internal=T)
  NONFOU = nonfounders(x, internal=T)

  # Adjust for "extra" founders
  treat_as_fou = attr(x, "treat_as_founder")
  if(length(treat_as_fou) > 0) {
    FOU = c(FOU, treat_as_fou)
    NONFOU = setdiff(NONFOU, treat_as_fou)
  }

  xsize = pedsize(x)
  males = males(x, internal=T)
  females = females(x, internal=T)
  fem_fou = .myintersect(females, FOU)
  fem_nonfou = .myintersect(females, NONFOU)

  is_nonfou = (1:xsize) %in% NONFOU

  offs = lapply(1:xsize, function(i) children(x, i, internal=TRUE))

  informative = logical(xsize)
  ncols_ny = lengths(genolist)/SEX  #males are vectors, females matrices w/ 2 rows
  for (k in seq_len(repeats)) {
    ncols = ncols_ny
    informative[males] = (ncols[males] < nall)
    informative[fem_fou] = (ncols[fem_fou] < nall * (nall + 1)/2)
    informative[fem_nonfou] = (ncols[fem_nonfou] < nall^2)
    for (i in males) {
      if (ncols[i] == 1) next
      g = genolist[[i]]
      if (is_nonfou[i] && informative[mor <- MIDX[i]])
        g = g[g %in% genolist[[mor]][1, ] | g %in% genolist[[mor]][2, ]]
      barn = offs[[i]]
      for (b in barn[informative[barn] & SEX[barn] == 2])
        g = g[g %in% genolist[[b]][1, ]]
      genolist[[i]] = g
    }
    for (i in females) {
      if (ncols[i] == 1) next
      g = genolist[[i]]
      if (is_nonfou[i] && informative[far <- FIDX[i]])
        g = g[, g[1, ] %in% genolist[[far]], drop = F]
      if (is_nonfou[i] && informative[mor <- MIDX[i]])
        g = g[, g[2, ] %in% genolist[[mor]][1, ] | g[2, ] %in% genolist[[mor]][2, ],
              drop = F]
      barn = offs[[i]]
      for (b in barn[informative[barn]]) {
        if (SEX[b] == 1)
          g = g[, g[1, ] %in% genolist[[b]] | g[2, ] %in% genolist[[b]], drop = F]
        else
          g = g[, g[1, ] %in% genolist[[b]][2, ] | g[2, ] %in% genolist[[b]][2, ], drop = F]
      }
      genolist[[i]] = g
    }
    ncols_ny = lengths(genolist)/SEX
    if (any(ncols_ny == 0)) {
      attr(genolist, "impossible") = TRUE
      return(genolist)
    }
    if (sum(ncols_ny) == sum(ncols))
      return(genolist)
  }
  genolist
}

