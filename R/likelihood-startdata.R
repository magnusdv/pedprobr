#### FUNCTIONS FOR CREATING THE INTITIAL HAPLOTYPE COMBINATIONS W/PROBABILITIES.

startdata_M = function(x, marker, eliminate = 0, treatAsFounder = NULL) {
  if(isXmarker(marker))
    startdata_M_X(x, marker, eliminate, treatAsFounder)
  else
    startdata_M_AUT(x, marker, eliminate, treatAsFounder)
}


startdata_M_AUT_new = function(x, marker, eliminate = 0, treatAsFounder = NULL) {#print("new startdata")

  # Sort each genotype
  swtch = marker[,1] > marker[,2]
  if(any(swtch))
    marker[swtch, 1:2] = marker[swtch, 2:1]

  nInd = length(x$ID)
  fouInt = founders(x, internal = TRUE)

  # Founders (except LB-copies): Genotypes need not be sorted
  unsortedFounder = logical(nInd)
  unsortedFounder[fouInt] = TRUE
  unsortedFounder[treatAsFounder] = TRUE
  unsortedFounder[x$LOOP_BREAKERS[, 2]] = FALSE

  # Founders with 1 child: Skip genotypes, draw alleles directly (assuming HWE)
  simpleFou = unsortedFounder & tabulate(c(x$FIDX, x$MIDX), nInd) == 1

  # Founder inbreeding lookup vector; 0's for nonfounders
  fi = numeric(nInd)
  if(!is.null(x$FOUNDER_INBREEDING$autosomal))
    fi[fouInt] = x$FOUNDER_INBREEDING$autosomal

  afr = attr(marker, "afreq")
  n = length(afr)
  nseq = seq_len(n)

  # Complete set of ordered genotypes (used below)
  COMPLETE = list(pat = rep(nseq, each = n), mat = rep.int(nseq, times = n))

  impossible = FALSE

  # Loop through all individuals
  glist = lapply(seq_len(nInd), function(i) {

    a = marker[i, 1]
    b = marker[i, 2]

    ### Founder with 1 child: Draw alleles directly (without genotype)
    if(simpleFou[i]) {
      if(a == 0 && b == 0) { # untyped founder
        allele = nseq
        prob = afr
      }
      else if (a == 0) { # partial founder genotype: -/b
        allele = nseq
        prob = afr/2
        prob[b] = prob[b] + 0.5
        if(fi[i] > 0) {
          prob = prob * (1-fi[i])
          prob[b] = prob[b] + fi[i]
        }
      }
      else if (a == b) { # homozygous founder
        allele = a
        prob = afr[a]^2
        if(fi[i] > 0)
          prob = prob * (1-fi[i]) + afr[a] * fi[i]
      }
      else { # heterozygous founder
        allele = c(a,b)
        prob = rep(afr[a]*afr[b]*(1-fi[i]), 2)   # 2*a*b * 0.5
      }

      return(list(allele = allele, prob = prob))
    }

    ### Otherwise

    if(a == 0 && b == 0)
      g = COMPLETE
    else if (a == 0)
      g = list(pat = c(nseq, rep(b, n - 1)), mat = c(rep(b, n), nseq[-b]))
    else if (a == b)
      g = list(pat = a, mat = b)
    else
      g = list(pat = c(a, b), mat = c(b, a))

    # No phasing for founders (except loop breakers)
    if(unsortedFounder[i]) {
      keep = g$pat <= g$mat
      g$pat = g$pat[keep]
      g$mat = g$mat[keep]
      g$prob = HWprob(g$pat, g$mat, afr, fi[i])
    }
    else
      g$prob = rep.int(1, length(g$mat))

    g
  })

  names(glist) = x$ID
  attr(glist, "impossible") = FALSE # TODO?

  # If mutations, don't eliminate any genotypes
  if (allowsMutations(marker))
    eliminate = 0

  if(eliminate)
    genolist = eliminate()

  if (attr(glist, "impossible"))
    return(structure(list(), impossible = TRUE))

  glist
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
  for(i in seq_along(glist)) {

    g = glist[[i]]
    if (i %in% FOU) {
      prob = HWprob(g$pat, g$mat, afr, f = FOU_INB[i])

      if (sum(prob) == 0){
        impossible = TRUE
        break
      }
    }
    else
      prob = rep.int(1, length(g$mat))

    glist[[i]]$prob = as.numeric(prob)
  }

  attr(glist, "impossible") = impossible
  glist
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

startdata_MM = function(x, marker1, marker2, eliminate = 0, treatAsFounder = NULL) {
  if(isXmarker(marker1))
    startdata_MM_X(x, marker1, marker2, eliminate = eliminate, treatAsFounder = treatAsFounder)
  else
    startdata_MM_AUT(x, marker1, marker2, eliminate = eliminate, treatAsFounder = treatAsFounder)
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
