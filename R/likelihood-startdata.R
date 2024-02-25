#### FUNCTIONS FOR CREATING THE INTITIAL HAPLOTYPE COMBINATIONS W/PROBABILITIES.

startdata_M = function(x, marker, eliminate = 0, treatAsFounder = NULL) {#print("new startdata")

  nInd = length(x$ID)
  Xchrom = isXmarker(marker)
  SEX = x$SEX

  # Founders (except LB-copies): Genotypes need not be phased
  fouInt = founders(x, internal = TRUE)
  isFounder = logical(nInd)
  isFounder[fouInt] = TRUE
  isFounder[treatAsFounder] = TRUE
  isFounder[x$LOOP_BREAKERS[, 2]] = FALSE

  # Founders with 1 child: Skip genotypes, draw alleles directly (assuming HWE)
  simpleFou = isFounder & tabulate(c(x$FIDX, x$MIDX), nInd) == 1

  # Founder inbreeding lookup vector; 0's for nonfounders
  fi = numeric(nInd)
  fval = x$FOUNDER_INBREEDING[[if(Xchrom) "x" else "autosomal"]]
  if(!is.null(fval))
    fi[fouInt] = fval

  marker = .sortGeno(marker)
  afr = attr(marker, "afreq")
  n = length(afr)

  # Complete set of genotypes - phased and unphased (for founders) - used below
  nseq = seq_len(n)
  allPhased   = list(pat = rep(nseq, each = n), mat = rep.int(nseq, times = n))
  allUnphased = list(pat = rep(nseq, n:1), mat = sequence.default(n:1, from = nseq))

  # Loop through all individuals
  glist = vector(nInd, mode = "list")
  imp = FALSE

  for(i in seq_along(glist)) {
    a = marker[i, 1]
    b = marker[i, 2]

    if(Xchrom && SEX[i] == 1) { # Hemizygous male
      mat = if(b == 0) nseq else b
      prob = if(isFounder[i]) afr[mat] else rep(1, length(mat))
      g = list(mat = mat, prob = prob)
    }
    else if(simpleFou[i]) # Simple founder: alleles directly (skip genotypes)
      g = .alleleDistrib(a, b, afr, f = fi[i])
    else if(isFounder[i]) # General founder: unphased genotypes with HW probs
      g = .genoDistribFounder(a, b, afr, f = fi[i], COMPLETE = allUnphased)
    else # Nonfounder: Phased genos with 1's as prob
      g = .genoDistribNonfounder(a, b, COMPLETE = allPhased)

    # Remove entries with prob = 0
    g = .reduce(g)
    glist[[i]] = g

    # If all zero: impossible!
    if(!length(g$prob)) {imp = TRUE; break}
  }

  names(glist) = x$ID
  attr(glist, "impossible") = imp

  ### Eliminate # TODO!! Previous implement was applied before probs
  # if (eliminate > 0 && !allowsMutations(marker))
    # glist = eliminate(glist, ...)

  if (attr(glist, "impossible"))
    glist = structure(list(), impossible = TRUE)

  glist
}


.alleleDistrib = function(a, b, afr, f = 0) {
  nseq = seq_along(afr)

  if(a == 0 && b == 0) { # untyped founder
    allele = nseq
    prob = afr
  }
  else if (a == 0) { # partial founder genotype: -/b
    allele = nseq
    prob = afr/2
    prob[b] = prob[b] + 0.5
    if(f > 0) {
      prob = prob * (1-f)
      prob[b] = prob[b] + f
    }
  }
  else if (a == b) { # homozygous founder
    allele = a
    prob = afr[a]^2
    if(f > 0)
      prob = prob * (1-f) + afr[a] * f
  }
  else { # heterozygous founder
    allele = c(a,b)
    prob = rep(afr[a] * afr[b] * (1-f), 2)   # 2*a*b * 0.5
  }

  list(allele = allele, prob = prob)
}

.genoDistribFounder = function(a, b, afr, f = 0, COMPLETE) {
  # Produce list of vectors pat, mat and prob
  # Unphased genotypes; use HW for probs
  # Assumes a <= b (single integers)

  if(a == 0 && b == 0)
    g = COMPLETE
  else if (a == 0) {
    n = max(COMPLETE$pat)
    g = list(pat = seq_len(n), mat = rep(b, n))
  }
  else
    g = list(pat = a, mat = b)

  g$prob = HWprob(g$pat, g$mat, afr, f)
  g
}

.genoDistribNonfounder = function(a, b, COMPLETE) {
  # Produce list of vectors pat, mat and prob (all 1's)
  # Includes both phases of heterozygous genotypes
  # Assumes a <= b (single integers)

  if(a == 0 && b == 0)
    g = COMPLETE
  else if (a == 0) {
    n = max(COMPLETE$pat)
    g = list(pat = c(seq_len(n), rep(b, n - 1)), mat = c(rep(b, n), seq_len(n)[-b]))
  }
  else if (a == b)
    g = list(pat = a, mat = b)
  else
    g = list(pat = c(a, b), mat = c(b, a))

  g$prob = rep(1, length(g$pat))
  g
}

# TODO: Old version, to be deleted
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

# TODO: Delete
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

  glist1 = startdata_M(x, marker1, eliminate = eliminate, treatAsFounder = treatAsFounder)
  glist2 = startdata_M(x, marker2, eliminate = eliminate, treatAsFounder = treatAsFounder)

  if (attr(glist1, "impossible") || attr(glist2, "impossible"))
    return(structure(list(), impossible = TRUE))

  nInd = length(x$ID)
  Xchrom = isXmarker(marker1)
  SEX = x$SEX

  # Founders (except LB-copies)
  fouInt = founders(x, internal = TRUE)
  isFounder = logical(nInd)
  isFounder[fouInt] = TRUE
  isFounder[treatAsFounder] = TRUE
  isFounder[x$LOOP_BREAKERS[, 2]] = FALSE

  # Loop through all individuals
  glist = vector(nInd, mode = "list")
  imp = FALSE

  for(i in seq_along(glist)) {
    g1 = glist1[[i]]
    g2 = glist2[[i]]
    len1 = length(g1$prob)
    len2 = length(g2$prob)
    idx1 = rep(seq_len(len1), each = len2)
    idx2 = rep(seq_len(len2), times = len1)

    # Same in all cases
    prob = as.numeric(g1$prob[idx1] * g2$prob[idx2])

    if(sum(prob == 0)) {
      imp = TRUE
      break
    }

    # Case 1: Hemizygous
    if(Xchrom && SEX[i] == 1) {
      g = list(mat1 = g1$mat[idx1], mat2 = g2$mat[idx2], prob = prob)
      glist[[i]] = .reduce(g)
      next
    }

    # Case 2: Simple founder
    if(!is.null(g1$allele)) {
      g = list(allele1 = g1$allele[idx1], allele2 = g2$allele[idx2], prob = prob)
      glist[[i]] = .reduce(g)
      next
    }

    # Main case
    pat1 = g1$pat[idx1]
    mat1 = g1$mat[idx1]
    pat2 = g2$pat[idx2]
    mat2 = g2$mat[idx2]

    # Doubly heterozygous founders: Include the other phase as well
    if (isFounder[i]) {
      doublyhet = pat1 != mat1 & pat2 != mat2
      if (any(doublyhet)) {
        pat1 = c(pat1, pat1[doublyhet])
        mat1 = c(mat1, mat1[doublyhet])
        p2 = pat2; m2 = mat2
        pat2 = c(p2, m2[doublyhet])  # note switch
        mat2 = c(m2, p2[doublyhet])

        # Split probabilities also
        prob[doublyhet] = prob[doublyhet]/2
        prob = c(prob, prob[doublyhet])
      }
    }

    g = list(pat1 = pat1, mat1 = mat1, pat2 = pat2, mat2 = mat2, prob = prob)
    glist[[i]] = .reduce(g)
  }

  names(glist) = x$ID
  attr(glist, "impossible") = imp
  glist
}


# TODO: Old version, to be deleted
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

# TODO: Delete
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

# TODO: No longer needed?
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

.reduce = function(g) {
  keep = g$prob > 0
  if(all(keep))
    return(g)
  lapply(g, function(vec) vec[keep])
}
