#### FUNCTIONS FOR CREATING THE INTITIAL HAPLOTYPE COMBINATIONS W/PROBABILITIES.

startdata_M = function(x, marker, eliminate = 0, treatAsFounder = NULL) {
  if(isXmarker(marker))
    startdata_M_X(x, marker, eliminate, treatAsFounder)
  else
    startdata_M_AUT(x, marker, eliminate, treatAsFounder)
}

### FOR DEBUGGING
.printG = function(x, ...) {
  als = attr(x, "als") %||% 1:max(unlist(x, use.names = FALSE))
  hasprob = length(x) && !is.null(x[[1]]$prob)
  pretty = lapply(x, function(y) {
    g = paste(als[y$pat], als[y$mat], sep = "/")
    if(hasprob) `names<-`(y$prob, g) else g
  })
  print(pretty)
}

startdata_M_AUT = function(x, marker, eliminate = 0, treatAsFounder = NULL,
                           dropout = 0) {
  nInd = length(x$ID)

  # Dropout loopup: TRUE if homozygous and nonzero factor
  if(length(dropout) == 1)
    dropout = rep_len(dropout, nInd)
  dropoutTF = marker[,1] > 0 & marker[,1] == marker[,2] & dropout > 0

  glist = .buildGenolist(x, marker, eliminate, treatAsFounder, dropoutTF = dropoutTF)

  if(attr(glist, "impossible"))
    return(structure(list(), impossible = TRUE))

  # Founder lookup
  isFOU = x$FID == 0

  # Founder inbreeding lookup, with 0 at nonfounders
  FOU_INB = rep_len(0, nInd)
  FOU_INB[isFOU] = x$FOUNDER_INBREEDING$autosomal %||% 0

  # Add 'extra' founders (NB: after FOU_INB)
  isFOU[treatAsFounder] = TRUE

  afr = attr(marker, "afreq")
  impossible = FALSE

  # Add probabilities to each genotype
  for(i in seq_along(glist)) {
    g = glist[[i]]
    if(isFOU[i])
      prob = HWprob(g$pat, g$mat, afr, f = FOU_INB[i]) # NB: No dropout here!
    else
      prob = rep_len(1, length(g$mat))

    # Dropout
    if(dropoutTF[i]) {
      di = dropout[i]
      prob[g$pat == g$mat] = prob[g$pat == g$mat] * (1 - di^2)
      prob[g$pat != g$mat] = prob[g$pat != g$mat] * di*(1-di)
    }

    if (sum(prob) == 0){
      impossible = TRUE
      break
    }

    glist[[i]]$prob = prob
  }

  attributes(glist) = list(names = x$ID, impossible = impossible,
                           als = attr(marker, "alleles"))
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
