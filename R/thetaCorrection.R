# v: vector of (integer) alleles; typically all founder alleles
# afr: output of afreq(m)
# theta: correction
thetaCorr = function(v, afr, theta) {
  p = sapply(seq_along(v), function(j) {
    a = v[j]
    bj = sum(v[seq_len(j - 1)] == a)
    (bj * theta + (1 - theta) * afr[a]) / (1 + (j - 2)*theta)
  })
  p
}


likTheta = function(x, m, theta, peeler, peelOrder) {

  #---------------------
  # Mostly for debug (peeler and peelOrder will not be missing in practice)
  if(missing(peeler)) {
    if(!isXmarker(m))
      peeler = function(dat, sub) .peel_M_AUT(dat, sub, mutmat = mutmod(m))
    else
      peeler = function(dat, sub) .peel_M_X(dat, sub, SEX = x$SEX, mutmat = mutmod(m))
  }

  if(missing(peelOrder))
    peelOrder = informativeSubnucs(x, m, peelOrder = peelingOrder(x))
  #----------------------

  # Freqs
  afr = afreq(m)

  # Starting point: all possible genotypes for each
  # NB: no symmetry reduction for founders
  glist = .buildGenolist(x, m, eliminate = 10, foundersUnordered = FALSE)

  # Number of possible genotypes for each
  nGeno = sapply(glist, function(g) length(g$pat))

  # Founders (internal ID)
  FOU = c(founders(x, internal = TRUE), attr(peelOrder, "treatAsFounder"))
  NONFOU = seq_len(pedsize(x))[-FOU]

  # List of all combinations of founder genotypes
  fouGrid = fastGrid(lapply(nGeno[FOU], seq_len), as.list = TRUE)

  # Useful later
  glistFOU = glist[FOU]

  # Initialise likelihood
  likel = 0

  # Loop over combos
  for (gIdx in fouGrid) { # gIdx is a vector of indices, of length same as FOU

    # Collect all founder alleles in that combo
    als = unlist(lapply(seq_along(FOU), function(i) {
      g = glistFOU[[i]]
      idx = gIdx[i]
      c(g$pat[idx], g$mat[idx])
    }))

    # Theta-corrected probs (when als are sampled sequentially)
    probs = thetaCorr(als, afr, theta)

    # Create startdata with single geno for each founder
    # TODO: better to use prod(probs) in first and 1,1,.. for rest?!
    gli = glist
    for(i in FOU)
      gli[[i]] = list(pat = als[2*i-1], mat = als[2*i],
                      prob = probs[2*i-1]*probs[2*i])
    for(i in NONFOU)
      gli[[i]]$prob = rep(1, nGeno[i])

    attr(gli, "impossible") = any(probs == 0)

    li = peelingProcess(x, m, startdata=gli, peeler, peelOrder)
    #print(lapply(gli, pasteGenoProb)); print(li)
    likel = likel + li
  }

  likel
}
