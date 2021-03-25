#### PEELING FUNCTIONS

choosePeeler = function(twolocus, rho, Xchrom, SEX, mutmat) {
  if(!twolocus && !Xchrom)
    .f = function(dat, sub) .peel_M_AUT(dat, sub, mutmat)
  else if(!twolocus && Xchrom)
    .f = function(dat, sub) .peel_M_X(dat, sub, SEX, mutmat)
  else if(twolocus && !Xchrom)
    .f = function(dat, sub) .peel_MM_AUT(dat, sub, rho)
  else if(twolocus && Xchrom)
    .f = function(dat, sub) .peel_MM_X(dat, sub, rho, SEX)
  .f
}

.peel_M_AUT = function(dat, sub, mutmat = NULL) {
  fa = sub$father
  mo = sub$mother
  link = sub$link

  faDat = dat[[fa]]
  moDat = dat[[mo]]
  likel = tcrossprod(faDat$prob, moDat$prob) # faDat$prob %*% t.default(moDat$prob)
  faLen = nrow(likel)
  moLen = ncol(likel)

  # Loop over the children, except the link if this is a child.
  for (ch in .mysetdiff(sub$children, link)) {
    chDat = dat[[ch]]
    chprob = chDat$prob
    chLen = length(chprob)
    transPat = .transProbM(faDat, chDat$pat, mutmat = mutmat$male)
    transMat = .transProbM(moDat, chDat$mat, mutmat = mutmat$female)
    dim(transMat) = c(chLen, moLen)
    transMat_rep = transMat[rep(seq_len(chLen), faLen), ] #as.numeric(do.call(rbind, rep(list(transMat), faLen)))
    mm = .colSums((transPat * chprob) * transMat_rep, chLen, faLen * moLen)
    likel = likel * mm
  }

  if (link == 0)
    return(sum(likel))

  # Goal is to update the probabilities of the link individual ("pivot"):
  pivDat = dat[[link]]

  if (link == fa)
    res = .rowSums(likel, faLen, moLen)
  else if (link == mo)
    res = .colSums(likel, faLen, moLen)
  else { # link is a child
    pivp = pivDat$prob
    pivLen = length(pivp)

    TRarray = array(0, dim = c(faLen, moLen, pivLen))
    transPat = .transProbM(faDat, pivDat$pat, mutmat = mutmat$male)
    dim(transPat) = c(pivLen, faLen)
    transMat = .transProbM(moDat, pivDat$mat, mutmat = mutmat$female)
    dim(transMat) = c(pivLen, moLen)
    for (i in seq_len(faLen)) {
      transpat = transPat[, i]
      for (j in seq_len(moLen))
        TRarray[i, j, ] = transpat * transMat[, j]
    }

    # DEBUG
    # print(array(TRarray, dim = c(faLen, moLen, pivLen), dimnames = lapply(list(faDat, moDat, pivh), pasteHap)))

    arr = as.vector(TRarray) * as.vector(likel)
    dim(arr) = dim(TRarray)
    res = .colSums(arr, faLen * moLen, pivLen)  # sum for each entry of haps[[link]]
    res = res * pivp
  }

  # Update the probabilities
  pivDat$prob = res

  # Reduce if possible
  if(any(res == 0)) {
    pivDat$pat = pivDat$pat[res > 0]
    pivDat$mat = pivDat$mat[res > 0]
    pivDat$prob = pivDat$prob[res > 0]
  }

  dat[[link]] = pivDat

  if (sum(res) == 0)
    attr(dat, "impossible") = TRUE

  dat
}


.peel_M_X = function(dat, sub, SEX, mutmat = NULL) {
  fa = sub$father
  mo = sub$mother
  link = sub$link

  faDat = dat[[fa]]
  moDat = dat[[mo]]
  likel = tcrossprod(faDat$prob, moDat$prob) # faDat$prob %*% t.default(moDat$prob)
  dims = dim(likel)
  faLen = dims[1L]
  moLen = dims[2L]

  # Loop over the children, except the link if this is a child.
  for (ch in .mysetdiff(sub$children, link)) {
    chDat = dat[[ch]]
    chprob = chDat$prob
    chLen = length(chprob)

    if (SEX[ch] == 1) {
      transMat = .transProbM(moDat, chDat$mat, mutmat = mutmat$female)
      mm = rep(.colSums(transMat * chprob, chLen, moLen), each = faLen)
    }
    else {
      transPat =
        if (is.null(mutmat))
          unlist(lapply(faDat$mat, function(fal) as.numeric(fal == chDat$pat)))
        else
          unlist(lapply(faDat$mat, function(fal) mutmat$male[fal, chDat$pat]))

      transMat = .transProbM(moDat, chDat$mat, mutmat = mutmat$female)
      dim(transMat) = c(chLen, moLen)
      transMat_rep = transMat[rep(seq_len(chLen), faLen), ] # as.numeric(do.call(rbind, rep(list(transMat), faLen)))
      mm = .colSums((transPat * chprob) * transMat_rep, chLen, faLen * moLen)
    }
    likel = likel * mm
  }


  if (link == 0)
    return(sum(likel))

  # Goal is to update the probabilities of the link individual ("pivot"):
  pivDat = dat[[link]]

  if (link == fa)
    res = .rowSums(likel, faLen, moLen)
  else if (link == mo)
    res = .colSums(likel, faLen, moLen)
  else { # link is a child
    pivp = pivDat$prob
    pivLen = length(pivp)

    if (SEX[link] == 1) {
     transMat = .transProbM(moDat, pivDat$mat, mutmat = mutmat$female)
     dim(transMat) = c(pivLen, moLen)
     TRarray = rep(t.default(transMat), each = faLen)
    }
    else {
     TRarray = array(0, dim = c(faLen, moLen, pivLen))
     transMat = .transProbM(moDat, pivDat$mat, mutmat = mutmat$female)
     dim(transMat) = c(pivLen, moLen)
     for (i in seq_len(faLen)) {
       transPat =
         if (is.null(mutmat)) as.numeric(faDat$mat[i] == pivDat$pat)
       else mutmat$male[faDat$mat[i], pivDat$pat]
       TRarray[i, , ] = t.default(transPat * transMat)  #TODO:make faster?
     }
    }

    # DEBUG
    # print(array(TRarray, dim = c(faLen, moLen, pivLen), dimnames = lapply(list(faDat, moDat, pivh), pasteHap)))

    arr = as.vector(TRarray) * as.vector(likel)
    dim(arr) = dim(TRarray)
    res = .colSums(arr, faLen * moLen, pivLen)  #sum for each entry of haps[[link]]
    res = res * pivp
  }

  # Update the probabilites
  pivDat$prob = res

  # Reduce if possible
  if(any(res == 0)) {
    if(SEX[link] == 2)
      pivDat$pat = pivDat$pat[res > 0]
    pivDat$mat = pivDat$mat[res > 0]
    pivDat$prob = pivDat$prob[res > 0]
  }

  dat[[link]] = pivDat

  if (sum(res) == 0)
    attr(dat, "impossible") = TRUE

  dat
}


.transProbM = function(parent, childhap, mutmat = NULL) {
  # parent = list of vectors pat, mat
  # childhap = vector of any length (parental allele);
  # mutmat = mutation matrix
  # output: Vector of probs, length ncol(par)*child. See debug below
  np = length(parent$pat)
  nc = length(childhap)
  hap1 = rep(parent$pat, each = nc)
  hap2 = rep(parent$mat, each = nc)
  chHap = rep(childhap, np)

  if (is.null(mutmat))
    prob = ((hap1 == chHap) + (hap2 == chHap))/2
  else
    prob = (mutmat[cbind(hap1, chHap)] + mutmat[cbind(hap2, chHap)])/2

  # DEBUG:
  # print(matrix(prob, ncol = ncol(parent), dimnames = list(childhap, pasteHap(parent))))

  prob
}


#####################
## TWO LINKED MARKERS
#####################

#TODO: add mutations!

.peel_MM_AUT = function(dat, sub, rho) {
  fa = sub$father
  mo = sub$mother
  link = sub$link

  faDat = dat[[fa]]
  moDat = dat[[mo]]
  likel = tcrossprod(faDat$prob, moDat$prob) # faDat$prob %*% t.default(moDat$prob)
  faLen = nrow(likel)
  moLen = ncol(likel)

  # Loop over the children, except the link if this is a child.
  for (ch in .mysetdiff(sub$children, link)) {
    chDat = dat[[ch]]
    chprob = chDat$prob
    chLen = length(chprob)

    transPat = .transProbMM(faDat, chDat[c('pat1', 'pat2')], rho = rho)
    transMat = .transProbMM(moDat, chDat[c('mat1', 'mat2')], rho = rho)
    dim(transMat) = c(chLen, moLen)
    transMat_rep = transMat[rep(seq_len(chLen), faLen), ] # as.numeric(do.call(rbind, rep(list(transMat), faLen)))
    mm = .colSums((transPat * chprob) * transMat_rep, chLen, faLen * moLen)
    likel = likel * mm
  }

  if (link == 0)
    return(sum(likel))

  # Goal is to update the probabilities of the link individual ("pivot"):
  pivDat = dat[[link]]

  if (link == fa)
    res = .rowSums(likel, faLen, moLen)
  else if (link == mo)
    res = .colSums(likel, faLen, moLen)
  else { # link is a child
    pivp = pivDat$prob
    pivLen = length(pivp)

    TRarray = array(0, dim = c(faLen, moLen, pivLen))
    transPat = .transProbMM(faDat, pivDat[c('pat1', 'pat2')], rho = rho)
    transMat = .transProbMM(moDat, pivDat[c('mat1', 'mat2')], rho = rho)
    dim(transPat) = c(pivLen, faLen)
    dim(transMat) = c(pivLen, moLen)

    for (i in seq_len(faLen)) {
      transpat = transPat[, i]
      for (j in seq_len(moLen))
        TRarray[i, j, ] = transpat * transMat[, j]
    }

    arr = as.vector(TRarray) * as.vector(likel)
    dim(arr) = dim(TRarray)
    res = .colSums(arr, faLen * moLen, pivLen)  #sum for each entry of haps[[link]]
    res = res * pivp
  }

  # Update the probabilities
  pivDat$prob = res

  # Reduce if possible
  if(any(res == 0)) {
    pivDat$pat1 = pivDat$pat1[res > 0]
    pivDat$mat1 = pivDat$mat1[res > 0]
    pivDat$pat2 = pivDat$pat2[res > 0]
    pivDat$mat2 = pivDat$mat2[res > 0]
    pivDat$prob = pivDat$prob[res > 0]
  }

  dat[[link]] = pivDat

  if (sum(res) == 0)
    attr(dat, "impossible") = TRUE

  dat
}

.peel_MM_X = function(dat, sub, rho, SEX) {
  fa = sub$father
  mo = sub$mother
  link = sub$link

  faDat = dat[[fa]]
  moDat = dat[[mo]]
  likel = tcrossprod(faDat$prob, moDat$prob) # faDat$prob %*% t.default(moDat$prob)
  faLen = nrow(likel)
  moLen = ncol(likel)

  # Loop over the children, except the link if this is a child.
  for (ch in .mysetdiff(sub$children, link)) {
    chDat = dat[[ch]]
    chprob = chDat$prob
    chLen = length(chprob)

    if (SEX[ch] == 1) {
      transMat = .transProbMM(moDat, chDat[c('mat1', 'mat2')], rho = rho)
      mm = rep(.colSums(transMat * chprob, chLen, moLen), each = faLen)
    }
    else {
      # trans paternal (could make separate function for this)
      fa1 = faDat$mat1
      fa2 = faDat$mat2
      chp1 = chDat$pat1
      chp2 = chDat$pat2
      fa1rep = rep(fa1, each = length(chp1))
      fa2rep = rep(fa2, each = length(chp1))
      chp1rep = rep(chp1, length(fa1))
      chp2rep = rep(chp2, length(fa1))

      transPat = as.numeric((fa1rep == chp1rep) & (fa2rep == chp2rep))

      transMat = .transProbMM(moDat, chDat[c('mat1', 'mat2')], rho = rho)
      dim(transMat) = c(chLen, moLen)
      transMat_rep = transMat[rep(seq_len(chLen), faLen), ] # as.numeric(do.call(rbind, rep(list(transMat), faLen)))
      mm = .colSums((transPat * chprob) * transMat_rep, chLen, faLen * moLen)
    }

    likel = likel * mm
  }

  if (link == 0)
    return(sum(likel))

  # Goal is to update the probabilities of the link individual ("pivot"):
  pivDat = dat[[link]]

  if (link == fa)
    res = .rowSums(likel, faLen, moLen)
  else if (link == mo)
    res = .colSums(likel, faLen, moLen)
  else { # link is a child
    pivp = pivDat$prob
    pivLen = length(pivp)

    if (SEX[link] == 1) {
      transMat = .transProbMM(moDat, pivDat[c('mat1', 'mat2')], rho = rho)
      dim(transMat) = c(pivLen, moLen)
      TRarray = rep(t.default(transMat), each = faLen)
    }
    else {
      TRarray = array(0, dim = c(faLen, moLen, pivLen))
      transMat = .transProbMM(moDat, pivDat[c('mat1', 'mat2')], rho = rho)
      dim(transMat) = c(pivLen, moLen)
      for (i in seq_len(faLen)) {
        transPat = as.numeric(faDat$mat1[i] == pivDat$pat1 & faDat$mat2[i] == pivDat$pat2)
        TRarray[i, , ] = t.default(transPat * transMat)  #TODO:make faster!
      }
    }
    arr = as.vector(TRarray) * as.vector(likel)
    dim(arr) = dim(TRarray)
    res = .colSums(arr, faLen * moLen, pivLen)  #sum for each entry of haps[[link]]
    res * pivp
  }

  # Update the probabilites
  pivDat$prob = res

  # Reduce if possible
  if(any(res == 0)) {
    if(SEX[link] == 2) {
      pivDat$pat1 = pivDat$pat1[res > 0]
      pivDat$pat2 = pivDat$pat2[res > 0]
    }
    pivDat$mat1 = pivDat$mat1[res > 0]
    pivDat$mat2 = pivDat$mat2[res > 0]
    pivDat$prob = pivDat$prob[res > 0]
  }

  dat[[link]] = pivDat

  if (sum(res) == 0)
    attr(dat, "impossible") = TRUE

  dat
}


.transProbMM = function(par, gamete, rho) {
  # parent = list(pat1, mat1, pat2, mat2); vectors of same length
  # gamete = list(pat1, pat2) (or mat); vecs of same length
  lenPar = length(par[[1]])
  lenGam = length(gamete[[1]])

  p1 = rep(par$pat1, each = lenGam)
  m1 = rep(par$mat1, each = lenGam)
  p2 = rep(par$pat2, each = lenGam)
  m2 = rep(par$mat2, each = lenGam)

  gam1 = rep(gamete[[1]], lenPar)
  gam2 = rep(gamete[[2]], lenPar)

  res =
    (p1 == gam1 & p2 == gam2) * (1 - rho) +
    (m1 == gam1 & m2 == gam2) * (1 - rho) +
    (p1 == gam1 & m2 == gam2) * rho +
    (m1 == gam1 & p2 == gam2) * rho

  res/2
}
