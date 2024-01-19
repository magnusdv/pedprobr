#### PEELING FUNCTIONS

# Not currently used (?)
choosePeeler = function(twolocus, rho, Xchrom, SEX, mutmat, mutmat2 = NULL) {
  if(!twolocus && !Xchrom)
    .f = function(dat, sub) .peel_M_AUT(dat, sub, mutmat)
  else if(!twolocus && Xchrom)
    .f = function(dat, sub) .peel_M_X(dat, sub, SEX, mutmat)
  else if(twolocus && !Xchrom)
    .f = function(dat, sub) .peel_MM_AUT(dat, sub, rho, mutmat, mutmat2)
  else if(twolocus && Xchrom)
    .f = function(dat, sub) .peel_MM_X(dat, sub, rho, SEX, mutmat, mutmat2)
  .f
}

.peel_M_AUT = function(dat, sub, mutmat = NULL, newalg = F) {
  fa = sub$father
  mo = sub$mother
  link = sub$link

  faDat = dat[[fa]]
  moDat = dat[[mo]]
  likel = tcrossprod(faDat$prob, moDat$prob) # faDat$prob %*% t.default(moDat$prob)
  faLen = nrow(likel)
  moLen = ncol(likel)

  for (ch in .mysetdiff(sub$children, link)) {
    chDat = dat[[ch]]
    transPat = .transProbM(faDat, chDat$pat, mutmat = mutmat$male)
    transMat = .transProbM(moDat, chDat$mat, mutmat = mutmat$female)

    if(newalg) {
      mm = crossprod(chDat$prob * transPat, transMat)
    }
    else {
      chprob = chDat$prob
      chLen = length(chprob)
      dim(transPat) = NULL
      transMat_rep = transMat[rep(seq_len(chLen), faLen), ]
      mm = .colSums((transPat * chprob) * transMat_rep, chLen, faLen * moLen)
    }
    likel = likel * mm
  }

  if (link == 0)
    return(sum(likel))

  # Update the probabilities of link (pivot) individual:
  pivDat = dat[[link]]

  if (link == fa)
    res = .rowSums(likel, faLen, moLen)
  else if (link == mo)
    res = .colSums(likel, faLen, moLen)
  else { # link is a child
    pivLen = length(pivDat$prob)
    transPat = .transProbM(faDat, pivDat$pat, mutmat = mutmat$male)
    transMat = .transProbM(moDat, pivDat$mat, mutmat = mutmat$female)

    if(newalg) {
      faPr = .rowSums(likel, faLen, moLen)
      moPr = .colSums(likel, faLen, moLen)
      chProbPat = .colSums(faPr * t.default(transPat), faLen, pivLen)
      chProbMat = .colSums(moPr * t.default(transMat), moLen, pivLen)
      res = pivDat$prob * chProbPat * chProbMat / sum(likel)
    }
    else {
      TRarray = array(0, dim = c(faLen, moLen, pivLen))
      for (i in seq_len(faLen)) {
        transpat = transPat[, i]
        for (j in seq_len(moLen))
          TRarray[i, j, ] = transpat * transMat[, j]
      }
      arr = as.vector(TRarray) * as.vector(likel)
      dim(arr) = dim(TRarray)
      res = .colSums(arr, faLen * moLen, pivLen)  # sum for each entry of haps[[link]]
      res = res * pivDat$prob
    }

  }

  # Update the probabilities
  pivDat$prob = res

  # Reduce if possible
  if(any(res == 0)) {
    pivDat$allele = pivDat$allele[res > 0] # or NULL
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
          unlist(lapply(faDat$mat, function(fal) as.numeric(fal == chDat$pat)), recursive = FALSE, use.names = FALSE)
        else
          unlist(lapply(faDat$mat, function(fal) mutmat$male[fal, chDat$pat]), recursive = FALSE, use.names = FALSE)

      transMat = .transProbM(moDat, chDat$mat, mutmat = mutmat$female)
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
  else {
    pivp = pivDat$prob
    pivLen = length(pivp)

    if (SEX[link] == 1) {
     transMat = .transProbM(moDat, pivDat$mat, mutmat = mutmat$female)
     TRarray = rep(t.default(transMat), each = faLen)
    }
    else {
     TRarray = array(0, dim = c(faLen, moLen, pivLen))
     transMat = .transProbM(moDat, pivDat$mat, mutmat = mutmat$female)
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
  # parent = list of vectors, either both 'pat' and 'mat', or 'allele'
  # childhap = vector of any length
  # mutmat = mutation matrix
  # output: Vector of nc * np. See debug below
  np = length(parent$prob)
  nc = length(childhap)

  # New: If parent has `allele`, use this directly
  if(!is.null(pal <- parent$allele)) {
    if(nc > 1)
      pal = rep(pal, each = nc)
    res = if(is.null(mutmat)) as.numeric(pal == childhap) else mutmat[cbind(pal, childhap, deparse.level = 0)]
    dim(res) = c(nc, np)
    return(res)
  }

  pat = if(nc > 1) rep(parent$pat, each = nc) else parent$pat
  mat = if(nc > 1) rep(parent$mat, each = nc) else parent$mat

  if (is.null(mutmat))
    res = ((pat == childhap) + (mat == childhap))/2
  else
    res = (mutmat[cbind(pat, childhap)] + mutmat[cbind(mat, childhap)])/2

  # DEBUG:
  # print(matrix(prob, nrow = nc, dimnames = list(childhap, pasteHap(parent))))

  dim(res) = c(nc, np)
  res
}


#####################
## TWO LINKED MARKERS
#####################

.peel_MM_AUT = function(dat, sub, rho, mut1 = NULL, mut2 = NULL, newalg = F) {
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
    transPat = .transProbMM(faDat, chDat[c('pat1', 'pat2')], rho = rho, mutmat1 = mut1$male, mutmat2 = mut2$male)
    transMat = .transProbMM(moDat, chDat[c('mat1', 'mat2')], rho = rho, mutmat1 = mut1$female, mutmat2 = mut2$female)
    if(newalg) {
      mm = crossprod(chDat$prob * transPat, transMat)
    }
    else {
      dim(transPat) = NULL
      chprob = chDat$prob
      chLen = length(chprob)
      transMat_rep = transMat[rep(seq_len(chLen), faLen), ]
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
    pivLen = length(pivDat$prob)
    transPat = .transProbMM(faDat, pivDat[c('pat1', 'pat2')], rho = rho, mutmat1 = mut1$male, mutmat2 = mut2$male)
    transMat = .transProbMM(moDat, pivDat[c('mat1', 'mat2')], rho = rho, mutmat1 = mut1$female, mutmat2 = mut2$female)

    if(newalg) {
      faPr = .rowSums(likel, faLen, moLen)
      moPr = .colSums(likel, faLen, moLen)
      chProbPat = .colSums(faPr * t.default(transPat), faLen, pivLen)
      chProbMat = .colSums(moPr * t.default(transMat), moLen, pivLen)
      res = pivDat$prob * chProbPat * chProbMat / sum(likel)
    }
    else {
      TRarray = array(0, dim = c(faLen, moLen, pivLen))
      for (i in seq_len(faLen)) {
        transpat = transPat[, i]
        for (j in seq_len(moLen))
          TRarray[i, j, ] = transpat * transMat[, j]
      }

      arr = as.vector(TRarray) * as.vector(likel)
      dim(arr) = dim(TRarray)
      res = .colSums(arr, faLen * moLen, pivLen)  #sum for each entry of haps[[link]]
      res = res * pivDat$prob
    }
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
    if(sum(res) == 0)
      attr(dat, "impossible") = TRUE
  }

  dat[[link]] = pivDat
  dat
}

.peel_MM_X = function(dat, sub, rho, SEX, mut1 = NULL, mut2 = NULL) {
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
      transMat = .transProbMM(moDat, chDat[c('mat1', 'mat2')], rho = rho, mutmat1 = mut1$female, mutmat2 = mut2$female)
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

      loc1 = if(is.null(mut1)) as.numeric(fa1rep == chp1rep) else mut1$male[cbind(fa1rep, chp1rep, deparse.level = 0)]
      loc2 = if(is.null(mut2)) as.numeric(fa2rep == chp2rep) else mut2$male[cbind(fa2rep, chp2rep, deparse.level = 0)]
      transPat = loc1 * loc2

      transMat = .transProbMM(moDat, chDat[c('mat1', 'mat2')], rho = rho, mutmat1 = mut1$female, mutmat2 = mut2$female)
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
      transMat = .transProbMM(moDat, pivDat[c('mat1', 'mat2')], rho = rho, mutmat1 = mut1$female, mutmat2 = mut2$female)
      TRarray = rep(t.default(transMat), each = faLen)
    }
    else {
      TRarray = array(0, dim = c(faLen, moLen, pivLen))
      transMat = .transProbMM(moDat, pivDat[c('mat1', 'mat2')], rho = rho, mutmat1 = mut1$female, mutmat2 = mut2$female)
      for (i in seq_len(faLen)) {
        loc1 = if(is.null(mut1)) as.numeric(faDat$mat1[i] == pivDat$pat1) else mut1$male[faDat$mat1[i], pivDat$pat1]
        loc2 = if(is.null(mut2)) as.numeric(faDat$mat2[i] == pivDat$pat2) else mut2$male[faDat$mat2[i], pivDat$pat2]
        transPat = loc1 * loc2
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


.transProbMM = function(par, gamete, rho, mutmat1 = NULL, mutmat2 = NULL, debug = FALSE) {
  # parent = list(pat1, mat1, pat2, mat2); vectors of same length
  # gamete = list(pat1, pat2) (or mat); vecs of same length
  gam1 = gamete[[1]]
  gam2 = gamete[[2]]
  np = length(par$prob)
  nc = length(gam1)

  # If parent has `allele1`, `allele2`, use this directly
  if(!is.null(par$allele1)) {
    pal1 = if(nc > 1) rep(par$allele1, each = nc) else par$allele1
    pal2 = if(nc > 1) rep(par$allele2, each = nc) else par$allele2
    prob1 = if(is.null(mutmat1)) as.numeric(pal1 == gam1) else mutmat1[cbind(pal1, gam1, deparse.level = 0)]
    prob2 = if(is.null(mutmat2)) as.numeric(pal2 == gam2) else mutmat2[cbind(pal2, gam2, deparse.level = 0)]
    res = prob1 * prob2
    dim(res) = c(nc, np)
    return(res)
  }

  # Suffixes below refer to locus (1 or 2)
  p1 = if(nc > 1) rep(par$pat1, each = nc) else par$pat1
  m1 = if(nc > 1) rep(par$mat1, each = nc) else par$mat1
  p2 = if(nc > 1) rep(par$pat2, each = nc) else par$pat2
  m2 = if(nc > 1) rep(par$mat2, each = nc) else par$mat2

  loc1pat = if(is.null(mutmat1)) as.integer(p1 == gam1) else mutmat1[cbind(p1, gam1, deparse.level = 0)]
  loc1mat = if(is.null(mutmat1)) as.integer(m1 == gam1) else mutmat1[cbind(m1, gam1, deparse.level = 0)]
  loc2pat = if(is.null(mutmat2)) as.integer(p2 == gam2) else mutmat2[cbind(p2, gam2, deparse.level = 0)]
  loc2mat = if(is.null(mutmat2)) as.integer(m2 == gam2) else mutmat2[cbind(m2, gam2, deparse.level = 0)]

  res = (loc1pat * loc2pat * (1 - rho) + loc1mat * loc2mat * (1 - rho) + loc1pat * loc2mat * rho + loc1mat * loc2pat * rho)/2

  # DEBUG:
  if(debug) {
    parHap = paste(pasteHap(par[c(1,3)], ""), pasteHap(par[c(2,4)], ""), sep = "|")
    print(matrix(res, ncol = np, dimnames = list(pasteHap(gamete, ""), parHap)))
  }

  dim(res) = c(nc, np)
  res
}

pasteHap = function(lst, sep = "/") paste(lst[[1]], lst[[2]], sep = sep)
