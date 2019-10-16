#### PEELING FUNCTIONS

choosePeeler = function(twolocus, theta, Xchrom, SEX, mutmat) {
  if(!twolocus && !Xchrom)
    .f = function(dat, sub) .peel_M_AUT(dat, sub, mutmat)
  else if(!twolocus && Xchrom)
    .f = function(dat, sub) .peel_M_X(dat, sub, SEX, mutmat)
  else if(twolocus && !Xchrom)
    .f = function(dat, sub) .peel_MM_AUT(dat, sub, theta)
  else if(twolocus && Xchrom)
    .f = function(dat, sub) .peel_MM_X(dat, sub, theta, SEX)
  .f
}

.peel_M_AUT = function(dat, sub, mutmat = NULL) {
  fa = sub$father
  mo = sub$mother
  link = sub$link

  farh = dat[[c(fa, 1)]]
  morh = dat[[c(mo, 1)]]
  likel = dat[[c(fa, 2)]] %*% t.default(dat[[c(mo, 2)]])
  fa_len = nrow(likel)
  mo_len = ncol(likel)

  # Loop over the children, except the link if this is a child.
  for (b in .mysetdiff(sub$children, link)) {
    bh = dat[[c(b,1)]]
    bp = dat[[c(b,2)]]
    bl = length(bp)
    trans_pats = .trans_M(farh, bh[1, ], mutmat = mutmat$male)
    trans_mats = .trans_M(morh, bh[2, ], mutmat = mutmat$female)
    dim(trans_mats) = c(bl, mo_len)
    trans_mats_rep = as.numeric(do.call(rbind, rep(list(trans_mats), fa_len)))
    mm = .colSums((trans_pats * bp) * trans_mats_rep, bl, fa_len * mo_len)
    likel = likel * mm
  }

  if (link == 0) return(sum(likel))

  if (link == fa)
    res = .rowSums(likel, fa_len, mo_len)
  else if (link == mo)
    res = .colSums(likel, fa_len, mo_len)
  else { # link is a child
    pivh = dat[[c(link, 1)]]
    pivp = dat[[c(link, 2)]]
    pi_len = length(pivp)

    TRarray = array(0, dim = c(fa_len, mo_len, pi_len))
    trans_pats = .trans_M(farh, pivh[1, ], mutmat = mutmat$male)
    dim(trans_pats) = c(pi_len, fa_len)
    trans_mats = .trans_M(morh, pivh[2, ], mutmat = mutmat$female)
    dim(trans_mats) = c(pi_len, mo_len)
    for (i in seq_len(fa_len)) {
      transpat = trans_pats[, i]
      for (j in seq_len(mo_len))
        TRarray[i, j, ] = transpat * trans_mats[, j]
    }

    # DEBUG
    # print(array(TRarray, dim = c(fa_len, mo_len, pi_len), dimnames = lapply(list(farh, morh, pivh), pasteHap)))

    arr = as.vector(TRarray) * as.vector(likel)
    dim(arr) = dim(TRarray)
    res = .colSums(arr, fa_len * mo_len, pi_len)  # sum for each entry of haps[[link]]
    res = res * pivp
  }

  if (sum(res) == 0)
    attr(dat, "impossible") = TRUE

  pivhapUpdate = dat[[c(link, 1)]][, res > 0, drop = F]
  dat[[link]] = list(hap = pivhapUpdate, prob = res[res > 0])
  dat
}


.peel_M_X = function(dat, sub, SEX, mutmat = NULL) {
  fa = sub$father
  mo = sub$mother
  link = sub$link

  farh = dat[[c(fa, 1)]]
  morh = dat[[c(mo, 1)]]
  likel = dat[[c(fa, 2)]] %*% t.default(dat[[c(mo, 2)]])
  dims = dim(likel)
  fa_len = dims[1L]
  mo_len = dims[2L]

  # Loop over the children, except the link if this is a child.
  for (b in .mysetdiff(sub$children, link)) {
    bh = dat[[c(b,1)]]
    bp = dat[[c(b,2)]]
    bl = length(bp)

    if (SEX[b] == 1) {
      trans_mats = .trans_M(morh, bh, mutmat = mutmat$female)
      mm = rep(.colSums(trans_mats * bp, bl, mo_len), each = fa_len)
    }
    else {
      trans_pats =
        if (is.null(mutmat)) unlist(lapply(farh, function(fh) as.numeric(fh == bh[1, ])))
        else unlist(lapply(farh, function(fh) mutmat$male[fh, bh[1, ]]))

      trans_mats = .trans_M(morh, bh[2, ], mutmat = mutmat$female)
      dim(trans_mats) = c(bl, mo_len)
      trans_mats_rep = as.numeric(do.call(rbind, rep(list(trans_mats), fa_len)))
      mm = .colSums((trans_pats * bp) * trans_mats_rep, bl, fa_len * mo_len)
    }
    likel = likel * mm
  }


  if (link == 0)
    return(sum(likel))

  if (link == fa)
    res = .rowSums(likel, fa_len, mo_len)
  else if (link == mo)
    res = .colSums(likel, fa_len, mo_len)
  else { # link is a child
    pivh = dat[[c(link, 1)]]
    pivp = dat[[c(link, 2)]]
    pi_len = length(pivp)

    if (SEX[link] == 1) {
     trans_mats = .trans_M(morh, pivh, mutmat = mutmat$female)
     dim(trans_mats) = c(pi_len, mo_len)
     TRarray = rep(t.default(trans_mats), each = fa_len)
    }
    else {
     TRarray = array(0, dim = c(fa_len, mo_len, pi_len))
     trans_mats = .trans_M(morh, pivh[2, ], mutmat = mutmat$female)
     dim(trans_mats) = c(pi_len, mo_len)
     for (i in seq_len(fa_len)) {
       trans_pats =
         if (is.null(mutmat)) as.numeric(farh[i] == pivh[1, ])
       else mutmat$male[farh[i], pivh[1, ]]
       TRarray[i, , ] = t.default(trans_pats * trans_mats)  #TODO:make faster?
     }
    }

    # DEBUG
    # print(array(TRarray, dim = c(fa_len, mo_len, pi_len), dimnames = lapply(list(farh, morh, pivh), pasteHap)))

    arr = as.vector(TRarray) * as.vector(likel)
    dim(arr) = dim(TRarray)
    res = .colSums(arr, fa_len * mo_len, pi_len)  #sum for each entry of haps[[link]]
    res = res * pivp
  }

  if (sum(res) == 0) attr(dat, "impossible") = TRUE

  if (SEX[link] == 1)
    pivhapUpdate = dat[[c(link, 1)]][res > 0]
  else
    pivhapUpdate = dat[[c(link, 1)]][, res > 0, drop = F]

  dat[[link]] = list(hap = pivhapUpdate, prob = res[res > 0])
  dat
}


.trans_M = function(parent, childhap, mutmat = NULL) {
  # parent = matrix with 2 rows; each column a genotype
  # childhap = vector of any length (parental allele);
  # mutmat = mutation matrix
  # output: Vector of probs, length ncol(par)*child. See debug below
  sq = seq_len(ncol(parent))
  if (is.null(mutmat))
    prob = unlist(lapply(sq, function(i)
      ((parent[1, i] == childhap) + (parent[2, i] == childhap))/2))
  else
    prob = unlist(lapply(sq, function(i)
      (mutmat[parent[1, i], childhap] + mutmat[parent[2, i], childhap])/2))

  # DEBUG:
  # print(matrix(prob, ncol = ncol(parent), dimnames = list(childhap, pasteHap(parent))))

  prob
}

#####################
## TWO LINKED MARKERS
#####################

#TODO: add mutations!

.peel_MM_AUT = function(dat, sub, theta) {
  fa = sub$father
  mo = sub$mother
  link = sub$link

  farh = dat[[c(fa, 1)]]
  morh = dat[[c(mo, 1)]]
  likel = dat[[c(fa, 2)]] %*% t.default(dat[[c(mo, 2)]])
  fa_len = nrow(likel)
  mo_len = ncol(likel)

  mm_init = numeric(length(likel))
  dim(mm_init) = dim(likel)

  # Loop over the children, except the link if this is a child.
  for (b in .mysetdiff(sub$children, link)) {
    mm = mm_init
    bh = dat[[c(b, 1)]]
    bp = dat[[c(b, 2)]]
    for (i in seq_len(fa_len)) {
      transfather = .trans_MM(farh[, i], bh[c(1, 3), , drop = F], theta)
      for (j in seq_len(mo_len))
        mm[i, j] = (transfather * .trans_MM(morh[, j], bh[c(2, 4), , drop = F], theta)) %*% bp
    }
    likel = likel * mm
  }

  if (link == 0)
    return(sum(likel))

  if (link == fa)
    res = .rowSums(likel, fa_len, mo_len)
  else if (link == mo)
    res = .colSums(likel, fa_len, mo_len)
  else { # link is a child
    pivh = dat[[c(link, 1)]]
    pivp = dat[[c(link, 2)]]
    pi_len = length(pivp)

    TRarray = array(0, dim = c(fa_len, mo_len, pi_len))
    for (i in seq_len(fa_len)) {
      transfather = .trans_MM(farh[, i], pivh[c(1, 3), , drop = F], theta)
      for (j in seq_len(mo_len))
        TRarray[i, j, ] = transfather * .trans_MM(morh[, j], pivh[c(2, 4), , drop = F], theta)
    }
    arr = as.vector(TRarray) * as.vector(likel)
    dim(arr) = dim(TRarray)
    res = .colSums(arr, fa_len * mo_len, pi_len)  #sum for each entry of haps[[link]]
    res = res * pivp
  }

  if (sum(res) == 0) attr(dat, "impossible") = TRUE

  pivhapUpdate = dat[[c(link, 1)]][, res > 0, drop = F]
  dat[[link]] = list(hap = pivhapUpdate, prob = res[res > 0])
  dat
}

.peel_MM_X = function(dat, sub, theta, SEX) {
  fa = sub$father
  mo = sub$mother
  link = sub$link

  farh = dat[[c(fa, 1)]]
  morh = dat[[c(mo, 1)]]
  likel = dat[[c(fa, 2)]] %*% t.default(dat[[c(mo, 2)]])
  fa_len = nrow(likel)
  mo_len = ncol(likel)

  mm_init = numeric(length(likel))
  dim(mm_init) = dim(likel)

  # Loop over the children, except the link if this is a child.
  for (b in .mysetdiff(sub$children, link)) {
    mm = mm_init
    bh = dat[[c(b, 1)]]
    bp = dat[[c(b, 2)]]
    if (SEX[b] == 1) {
      for (j in seq_len(mo_len))
        mm[, j] = .trans_MM(morh[, j], bh, theta) %*% bp
    }
    else {
      for (j in seq_len(mo_len)) {
        transmother = .trans_MM(morh[, j], bh[c(2, 4), , drop = F], theta)
        for (i in seq_len(fa_len))
          mm[i, j] = (as.numeric(farh[1, i] == bh[1, ] & farh[2, i] == bh[3, ]) * transmother) %*% bp
      }
    }
    likel = likel * mm
  }

  if (link == 0)
    return(sum(likel))

  if (link == fa)
    res = .rowSums(likel, fa_len, mo_len)
  else if (link == mo)
    res = .colSums(likel, fa_len, mo_len)
  else { # link is a child
    pivh = dat[[c(link, 1)]]
    pivp = dat[[c(link, 2)]]
    pi_len = length(pivp)

    TRarray = array(0, dim = c(fa_len, mo_len, pi_len))
    if (SEX[link] == 1) {
      for (j in seq_len(mo_len)) {
        transmother = .trans_MM(morh[, j], pivh, theta)
        for (i in seq_len(fa_len))
          TRarray[i, j, ] = transmother
      }
    } else {
      for (j in seq_len(mo_len)) {
        transmother = .trans_MM(morh[, j], pivh[c(2, 4), ], theta)
        for (i in seq_len(fa_len))
          TRarray[i, j, ] = as.numeric(farh[1, i] == pivh[1, ] & farh[2, i] == pivh[3, ]) * transmother
      }
    }
    arr = as.vector(TRarray) * as.vector(likel)
    dim(arr) = dim(TRarray)
    res = .colSums(arr, fa_len * mo_len, pi_len)  #sum for each entry of haps[[link]]
    res * pivp
  }

  if (sum(res) == 0)  attr(dat, "impossible") = TRUE

  pivhapUpdate = dat[[c(link, 1)]][, res > 0, drop = F]
  dat[[link]] = list(hap = pivhapUpdate, prob = res[res > 0])
  dat
}


.trans_MM <- function(parent.haps, gamete.hap, theta) {
  # parent.haps = c(M1_1, M1_2, M2_1, M2_2)
  if (is.matrix(gamete.hap))
    vapply(seq_len(ncol(gamete.hap)), function(kol)
      .trans_MM(parent.haps, gamete.hap[, kol], theta), 1)
  else sum(c(
    all(parent.haps[c(1, 3)] == gamete.hap),
    all(parent.haps[c(2, 4)] == gamete.hap),
    all(parent.haps[c(1, 4)] == gamete.hap),
    all(parent.haps[c(2, 3)] == gamete.hap)) *
      c(1 - theta, 1 - theta, theta, theta))/2
}
