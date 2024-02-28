#### FUNCTIONS FOR CREATING THE INTITIAL HAPLOTYPE COMBINATIONS W/PROBABILITIES.

startdata_M = function(x, marker, pedInfo = NULL) {#print("new startdata")

  if(is.null(pedInfo))
    pedInfo = .pedInfo(x, Xchrom = isXmarker(marker))

  nInd = length(x$ID)
  SEX = x$SEX
  FIDX = x$FIDX
  MIDX = x$MIDX
  Xchrom = pedInfo$Xchrom
  isFounder = pedInfo$isFounder
  simpleFou = pedInfo$simpleFou
  fi = pedInfo$fouInb
  offs = pedInfo$offs

  marker = .sortGeno(marker)
  afr = attr(marker, "afreq")
  n = length(afr)

  # Complete set of genotypes - phased and unphased (for founders) - used below
  nseq = seq_len(n)
  allPhased   = list(pat = rep(nseq, each = n), mat = rep.int(nseq, times = n))
  allUnphased = list(pat = rep(nseq, n:1), mat = sequence.default(n:1, from = nseq))

  # Initialise output
  glist = vector(nInd, mode = "list")
  indOrder = seq_along(glist)
  imp = FALSE

  # Elimination: not for SNPs, and only if no mutation model
  eliminate = n > 2 && any(informative <- marker[,2] > 0) && is.null(attr(marker, "mutmod"))

  indOrder = if(eliminate) order(!informative) else seq_len(nInd)

  # Loop through all individuals
  for(i in indOrder) {
    a = marker[i, 1]
    b = marker[i, 2]

    Xmale = Xchrom && SEX[i] == 1 # Hemizygous male
    if(Xmale) {
      mat = if(b == 0) nseq else b
      if(isFounder[i])
        g = list(mat = mat, prob = afr[mat]) |> .reduce()
      else
        g = list(mat = mat, prob = rep(1, length(mat)))
    }
    else if(simpleFou[i]) # Simple founder: alleles directly (skip genotypes)
      g = .alleleDistrib(a, b, afr, f = fi[i])
    else if(isFounder[i]) # General founder: unphased genotypes with HW probs
      g = .genoDistribFounder(a, b, afr, f = fi[i], COMPLETE = allUnphased)
    else # Nonfounder: Phased genos with 1's as prob
      g = .genoDistribNonfounder(a, b, COMPLETE = allPhased)

    # Eliminate genotypes based on parents/children
    if(eliminate && length(g$prob) > 1) {   #.bef = length(g$prob)
      g = .elim(g, glist, FIDX[i], MIDX[i], offs[[i]], nall = n, sex = SEX[i], Xmale = Xmale)
      #cat(sprintf("Elim of '%s': %d -> %d\n", x$ID[i], .bef, length(g$prob)))
    }

    glist[[i]] = g

    # If all zero: impossible!
    if(!length(g$prob)) {imp = TRUE; break}
  }

  names(glist) = x$ID
  attr(glist, "impossible") = imp
  #class(glist) = c("glist", class(glist))
  glist
}


# print.glist = function(x, sep = "/") {
#   y = lapply(x, function(g) {
#     nms = g$allele %||% if(is.null(g$pat)) g$mat else paste(g$pat, g$mat, sep = sep)
#     `names<-`(g$prob, nms)
#   })
#   print(y)
# }

# Various info used repeatedly
.pedInfo = function(x, treatAsFounder = NULL, Xchrom = FALSE) {
  nInd = length(x$ID)

  # Founders (except LB-copies): Genotypes need not be phased
  fouInt = founders(x, internal = TRUE)
  isFounder = logical(nInd)
  isFounder[fouInt] = TRUE
  isFounder[treatAsFounder] = TRUE
  isFounder[x$LOOP_BREAKERS[, 2]] = FALSE

  # Founders with 1 child
  simpleFou = isFounder & tabulate(c(x$FIDX, x$MIDX), nInd) == 1

  # Founder inbreeding lookup vector; 0's for nonfounders
  fi = numeric(nInd)
  fval = x$FOUNDER_INBREEDING[[if(Xchrom) "x" else "autosomal"]]
  if(!is.null(fval))
    fi[fouInt] = fval

  # List of offspring (internal)
  offs = .allchildren(x)

  list(nInd = nInd, isFounder = isFounder, simpleFou = simpleFou, fouInb = fi,
       offs = offs, Xchrom = Xchrom)
}


.elim = function(g, glist, fa, mo, ch, nall, sex, Xmale = FALSE) {

  # Case 1: Simple founder
  if(!is.null(g$allele)) {
    for(b in ch) {
      bhap = glist[[b]][[if(sex == 1) "pat" else "mat"]]
      if(!is.null(bhap) && length(bhap) <= nall) {
        kp = g$allele %in% bhap
        g$allele = g$allele[kp]
        g$prob = g$prob[kp]
      }
    }
    return(g)
  }

  # Case 2: ordinary, i.e., with pat & mat
  fag = if(!Xmale && fa > 0) glist[[fa]] else NULL
  if(!is.null(fag) && length(fag$prob) <= nall) {
    kp = g$pat %in% c(fag$allele, fag$pat, fag$mat)
    g$pat = g$pat[kp]
    g$mat = g$mat[kp]
    g$prob = g$prob[kp]
  }
  mag = if(mo > 0) glist[[mo]] else NULL
  if(!is.null(mag) && length(mag$prob) <= nall) {
    kp = g$mat %in% c(mag$allele, mag$pat, mag$mat)
    g$pat = g$pat[kp]
    g$mat = g$mat[kp]
    g$prob = g$prob[kp]
  }
  for(b in ch) {
    bhap = glist[[b]][[if(sex == 1) "pat" else "mat"]]
    if(!is.null(bhap) && length(bhap) <= nall) {
      kp = if(Xmale) g$mat %in% bhap else g$pat %in% bhap | g$mat %in% bhap
      g$pat = g$pat[kp]
      g$mat = g$mat[kp]
      g$prob = g$prob[kp]
    }
  }
  g
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

  g = list(allele = allele, prob = prob)
  .reduce(g)
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
  .reduce(g)
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



##################################
# STARTDATA FOR TWO LINKED MARKERS
##################################


startdata_MM = function(x, marker1, marker2, pedInfo = NULL) {

  if(is.null(pedInfo))
    pedInfo = .pedInfo(x, Xchrom = isXmarker(marker1))

  glist1 = startdata_M(x, marker1, pedInfo = pedInfo)
  glist2 = startdata_M(x, marker2, pedInfo = pedInfo)

  if (attr(glist1, "impossible") || attr(glist2, "impossible"))
    return(structure(list(), impossible = TRUE))

  nInd = length(x$ID)
  SEX = x$SEX
  Xchrom = pedInfo$Xchrom
  isFounder = pedInfo$isFounder

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


.reduce = function(g) {
  keep = g$prob > 0
  if(all(keep))
    return(g)
  lapply(g, function(vec) vec[keep])
}

.allchildren = function(x) {
  FIDX = x$FIDX
  MIDX = x$MIDX
  SEX = x$SEX
  nInd = length(FIDX)

  nonf = which(FIDX > 0)

  res = split.default(c(nonf, nonf), factor(c(FIDX[nonf], MIDX[nonf]), levels = seq_len(nInd)))
  names(res) = NULL

  ### Slower in almost all cases:
  # res = lapply(seq_len(nInd), function(i)
  #  if(SEX[i] == 1) which(FIDX == i) else which(MIDX == i))

  res
}
