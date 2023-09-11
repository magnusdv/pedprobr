#### .buildGenolist and .eliminate

.genotypeHaploList = function(gt, n, unordered, complete = NULL, dropoutTF = FALSE) {
  nseq = seq_len(n)

  # The complete matrix can (should!) be supplied to avoid making it each time
  if(is.null(complete))
    complete = list(pat = rep(nseq, each = n), mat = rep.int(nseq, times = n))

  a = gt[1]
  b = gt[2]

  if(a == 0 && b == 0)
    g = complete
  else if (a == 0 && b > 0)
    g = list(pat = c(nseq, rep(b, n - 1)), mat = c(rep(b, n), nseq[-b]))
  else if (a > 0 && b == 0)
    g = list(pat = c(nseq, rep(a, n - 1)), mat = c(rep(a, n), nseq[-a]))
  else if (a == b && dropoutTF)
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

.buildGenolist = function(x, marker, eliminate = 0, treatAsFounder = NULL,
                          foundersUnordered = TRUE, dropoutTF = NULL) {
  nInd = length(x$ID)
  nall = nAlleles(marker)

  # Founders (except loop breaker copies) need only *unordered* genotypes
  unordered = logical(nInd)
  if(foundersUnordered) {
    unordered[founders(x, internal = TRUE)] = TRUE
    unordered[treatAsFounder] = TRUE
    unordered[x$LOOP_BREAKERS[, 2]] = FALSE
  }

  # A matrix containing a complete set of ordered genotypes
  nseq = seq_len(nall)
  COMPLETE = list(pat = rep(nseq, each = nall),
                  mat = rep.int(nseq, times = nall))

  if(is.null(dropoutTF))
    dropoutTF = logical(nInd)

  # Building a list of genotypes for each indiv.
  genolist = lapply(1:nInd, function(i)
    .genotypeHaploList(marker[i, ], nall, unordered[i], COMPLETE, dropoutTF = dropoutTF[i]))

  #.printG(genolist)

  attr(genolist, "impossible") = FALSE

  # If mutations, don't eliminate any genotypes
  if (allowsMutations(marker))
    return(genolist)

  .eliminate(x, genolist, nall, repeats = eliminate, treatAsFounder = treatAsFounder)
}



.eliminate = function(x, genolist, nall, repeats = 0, treatAsFounder = NULL) {
  if (repeats == 0 || attr(genolist, "impossible"))
    return(genolist)
  N = pedsize(x)

  FOU = founders(x, internal = TRUE)
  NONFOU = nonfounders(x, internal = TRUE)
  FIDX = x$FIDX
  MIDX = x$MIDX
  SEX = x$SEX

  # Adjust for "extra" founders
  if(length(treatAsFounder) > 0) {
    FOU = c(FOU, treatAsFounder)
    NONFOU = .mysetdiff(NONFOU, treatAsFounder)
  }

  offs = lapply(1:N, function(i) children(x, i, internal = TRUE))

  # Current number of genotypes of each
  Ngt = unlist(lapply(genolist, function(g) length(g$pat)),
               recursive = FALSE, use.names = FALSE)

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
      sx = SEX[i]

      if (i %in% NONFOU) {
        if(informative[fa <- FIDX[i]]) {
          kp = g$pat %in% genolist[[fa]]$pat | g$pat %in% genolist[[fa]]$mat
          if(!all(kp)) {
            g$pat = g$pat[kp]
            g$mat = g$mat[kp]
          }
        }
        if(informative[mo <- MIDX[i]]) {
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

    Ngt = unlist(lapply(genolist, function(g) length(g$pat)),
                 recursive = FALSE, use.names = FALSE)
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

  Ngt = unlist(lapply(genolist, function(g) length(g$mat)),
               recursive = FALSE, use.names = FALSE)

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
        if(informative[fa <- FIDX[i]]) {
          kp = g$pat %in% genolist[[fa]]$mat
          if(!all(kp)) {
            g$pat = g$pat[kp]
            g$mat = g$mat[kp]
          }
        }
        if(informative[mo <- MIDX[i]]) {
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

    Ngt = unlist(lapply(genolist, function(g) length(g$mat)),
                 recursive = FALSE, use.names = FALSE)

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
