stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

# Test that input is a positive (or similar) integer.
isCount = function(x, minimum = 1) {
  isTRUE(length(x) == 1 &&
         (is.integer(x) || (is.numeric(x) && x == as.integer(x))) &&
         x >= minimum)
}

`%||%` = function(x, y) {
  if(is.null(x)) y else x
}

.mysetdiff = function(x, y) {
  unique.default(x[match(x, y, 0L) == 0L])
}

# Fast intersection. NB: assumes no duplicates!
.myintersect = function (x, y) {
  y[match(x, y, 0L)]
}

# Stripped version of expand.grid
fastGrid = function(argslist, as.list = FALSE) {
  nargs = length(argslist)
  orep = nr = prod(lengths(argslist))
  if (nargs == 0L || nr == 0L)
    return(matrix(ncol = 0, nrow = 0))

  rep.fac = 1L
  res = NULL
  for (x in argslist) {
    nx = length(x)
    orep = orep/nx
    res = c(res, x[rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep)])  #this is res[, i]
    rep.fac = rep.fac * nx
  }
  dim(res) = c(nr, nargs)
  if (as.list)
    res = lapply(seq_len(nr), function(r) res[r, ])
  res
}


log_or_not = function(x, logbase) {
  if (is.numeric(logbase))
    log(x, logbase)
  else x
}


# Equivalent to t.default(combn(n, 2)), but ~6 times faster.
.comb2 = function(n) {
  if (n < 2)
    return(matrix(nrow = 0, ncol = 2))
  v1 = rep.int(seq_len(n - 1), (n - 1):1)
  v2 = NULL
  for (i in 2:n) v2 = c(v2, i:n)
  cbind(v1, v2, deparse.level = 0)
}


#' Hardy-Weinberg probabilities
#'
#' @param allele1,allele2 Vectors of equal length, containing alleles in the
#'   form of indices of `afreq`
#' @param afreq A numeric vector with allele frequencies
#' @param f A single number in `[0, 1]`; the inbreeding coefficient
#'
#' @return A numeric vector of the same length as `allele1` and  `allele2`
#'
#' @examples
#' p = 0.1; q = 1-p
#' hw = HWprob(c(1,1,2), c(1,2,2), c(p, q))
#' stopifnot(all.equal(hw, c(p^2, 2*p*q, q^2)))
#'
#' @export
HWprob = function(allele1, allele2, afreq, f = 0) {
  homoz = allele1 == allele2
  hw = afreq[allele1] * afreq[allele2] * (2 - homoz)

  if(!is.na(f) && f > 0)
    hw = afreq[allele1] * homoz * f + hw * (1 - f)

  as.numeric(hw) # remove names; slightly faster than unname
}


#' Genotype matrix
#'
#' An autosomal marker with `n` alleles has `choose(n+1, 2)` possible
#' unordered genotypes. This function returns these as rows in a
#' matrix.
#'
#' @param n A positive integer.
#' @return An integer matrix with two columns and `choose(n+1, 2)` rows.
#'
#' @examples
#' allGenotypes(3)
#'
#' @export
allGenotypes = function(n) {
  # rbind(cbind(seq_len(n), seq_len(n)), .comb2(n))
  if (n < 1)
    return(matrix(integer(0), ncol = 2))
  nseq = seq_len(n)
  cbind(
    rep.int(nseq, times = n:1),
    unlist(lapply(nseq, function(i) i:n), recursive = FALSE, use.names = FALSE)
  )
}

# Debug tool, pretty-print "startdata" probs
pasteGenoProb = function(g) {
  p = g$prob
  names(p) = paste(g$pat, g$mat, sep="/")
  p
}


#' @importFrom pedmut isStationary sexEqual
hasStationaryModel = function(m) {
  mut = attr(m, 'mutmod')
  if(is.null(mut)) return(TRUE)

  afr = afreq(m)

  isStationary(mut$male, afr) &&
    (sexEqual(mut) || isStationary(mut$female, afr))
}


fixMerlinLog = function(a, logbase = NULL) {
  if(!is.null(logbase)) {
    if(length(logbase) != 1 || !is.numeric(logbase) || logbase <= 0)
      stop2("`logbase` must be a positive number: ", logbase)

    if(logbase == exp(1))
      res = a
    else
      res = round(a/log(logbase),3)
  }
  else {
    res = signif(exp(a), digits = 3)
    uflow = res == 0 & a > -Inf
    if(any(uflow))
      warning("Underflow!\nSome lnLikelihoods reported by MERLIN are too small to be exp'ed.",
              immediate. = TRUE, call. = FALSE)
    oflow = res == Inf & a < Inf
    if(any(oflow))
      warning("Overflow!\nSome lnLikelihoods reported by MERLIN are too large to be exp'ed.",
              immediate. = TRUE, call. = FALSE)
  }
  res
}


checkRho = function(rho, max = 0.5) {
  if(is.null(rho))
    stop2("Argument `rho` is missing")
  if(!is.numeric(rho))
    stop2("Argument `rho` should be a number, not: ", class(rho))
  if(length(rho) != 1)
    stop2("Argument `rho` must have length 1: ", rho)
  if(rho < 0)
    stop2("Argument `rho` cannot be negative: ", rho)
  if(rho > max)
    stop2(sprintf("Argument `rho` cannot exceed %f: %f", max, rho))
}

