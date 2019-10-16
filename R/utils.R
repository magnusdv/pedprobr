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
#' An autosomal marker with `n` alleles has `G = choose(n+1, 2)` possible
#' genotypes (when the allele order within a genotype does not matter). This
#' function returns these in the form of a `G * 2` matrix where each row is a
#' genotype. The computation is trivial; the main point is to provide a
#' consistent ordering of the genotypes. The homozygous genotypes come first,
#' followed by the heterozygous.
#'
#' @param n A positive integer.
#'
#' @examples
#' allGenotypes(3)
#'
#' @export
allGenotypes = function(n) {
  rbind(cbind(seq_len(n), seq_len(n)), .comb2(n))
}

# Debug tools: Paste genotypes given as 2*k matrices
pasteHap = function(hapmat) {
  if(is.matrix(hapmat)) {
    stopifnot(nrow(hapmat) == 2)
    return(paste(hapmat[1, ], hapmat[2, ], sep = "/"))
  }
  stopifnot(is.numeric(hapmat))
  as.character(hapmat)
}


#' @importFrom pedmut isStationary
hasStationaryModel = function(m) {
  mut = attr(m, 'mutmod')
  if(is.null(mut)) return(TRUE)

  sexEq = attr(mut, 'sexEqual')
  afr = afreq(m)

  isStationary(mut$male, afr) &&
    (isTRUE(sexEq) || isStationary(mut$female, afr))
}
