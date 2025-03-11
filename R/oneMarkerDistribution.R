#' Genotype distribution for a single marker
#'
#' Computes the genotype probability distribution of one or several pedigree
#' members, possibly conditional on known genotypes for the marker.
#'
#' @param x A `ped` object or a list of such.
#' @param ids A vector of ID labels of one or more members of `x`.
#' @param marker Either a `marker` object or the name (or index) of a marker
#'   attached to `x`. If `x` has multiple components, only the latter is
#'   allowed.
#' @param partialmarker (Deprecated) An alias for `marker`.
#' @param loopBreakers (Only relevant if the pedigree has loops). A vector with
#'   ID labels of individuals to be used as loop breakers. If NULL (default)
#'   loop breakers are selected automatically. See [breakLoops()].
#' @param grid.subset (Optional; not relevant for most users.) A numeric matrix
#'   describing a subset of all marker genotype combinations for the `ids`
#'   individuals. The matrix should have one column for each of the `ids`
#'   individuals, and one row for each combination: The genotypes are described
#'   in terms of the matrix `M = allGenotypes(n)`, where `n` is the number of
#'   alleles for the marker. If the entry in column `j` is the integer `k`, this
#'   means that the genotype of individual `ids[j]` is row `k` of `M`.
#' @param output A character string, either `"array"` (default), "table" or
#'   "sparse". See Value.
#' @param verbose A logical.
#'
#' @return The output format depends on the `output` argument:
#'
#' * "array": A named `k`-dimensional array, where `k = length(ids)`, with the
#'   joint genotype distribution for the `ids` individuals, conditional on the
#'   known genotypes if present.
#' * "table": A data frame with `k+1` columns, where each row corresponds
#'   to a genotype combination, and the last column `prob` gives the
#'   probability.
#' * "sparse": A data frame with the same structure as the "table" output,
#'   but only combinations with non-zero probability are included.
#'
#' @seealso [twoMarkerDistribution()]
#'
#' @examples
#'
#' # Trivial example: Hardy-Weinberg probabilities for an equifrequent SNP
#' s = singleton(id = 1) |> addMarker(alleles = 1:2, afreq = c(0.5, 0.5))
#' oneMarkerDistribution(s, ids = 1)
#'
#' # Conditioning on a partial genotype
#' s = setGenotype(s, ids = 1, geno = "1/-")
#' oneMarkerDistribution(s, ids = 1)
#'
#' # Genotype distribution for a child of heterozygous parents
#' trio = nuclearPed(father = "fa", mother = "mo", child = "ch") |>
#'   addMarker(fa = "1/2", mo = "1/2")
#' oneMarkerDistribution(trio, ids = "ch")
#'
#' # Joint distribution of the parents, given that the child is heterozygous
#' trio = addMarker(trio, ch = "1/2")
#' ids = c("fa", "mo")
#' oneMarkerDistribution(trio, ids = ids, marker = 2)
#'
#' # Table output of the previous example
#' oneMarkerDistribution(trio, ids = ids, marker = 2, output = "table")
#' oneMarkerDistribution(trio, ids = ids, marker = 2, output = "sparse")
#'
#' # A different example: The genotype distribution of an individual (id = 8)
#' # whose half cousin (id = 9) is homozygous for a rare allele.
#' y = halfCousinPed(degree = 1) |>
#'   addMarker("9" = "a/a", afreq = c(a = 0.01, b = 0.99))
#'
#' oneMarkerDistribution(y, ids = 8)
#'
#' # Multi-component (trivial) example
#' z = singletons(1:2) |> addMarker(`1` = "1/2", `2` = "1/2", alleles = 1:2)
#' oneMarkerDistribution(z, 1:2)
#' oneMarkerDistribution(z, 1:2, output = "sparse")
#'
#' @export
oneMarkerDistribution = function(x, ids, marker = 1, loopBreakers = NULL,
                                 grid.subset = NULL, partialmarker = NULL,
                                 output = c("array", "table", "sparse"),
                                 verbose = TRUE) {

  if(!is.null(partialmarker)) {
    cat("The argument `partialmarker` has been renamed to `marker` and will be removed in a future version.\n")
    marker = partialmarker
  }

  ids = as.character(ids)

  output = match.arg(output)
  formatResult = function(res) {
    switch(output, array = res,
           table = omd2df(res, ids, sparse = FALSE),
           sparse = omd2df(res, ids, sparse = TRUE))
  }

  if(is.pedList(x)) {
    if(is.marker(marker))
      stop2("When `x` has multiple components, `marker` cannot be an unattached marker object")

    pednr = getComponent(x, ids, checkUnique = TRUE, errorIfUnknown = TRUE)
    if(all(pednr == pednr[1]))
      x = x[[pednr[1]]]
    else {
      compRes = lapply(unique.default(pednr), function(i) {
        idsC = ids[pednr == i]
        lb = if(is.null(loopBreakers)) NULL else .myintersect(loopBreakers, x[[i]]$ID)
        gs = if(is.null(grid.subset)) NULL else unique.matrix(grid.subset[, pednr == i, drop = FALSE])
        oneMarkerDistribution(x[[i]], idsC, marker = marker, loopBreakers = lb, grid.subset = gs, verbose = FALSE)
      })
      res = Reduce(`%o%`, compRes)
      return(formatResult(res))
    }
  }

  if(!is.ped(x))
    stop2("Input is not a pedigree")

  m = marker

  if (!is.marker(m)) {
    if(length(m) != 1)
      stop2("`marker` must have length 1")
    m = getMarkers(x, markers = m)[[1]]
  }

  # Special case: Unconditional & all founders
  # TODO: Generalisation to unrelated clusters - requires ribd!
  if(length(ids) > 1 && all(m == 0) & all(ids %in% founders(x))) {
    compRes = lapply(seq_along(ids), function(i) {
      gs = if(is.null(grid.subset)) NULL else unique.matrix(grid.subset[, i, drop = FALSE])
      oneMarkerDistribution(x, ids[i], marker = m, loopBreakers = loopBreakers,
                            grid.subset = gs, output = "array", verbose = FALSE)
    })
    res = Reduce(`%o%`, compRes)
    return(formatResult(res))
  }

  # Reorder if needed
   if(!hasParentsBeforeChildren(x)) {
    x = parentsBeforeChildren(setMarkers(x, m))
    m = x$MARKERS[[1]]
   }

  if (!is.null(x$LOOP_BREAKERS))
    stop2("`ped` objects with pre-broken loops are not allowed as input to `oneMarkerDistribution()`")

  alleles = alleles(m)
  Xchrom = isXmarker(m)

  if (verbose) {
    cat("Known genotypes:\n")
    print(m)
    cat("\nChromosome type    :", ifelse(Xchrom, "X-linked", "autosomal"))
    cat("\nTarget individuals :", toString(ids), "\n")
  }

  st = Sys.time()

  # Compute grid before loop breaking (works better with eliminate2)
  if (is.null(grid.subset))
    grid.subset = genoCombinations(x, m, ids, make.grid = TRUE)
  else
    grid.subset = as.matrix(grid.subset)

  if (x$UNBROKEN_LOOPS) {
    x = breakLoops(setMarkers(x, m), loopBreakers = loopBreakers, verbose = verbose)
    m = x$MARKERS[[1]]
  }

  int.ids = internalID(x, ids)
  allgenos = allGenotypes(nAlleles(m))

  # Character with genotype labels
  gt.strings = paste(alleles[allgenos[, 1]], alleles[allgenos[, 2]], sep = "/")
  if(Xchrom) {
    sx = getSex(x, ids)
    geno.names =  list(alleles, gt.strings)[sx]
  }
  else {
    geno.names = rep(list(gt.strings), length(ids))
  }

  # Create output array. Will hold likelihood of each genotype combo
  probs = array(0, dim = lengths(geno.names, use.names = FALSE), dimnames = geno.names)

  # Subset of `probs` that is affected by grid.subset
  probs.subset = grid.subset

  # Needs adjustment for X (in male columns)
  if(Xchrom) {
    homoz = which(allgenos[,1] == allgenos[,2])
    probs.subset[, sx == 1] = match(probs.subset[, sx == 1], homoz)
  }

  nr = nrow(grid.subset)

  # Compute marginal
  marginal = likelihood(x, markers = m)
  if (marginal == 0)
      stop2("Partial marker is impossible")

  if(verbose) {
    cat("Marginal likelihood:", marginal, "\n")
    cat("Calculations needed:", nr, "\n")
  }

  # Create list of all markers
  mlist = lapply(seq_len(nr), function(i) {
    r = grid.subset[i,]
    m[int.ids, ] = allgenos[r, ]; m})

  # Calculate likelihoods and insert in result array
  probs[probs.subset] = likelihood(x, mlist, allX = Xchrom)

  # Timing
  if(verbose)
    cat("\nAnalysis finished in", format(Sys.time() - st, digits = 3), "\n")

  res = probs/marginal
  formatResult(res)
}


# Convert output array to data frame
omd2df = function(arr, ids, sparse = FALSE) {
  dn = dimnames(arr)
  args = c(dn, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  df = do.call(expand.grid, args)
  names(df) = ids
  df$prob = as.vector(arr)

  if(sparse) {
    df = df[df$prob > 0, , drop = FALSE]
    rownames(df) = NULL
  }

  df
}
