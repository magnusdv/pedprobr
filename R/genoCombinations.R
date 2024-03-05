#' Genotype combinations
#'
#' Returns the possible genotype combinations in a pedigree, given partial
#' marker data. This function is mainly for internal use.
#'
#' @param x a [ped()] object.
#' @param partialmarker a [marker()] object compatible with `x`.
#' @param ids a vector with ID labels of one or more pedigree members.
#' @param make.grid a logical indicating if the result should be simplified to a
#'   matrix.
#'
#' @return If `make.grid = FALSE` (the default) the function returns a list of
#'   integer vectors, one vector for each element of `ids`. Each integer
#'   represents a genotype, in the form of a row number of the matrix
#'   `allGenotypes(n)`, where `n` is the number of alleles of the marker.
#'
#'   If `make.grid = TRUE`, the Cartesian product of the vectors is taken,
#'   resulting in a matrix with one column for each element of `ids`.
#'
#' @export
genoCombinations = function(x, partialmarker = x$MARKERS[[1]], ids, make.grid = TRUE) {
  ids = as.character(ids)
  int.ids = internalID(x, ids)

  allg = allGenotypes(nAlleles(partialmarker))
  homoz = which(allg[,1] == allg[,2])
  allgRef = 1000 * (allg[, 1] + allg[, 2]) + abs(allg[, 1] - allg[, 2])

  matchRefRows = function(g) { # get corresponding row numbers of 'allg'
    sort.int(unique.default(match(1000 * (g$pat + g$mat) + abs(g$pat - g$mat), allgRef)))
  }

  Xchrom = isXmarker(partialmarker)

  if (Xchrom) {
    SEX = x$SEX
    glist = .buildGenolistX(x, partialmarker, eliminate = 10)
    if (attr(glist, "impossible"))
      stop2("Impossible partial marker")
    rows = lapply(int.ids, function(i)
      switch(SEX[i],
        homoz[glist[[i]]$mat],
        matchRefRows(glist[[i]])))
  }
  else {
    glist = .buildGenolist(x, partialmarker, eliminate = 10)
    if (attr(glist, "impossible"))
      stop2("Impossible partial marker")
    rows = lapply(glist[int.ids], matchRefRows)
  }

  if(!make.grid)
    return(rows)

  if(allowsMutations(partialmarker)) {
    grid = fastGrid(rows)
    colnames(grid) = ids
    return(grid)
  }

  ### If no mutations: Restrict combinations using PO pairs

  if(!hasParentsBeforeChildren(x))
    stop2("Pedigree is not sorted with parents before children")

  idsSorted = ids[order(int.ids)]
  rows = rows[order(int.ids)]

  # Parents
  if(Xchrom)
    pars = lapply(idsSorted, function(i) mother(x, i) |> match(idsSorted, nomatch = 0) |> .mysetdiff(0))
  else
    pars = lapply(idsSorted, function(i) parents(x, i) |> match(idsSorted, nomatch = 0) |> .mysetdiff(0))

  # Overlapping genotypes (as row numbers in `allg`)
  compat = lapply(1:nrow(allg), function(i)
    which(allg[i,1] == allg[,1] | allg[i,1] == allg[,2] | allg[i,2] == allg[,1] | allg[i,2] == allg[,2]))

  # Make grid restricting to compatible genotypes for PO-pairs
  grid = fastGridRestricted(rows, linkedWith = pars, compatible = compat)
  colnames(grid) = idsSorted

  # Original order
  grid = grid[, ids, drop = FALSE]

  grid
}
