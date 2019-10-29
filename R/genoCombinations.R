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
#'   If `make.grid = TRUE`, the cartesian product of the vectors is taken,
#'   resulting in a matrix with one column for each element of `ids`.
#'
#' @export
genoCombinations = function(x, partialmarker, ids, make.grid = TRUE) {
    int.ids = internalID(x, ids)
    mutations = allowsMutations(partialmarker)

    allg = allGenotypes(nAlleles(partialmarker))
    homoz = which(allg[,1] == allg[,2])
    allgRef = 1000 * (allg[, 1] + allg[, 2]) + abs(allg[, 1] - allg[, 2])

    matchRefRows = function(genomatr) {
      # In: matrix with 2 rows (each column a genotype). Out: vector of 'allg' row numbers
      row1 = genomatr[1, ]
      row2 = genomatr[2, ]
      sort.int(unique.default(match(1000 * (row1 + row2) + abs(row1 - row2), allgRef)))
    }

    if (isXmarker(partialmarker)) {
      SEX = x$SEX
      glist = .buildGenolistX(x, partialmarker, eliminate = ifelse(mutations, 0, 100))
      if (attr(glist, "impossible")) stop2("Impossible partial marker")
      rows = lapply(int.ids, function(i) switch(SEX[i], homoz[glist[[i]]], matchRefRows(glist[[i]])))
    }
    else {
      glist = .buildGenolist(x, partialmarker, eliminate = ifelse(mutations, 0, 100))
      if (attr(glist, "impossible"))
        stop2("Impossible partial marker")
      rows = lapply(glist[int.ids], matchRefRows)
    }

    if (make.grid)
        fastGrid(rows) else rows
}
