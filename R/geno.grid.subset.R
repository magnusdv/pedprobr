#' Genotype combinations
#'
#' Auxiliary functions computing possible genotype combinations in a pedigree.
#' These are not normally intended for end users.
#'
#'
#' @param x a [ped()] object.
#' @param partialmarker a [marker()] object compatible with `x`.
#' @param ids a vector with ID labels of one or more pedigree members.
#' @param make.grid a logical. If FALSE, a list is returned, otherwise
#' `fast.grid` is applied to the list before returning it.
#'
#' @export
geno.grid.subset = function(x, partialmarker, ids, make.grid = T) {
    int.ids = internalID(x, ids)
    nall = nAlleles(partialmarker)
    mutations = !is.null(attr(partialmarker, "mutmat"))
    chrom = if (is_Xmarker(partialmarker)) "X" else "AUTOSOMAL"

    allg = allGenotypes(nall)
    allg_ref = 1000 * (allg[, 1] + allg[, 2]) + abs(allg[, 1] - allg[, 2])

    match_ref_rows = function(genomatr) {
        # In: matrix with 2 rows (each column a genotype). Out: vector of 'allg' row numbers
      row1 = genomatr[1, ]
      row2 = genomatr[2, ]
      sort.int(unique.default(match(1000 * (row1 + row2) + abs(row1 - row2), allg_ref)))
    }

    switch(chrom, AUTOSOMAL = {
        glist = .build_genolist(x, partialmarker, eliminate = ifelse(mutations, 0, 100))
        if (attr(glist, "impossible")) stop2("Impossible partial marker")
        rows = lapply(glist[int.ids], match_ref_rows)
    }, X = {
        SEX = x$SEX
        glist = .build_genolist_X(x, partialmarker, eliminate = ifelse(mutations, 0, 100))
        if (attr(glist, "impossible")) stop2("Impossible partial marker")
        rows = lapply(int.ids, function(i) switch(SEX[i], glist[[i]], match_ref_rows(glist[[i]])))
    })
    if (make.grid)
        fast.grid(rows) else rows
}
