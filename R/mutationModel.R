# Re-export helper function for creating mutation models

#' @importFrom pedmut mutationModel
#' @export
pedmut::mutationModel



#' @importFrom pedmut isStationary
hasStationaryModel = function(m) {
  mut = mutmod(m)
  if(is.null(mut)) return(TRUE)

  sexEq = attr(mut, 'sexEqual')
  afr = afreq(m)

  isStationary(mut$male, afr) &&
    (isTRUE(sexEq) || isStationary(mut$female, afr))
}
