#' Genetic map functions
#'
#' Simple implementations of the classical map functions of Haldane and Kosambi,
#' relating the genetic distance and the recombination rate between two linked loci.
#'
#' @param cM A numeric vector with genetic distances in centiMorgan, or NULL.
#' @param rho A numeric vector with recombination rates, or NULL.
#'
#' @returns A numeric of the same length as the input.
#'
#' @examples
#' cM = 0:200
#' dat = cbind(Haldane = haldane(cM = cM),
#'             Kosambi = kosambi(cM = cM))
#' matplot(cM, dat, ylab = "Recombination rate", type = "l")
#' legend("topleft", legend = colnames(dat), col = 1:2, lty = 1:2)
#'
#'
#' rho = seq(0, 0.49, length = 50)
#' dat2 = cbind(Haldane = haldane(rho = rho),
#'              Kosambi = kosambi(rho = rho))
#' matplot(rho, dat2, xlab = "Recombination rate", ylab = "cM", type = "l")
#' legend("topleft", legend = colnames(dat), col = 1:2, lty = 1:2)
#'
#' @export
haldane = function(cM = NULL, rho = NULL) {
  if(is.null(cM) + is.null(rho) != 1)
    stop2("Exactly one of `cM` and `rho` must be NULL")

  if(is.null(cM))
    -50 * log(1 - 2 * rho)
  else
    .5 * (1 - exp(-cM/50))
}


#' @rdname haldane
#' @export
kosambi = function(cM = NULL, rho = NULL) {
  if(is.null(cM) + is.null(rho) != 1)
    stop2("Exactly one of `cM` and `rho` must be NULL")

  if(is.null(cM))
    25 * log((1 + 2*rho)/(1-2*rho))
  else
    .5 * (exp(cM/25) - 1)/(exp(cM/25) + 1)
}
