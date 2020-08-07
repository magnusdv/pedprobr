#' Set a mutation model
#'
#' This function offers a convenient way to attach mutation models to a pedigree
#' with marker data. It wraps [pedmut::mutationModel()], which does the main
#' work of creating the models, but relieves the user from having to loop
#' through the markers in order to supply the correct alleles and frequencies
#' for each marker.
#'
#' Currently, the following models are implemented in the `pedmut` package:
#'
#' * `equal` :  All mutations equally likely; probability \eqn{1-rate} of no
#' mutation
#'
#' * `proportional` : Mutation probabilities are proportional to the target
#' allele frequencies
#'
#' * `onestep`: A mutation model for microsatellite markers, allowing mutations
#' only to the nearest neighbours in the allelic ladder. For example, '10' may
#' mutate to either '9' or '11', unless '10' is the lowest allele, in which case
#' '11' is the only option. This model is not applicable to loci with
#' non-integral microvariants.
#'
#' * `stepwise`: A common model in forensic genetics, allowing different
#' mutation rates between integer alleles (like '16') and non-integer
#' "microvariants" like '9.3'). Mutations also depend on the size of the
#' mutation if the parameter 'range' differs from 1.
#'
#' * `custom` : Allows any mutation matrix to be provided by the user, in the
#' `matrix` parameter
#'
#' * `random` : This produces a matrix of random numbers, where each row is
#' normalised so that it sums to 1
#'
#' * `trivial` : The identity matrix; i.e. no mutations are possible.
#'
#' @param x A `ped` object or a list of such.
#' @param markers A vector of names or indices referring to markers attached to
#'   `x`.
#' @param model A model name implemented by [pedmut::mutationModel()]. See
#'   Details.
#' @param ... Arguments forwarded to [pedmut::mutationModel()], e.g., `rate`.
#'
#' @return An object similar to `x`.
#'
#' @examples
#' ### Example requires the pedmut package ###
#' if (requireNamespace("pedmut", quietly = TRUE)){
#'
#' # A pedigree with data from a single marker
#' x = nuclearPed(1)
#' x = setMarkers(x, marker(x, geno = c("a/a", NA, "b/b"))) # mutation!
#'
#' # Set `equal` model
#' x = setMutationModel(x, marker = 1, model = "equal", rate = 0.01)
#'
#' # Inspect model
#' mutmod(x, 1)
#'
#' }
#'
#' @importFrom pedmut mutationModel
#' @export
setMutationModel = function(x, markers = NULL, model, ...) {
  if (!requireNamespace("pedmut", quietly = TRUE))
    stop2("Package `pedmut` must be installed in order to include mutation models")

  opts = list(...)

  mIdx = whichMarkers(x, markers)
  for(i in mIdx) {
    args = c(list(model = model, alleles = alleles(x, i), afreq = afreq(x, i)), opts)

    modi = do.call(pedmut::mutationModel, args)
    mutmod(x, i) = modi
  }

  x
}
