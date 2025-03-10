#' Genotype distribution for two linked markers
#'
#' Computes the joint genotype distribution of two markers for a specified
#' pedigree member, conditional on known genotypes and the recombination rate
#' between the markers.
#'
#' @param x A `ped` object or a list of such.
#' @param id A single ID label.
#' @param partialmarker1,partialmarker2 Either a `marker` object, or the name (or
#'   index) of a marker attached to `x`.
#' @param rho A single numeric in the interval `[0, 0.5]`: the recombination
#'   fraction between the two markers.
#' @param loopBreakers (Only relevant if the pedigree has loops). A vector with
#'   ID labels of individuals to be used as loop breakers. If NULL (default)
#'   loop breakers are selected automatically. See [breakLoops()].
#' @param verbose A logical.
#'
#' @return A named matrix giving the joint genotype distribution.
#'
#' @author Magnus Dehli Vigeland
#' @seealso [oneMarkerDistribution()]
#'
#' @examples
#'
#' # A sib-pair pedigree
#' x = nuclearPed(children = c("bro1", "bro2"))
#'
#' # Two SNP markers; first brother homozygous for the `1` allele
#' SNP1 = SNP2 = marker(x, bro1 = "1/1", afreq = c("1" = 0.5, "2" = 0.5))
#'
#' plot(x, marker = list(SNP1, SNP2))
#'
#' # Genotype distribution for the brother depends on linkage
#' twoMarkerDistribution(x, id = "bro2", SNP1, SNP2, rho = 0)
#' twoMarkerDistribution(x, id = "bro2", SNP1, SNP2, rho = 0.5)
#'
#' # X-linked
#' chrom(SNP1) = chrom(SNP2) = "X"
#'
#' plot(x, marker = list(SNP1, SNP2))
#'
#' twoMarkerDistribution(x, id = "bro2", SNP1, SNP2, rho = 0)
#' twoMarkerDistribution(x, id = "bro2", SNP1, SNP2, rho = 0.5)
#'
#' @export
twoMarkerDistribution <- function(x, id, partialmarker1, partialmarker2, rho = NULL,
                                  loopBreakers = NULL, verbose = TRUE) {

  if(length(id) != 1)
    stop2("Argument `id` must have length 1: ", id)

  if(is.pedList(x)) {
    if(is.marker(partialmarker1) || is.marker(partialmarker2))
      stop2("When `x` has multiple components, the partial markers must be attached")

    pednr = getComponent(x, id, checkUnique = TRUE)
    x = x[[pednr]]
  }

  if(!is.ped(x))
    stop2("Input is not a pedigree")

  m1 = partialmarker1
  if (!is.marker(m1)) {
    if(length(m1) != 1)
      stop2("Argument `partialmarker1` must have length 1: ", partialmarker1)
    m1 = getMarkers(x, markers = m1)[[1]]
  }

  m2 = partialmarker2
  if (!is.marker(m2)) {
    if(length(m2) != 1)
      stop2("Argument `partialmarker2` must have length 1: ", partialmarker2)
    m2 = getMarkers(x, markers = m2)[[1]]
  }
  if (!is.null(x$LOOP_BREAKERS))
    stop2("Pedigrees with pre-broken loops are not allowed in this function")

  if (!identical(chrom(m1), chrom(m2)))
    stop2("Partial markers are on different chromosomes: ", toString(c(chrom(m1), chrom(m2))))

  checkRho(rho)

  onX = isXmarker(m1)
  XandMale = onX && getSex(x, id) == 1

  if (verbose) {
    cat("Known genotypes:\n")

    df = as.data.frame(setMarkers(x, list(m1, m2)))[-(2:4)]

    # Add arrow in fourth column
    df = cbind(df, arrow = "", stringsAsFactors = FALSE)
    df$arrow[internalID(x, id)] = " <==="
    names(df)[4] = ""

    print(df, row.names = FALSE)

    cat("\nAllele frequencies, marker 1:\n")
    print(data.frame(as.list(afreq(m1)), check.names = FALSE), row.names = FALSE)
    cat("\nAllele frequencies, marker 2:\n")
    print(data.frame(as.list(afreq(m2)), check.names = FALSE), row.names = FALSE)

    cat("\nRecombination rate :", rho)
    cat("\nChromosome type    :", ifelse(onX, "X-linked", "autosomal"))
    cat("\nTarget individual  :", id, "\n")
  }

  # Start timer
  starttime = Sys.time()

  # Do this before loop breaking, since eliminate2 works better WITH the loops.
  grid.subset = fastGrid(c(genoCombinations(x, m1, id, make.grid = FALSE),
                           genoCombinations(x, m2, id, make.grid = FALSE)))

  if (x$UNBROKEN_LOOPS) {
    x = breakLoops(setMarkers(x, list(m1, m2)), loopBreakers = loopBreakers, verbose = verbose)
    m1 = x$MARKERS[[1]]
    m2 = x$MARKERS[[2]]
  }

  int.id = internalID(x, id)
  allgenos1 = allGenotypes(nAlleles(m1))
  allgenos2 = allGenotypes(nAlleles(m2))
  alleles1 = alleles(m1)
  alleles2 = alleles(m2)

  # Character with genotype labels
  if (XandMale) {
    geno.names = list(alleles1, alleles2)
  }
  else {
    gt1.strings = paste(alleles1[allgenos1[, 1]], alleles1[allgenos1[, 2]], sep = "/")
    gt2.strings = paste(alleles2[allgenos2[, 1]], alleles2[allgenos2[, 2]], sep = "/")
    geno.names = list(gt1.strings, gt2.strings)
  }

  # Create output array. Will hold likelihood of each genotype combo
  probs = array(0, dim = lengths(geno.names, use.names = FALSE), dimnames = geno.names)

  # Subset of `probs` that is affected by grid.subset
  probs.subset = grid.subset

  # Needs adjustment for X if `id` is male
  if(XandMale) {
    homoz1 = which(allgenos1[,1] == allgenos1[,2])
    homoz2 = which(allgenos2[,1] == allgenos2[,2])
    probs.subset[, 1] = match(probs.subset[, 1], homoz1)
    probs.subset[, 2] = match(probs.subset[, 2], homoz2)
  }

  marginal = likelihood2(x, marker1 = m1, marker2 = m2, rho = rho)
  if (marginal == 0)
    stop2("Partial marker data is impossible")

  if(verbose) {
    cat("Marginal likelihood:", marginal, "\n")
    cat("Calculations needed:", nrow(grid.subset), "\n")
  }

  probs[probs.subset] = apply(grid.subset, 1, function(allg_rows) {
    m1[int.id, ] = allgenos1[allg_rows[1], ]
    m2[int.id, ] = allgenos2[allg_rows[2], ]
    likelihood2(x, marker1 = m1, marker2 = m2, rho = rho)
  })

  # Timing
  totalTime = format(Sys.time() - starttime, digits = 3)
  if(verbose)
    cat("\nAnalysis finished in", totalTime, "\n")

  probs/marginal
}
