#' Pedigree likelihood computed by MERLIN
#'
#' For these functions to work, the program MERLIN (see References below) must be
#' installed and correctly pointed to in the PATH variable. The `merlin()`
#' function is a general wrapper which runs MERLIN with the indicated options,
#' after creating the appropriate input files. For convenience, MERLIN's
#' "--likelihood" functionality is wrapped in a separate function.
#'
#' The `merlin()` function creates input files "_merlin.ped", "_merlin.dat",
#' "_merlin.map" and "_merlin.freq" in the `dir` directory, and then runs the
#' following command through a call to [system()]:
#'
#' \preformatted{merlin -p _merlin.ped -d _merlin.dat -m _merlin.map -f
#' _merlin.freq  <options> }
#'
#' `likelihoodMerlin()` first runs `merlin()` with `options = "--likelihood
#' --bits:100 --megabytes:4000 --quiet"`, and then extracts the likelihood
#' values from the MERLIN output. Note that the output is the *total* likelihood
#' including all markers.
#'
#' For likelihood computations with linked markers, the argument `rho` should
#' indicate the recombination fractions between each consecutive pair of markers
#' (i.e., `rho[i]` is the recombination rate between markers `i-1` and `i`).
#' These will be converted to centiMorgan distances using Haldane's map
#' function, and used to create genetic marker map in a MERLIN-friendly format.
#'
#' @param x a [`ped`] object.
#' @param options a single string containing all arguments to merlin except for
#'   the input file indications.
#' @param markers a vector of names or indices of markers attached to `x`.
#'   (Default: all markers).
#' @param verbose a logical.
#' @param generateFiles a logical. If TRUE (default), input files to MERLIN
#'   named '_merlin.ped', '_merlin.dat', '_merlin.map', and '_merlin.freq' are
#'   created in the directory indicated by `dir`. If FALSE, no files are
#'   created.
#' @param cleanup a logical. If TRUE (default), the MERLIN input files are
#'   deleted after the call to MERLIN.
#' @param dir the name of the directory where input files should be written.
#' @param logfile a character. If this is given, the MERLIN screen output will
#'   be dumped to a file with this name.
#' @param merlinpath the path to the folder containing the merlin executables.
#'   If the executables are on the system's search path, this can be left as
#'   NULL (default).
#' @param rho A vector of length one less than the number of markers, specifying
#'   the recombination rate between each consecutive pair.
#'
#' @return `merlin()` returns the screen output of MERLIN invisibly.
#'
#'  `likelihoodMerlin()` returns a single number; the total likelihood using
#'   all indicated markers.
#'
#'  `checkMerlin()` returns TRUE if MERLIN is installed and available on the
#' system path, and FALSE otherwise.
#'
#' @author Magnus Dehli Vigeland
#' @references <http://csg.sph.umich.edu/abecasis/Merlin/>
#'
#' @examples
#'
#' merlinInstalled = checkMerlin()
#'
#' if(merlinInstalled) {
#'
#' ### Trivial example for validation
#' x = nuclearPed(1)
#' m1 = marker(x, "1" = 1:2)           # likelihood = 1/2
#' m2 = marker(x, "1" = 1, "3" = 1:2)  # likelihood = 1/8
#' x = setMarkers(x, list(m1, m2))
#'
#' # MERLIN likelihoods
#' lik1 = likelihoodMerlin(x, markers = 1, verbose = FALSE)
#' lik2 = likelihoodMerlin(x, markers = 2, verbose = FALSE)
#' likTot = likelihoodMerlin(x, verbose = FALSE)
#' stopifnot(all.equal(
#'   round(c(lik1, lik2, likTot), c(3,3,4)), c(1/2, 1/8, 1/16)))
#'
#' # Example with ped lists
#' y = list(singleton(1), singleton(2))
#' y = setMarkers(y, locus = list(alleles = 1:2))
#' genotype(y[[1]], marker = 1, id = '1') = 1:2
#' genotype(y[[2]], marker = 1, id = '2') = 1
#' lik = likelihoodMerlin(y, verbose = FALSE)
#' stopifnot(all.equal(round(lik, 3), 1/8))
#'
#' ### Linked markers
#' z = nuclearPed(2)
#' m = marker(z, geno = c("1/1", "1/2", "1/2", "1/2"))
#' z = setMarkers(z, list(m, m))
#'
#' # By MERLIN...
#' L1 = likelihoodMerlin(z, markers = 1:2, rho = 0.25, verbose = FALSE)
#'
#' # ...and by pedprobr
#' L2 = likelihood2(z, marker1 = 1, marker2 = 2, rho = 0.25)
#'
#' stopifnot(all.equal(round(L1, 6), round(L2, 6)))
#' }
#'
#' @export
merlin = function(x, options, markers = NULL, verbose = TRUE,
                  generateFiles = TRUE, cleanup = TRUE, dir = tempdir(),
                  logfile = NULL, merlinpath = NULL) {

  # Select markers
  if (!hasMarkers(x))
    stop2("Pedigree has no attached markers")
  if(is.null(markers))
    markers = seq_len(nMarkers(x))
  x = selectMarkers(x, markers)

  # MERLIN or MINX?
  xchrom = isXmarker(x)
  if(all(xchrom)) {
    if(verbose) cat("All markers are X-linked; calling MINX\n")
    program = "minx.exe"
  }
  else if(all(!xchrom)) {
    program = "merlin.exe"
  }
  else
    stop2("Both autosomal and X-linked markers are selected\n",
          "Please use the `markers` argument to run these in separate calls")

  if(!is.null(merlinpath))
    program = file.path(merlinpath, program)

  prefix = file.path(dir, "_merlin")
  # Generate input files to MERLIN/MINX
  if (generateFiles) {
    files = writePed(x, prefix = prefix, merlin = TRUE, verbose = verbose)

    if(cleanup)
      on.exit({unlink(files); if (verbose) cat("MERLIN input files removed\n")})
  }

  commandArgs = c(program,
                  sprintf("-p %s.ped", prefix),
                  sprintf("-d %s.dat", prefix),
                  sprintf("-m %s.map", prefix),
                  sprintf("-f %s.freq", prefix),
                  options)
  command = paste(commandArgs, collapse = " ")

  if (verbose)
    cat("\nExecuting the following command:\n", paste0(commandArgs, collapse = "\n "), "\n", sep = "")

  # Run MERLIN and store output
  mout = suppressWarnings(system(command, intern = TRUE))

  # Write logfile if indicated
  if (!is.null(logfile))
    write(mout, logfile)

  err = NULL
  if (any(fatal <- substr(mout, 1, 11) == "FATAL ERROR"))
    err = mout[which(fatal)[1]:length(mout)]
  else if (any(warn <- substr(mout, 2, 8) == "WARNING"))
    err = mout[which(warn)[1] + 0:5]

  if(!is.null(err))
    warning(paste0(err, collapse = "\n"), call. = FALSE)
  else if (verbose)
    cat("\nMERLIN run completed\n")

  invisible(mout)
}

#' @param ... Further arguments passed on to `merlin()`.
#'
#' @rdname merlin
#' @export
likelihoodMerlin = function(x, markers = NULL, rho = NULL,
                            options = "--likelihood --bits:100 --megabytes:4000 --quiet",
                            ...) {

  # Select markers
  x = selectMarkers(x, markers %||% seq_len(nMarkers(x)))

  # If rho given, replace map
  if(!is.null(rho)) {
    if(length(rho) != nMarkers(x) - 1)
      stop2("Argument `rho` must have length one less than the number of markers")

    # Avoid infinities
    rho[rho == 0.5] = haldane(cM = 500)

    # If no chromosome info given, place all markers on chrom 1
    if(all(is.na(chrom(x))))
      chrom(x) = 1

    # Set centiMorgan positions (using the posMb slot, but this is interpreted as cM by Merlin)
    posMb(x) = c(0, haldane(rho = rho))
  }

  # Run MERLIN
  mout = merlin(x, options = options, ...)

  # Catch possible errors
  if (any(skipped <- substr(mout, 3, 9) == "SKIPPED"))
    stop2(paste(mout[sort(c(which(skipped)-1, which(skipped)))], collapse = "\n"))

  # Bad inheritance? Return 0
  if (any(grepl("Skipping Marker .* [BAD INHERITANCE]", mout)))
    return(0)

  # Different chromosomes?
  chromLines = which(substr(mout, 1, 20) == "Analysing Chromosome")

  # Lines with loglik results
  nFam = if(is.pedList(x)) length(x) else 1
  likLines = which(substr(mout, 1, 27) == sprintf("lnLikelihood for %d families", nFam))

  if(length(likLines) == 0)
    return(0)

  # If single output value: Return likelihood
  if(length(chromLines) == 0 && length(likLines) == 1) {
    lnlik = as.numeric(strsplit(mout[likLines]," = ")[[1]][2])
    return(exp(lnlik))
  }

  #---------------
  # Otherwise: Return total likelihood
  #---------------

  if(length(chromLines) != length(likLines))
    stop2(mout)

  # Extract log-likelihoods
  lnliks = as.numeric(unlist(lapply(strsplit(mout[likLines]," = "), '[', 2)))

  # Return total
  total = sum(lnliks)
  return(exp(total))
}



#' @rdname merlin
#' @export
checkMerlin = function() {
  Sys.which("merlin.exe") != ""
}
