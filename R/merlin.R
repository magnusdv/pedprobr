#' Pedigree likelihood computed by MERLIN
#'
#' For these functions to work, the program MERLIN (see References below) must
#' be installed and correctly pointed to in the PATH variable. The `merlin()`
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
#' @param linkageMap a data frame with three columns (chromosome; marker name;
#'   centiMorgan position) to be used as the marker map by MERLIN.
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
#' @param logbase Either NULL (default) or a positive number indicating the
#'   basis for logarithmic output. Typical values are `exp(1)` and 10.
#' @param program A character containing "merlin", "minx" or both (default),
#'   optionally including full paths.
#' @param version A logical. If TRUE (default), it is checked that running
#'   `program` produces a printout starting with "MERLIN 1.1.2".
#' @param error A logical, indicating if an error should be raised if `program`
#'   is not found. Default: FALSE.
#'
#' @return `merlin()` returns the screen output of MERLIN invisibly.
#'
#'   `likelihoodMerlin()` returns a single number; the total likelihood using
#'   all indicated markers.
#'
#'   `checkMerlin()` returns TRUE if the MERLIN executable indicated by
#'   `program` is found on the system. Otherwise FALSE, or (if `error = TRUE`)
#'   an error is raised.
#'
#' @author Magnus Dehli Vigeland
#' @references <https://csg.sph.umich.edu/abecasis/merlin/>
#'
#' @examples
#'
#' if(checkMerlin()) {
#'
#' ### Trivial example for validation
#' x = nuclearPed(1) |>
#'   addMarker("1" = "1/2") |>            # likelihood = 1/2
#'   addMarker("1" = "1/1", "3" = "1/2")  # likelihood = 1/8
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
#' genotype(y[[1]], marker = 1, id = "1") = "1/2"
#' genotype(y[[2]], marker = 1, id = "2") = "1/1"
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
#' # stopifnot(all.equal(signif(L1, 3), signif(L2, 3)))
#' }
#'
#' @export
merlin = function(x, options, markers = NULL, linkageMap = NULL, verbose = TRUE,
                  generateFiles = TRUE, cleanup = TRUE, dir = tempdir(),
                  logfile = NULL, merlinpath = NULL) {

  # Select markers
  if (!hasMarkers(x))
    stop2("Pedigree has no attached markers")

  # Set linkage map if provided
  if(!is.null(linkageMap))
    x = setMap(x, map = linkageMap, matchNames = NA)

  if(is.null(markers))
    markers = seq_len(nMarkers(x))
  x = selectMarkers(x, markers)

  # MERLIN or MINX?
  xchrom = isXmarker(x)
  if(all(xchrom)) {
    if(verbose) cat("All markers are X-linked; calling MINX\n")
    program = "minx"
  }
  else if(all(!xchrom)) {
    program = "merlin"
  }
  else
    stop2("Both autosomal and X-linked markers are selected\n",
          "Please use the `markers` argument to run these in separate calls")

  if(!is.null(merlinpath))
    program = file.path(merlinpath, program)

  checkMerlin(program, version = FALSE, error = TRUE)

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
likelihoodMerlin = function(x, markers = NULL, linkageMap = NULL, rho = NULL, logbase = NULL,
                            options = "--likelihood --bits:100 --megabytes:4000 --quiet",
                            ...) {

  # Select markers
  x = selectMarkers(x, markers %||% seq_len(nMarkers(x)))

  # If rho given, replace map
  if(!is.null(rho)) {

    if(!is.null(linkageMap))
      stop2("At least one of `rho` and `linkageMap` must be NULL")

    if(length(rho) != nMarkers(x) - 1)
      stop2("Argument `rho` must have length one less than the number of markers")

    # Avoid infinities
    rho[rho == 0.5] = haldane(cM = 500)

    # If no chromosome info given, place all markers on chrom 1
    chr = chrom(x)
    if(all(is.na(chr)))
      chr = rep_len(1, nMarkers(x))

    # Convert to centiMorgan positions
    cm = c(0, haldane(rho = rho))

    linkageMap = data.frame(CHROM = chr, MARKER = name(x), CM = cm)
  }

  # Run MERLIN
  mout = merlin(x, linkageMap = linkageMap, options = options, ...)

  # Catch possible errors
  if (any(skipped <- substr(mout, 3, 9) == "SKIPPED"))
    stop2(paste(mout[sort(c(which(skipped)-1, which(skipped)))], collapse = "\n"))

  # Bad inheritance? Return lnLik = -Inf (i.e. L = 0)
  if (any(grepl("Skipping Marker .* [BAD INHERITANCE]", mout)))
    return(fixMerlinLog(-Inf, logbase = logbase))

  # Different chromosomes?
  chromLines = which(substr(mout, 1, 20) == "Analysing Chromosome")

  # Lines with loglik results
  nFam = if(is.pedList(x)) length(x) else 1
  likLines = which(substr(mout, 1, 27) == sprintf("lnLikelihood for %d families", nFam))

  if(length(likLines) == 0)
    return(fixMerlinLog(-Inf, logbase = logbase))

  # If single output value: Return likelihood
  if(length(chromLines) == 0 && length(likLines) == 1) {
    lnlik = as.numeric(strsplit(mout[likLines]," = ")[[1]][2])
    return(fixMerlinLog(lnlik, logbase = logbase))
  }

  #---------------
  # Otherwise: Return total likelihood
  #---------------

  if(length(chromLines) != length(likLines))
    stop2(mout)

  # Extract log-likelihoods
  lnliks = as.numeric(unlist(lapply(strsplit(mout[likLines]," = "), '[', 2)))

  # Return total
  totalLnLik = sum(lnliks)

  fixMerlinLog(totalLnLik, logbase = logbase)
}



#' @rdname merlin
#' @export
checkMerlin = function(program = NULL, version = TRUE, error = FALSE) {

  # Default NULL for back compatibility
  if(is.null(program))
    return(checkMerlin("merlin", version = version, error = error) &&
             checkMerlin("minx", version = version, error = error))

  if(!is.character(program) || length(program) != 1)
    stop2("`program` must be a character of length 1")

  if(!nzchar(pth <- Sys.which(program))) {
    if(error)
      stop2("Executable not found. Use `merlinpath` to supply the path to the MERLIN folder")
    return(FALSE)
  }

  if(version) {
    out = tryCatch(
      suppressWarnings(system(program, intern = TRUE)),
      error = function(e) if(error) stop2(e) else return(FALSE))

    if(!isTRUE(startsWith(out[1], "MERLIN 1.1.2"))) {
      if(error)
        stop2(pth, "\nExpected printout to start with 'MERLIN 1.1.2', but received:\n", paste(out, collapse = "\n"))
      return(FALSE)
    }
  }

  TRUE
}
