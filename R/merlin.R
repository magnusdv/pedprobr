#' Pedigree likelihood computed by MERLIN
#'
#' For this functions to work, the program MERLIN (see References below) must be
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
#'
#' @return `merlin()` returns the screen output of MERLIN invisibly.
#'
#'   `likelihoodMerlin()` returns a single number; the total likelihood using
#'   all indicated markers.
#'
#' @author Magnus Dehli Vigeland
#' @references <http://csg.sph.umich.edu/abecasis/Merlin/>
#'
#' @examples
#'
#' \donttest{
#' ### Requires MERLIN to be installed ###
#'
#' x = nuclearPed(1)
#' m1 = marker(x, "1" = 1:2)           # likelihood = 1/2
#' m2 = marker(x, "1" = 1, "3" = 1:2)    # likelihood = 1/8
#' x = setMarkers(x, list(m1,m2))
#'
#' # Likelihood computation by MERLIN:
#' likelihoodMerlin(x, verbose = FALSE)
#' likelihoodMerlin(x, markers = 1, verbose = FALSE)
#' }
#'
#' @export
merlin = function(x, options, markers = NULL, verbose = TRUE,
                  generateFiles = TRUE, cleanup = TRUE, dir = tempdir(),
                  logfile = NULL) {
  #if(!is.ped(x)) stop2("Input is not a `ped` object")

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
    program = "minx"
  }
  else if(all(!xchrom)) {
    program = "merlin"
  }
  else
    stop2("Both autosomal and X-linked markers are selected\n",
          "Please use the `markers` argument to run these in separate calls")

  prefix = file.path(dir, "_merlin")
  # Generate input files to MERLIN/MINX
  if (generateFiles) {
    files = writePed(x, prefix = prefix, merlin = TRUE,
                     what = c("ped", "dat", "map", "freq"), verbose = verbose)

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

  if (any(substr(mout, 1, 11) == "FATAL ERROR")) {
    warning(paste0(mout, collapse = "\n"), "\nFATAL ERROR reported by merlin")
  }
  else if (verbose) cat("\nMERLIN run completed\n")

  invisible(mout)
}

#' @param ... Further arguments passed on to `merlin`
#'
#' @rdname merlin
#' @export
likelihoodMerlin = function(x, ...) {

  # Run MERLIN
  args = "--likelihood --bits:100 --megabytes:4000 --quiet"
  mout = merlin(x, args, ...)

  # Catch possible error
  if (!is.na(skipped <- which(substr(mout, 3, 9) == "SKIPPED")[1]))
    stop2(paste(mout[c(skipped - 1, skipped)], collapse = "\n"))

  # Different chromosomes?
  chromLines = which(substr(mout, 1, 20) == "Analysing Chromosome")

  # Lines with loglik results
  nFam = if(is.pedList(x)) length(x) else 1
  likLines = which(substr(mout, 1, 27) == sprintf("lnLikelihood for %d families", nFam))

  # If single output value: Return likelihood
  if(length(chromLines) == 0 && length(likLines) == 1) {
    lnlik = as.numeric(strsplit(mout[likLines]," = ")[[1]][2])
    return(exp(lnlik))
  }

  #---------------
  # Otherwise: Return total likelihood
  #---------------
  message("NB: Several chromosomes - output is total likelihood.")

  if(length(chromLines) != length(likLines))
    stop2(mout)

  # Extract log-likelihoods
  lnliks = as.numeric(unlist(lapply(strsplit(mout[likLines]," = "), '[', 2)))

  # Return total
  total = sum(lnliks)
  return(round(exp(total), 3))
}

