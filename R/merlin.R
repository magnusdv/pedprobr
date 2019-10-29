#' Pedigree likelihood computed by MERLIN
#'
#' This function is a wrapper of the "--likelihood" functionality of the MERLIN
#' software. It computes the total likelihood of the pedigree given the
#' indicated marker data. For this function to work, MERLIN must be installed
#' and correctly pointed to in the PATH variable.
#'
#' By default the following MERLIN command is run via [system()], after creating
#' appropriate files in the current working directory:
#'
#' \preformatted{% merlin -p _merlin.ped -d _merlin.dat -m _merlin.map -f
#' _merlin.freq --likelihood --bits:100 --megabytes:4000 --quiet }
#'
#' @param x a [`ped`] object
#' @param markers a vector of names or indices of markers attached to `x`.
#'   (Default: all markers).
#' @param logbase a positive number, or NULL. If numeric, the log-likelihood is
#'   returned, with `logbase` as basis for the logarithm.
#' @param verbose a logical.
#' @param generateFiles a logical. If TRUE, input files to MERLIN named
#'   '_merlin.ped', '_merlin.dat', '_merlin.map', and '_merlin.freq' are created
#'   in the current directory. If FALSE, no files are created.
#' @param cleanup a logical. If TRUE, the MERLIN input files are deleted after
#'   the call to MERLIN.
#' @param logfile a character. If this is given, the MERLIN screen output will
#'   be written to a file with this name.
#'
#' @return A number.
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
#' m = marker(x, "3" = 1:2)
#' x = setMarkers(x, m)
#'
#' # Likelihood computation by MERLIN:
#' likelihoodMerlin(x)
#' }
#'
#' @export
likelihoodMerlin = function(x, markers = seq_len(nMarkers(x)), logbase = NULL,
                            verbose = FALSE, generateFiles = TRUE,
                            cleanup = generateFiles, logfile = "") {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(!is.null(logbase) && (!is.numeric(logbase) || length(logbase) != 1 || logbase <= 0))
    stop2("`logbase` must be a single positive number: ", logbase)

  # Select markers
  if (nMarkers(x) == 0)
    stop2("Pedigree has no attached markers")
  x = selectMarkers(x, markers)

  # MERLIN or MINX?
  xchrom = vapply(x$MARKERS, isXmarker, FUN.VALUE = FALSE)
  if(all(xchrom)) {
    message("All markers are X-linked; calling MINX")
    program = "minx"
  }
  else if(all(!xchrom)) {
    message("All markers are autosomal; calling MERLIN")
    program = "merlin"
  }
  else
    stop2("Both autosomal and X-linked markers are selected\n",
          "Please use the `markers` argument to run these in separate calls")

  # Generate input files to MERLIN/MINX
  if (generateFiles) {
    files = writePed(x, prefix = "_merlin", merlin = TRUE,
                      what = c("ped", "dat", "map", "freq"), verbose = verbose)
  }

  # Utility for cleaning up afterwards
  clean = function(cleanup, verbose, files)
    if (cleanup) {
      unlink(files)
      if (verbose)  message("Temporary files removed")
    }

  command = paste(program,
                  "-p _merlin.ped",
                  "-d _merlin.dat",
                  "-m _merlin.map",
                  "-f _merlin.freq",
                  "--likelihood",
                  "--bits:100",
                  "--megabytes:4000",
                  "--quiet")

  if (verbose)
    message("\nExecuting the following command:\n", command, "\n")

  # Run MERLIN and store output
  mout = suppressWarnings(system(command, intern = TRUE))

  # Clean up
  clean(cleanup, verbose, files)

  # Write logfile if indicated
  if (nzchar(logfile))
    write(mout, logfile)

  if (any(substr(mout, 1, 11) == "FATAL ERROR")) {
    warning(paste0(mout, collapse = "\n"), "\nFATAL ERROR reported by merlin")
    return(invisible())
  }

  if (verbose) message("MERLIN run completed")

  if (!is.na(skipped <- which(substr(mout, 3, 9) == "SKIPPED")[1]))
    stop2(paste(mout[c(skipped - 1, skipped)], collapse = "\n"))

  # Extract likelihood value
  chromLines = which(substr(mout, 1, 20) == "Analysing Chromosome")
  likLines = which(substr(mout, 1, 27) == "lnLikelihood for 1 families")

  if(length(chromLines) == 0 && length(likLines) == 1) {
    lnlik = as.numeric(strsplit(mout[likLines]," = ")[[1]][2])
    if(!is.null(logbase))
      return(lnlik * log(exp(1), logbase))
    else
      return(exp(lnlik))
  }

  if(!length(chromLines) == length(likLines))
    stop2(mout)
  message("NB: Several chromosomes - output is total likelihood.")

  # Not used yet
  # nchars = nchar(mout)
  # chroms = as.numeric(substr(mout[chromLines], 22, nchars[chromLines]))

  # Extract log-likelihoods
  lnliks = as.numeric(unlist(lapply(strsplit(mout[likLines]," = "), '[', 2)))

  total = sum(lnliks)
  return(round(exp(total), 3))
}


