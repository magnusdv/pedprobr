#' Pedigree likelihood computed by MERLIN
#'
#' For this function to work, MERLIN must be installed and the path to
#' merlin.exe included in the PATH variable. The `likelihood_merlin` function is
#' a wrapper for one particular functionality of MERLIN.
#'
#' By default the following MERLIN command is run via [system()],
#' after creating appropriate files in the current working directory:
#'
#' \preformatted{%
#'   merlin -p _merlin.ped -d _merlin.dat -m _merlin.map -f _merlin.freq
#'          --lik --bits:100 --megabytes:4000 --quiet
#' }
#'
#' @param x a [`ped`] object
#' @param markers an integer vector indicating which markers to use (default:
#'   all).
#' @param logbase a positive number, or NULL. If numeric, the log-likelihood is
#'   returned, with `logbase` as basis for the logarithm.
#' @param verbose a logical: Show MERLIN output and other information, or not.
#' @param generate.files a logical. If TRUE, the files 'merlin.ped',
#'   'merlin.dat', 'merlin.map', and 'merlin.freq'.
#' @param cleanup a logical: Should the MERLIN files be deleted automatically?
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
#' \dontrun{
#' x = pedtools::nuclearPed(1)
#' m = pedtools::marker(x, "3"=1:2)
#' x = pedtools::setMarkers(x, m)
#'
#' # Likelihood computation by MERLIN:
#' likelihood_merlin(x)
#' }
#'
#' @export
likelihood_merlin = function(x, markers = seq_len(nMarkers(x)), logbase=NULL, verbose = FALSE,
                             generate.files = TRUE, cleanup = generate.files, logfile = "") {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(!is.null(logbase))
    assert_that(assertthat::is.number(logbase), logbase>0)

  clean = function(cleanup, verbose, files)
    if (cleanup) {
      unlink(files)
      if (verbose)  cat("Files successfully removed\n")
  }

  if (nMarkers(x) == 0) stop2("Pedigree has no attached markers")
  x = selectMarkers(x, markers)

  if (generate.files) {
    files = writePed(x, prefix = "_merlin", merlin = TRUE,
                      what = c("ped", "dat", "map", "freq"), verbose=verbose)
  }
  # TODO: fix
  program = if (is_Xmarker(x$markerdata[[1]])) "minx" else "merlin"

  command = paste(program,
                  "-p _merlin.ped",
                  "-d _merlin.dat",
                  "-m _merlin.map",
                  "-f _merlin.freq",
                  "--lik",
                  "--bits:100",
                  "--megabytes:4000",
                  "--quiet")

  if (verbose)
    cat("\nExecuting the following command:\n", command, "\n\n", sep = "")

  # Run MERLIN and store output
  mout = suppressWarnings(system(command, intern = T))

  clean(cleanup, verbose, files)
  if (nzchar(logfile))
    write(mout, logfile)
  if (any(substr(mout, 1, 11) == "FATAL ERROR")) {
    cat("\n====================================\n",
        paste(mout[-(2:10)], collapse = "\n"),
        "====================================\n\n")
    return(invisible())
  }

  if (verbose) cat("Merlin run completed\n")

  if (!is.na(skipped <- which(substr(mout, 3, 9) == "SKIPPED")[1]))
    stop2(paste(mout[c(skipped - 1, skipped)], collapse = "\n"))

  # Extract likelihood value
  nchars = nchar(mout)
  chromLines = which(substr(mout, 1, 20) == "Analysing Chromosome")
  likLines = which(substr(mout, 1, 27) == "lnLikelihood for 1 families")

  if(length(chromLines)==0 && length(likLines)==1) {
    lnlik = as.numeric(strsplit(mout[likLines]," = ")[[1]][2])
    if(!is.null(logbase))
      return(lnlik * log(exp(1), logbase))
    else
      return(exp(lnlik))
  }

  if(!length(chromLines)==length(likLines))
    stop2(mout)
  message("NB: Several chromosomes - output is total likelihood.")
  chroms = as.numeric(substr(mout[chromLines], 22, nchars[chromLines]))
  lnliks = as.numeric(unlist(lapply(strsplit(mout[likLines]," = "), '[', 2)))

  total = sum(lnliks)
  return(round(exp(total), 3))
}


