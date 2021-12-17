#' Retrieve NCBI taxonomy file from public database
#'
#' This script downloads the full NCBI taxonomy from
#' <ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz> and
#' extracts it to the folder in `target.dir`.
#'
#' @param target.dir path to the folder where the files will be saved (
#' accepts relative and absolute paths)
#' @param method Method to be used for downloading files. Current download
#' methods are "internal", "wininet" (Windows only) "libcurl", "wget" and
#' "curl", and there is a value "auto": see _Details_ and _Note_ in the
#' documentation of \code{utils::download.file()}.
#' @param unzip The unzip method to be used. See the documentation of
#' \code{utils::unzip()} for details.
#' @param url URL of the full NCBI taxonomy file.
#' @param timeout maximum time allowed for the download
#' @param ... additional attributes (currently ignored)
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   TaxonSampling::retrieve_NCBI_taxonomy(target.dir = "data_files/taxdump")
#' }
#'
#' @return No return value, called for side effects (see Description).

retrieve_NCBI_taxonomy <- function(target.dir,
                                   method  = "auto",
                                   unzip   = getOption("unzip"),
                                   url     = "https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip",
                                   timeout = 600,
                                   ...){

  # ================== Sanity checks ==================
  assertthat::assert_that(is.character(target.dir),
                          is.character(url),
                          length(url) == 1,
                          is.numeric(timeout),
                          length(timeout) == 1,
                          timeout > 0)

  if(!dir.exists(target.dir)){
    dir.create(target.dir, recursive = TRUE)
  } else {
    filelist <- dir(target.dir, full.names = TRUE)
    unlink(filelist, recursive = TRUE, force = TRUE)
  }

  oldtimeout <- getOption("timeout")
  options(timeout = timeout)
  res1 <- utils::download.file(url,
                               quiet    = FALSE,
                               destfile = paste0(target.dir, "/tmpdata.zip"),
                               cacheOK  = FALSE,
                               method   = method)
  options(timeout = oldtimeout)

  if(res1 != 0) stop("Error downloading file \n", url)

  utils::unzip(paste0(target.dir, "/tmpdata.zip"),
               unzip = unzip,
               exdir = target.dir)
  unlink(paste0(target.dir, "/__MACOSX"), recursive = TRUE, force = TRUE)

  file.remove(paste0(target.dir, "/tmpdata.zip"))

  invisible(TRUE)
}
