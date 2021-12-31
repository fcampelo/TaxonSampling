#' Retrieve NCBI taxonomy file from public database
#'
#' This script downloads the full NCBI taxonomy from
#' <ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz> and
#' extracts it to the folder in `target.dir`.
#'
#' @param target.dir path to the folder where the files will be saved (
#' accepts relative and absolute paths)
#' @param spp_file path to file for saving the species count for each taxon ID
#' (saved as a tsv file)
#' two-column data frame with the input taxon IDs in the first
#' column, and the corresponding number of known species in the second column.
#' Ignored if `ids_file` is not `NULL`.
#' @param method Method to be used for downloading files. Current download
#' methods are "internal", "wininet" (Windows only) "libcurl", "wget" and
#' "curl", and there is a value "auto": see _Details_ and _Note_ in the
#' documentation of \code{utils::download.file()}.
#' @param unzip The unzip method to be used. See the documentation of
#' \code{utils::unzip()} for details.
#' @param url URL of the full NCBI taxonomy file.
#' @param timeout maximum time allowed for the download, in seconds.
#' Increase when under a slow connection.
#' @param ncpus number of cores to use for species counting.
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
                                   spp_file = NULL,
                                   method  = "auto",
                                   unzip   = getOption("unzip"),
                                   url     = "https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip",
                                   timeout = 600,
                                   ncpus   = 1,
                                   ...){

  # ================== Sanity checks ==================
  assertthat::assert_that(is.character(target.dir),
                          is.character(url),
                          length(url) == 1,
                          is.numeric(timeout),
                          length(timeout) == 1,
                          timeout > 0,
                          is.null(spp_file) ||
                            (is.character(spp_file) &&
                               length(spp_file) == 1))

  if(!dir.exists(target.dir)){
    dir.create(target.dir, recursive = TRUE)
  }
  if(file.exists(paste0(target.dir, "/tmp_NCBI_taxdump.zip"))){
    file.remove(paste0(target.dir, "/tmp_NCBI_taxdump.zip"))
  }

  oldtimeout <- getOption("timeout")
  options(timeout = timeout)
  res1 <- utils::download.file(url,
                               quiet    = FALSE,
                               destfile = paste0(target.dir, "/tmp_NCBI_taxdump.zip"),
                               cacheOK  = FALSE,
                               method   = method)
  options(timeout = oldtimeout)

  if(res1 != 0) stop("Error downloading file \n", url)

  utils::unzip(paste0(target.dir, "/tmp_NCBI_taxdump.zip"),
               unzip = unzip,
               exdir = target.dir)
  unlink(paste0(target.dir, "/__MACOSX"), recursive = TRUE, force = TRUE)

  file.remove(paste0(target.dir, "/tmp_NCBI_taxdump.zip"))

  if(!is.null(spp_file)){
    if(!dir.exists(dirname(spp_file))) {
      dir.create(dirname(spp_file), recursive = TRUE)
    }

    message("Extracting species counts. This may take a while...",
            "\n(Takes longer in Windows than Unix-based systems)")

    nodes <- as.data.frame(
      data.table::fread(paste(target.dir, "nodes.dmp", sep = "/"),
                        sep = "|", strip.white = TRUE,
                        colClasses = c("numeric", "numeric", "character", rep("NULL", 16)),
                        col.names = c("id", "parent", "level"),
                        verbose = FALSE))
    nodes$level <- gsub("\\t", "", nodes$level)

    cl <- parallel::makeCluster(min(parallel::detectCores() - 1, ncpus))
    parallel::clusterExport(cl, "nodes")

    x <- table(unlist(
      pbapply::pblapply(nodes$id[nodes$level == "species"],
                        CHNOSZ::allparents,
                        nodes = nodes,
                        cl    = cl)))

    x <- data.frame(id            = names(x),
                    species_count = unname(as.numeric(x)))

    parallel::stopCluster(cl)

    utils::write.table(x, file = spp_file, sep = "\t", row.names = FALSE)

  }

  invisible(TRUE)
}
