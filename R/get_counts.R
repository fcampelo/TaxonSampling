#' Return the relevant counts for running TaxonSampling
#'
#' This function evaluates the taxonomical structure from a set of input IDs
#' and returns the count of each taxon ID with non-zero occurrences in the
#' taxonomy graph, as well as the count of known species under each taxon.
#'
#' @param taxonomy_path path to folder containing the NCBI taxonomy files
#' (i.e., the extracted contents of _taxdump.zip_, which can be downloaded from
#' <ftp://ftp.ncbi.nih.gov/pub/taxonomy/> or retrieved using
#' [retrieve_NCBI_taxonomy()]).
#' @param ids_file path to a tab-separated file containng two two columns,
#' with the input taxon IDs in the first
#' column, and the corresponding sequence IDs (corresponding
#' to the ID strings int the multifasta input file, without
#' the ">" line starter) in the second column.
#' @param ids_df two-column data frame with the input taxon IDs in the first
#' column, and the corresponding sequence IDs (corresponding
#' to the ID strings int the multifasta input file, without
#' the ">" line starter) in the second column. Ignored if
#' `ids_file` is not `NULL`.
#' @param spp_df two-column data frame with the input taxon IDs in the first
#' column, and the corresponding number of known species in the second column.
#' Ignored if `spp_file` is not `NULL`.
#' @param spp_file path to a tab-separated file containng two two columns,
#' with the input taxon IDs in the first column, and the corresponding number of
#' known species in the second column. This can be (i) generated locally by
#' [get_taxID_spp_counts()] (**note**: this is a very time-consuming function,
#' and can take a few days to compute); or a relatively recent  downloaded as part of the
#' package data files using [retrieve_data_files()].
#' @param verbose logical: regulates function echoing to console.
#' @param nodes data.frame containing the pre-processed information about
#' the NCBI taxonomy structure. This is generated either by using
#' [CHNOSZ::getnodes()], or as a result of a previous call to
#' [get_taxonomy_counts()]. If `nodes` is not `NULL` then `taxonomy_path` is
#' ignored.
#'
#' @return list object of class _taxonsampling_, containing:
#' \itemize{
#'     \item `$ids_df`: data.frame with taxon IDs in column 1 and corresponding
#'     sequence IDs in column 2, as loaded from `ids_file` or passed directly as
#'     input. Filtered to maintain only IDs that exist in `$nodes` and to remove
#'     duplicated IDs.
#'     \item `$nodes`: data.frame containing the pre-processed information about
#'     the NCBI taxonomy structure, extracted from file _nodes.dmp_ of the
#'     taxonomy files. Filtered to keep only nodes with IDs present in `$ids_df`
#'     and with a non-zero total count.
#'     \item `$countIDs`: numeric vector with the counts of the number of
#'     taxonomy nodes (of all levels) under each taxon ID.
#'     \item `spp_df`: two-column data frame with the input taxon IDs in the
#'     first column, and the corresponding number of known species in the second
#'    column. Filtered to have only the IDs present in `$countIDs`.
#'    \item `$countSpp`, a numeric vector with the counts of (known) species for each
#'    taxon from `taxlist$countIDs`.
#' }
#'
#' @export

get_counts <- function(taxonomy_path = NULL,
                       ids_file      = NULL,
                       ids_df        = NULL,
                       spp_file      = NULL,
                       spp_df        = NULL,
                       nodes         = NULL,
                       verbose       = TRUE) {

  # ===========================================================================
  # Sanity checks
  assertthat::assert_that(is.null(ids_df) ||
                            (is.data.frame(ids_df) && nrow(ids_df) > 0),
                          is.null(ids_file) ||
                            (is.character(ids_file) && length(ids_file) == 1),
                          (is.null(ids_df) + is.null(ids_file)) < 2,
                          is.null(taxonomy_path) ||
                            (is.character(taxonomy_path) &&
                               length(taxonomy_path) == 1 &&
                               dir.exists(taxonomy_path)),
                          is.null(spp_df) || (
                            is.data.frame(spp_df) &&
                              ncol(spp_df) >= 2 &&
                              nrow(spp_df) > 0),
                          is.null(spp_file) ||
                            (is.character(spp_file) &&
                               length(spp_file) == 1 &&
                               file.exists(spp_file)),
                          is.null(spp_df) + is.null(spp_file) < 2,
                          is.logical(verbose),
                          length(verbose) == 1,
                          is.null(nodes) || is.data.frame(nodes))

  # ===========================================================================
  # Load required files

  odt <- options("datatable.showProgress")
  options(datatable.showProgress = FALSE)
  # Load ids from file if required
  if(!is.null(ids_file)) {
    if(file.exists(ids_file)){
      ids_df <- as.data.frame(
        data.table::fread(ids_file, sep = "\t",
                          col.names = c("taxID", "seqID"),
                          verbose = FALSE))
    } else {
      stop("File ", ids_file, " not found.")
    }
  } else if(!is.null(ids_df)){
    names(ids_df) <- c("taxID", "seqID")
  }

  # Load nodes from file if required
  if(is.null(nodes)){
    nodes <- as.data.frame(
      data.table::fread(paste(taxonomy_path, "nodes.dmp", sep = "/"),
                        sep = "|", strip.white = TRUE,
                        colClasses = c("numeric", "numeric", "character", rep("NULL", 16)),
                        col.names = c("id", "parent", "level"),
                        verbose = FALSE))
    nodes$level <- gsub("\\t", "", nodes$level)
  } else {
    names(nodes)[1:3] <- c("id", "parent", "level")
  }

  # Load ids from file if required
  if(!is.null(spp_file)){
    spp_df <- as.data.frame(
      data.table::fread(spp_file, sep = "\t",
                        col.names = c("taxID", "species_count"),
                        verbose = FALSE))
  } else {
    names(spp_df)[1:2] <- c("taxID", "species_count")
  }

  # Return data.table options to previous state
  options(odt)

  # ===========================================================================
  # Process taxonomy counts

  # Filter IDs that aren't part of NCBI notation.
  idx <- which(!(ids_df$taxID %in% nodes$id))
  if (length(idx) > 0) {
    warning("The following IDs are not found in the NCBI taxonomy files and will be ignored:\n",
            paste(ids_df$taxID[idx], collapse = "\n"))
    ids_df <- ids_df[-idx, ]
  }

  # Remove duplicates
  idx <- which(duplicated(ids_df$taxID))
  if (length(idx) > 0) {
    warning("Some IDs are duplicated, only the first occurrence will be used.")
    ids_df <- ids_df[-idx, ]
  }

  # Count how often each taxonomy ID occurs in the dataset
  countIDs            <- rep.int(0, nrow(nodes))
  names(countIDs)     <- nodes$id
  nodes_parent        <- nodes$parent
  names(nodes_parent) <- nodes$id
  searchIDs           <- ids_df$taxID
  countIDs[as.character(searchIDs)] <- 1

  if(verbose) message("Counting taxonomy IDs")
  while (length(searchIDs) > 0) {
    searchIDs <- nodes_parent[as.character(searchIDs)]
    parentage <- table(searchIDs)
    countIDs[names(parentage)] <- countIDs[names(parentage)] + parentage
    searchIDs <- searchIDs[searchIDs != 1]
    if(verbose) cat("\r", sprintf("--> %06d taxIDs still being counted...", length(searchIDs)))
  }
  if(verbose) cat("\r", paste(rep(" ", 50), collapse = ""))

  countIDs <- countIDs[countIDs > 0]

  # ===========================================================================
  # Process species counts

  if(verbose) message("\rCounting species")
  # Filter only the Spp counts for taxons listed in countIDs
  spp_df   <- spp_df[spp_df$taxID %in% names(countIDs), ]
  countSpp <- spp_df$species_count
  names(countSpp) <- spp_df$taxID

  # Species TaxIDs may have a count of zero (depending on how counting is done),
  # so increment by 1
  countSpp[countSpp == 0] <- 1

  # Reduces the search size of nodes to the relevant input taxIDs and their
  # related ancestors/offspring. Can greatly reduce search/running time for
  # the remaining functions.
  nodes <- nodes[nodes$id %in% as.numeric(names(countIDs)), ]

  taxlist <- list(nodes    = nodes,
                  ids_df   = ids_df,
                  countIDs = countIDs,
                  spp_df   = spp_df,
                  countSpp = countSpp)

  class(taxlist) <- c(class(taxlist), "taxonsampling")

  return(taxlist)
}