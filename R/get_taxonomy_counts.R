#' Return the count of each taxon in a taxonomy graph
#'
#' This function evaluates the taxonomy structure from a set of input IDs and
#' returns the count of each taxon ID with non-zero occurrences in the
#' taxonomy graph.
#'
#' @param taxonomy_path path to folder containing the NCBI taxonomy files
#' (i.e., the contents of _taxdump.zip_, which can be downloaded from
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
#' @param verbose logical: regulates function echoing to console.
#' @param nodes data.frame containing the pre-processed information about
#' the NCBI taxonomy structure. This is generated either by using
#' [CHNOSZ::getnodes()], or as a result of a previous call to this
#' function. If `nodes` is not `NULL` then `taxonomy_path` is ignored.
#'
#' @return list object containing:
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
#' }
#'
#' @export

get_taxonomy_counts <- function(taxonomy_path = NULL,
                                ids_file      = NULL,
                                ids_df        = NULL,
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
                          is.logical(verbose),
                          length(verbose) == 1,
                          is.null(nodes) || is.data.frame(nodes))

  # Load ids from file if required
  if(!is.null(ids_file)) {
    if(file.exists(ids_file)){
      ids_df <- data.table::fread(ids_file, sep = "\t",
                                  col.names = c("taxID", "seqID"))
    } else {
      stop("File ", ids_file, " not found.")
    }
  }

  # ===========================================================================

  # Extract all nodes from NCBI taxonomy
  if(is.null(nodes)){
    nodes <- as.data.frame(
      data.table::fread("./data_files/taxdump/nodes.dmp",
                        sep = "|", strip.white = TRUE,
                        colClasses = c("numeric", "numeric", rep("NULL", 17)),
                        col.names = c("id", "parent")))
  }

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
  if(verbose) cat("\n")

  countIDs <- countIDs[countIDs > 0]

  # Reduces the search size of nodes to the relevant input taxIDs and their
  # related ancestors/offsprings. Can greatly reduce search/running time for
  # the remaining functions.
  nodes <- nodes[nodes$id %in% as.numeric(names(countIDs)), ]

  taxlist <- list(nodes    = nodes,
                  ids_df   = ids_df,
                  countIDs = countIDs)
  class(taxlist) <- c(class(taxlist), "taxlist")

  return(taxlist)
}
