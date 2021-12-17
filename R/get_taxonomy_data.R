#' Return the count of each taxon in a taxonomy graph
#'
#' This function evaluates the taxonomy structure from a set of input IDs and
#' returns the count of each taxon ID with non-zero occurrences in the
#' taxonomy graph.
#'
#' @param ids_df two-column data frame with the input taxon IDs in the first
#' column, and the corresponding sequence IDs (corresponding
#' to the ID strings int the multifasta input file, without
#' the ">" line starter) in the second column. Ignored if
#' `ids_file` is not `NULL`.
#' @param ids_file path to a tab-separated file containng two two columns,
#' with the input taxon IDs in the first
#' column, and the corresponding sequence IDs (corresponding
#' to the ID strings int the multifasta input file, without
#' the ">" line starter) in the second column.
#' @param taxonomy_path path to folder containing the NCBI taxonomy files
#' (i.e., the contents of _taxdump.zip_, which can be downloaded from
#' <ftp://ftp.ncbi.nih.gov/pub/taxonomy/> or retrieved using
#' [retrieve_NCBI_taxonomy()]).
#' @param verbose logical: regulates function echoing to console.
#'
#' @return list object containing:
#' \itemize{
#'     \item `$nodes`: data.frame containing the pre-processed information about
#'     the NCBI taxonomy structure. This is extracted internally by calling
#'     [CHNOSZ::getnodes()].
#'     \item `$countIDs`: numeric vector with the counts of the number of
#'     taxnomoy IDs belonging to each taxon.
#' }
#'
#' @export

get_taxonomy_data <- function(ids_df = NULL, ids_file = NULL, taxonomy_path,
                              verbose = TRUE) {

  # ===========================================================================
  # Sanity checks
  assertthat::assert_that(is.null(ids_df) || is.data.frame(ids_df),
                          is.null(ids_file) ||
                            (is.character(ids_file) && length(ids_file) == 1),
                          (is.null(ids_df) + is.null(ids_file)) < 2,
                          is.character(taxonomy_path),
                          length(taxonomy_path) == 1,
                          dir.exists(taxonomy_path),
                          is.logical(verbose),
                          length(verbose) == 1)

  # Load ids from file if required
  if(!is.null(ids_file) && !file.exists(ids_file)){
    stop("File ", ids_file, " not found.")
  } else {
    ids_df <- utils::read.table(ids_file, sep = "\t", header = FALSE)
  }

  # Assert data frame size
  assertthat::assert_that(nrow(ids_df) > 0 && ncol(ids_df) == 2)

  # ===========================================================================

  # Extract all nodes from NCBI taxonomy
  if(verbose) message("Parsing NCBI Taxonomy data")
  nodes <- suppressMessages(CHNOSZ::getnodes(taxonomy_path))

  # Filter IDs that aren't part of NCBI notation.
  idx <- which(!(ids_df[, 1] %in% nodes$id))
  if (length(idx) > 0) {
    warning("The following IDs are not found in the NCBI taxonomy files and will be ignored:\n",
            paste(ids_df[idx, 1], collapse = "\n"))
    ids_df <- ids_df[-idx, ]
  }

  # Remove duplicates
  idx <- which(duplicated(ids_df[, 1]))
  if (length(idx) > 0) {
    warning("Some IDs are duplicated, only the first occurrence will be used.")
    ids_df <- ids_df[-idx, ]
  }

  # Count how often each taxonomy ID occurs in the dataset
  countIDs            <- rep.int(0, nrow(nodes))
  names(countIDs)     <- nodes$id
  nodes_parent        <- nodes$parent
  names(nodes_parent) <- nodes$id
  searchIDs           <- ids_df[, 1]
  countIDs[as.character(searchIDs)] <- 1

  if(verbose) message("Counting taxonomy IDs")
  while (length(searchIDs) > 0) {
    if(verbose) cat(".")
    searchIDs <- nodes_parent[as.character(searchIDs)]
    parentage <- table(searchIDs)
    countIDs[names(parentage)] <- countIDs[names(parentage)] + parentage
    searchIDs <- searchIDs[searchIDs != 1]
  }

  countIDs <- countIDs[countIDs > 0]
  return(list(nodes    = nodes,
              countIDs = countIDs))
}
