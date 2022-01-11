#' Generate a summary dataframe from the result of run_TS
#'
#' Returns a summary data frame containing relevant counts for each taxID
#' associated with the sequences passed to [get_counts()]
#' (parameters `ids_df` or `ids_file`) in the pre-processing stage before the
#' call to [run_TS()].
#'
#' @param taxlist list object of class _taxonsampling_, returned by
#' [run_TS()]
#'
#' @return Data frame containing fields `$taxID`, `$level` (taxonomic level),
#' `$name` (scientific name), `species_count` (total number of species or leaf
#' nodes under that taxonomy ID), `$sample_count` (number of sampled species or
#' leaf nodes under that taxonomy ID) and `$seq_count` (number of available
#' sequences under that taxonomy ID, passed in the call to [get_counts()]).
#'
#' @export
#'

get_summary_df <- function(taxlist){
  # Get ancestry list of (filtered) nodes
  ancestry <- get_taxID_spp_counts(nodes = taxlist$nodes,
                                   start_from_species = FALSE,
                                   what = "ancestry",
                                   verbose = FALSE)

  # Isolate only ancestries of taxIDs with sequences
  idx <- which(sapply(ancestry,
                      function(x, tax){(tax %in% x)},
                      tax = taxlist$ts.params$taxon))

  X <- table(unlist(ancestry[idx]))
  X <- data.table::data.table(taxID = as.integer(names(X)),
                              seq_count = as.integer(X))
  # Remove taxIDs above the root taxon
  torm <- CHNOSZ::allparents(taxlist$ts.params$taxon, nodes = taxlist$nodes)
  X <- X[!(X$taxID %in% torm), ]

  # Isolate only ancestries of sampled taxIDs
  idx <- which(sapply(ancestry,
                      function(x, ids){any(ids %in% x)},
                      ids = taxlist$outputIDs))

  Y <- table(unlist(ancestry[idx]))
  Y <- data.table::data.table(taxID = as.integer(names(Y)),
                              sample_count = as.integer(Y))

  # Prepare nodes for joining
  nodes <- taxlist$nodes[, -2]
  names(nodes)[1] <- "taxID"

  # Join all info
  X <- Y[X, on = "taxID"]
  X <- data.table::as.data.table(taxlist$spp_df)[X, on = "taxID"]
  X <- nodes[X, on = "taxID"]
  X$sample_count <- ifelse(is.na(X$sample_count), 0, X$sample_count)

  return(as.data.frame(X[, c(1:4,6,5)]))
}
