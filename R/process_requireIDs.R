process_requireIDs <- function(taxlist) {

  taxlist$ts.process$requireIDs.list <- NULL
  taxlist$ts.process$outputIDs       <- character()

  requireIDs <- taxlist$ts.params$requireIDs
  if (!is.null(requireIDs)) {
    # Filter IDs that aren't part of NCBI notation.
    idx <- which(!(requireIDs %in% taxlist$nodes$id))
    if (length(idx) > 0) {
      warning("The following IDs are not found in the NCBI taxonomy files and will be ignored:\n",
              paste(requireIDs[idx], collapse = "\n"))
      requireIDs <- requireIDs[-idx, ]
    }

    # Remove duplicates
    idx <- which(duplicated(requireIDs))
    if (length(idx) > 0) {
      warning("Some IDs are duplicated, only the first occurrence will be used.")
      requireIDs <- requireIDs[-idx]
    }

    if (length(requireIDs) > 0) {
      # TODO: Why generate this list?
      # taxlist$ts.process$requireIDs.list <- get_taxonomy_counts(
      #   ids_df  = data.frame(taxID = requireIDs,
      #                        seqID = NA),
      #   nodes   = taxlist$nodes,
      #   verbose = FALSE)
      taxlist$ts.process$outputIDs <- requireIDs
      taxlist$ts.process$m <- taxlist$ts.process$m - length(requireIDs)

    }
    taxlist$ts.params$requireIDs <- requireIDs
  }

  return(taxlist)
}
