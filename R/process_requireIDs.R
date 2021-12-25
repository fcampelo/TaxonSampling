process_requireIDs <- function(taxlist) {
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

    # Remove IDs that have any ignoreIDs as a parent.
    # (shouldn't happen, since we run process_ignoreIDs prior to calling
    # this function - but just to be on the safe side)
    # TODO: make this refer directly to taxlist$ts.params$ignoreIDs
    idx <- which(!(requireIDs %in% names(taxlist$countIDs)))
    if (length(idx) > 0) {
      warning("The following IDs are either children or part of your ignoreIDs and will be ignored:\n",
              paste(requireIDs[idx], collapse = "\n"))
      requireIDs <- requireIDs[-idx]
    }

    if (length(requireIDs) > 0) {
      # TODO: Why generate this list?
      taxlist$ts.process$requireIDs.list <- get_taxonomy_counts(
        ids_df  = data.frame(taxID = requireIDs,
                             seqID = NA),
        nodes   = taxlist$nodes,
        verbose = FALSE)
      taxlist$ts.process$outputIDs <- requireIDs
      taxlist$ts.process$m <- taxlist$ts.process$m - length(requireIDs)

    } else {
      taxlist$ts.process$requireIDs.list <- NULL
    }
    taxlist$ts.params$requireIDs <- requireIDs
  }

  return(taxlist)
}
