process_requireIDs <- function(taxlist, ignoreIDs) {
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
    idx <- which(!(requireIDs %in% names(taxlist$countIDs)))
    if (length(idx) > 0) {
      warning("The following IDs are either children or part of your ignoreIDs and will be ignored:\n",
              paste(requireIDs[idx], collapse = "\n"))
      requireIDs <- requireIDs[-idx]
    }

    if (length(requireIDs) > 0) {
      requireIDs <- get_taxonomy_counts(requireIDs, nodes)
    } else {
      requireIDs <- NULL
      #      print("Here!")
    }

    return(requiredIDs)
  }
}
