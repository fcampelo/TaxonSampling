process_ignoreNonLeafIDs <- function(taxlist) {
  ignoreNonLeafIDs <- taxlist$ts.params$ignoreNonLeafIDs

  if (!is.null(ignoreNonLeafIDs)) {
    # Filter IDs that aren't part of NCBI notation.
    idx <- which(!(ignoreNonLeafIDs %in% taxlist$nodes$id))
    if (length(idx) > 0) {
      warning("The following IDs are not found in the NCBI taxonomy files and will be ignored:\n",
              paste(ignoreNonLeafIDs[idx], collapse = "\n"))
      ignoreNonLeafIDs <- ignoreNonLeafIDs[-idx, ]
    }

    # Remove duplicates
    idx <- which(duplicated(ignoreNonLeafIDs))
    if (length(idx) > 0) {
      warning("Some IDs are duplicated, only the first occurrence will be used.")
      ignoreNonLeafIDs <- ignoreNonLeafIDs[-idx]
    }

    # Make sure every ignoreNonLeafID is not a leaf node
    # (i.e. not just a consequence of its children).
    for (id in ignoreNonLeafIDs) {
      children <- taxlist$nodes$id[taxlist$nodes$parent == id & taxlist$nodes$id != id]
      children <- intersect(children, names(taxlist$countIDs))
      childSum <- sum(taxlist$countIDs[as.character(children)])
      if (childSum > 0 & childSum < taxlist$countIDs[as.character(id)]) {
        taxlist$countIDs[as.character(id)] <- childSum
      }
    }
    taxlist$ts.params$ignoreNonLeafIDs <- ignoreNonLeafIDs
  }

  return(taxlist)
}
