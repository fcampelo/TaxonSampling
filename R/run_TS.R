# Wrapper for the Taxon Sampling. Provides support for performance, for
# ignoring specific taxon IDs, and for requiring specific taxon IDs.
#
# Args:
#   taxon: (char/integer) Taxon from which to start sampling children taxa.
#   m: (integer) size of the sample to generate.
#   nodes: (data.frame) pre-processed information about the NCBI taxonomy
#                       structure. Created by getnodes() from the CHNOSZ
#                       package.
#   countIDs: (vector) count of how many taxnomoy IDs belong to each taxon,
#                      created by TS_TaxonomyData().
#   replacement: (char) whether the algorithm allows to repeat IDs in order
#                       to maximize taxonomy diversity and to reach m IDs
#                       in the output.
#   randomize: (char) whether the algorithm will choose IDs randomly or
#                  maintaining a balanced allocation (m_i differing by no
#                  more than 1 if the maximum possible value wasn't reached).
#   method: (char) whether it favors balanced taxa representation
#                  (method = "balance") or maximized taxa representation
#                  (method = "diversity").
#   ignoreIDs: (char) IDs that mustn't appear in the output.
#   requireIDs: (char) IDs that must appear in the output.
#   ignoreNonLeafID: (char) (testing) a non-leaf ID to ignore; won't apply
#                           to leaf nodes and won't exclude its children from
#                           the analysis, unlike ignoreIDs.
#   sampling: (char) whether to sample species in an agnostic manner ("agnostic")
#                    or based on known species diversity ("known_species")
#
# Returns:
#   outputIDs: (vector) vector of IDs with maximized taxonomy balance
#                       or diversity.

run_TS <- function(taxlist, taxon, m, method = "diversity",
                   randomize = "no", replacement = FALSE,
                   ignoreIDs = NULL, requireIDs = NULL,
                   ignoreNonLeafIDs = NULL, sampling = "agnostic") {

  # ===========================================================================
  # Sanity checks
  assertthat::assert_that(is.list(taxlist),
                          all(c("countIDs", "nodes") %in% names(taxlist)),
                          is.character(taxon) || is.numeric(taxon),
                          assertthat::is.count(m),
                          is.character(randomize),
                          randomize %in% c("yes", "no", "after_first_round"),
                          is.logical(replacement),
                          length(replacement) == 1,
                          is.null(ignoreIDs) || is.numeric(ignoreIDs) || is.character(ignoreIDs),
                          is.null(requireIDs) || is.numeric(requireIDs) || is.character(requireIDs),
                          is.null(ignoreNonLeafIDs) || is.numeric(ignoreNonLeafIDs) || is.character(ignoreNonLeafIDs),
                          is.character(sampling),
                          length(sampling) == 1,
                          sampling %in% c("agnostic", "known_species"))

  # Force IDs to integer vectors
  if(is.character(taxon)) taxon <- as.integer(taxon)
  if(is.character(ignoreIDs)) ignoreIDs <- as.integer(ignoreIDs)
  if(is.character(requireIDs)) requireIDs <- as.integer(requireIDs)
  if(is.character(ignoreNonLeafIDs)) ignoreNonLeafIDs <- as.integer(ignoreNonLeafIDs)

  # ===========================================================================

  # Process ignoreIDs and ignoreNonLeafID
  taxlist <- process_ignoreIDs(taxlist, ignoreIDs)
  taxlist <- process_ignoreNonLeafIDs(taxlist, ignoreNonLeafIDs)


  if (!is.null(requireIDs)) {
    requireIDs <- as.integer(requireIDs)

    # Sanity check: only IDs that are present in the input.
    if (!all(is.element(requireIDs, nodes$id))) {
      cat("Warning: the following required IDs are not part of your input",
          "and will be ignored.\n",
          requireIDs[!is.element(requireIDs, nodes$id)], "\n")
      requireIDs <- requireIDs[is.element(requireIDs, nodes$id)]
    }

    # Sanity check: remove IDs that has any ignoreIDs as its parent.
    if (!all(is.element(requireIDs, names(countIDs)))) {
      cat("Warning: the following required IDs are children or part of your",
          "ignored IDs and will be ignored.\n",
          requireIDs[!is.element(requireIDs, names(countIDs))], "\n")
      requireIDs <- requireIDs[is.element(requireIDs, names(countIDs))]
    }

    # Sanity check: unique ID inputs, remove duplicates.
    if (any(duplicated(requireIDs))) {
      cat("Warning: some required IDs are repeated,",
          "using only one instance of each.\n",
          requireIDs[duplicated(requireIDs)], "\n")
      requireIDs <- unique(requireIDs)
    }

    if (length(requireIDs) > 0) {
      requireIDs <- TS_TaxonomyData(requireIDs, nodes)
    } else {
      requireIDs <- NULL
      #      print("Here!")
    }
  }

  # Reduce, once again, the node information to the necessary only,
  # reduces search time.
  Simplify_Nodes(nodes, countIDs)

  # Ensure m <= number of valid ids.
  if (!is.null(idsFile) & isTRUE(file.exists(as.character(idsFile)))) { #parsing variable if file
    ids <- read.table(idsFile, sep = "\t", header = FALSE, comment.char = "")
    ids <- ids[, 1]
  } else { #getting ids if list of IDs
    ids <- as.integer(idsFile)  # assume it's a test.name or back.name, for now
  }

  if (m > length(intersect(ids, names(countIDs)))) {
    m <- length(intersect(ids, names(countIDs)))
  }



  # Call the TS algorithm itself.
  outputIDs <- TS_Algorithm_Recursion(taxon, m, nodes, countIDs, method,
                                      randomize, replacement, ignoreIDs, requireIDs, ignoreNonLeafID, sampling)

  return(outputIDs)
}
