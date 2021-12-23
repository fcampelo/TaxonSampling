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

  # Force all IDs to integer
  if(is.character(taxon)) taxon <- as.integer(taxon)
  if(is.character(ignoreIDs)) ignoreIDs <- as.integer(ignoreIDs)
  if(is.character(requireIDs)) requireIDs <- as.integer(requireIDs)
  if(is.character(ignoreNonLeafIDs)) ignoreNonLeafIDs <- as.integer(ignoreNonLeafIDs)

  # Add ID lists to taxlist
  taxlist$ts.params$ignoreIDs        <- ignoreIDs
  taxlist$ts.params$requireIDs       <- requireIDs
  taxlist$ts.params$ignoreNonLeafIDs <- ignoreNonLeafIDs

  # ===========================================================================

  # Process ignoreIDs, ignoreNonLeafIDs and requireIDs
  # TODO: Check process_ignoreNonLeafIDs
  taxlist <- process_ignoreIDs(taxlist)
  taxlist <- process_ignoreNonLeafIDs(taxlist)
  taxlist <- process_requireIDs(taxlist)

  # Reduce node information to the necessary only, reduces search time.
  taxlist$nodes <- taxlist$nodes[taxlist$nodes$id %in% as.numeric(names(taxlist$countIDs)), 1:2]

  # Ensure m <= number of valid ids.
  m <- min(m, length(intersect(taxlist$ids_df[, 1], names(taxlist$countIDs))))

  # Call the TS algorithm itself.
  outputIDs <- ts_recursive(taxon = taxon, m = m, taxlist = taxlist,
                            method = method, randomize = randomize,
                            replacement = replacement,
                            ignoreIDs = ignoreIDs,
                            requireIDs = requireIDs,
                            ignoreNonLeafIDs = ignoreNonLeafIDs,
                            sampling = sampling)


  # Assemble output list
  taxlist$run_TS.params <- list(taxon = taxon, m = m, method = method,
                                randomize = randomize,
                                replacement = replacement,
                                ignoreIDs = ignoreIDs, requireIDs = requireIDs,
                                ignoreNonLeafIDs = ignoreNonLeafIDs,
                                sampling = sampling)

  taxlist$outputIDs <- outputIDs

  return(taxlist)
}
