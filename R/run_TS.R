#' Run Taxon Sampling
#'
#' Run the TaxonSampling method to return a sample of taxonomic IDs according
#' to the desired balance / diversity.
#'
#' @param taxlist list object returned by [get_species_count()]
#' @param taxon Taxon ID from which to start sampling children taxa (single
#' character or integer value)
#' @param m desired sample size
#' @param method sampling method to use. Accepts "balance" (favors balanced taxa
#' representation) or "diversity" (favors maximized taxa representation)
#' @param randomize randomization strategy: should the algorithm choose IDs
#' randomly ("yes"), maintaining a balanced allocation ("no"), or with a
#' balanced allocation at the top taxonomic level and randomized afterwards
#' ("after_first_round")?
#' @param replacement logical flag: should the algorithm allow repeated IDs in
#' the output (if needed to reach m IDs in the output with maximized taxonomy
#' diversity).
#' @param ignoreIDs vector (character or integer) of IDs that must not appear in
#' the output.
#' @param requireIDs vector (character or integer) of IDs that must appear
#' in the output. Notice that `ignoreIDs` has precedence over `requireIDs`,
#' i.e., IDs that occur in both will be ignored. `requireIDs` that are children
#' of any `ignoreIDs` will also be ignored.
#' @param ignoreNonLeafIDs non-leaf IDs to ignore; won't apply to leaf nodes
#' and won't exclude children from the sampling (unlike ignoreIDs).
#' @param sampling sampling mode. Accepts "agnostic" (sample species in a
#' diversity-agnostic manner) or "known_species" (sample based on known species
#' diversity).
#'
#' @return
#' Updated `taxlist` containing vector `$outputIDs` of sampled IDs.
#'
#' @export

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
                          is.character(randomize), length(randomize) == 1,
                          is.character(sampling), length(sampling) == 1,
                          is.character(method), length(method) == 1,
                          sampling  %in% c("agnostic", "known_species"),
                          randomize %in% c("yes", "no", "after_first_round"),
                          method    %in% c("diversity", "balanced"),
                          is.logical(replacement), length(replacement) == 1,
                          is.null(ignoreIDs) ||
                            is.numeric(ignoreIDs) ||
                            is.character(ignoreIDs),
                          is.null(requireIDs) ||
                            is.numeric(requireIDs) ||
                            is.character(requireIDs),
                          is.null(ignoreNonLeafIDs) ||
                            is.numeric(ignoreNonLeafIDs) ||
                            is.character(ignoreNonLeafIDs))

  if (randomize == "after_first_round" && method == "balance"){
    stop('Combination of randomize == "after_first_round" and method == "balance" not possible')
  }

  # Add all input parameters to taxlist
  taxlist$ts.params         <- as.list(environment())
  taxlist$ts.params$taxlist <- NULL

  # Force all IDs to integer
  if(is.character(taxon)) taxlist$ts.params$taxon <- as.integer(taxon)
  if(is.character(ignoreIDs)) taxlist$ts.params$ignoreIDs <- as.integer(ignoreIDs)
  if(is.character(requireIDs)) taxlist$ts.params$requireIDs <- as.integer(requireIDs)
  if(is.character(ignoreNonLeafIDs)) taxlist$ts.params$ignoreNonLeafIDs <- as.integer(ignoreNonLeafIDs)

  # ===========================================================================

  # Process ignoreIDs, ignoreNonLeafIDs and requireIDs
  # TODO: Check process_ignoreNonLeafIDs
  taxlist <- process_ignoreIDs(taxlist)
  taxlist <- process_ignoreNonLeafIDs(taxlist)
  taxlist <- process_requireIDs(taxlist)

  # Reduce node information to the necessary only, reduces search time.
  taxlist$nodes <- taxlist$nodes[taxlist$nodes$id %in% as.numeric(names(taxlist$countIDs)), ]

  # Ensure m <= number of valid ids.
  m <- min(m, length(intersect(taxlist$ids_df$taxID, names(taxlist$countIDs))))
  taxlist$ts.process$m <- m

  # Call the TS algorithm itself.
  taxlist <- ts_recursive(taxlist)

  return(taxlist)
}
