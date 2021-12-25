# Receives a group of Taxonomy IDs and the desired sample size m.
# Returns a sampled vector with maximized taxonomy diversity.
# Assumes that every input id is unique.
# May return repeated IDs if replacement == TRUE, to increase
# diversity when method == "diversity".
#
# Args:
# list containing fields:
#   taxlist$ts.params$taxon:  (char/integer) Taxon from which to start sampling
#                             children taxa.
#   taxlist$ts.params$m:      (integer) desired sample size.
#   taxlist$nodes:            (data.frame) pre-processed information from NCBI
#                             taxonomy.
#   taxlist$ts.params$countIDs: (numeric) count of taxnomoy IDs belonging to
#                             each taxon, created by get_taxonomy_counts()
#   replacement:              (logical) whether the algorithm allows to repeat
#                             IDs in order to maximize taxonomy diversity and to
#                             reach m IDs in the output.
#   randomize:                (char) whether the algorithm will choose IDs
#                             randomly or maintaining a balanced allocation (m_i
#                             differing by at most 1 if the maximum possible
#                             value wasn't reached).
#   method:                   (char) whether it favors balanced taxa
#                             representation ("balance") or maximized taxa
#                             representation ("diversity").
#   requireIDs:               (char) IDs that must appear in the output.
#
# Returns:
#   updated taxlist containing field:
#   outputIDs:                (vector) vector of IDs with maximized taxonomy
#                             balance or diversity.

ts_recursive <- function(taxlist, verbose = TRUE) {

  taxlist$ts.process$recursion.counter <- taxlist$ts.process$recursion.counter + 1

  if(verbose) cat("\r--> Recursion counter:", taxlist$ts.process$recursion.counter)

  # extract relevant variables
  taxon <- taxlist$ts.process$taxon
  m     <- taxlist$ts.process$m

  # Sanity check
  if (m <= 0) stop("m less or equal than zero during recursion.")

  # First step: Find the sub-taxa (children nodes) of the current taxon
  # that has members in user-provided data
  children <- taxlist$nodes$id[taxlist$nodes$parent == taxon & taxlist$nodes$id != taxon]
  children <- intersect(children, names(taxlist$countIDs))

  # Condition to end recursion
  if (length(children) == 0) {
    if(taxlist$ts.params$replacement){
      return(rep(taxon, m))
    } else {
      return(taxon)
    }
  }

  childrenCount    <- taxlist$countIDs[as.character(children)]
  childrenCountSpp <- taxlist$countSpp[as.character(children)]

  # Sanity check
  # In a few cases, the number of sequences may be greater than the number of
  # species. We saw this happening only in terminal nodes where there's both a
  # species and a subspecies, as our mapping of taxID2species ends at the
  # species level.
  # If that's the case, replace the number of known species by the number of
  # child sequences.
  # TODO: check this with Francisco
  childrenCountSpp[childrenCount > childrenCountSpp] <- childrenCount[childrenCount > childrenCountSpp]


  # For cases when one taxon isn't a leaf node, but is an input ID that should
  # be available to be sampled along with its children.
  # Example of this case is an input with Homo sapiens (9606),
  # H. sapiens neanderthalensis (63221) and H. sapiens ssp. denisova (741158).
  if (sum(childrenCount) < taxlist$countIDs[as.character(taxon)]) {
    childrenCount    <- c(childrenCount, taxlist$countIDs[as.character(taxon)])
    childrenCountSpp <- c(childrenCountSpp, taxlist$countIDs[as.character(taxon)])
    childrenCountSpp[as.character(taxon)] <- 1
    childrenCount[as.character(taxon)]    <- 1
  }

  # Allocation function to use:
  fname <- paste("sample",
                 substr(taxlist$ts.params$sampling, 1, 1),
                 substr(taxlist$ts.params$method, 1, 1),
                 substr(taxlist$ts.params$randomize, 1, 1),
                 sep = "_")

  # Perform sampling at current recursion level
  m_i <- do.call(fname,
                 args = list(m                = m,
                             m_i              = 0 * childrenCount,
                             childrenCount    = childrenCount,
                             childrenCountSpp = childrenCountSpp))

  # TODO: check the whole requireIDs thing.
  # Require specific IDs (and their parents/ancestors) to be present.
  # If any required ID has a lower m_i allocated than needed (informed by
  # requireIDs), reallocate from another ID in m_i.
  #  if (!is.null(requireIDs)) {
  #    req <- is.element(names(m_i), names(requireIDs))
  #    toAdd <- m_i[req] < requireIDs[names(m_i[req])]
  #    if (randomize == "no") {
  #      while (any(toAdd)) {
  #        subtract <- m_i[setdiff(names(m_i), names(toAdd))]
  #        subtract <- sample(names(subtract[subtract == max(subtract)]), 1)
  #        m_i[subtract] <- m_i[subtract] - 1
  #        add <- sample(names(toAdd[toAdd == TRUE]), 1)
  #        m_i[add] <- m_i[add] + 1
  #        toAdd <- m_i[req] < requireIDs[names(m_i[req])]
  #      }
  #    } else if (randomize == "yes") {
  #      while (any(toAdd)) {
  #        subtract <- m_i[setdiff(names(m_i), names(toAdd))]
  #        subtract <- sample(names(subtract[subtract > 0]), 1)
  #        m_i[subtract] <- m_i[subtract] - 1
  #        add <- sample(names(toAdd[toAdd == TRUE]), 1)
  #        m_i[add] <- m_i[add] + 1
  #        toAdd <- m_i[req] < requireIDs[names(m_i[req])]
  #      }
  #    }
  #  }

  outputIDs <- character()
  for (id in names(m_i)) {
    if (m_i[id] > 0){
      if (id == as.character(taxon)){
        outputIDs <- c(outputIDs, id)
      } else {
        taxlist$ts.process$m     <- m_i[id]
        taxlist$ts.process$taxon <- id
        outputIDs <- c(outputIDs, ts_recursive(taxlist))
      }
    }
  }

  return(outputIDs)
}
