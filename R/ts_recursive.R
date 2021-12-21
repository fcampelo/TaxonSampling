# An algorithm that receives a group of Taxonomy IDs and the size m of the
# sample to obtain from them. Returns a vector with a maximized taxonomy
# diversity. Assumes that every input id is unique.
# If replacement = "yes", then it may return repeated IDs if it increases
# the taxonomy diversity when method = "diversity".
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
#   requireIDs: (char) IDs that must appear in the output.
# Returns:
#   outputIDs: (vector) vector of IDs with maximized taxonomy balance
#                       or diversity.
ts_recursive <- function(taxon, m, taxlist,
                         method, randomize, replacement,
                         ignoreIDs, requireIDs, ignoreNonLeafID,
                         sampling) {


  # Sanity check
  if (m <= 0) error("m less or equal than zero during recursion.")

  # First step: Find the sub-taxa (children nodes) of the current taxon
  # that has members in user-provided data
  children <- taxlist$nodes$id[taxlist$nodes$parent == taxon & taxlist$nodes$id != taxon]
  children <- intersect(children, names(taxlist$countIDs))

  # Condition to end recursion
  if (length(children) == 0) {
    if (replacement == "no") {
      return(taxon)
    } else {
      return(rep(taxon, m))
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
  childrenCountSpp[childrenCount > childrenCountSpp] <- childrenCount[childrenCount > childrenCountSpp]


  # For cases when one taxon isn't a leaf node, but is an input ID that should
  # be available to sample along its children.
  # Example of this case is an input with Homo sapiens (9606),
  # H. sapiens neanderthalensis (63221) and H. sapiens ssp. denisova (741158).
  if (sum(childrenCount) < taxlist$countIDs[as.character(taxon)]) {
    childrenCount    <- c(childrenCount, countIDs[as.character(taxon)])
    childrenCountSpp <- c(childrenCountSpp, taxlist$countIDs[as.character(taxon)])
    childrenCountSpp[as.character(taxon)] <- 1
    childrenCount[as.character(taxon)]    <- 1
  }

  # Allocating m_i to the children taxa at last.
  # Diversity mode:
  #   Normal: allocate ensuring that any two taxa differ by at most 1 if
  #           m_i < n_i (or m_i < childrenCount), but allow more if one taxon
  #           has m_i == n_i.
  #   Randomize: if we don't have the m fully distributed over m_i (children),
  #              choose a random child from childrenCount that still has taxa
  #              available to choose.
  # Balance mode:
  #   Normal: ensure two taxa allocation (m_i) differ at most by 1.
  #   Randomize: fully random child sampling, uniform distribution among taxa.
  m_i        <- rep.int(0, length(childrenCount))
  names(m_i) <- names(childrenCount)

  if (sampling == "agnostic") {
    if (method == "diversity"){
      if(randomize == "no") {
        while (m > 0 & length(childrenCount[childrenCount > m_i]) <= m) {
          child <- names(childrenCount[childrenCount > m_i])
          m_i[child] <- m_i[child] + 1
          m <- m - length(child)
        }
        child <- sample(names(childrenCount[childrenCount > m_i]), m)
        m_i[child] <- m_i[child] + 1

      } else if (randomize == "after_first_round") {
        first_round = 0
        # if there's enough to sample at least one entry per lineage,
        # do it to increase diversity
        # TODO: check if this while loop makes sense (it will only happen 0 or 1 times due to first_round update)
        while (m > 0 & length(childrenCount[childrenCount > m_i]) <= m & first_round == 0) {
          child       <- names(childrenCount[childrenCount > m_i])
          m_i[child]  <- m_i[child] + 1
          m           <- m - length(child)
          first_round <- 1
        }
        while (m > 0 & length(childrenCount[childrenCount > m_i]) > 0) {
          child      <- sample(names(childrenCount[childrenCount > m_i]), 1)
          m_i[child] <- m_i[child] + 1
          m          <- m - 1
        }

      } else if (randomize == "yes") {
        while (m > 0 & length(childrenCount[childrenCount > m_i]) > 0) {
          child      <- sample(names(childrenCount[childrenCount > m_i]), 1)
          m_i[child] <- m_i[child] + 1
          m          <- m - 1
        }
      }
    } else if (method == "balance"){
      if(randomize == "no") {
        m_i                  <- m_i + floor(m / length(names(childrenCount)))
        sampledChildren      <- sample(names(childrenCount), m - sum(m_i))
        m_i[sampledChildren] <- m_i[sampledChildren] + 1

      } else if (randomize == "yes") {
        m_i <- table(sample(names(childrenCount), m, replace = TRUE))
      }
      # TODO: is there a randomize = "after_first_round" option under method == "balance"?

    }

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
  } else if (sampling == "known_species") {
    if (method == "diversity"){
      if (randomize == "no") {
        while (m > 0 & length(childrenCount[childrenCount > m_i]) <= m) {
          child      <- names(childrenCount[childrenCount > m_i])
          m_i[child] <- m_i[child] + 1
          m          <- m - length(child)
        }
        # If there are still sequences to sample, get one using the taxon
        # diversity as probability values
        if (m > 0) {
          child      <- sample(names(childrenCount[childrenCount > m_i]),
                               size = m,
                               prob = childrenCountSpp[childrenCount > m_i] / sum(childrenCountSpp[childrenCount > m_i]))
          m_i[child] <- m_i[child] + 1
        }

      } else if (randomize == "after_first_round") {
        first_round = 0
        # TODO: check if this while loop makes sense (it will only happen 0 or 1 times due to first_round update)
        while (m > 0 & length(childrenCount[childrenCount > m_i]) <= m & first_round == 0) {
          child       <- names(childrenCount[childrenCount > m_i])
          m_i[child]  <- m_i[child] + 1
          m           <- m - length(child)
          first_round <-  first_round + 1
        }
        while (m > 0 & length(childrenCount[childrenCount > m_i]) > 0) {
          child      <- sample(names(childrenCount[childrenCount > m_i]),
                               size = 1,
                               prob = childrenCountSpp[childrenCount > m_i] / sum(childrenCountSpp[childrenCount > m_i]))
          m_i[child] <- m_i[child] + 1
          m          <- m - 1
        }

      } else if (randomize == "yes") {
        while (m > 0 & length(childrenCount[childrenCount > m_i]) > 0) {
          tmp        <- childrenCountSpp[childrenCount > m_i] / sum(childrenCountSpp[childrenCount > m_i])
          child      <- sample(names(childrenCount[childrenCount > m_i]),
                               size = 1,
                               prob = childrenCountSpp[childrenCount > m_i] / sum(childrenCountSpp[childrenCount > m_i]))
          m_i[child] <- m_i[child] + 1
          m <- m - 1
        }
      }

    } else if (method == "balance") {
      if(randomize == "no") {
        m_i                  <- m_i + floor(m / length(names(childrenCount)))
        sampledChildren      <- sample(names(childrenCount),
                                       size = m - sum(m_i),
                                       prob = childrenCountSpp/sum(childrenCountSpp))
        m_i[sampledChildren] <- m_i[sampledChildren] + 1

      } else if (randomize == "yes") {
        m_i <- table(sample(names(childrenCount),
                            size    = m,
                            replace = TRUE,
                            prob    = childrenCountSpp / sum(childrenCountSpp)))
      }
    }

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
  } else {
    stop("Parameter sampling must be either agnostic or known_species")
  }

  outputIDs <- character(0)
  for (id in names(m_i)) {
    if (m_i[id] == 0) {
      next
    } else if (id == as.character(taxon)) {
      outputIDs <- c(outputIDs, id)
    } else {
      outputIDs <- c(outputIDs,
                     ts_recursive(taxon = id, m = m_i[id], taxlist = taxlist,
                                  method = method, randomize = randomize,
                                  replacement = replacement,
                                  requireIDs = requireIDs,
                                  sampling = sampling))
    }
  }

  return(outputIDs)
}
