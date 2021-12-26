evaluate_TS <- function(taxlist) {
  # Collect metrics about the taxon diversity of a sample (outputIDs).
  #
  # Args:
  #   taxlist: _taxonsampling_ list object returned by [run_TS()]
  # Returns:
  #   listTaxon: (char) children taxa selected at least once per level.

  selectedIDs <-  get_taxonomy_counts(ids_df = data.frame(taxID = taxlist$outputIDs,
                                                          seqID = NA),
                                      nodes = taxlist$nodes)
  listTaxon   <- list()
  listBias    <- list()

  # Find the sub-taxa (children nodes) of the current taxon
  children <- 1
  while (length(children) > 0) {
    taxon    <- as.integer(children)
    children <- selectedIDs$nodes$id[(selectedIDs$nodes$parent %in% taxon) & !(selectedIDs$nodes$id %in% taxon)]
    children <- intersect(children, names(selectedIDs$countIDs))

    selectedChildren <- intersect(children, names(selectedIDs$countIDs))
    listTaxon        <- c(listTaxon, list(selectedChildren))
  }

  # Metric: bias
  # Description: how much the data differs from the uniform distribution among
  #              the children sharing the same parent taxon, over the
  #              current taxonomic level.

  #  children <- 1
  #  while (length(children) > 0) {
  #    taxon <- as.integer(children)
  #    children <- nodes$id[is.element(nodes$parent, taxon) &
  #                         !is.element(nodes$id, taxon)]
  #    children <- intersect(children, names(countIDs))
  #
  #    levelBias <- 0
  #    for (parent in taxon) {
  #      parent <- as.character(parent)
  #      countChildren <- countIDs[nodes$parent == parent & nodes$id != parent]
  #      levelBias <- levelBias + sd(countChildren)
  #    }
  #    listBias <- c(listBias, list(levelBias))
  #  }

  return(listTaxon)
}

# Parameters
rand <- c("yes", "no", "after_first_round")
meth <- c("balanced", "diversity")
samp <- c("agnostic", "known_species")
repl <- c(TRUE, FALSE)
pars <- expand.grid(rand=rand, meth=meth, samp=samp, repl=repl,
                    stringsAsFactors = FALSE)
pars <- pars[-which(pars$meth == "balanced" & pars$rand == "after_first_round"), ]

# preprocess stuff:
taxlist <- get_taxonomy_counts(taxonomy_path = "data_files/taxdump/",
                               ids_file      = "data_files/metadata/TaxID2SeqID.txt") %>%
  get_species_counts(spp_file = "data_files/metadata/TaxID2sppCounts.tsv")

# Run test - are there errors?
output <- vector("list", nrow(pars))
for (i in 1:nrow(pars)){
  cat(sprintf("\nTrying %02d/%02d: [%s]", i, nrow(pars), paste(pars[i, ], collapse = ",")))
  output[[i]] <- run_TS(taxlist          = taxlist,
                        taxon            = 1,
                        m                = 100,
                        seq_file         = NULL,
                        out_file         = NULL,
                        method           = pars$meth[i],
                        randomize        = pars$rand[i],
                        replacement      = pars$repl[i],
                        ignoreIDs        = 9598,
                        requireIDs       = 9606,
                        ignoreNonLeafIDs = NULL,
                        sampling         = pars$samp[i],
                        verbose          = FALSE)
}
