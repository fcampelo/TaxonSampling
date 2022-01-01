#' Return the count of known species in each taxon
#'
#' This function processes a list of NCBI Taxon IDs and known number of species
#' and determines the count of species under each taxon.
#'
#' @param taxlist a list object of class _taxonsampling_, returned by
#' [get_taxonomy_counts()]).
#' @param spp_df two-column data frame with the input taxon IDs in the first
#' column, and the corresponding number of known species in the second column.
#' Ignored if `ids_file` is not `NULL`.
#' @param spp_file path to a tab-separated file containng two two columns,
#' with the input taxon IDs in the first
#' column, and the corresponding number of known species in the second column.
#' @param verbose logical: regulates function echoing to console.
#'
#' @return Input object `taxlist` updated to contain the additional fields:
#' \itemize{
#'    \item `spp_df`: two-column data frame with the input taxon IDs in the
#'    first column, and the corresponding number of known species in the second
#'    column. Filtered to have only the IDs present in `$countIDs`.
#'    \item `$countSpp`, a numeric vector with the counts of (known) species for each
#' taxon from `taxlist$countIDs`.
#' }
#'
#'
#' @export

get_species_counts <- function(taxlist,
                               spp_df   = NULL,
                               spp_file = NULL,
                               verbose  = TRUE) {
  # ===========================================================================
  # Sanity checks
  assertthat::assert_that(is.list(taxlist),
                          "countIDs" %in% names(taxlist),
                          is.null(spp_df) || (
                            is.data.frame(spp_df) &&
                              ncol(spp_df) >= 2 &&
                              nrow(spp_df) > 0),
                          is.null(spp_file) ||
                            (is.character(spp_file) &&
                               length(spp_file) == 1 &&
                               file.exists(spp_file)),
                          is.null(spp_df) + is.null(spp_file) < 2,
                          is.logical(verbose),
                          length(verbose) == 1)

  # ===========================================================================

  if(!is.null(spp_file)){
    # Load ids from file if available
    spp_df <- as.data.frame(
      data.table::fread(spp_file, sep = "\t",
                        col.names = c("taxID", "species_count"),
                        verbose = FALSE))

  } else {
    # Get ids from df if passed
    names(spp_df)[1:2] <- c("taxID", "species_count")

  }

  # ===========================================================================

  # Filter only the Spp counts for taxons listed in taxlist$countIDs
  spp_df   <- spp_df[spp_df$taxID %in% names(taxlist$countIDs), ]
  countSpp <- spp_df$species_count
  names(countSpp) <- spp_df$taxID

  # Species TaxIDs have zero child nodes, so we increment by 1
  countSpp[countSpp == 0] <- 1

  # Update taxlist and return
  taxlist$spp_df   <- spp_df
  taxlist$countSpp <- countSpp

  # Reduces the search size of nodes to the relevant input taxIDs and their
  # related ancestors/offsprings. Can greatly reduce search/running time for
  # the remaining functions.
  taxlist$nodes <- taxlist$nodes[taxlist$nodes$id %in% as.numeric(names(taxlist$countIDs)), ]

  class(taxlist) <- unique(c(class(taxlist), "taxonsampling"))

  return(taxlist)

}
