#' Count number of species under each taxonomy ID
#'
#' Extract the species count for each taxonomy ID in the _nodes_ structure
#' of the NCBI taxonomy.
#'
#' @param taxonomy_path path to folder containing the NCBI taxonomy files
#' (i.e., the extracted contents of _taxdump.zip_, which can be downloaded from
#' <ftp://ftp.ncbi.nih.gov/pub/taxonomy/> or retrieved using
#' [retrieve_NCBI_taxonomy()]). Ignored if `nodes` is not `NULL`.
#' @param nodes data.frame containing the pre-processed information about
#' the NCBI taxonomy structure. This can be generated, e.g., by using
#' [CHNOSZ::getnodes()].
#' @param spp_file path to file for saving the species count for each taxon ID
#' (saved as a tsv file). Ignored if `NULL`.
#' @param start_from_species logical. If `TRUE` the counting starts at
#' the 'species' level - i.e., taxon IDs below species (e.g., 'subspecies')
#' receive a count of zero and are not included in the other counts. If
#' `FALSE`, counting starts at the terminal nodes of the taxonomy, regardless of
#' taxonomic level).
#' @param verbose logical: regulates function echoing to console.
#'
#' @return data frame with species counts per taxon ID.
#' @importFrom data.table := .N
#'
#' @export

get_taxID_spp_counts <- function(taxonomy_path = NULL,
                                 nodes    = NULL,
                                 spp_file = NULL,
                                 start_from_species = FALSE,
                                 verbose = TRUE){

  # ===========================================================================
  # Sanity checks
  assertthat::assert_that(is.null(taxonomy_path) ||
                            (is.character(taxonomy_path) &&
                               length(taxonomy_path) == 1 &&
                               dir.exists(taxonomy_path)),
                          is.null(nodes) || is.data.frame(nodes),
                          is.null(nodes) + is.null(taxonomy_path) < 2,
                          is.null(spp_file) ||
                            (is.character(spp_file) &&
                               length(spp_file) == 1),
                          is.logical(start_from_species),
                          length(start_from_species) == 1,
                          is.logical(verbose),
                          length(verbose) == 1)

  # ===========================================================================
  # Get nodes table
  if(is.null(nodes)){
    # read taxonomy file
    odt <- options("datatable.showProgress")
    options(datatable.showProgress = FALSE)
    nodes <- data.table::fread(paste(taxonomy_path, "nodes.dmp", sep = "/"),
                               sep = "|", strip.white = TRUE,
                               colClasses = c("numeric", "numeric", "character", rep("NULL", 16)),
                               col.names = c("id", "parent", "level"),
                               verbose = FALSE)
    options(odt)
  } else {
    nodes <- data.table::as.data.table(nodes[, 1:3])
    names(nodes) <- c("id", "parent", "level")
  }

  nodes$level <- gsub("\\t", "", nodes$level)

  # Add fictional super-root node "0" (makes it easier to remove later)
  nodes$parent[nodes$id == 1] <- 0
  nodes <- rbind(nodes[1, ], nodes)
  nodes[1, 1:2] <- 0

  # Initialise data.table variable names to prevent warnings
  id     <- NULL
  parent <- NULL

  if(verbose) message('Extracting parentage structure')

  # initialise parentage table
  ids <- nodes[-1, 1:2]
  if(start_from_species) ids <- ids[nodes$level[-1] == 'species', ]

  while(any(ids$parent != 0)){
    if(verbose) cat(sprintf('\r--> Remaining IDs: %10d', sum(ids$parent != 0)))
    # Rename columns to facilitate processing (our query is always called 'id')
    names(ids) <- c(paste0("o", (ncol(ids)-1):1), "id")
    # join the parents of our (updated) 'id' column
    ids[nodes, on = list(id), parent := list(parent)]
  }

  if(verbose) {
    cat("\r", paste(rep(" ", 50), collapse = ""), "\r")
    message('Counting species/leaf nodes per taxon')
  }

  # Count occurrences of each taxonID in the parentage table
  ids <- unname(unlist(ids[, -c("parent")]))
  ids <- table(ids[ids != 0])

  # Build output
  spp_df <- data.table::data.table(TaxID = names(ids),
                                   species_count = as.numeric(ids))

  # If start_from_species is TRUE, add remaining IDs with a count of zero.
  toadd <- nodes$id[!(nodes$id %in% c(0, spp_df$TaxID))]
  spp_df <- rbind(spp_df,
                  data.table::data.table(TaxID = toadd,
                                         species_count = numeric(length(toadd))))

  # Save to file
  if(!is.null(spp_file)){
    if(!dir.exists(dirname(spp_file))) {
      dir.create(dirname(spp_file), recursive = TRUE)
    }
    utils::write.table(spp_df, file = spp_file, sep = "\t", row.names = FALSE)
  }

  return(as.data.frame(spp_df))

}
