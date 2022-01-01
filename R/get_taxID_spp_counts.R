#' Count number of species under each taxonomy ID
#'
#' Extract the species count for each taxonomy ID in the _nodes_ structure
#' of the NCBI taxonomy.
#'
#' @param taxonomy_path path to folder containing the NCBI taxonomy files
#' (i.e., the extracted contents of _taxdump.zip_, which can be downloaded from
#' <ftp://ftp.ncbi.nih.gov/pub/taxonomy/> or retrieved using
#' [retrieve_NCBI_taxonomy()]).
#' @param spp_file path to file for saving the species count for each taxon ID
#' (saved as a tsv file)
#' @param ncpus number of cores to use for species counting.
#' @param verbose logical: regulates function echoing to console.
#'
#' @return data frame with species counts per taxon ID.
#'
#' @export

get_taxID_spp_counts <- function(taxonomy_path,
                                 spp_file,
                                 ncpus = 1,
                                 verbose = TRUE){

  # ===========================================================================
  # Sanity checks
  assertthat::assert_that(is.character(spp_file),
                          length(spp_file) == 1,
                          is.character(taxonomy_path),
                          length(taxonomy_path) == 1,
                          dir.exists(taxonomy_path),
                          assertthat::is.count(ncpus),
                          is.logical(verbose),
                          length(verbose) == 1)

  # Set up multicore processing.
  if(ncpus > parallel::detectCores() - 1){
    ncpus <- parallel::detectCores() - 1
    message("Note: Only ", parallel::detectCores(), " nodes available. Using ",
            ncpus, " for counting species.")
  }

  if (ncpus > 1){
    if(.Platform$OS.type == "windows"){
      cl <- parallel::makeCluster(ncpus, setup_strategy = "sequential")
      parallel::clusterExport(cl, "nodes")
    } else {
      cl <- ncpus
    }
  } else {
    cl <- 1
  }

  # ===========================================================================

  # read taxonomy file
  nodes <- as.data.frame(
    data.table::fread(paste(taxonomy_path, "nodes.dmp", sep = "/"),
                      sep = "|", strip.white = TRUE,
                      colClasses = c("numeric", "numeric", "character", rep("NULL", 16)),
                      col.names = c("id", "parent", "level"),
                      verbose = FALSE))
  nodes$level <- gsub("\\t", "", nodes$level)

  # Extract species counts
  if(verbose) {
    message("Extracting species counts. This may take a while...")
    x <- table(unlist(
      pbapply::pblapply(nodes$id[nodes$level == "species"],
                        CHNOSZ::allparents,
                        nodes = nodes,
                        cl    = cl)))
  } else {
    if (inherits(cl, "cluster")) {
      x <- table(unlist(
        parallel::parLapply(cl    = cl,
                            X     = nodes$id[nodes$level == "species"],
                            fun   = CHNOSZ::allparents,
                            nodes = nodes)))
    } else {
      x <- table(unlist(
        parallel::mclapply(mc.cores = cl,
                           X        = nodes$id[nodes$level == "species"],
                           FUN      = CHNOSZ::allparents,
                           nodes    = nodes)))
    }
  }

  if (inherits(cl, "cluster")) parallel::stopCluster(cl)

  spp_df <- data.frame(taxID         = as.character(names(x)),
                       species_count = unname(as.numeric(x)))

  # Save to file
  if(!dir.exists(dirname(spp_file))) {
    dir.create(dirname(spp_file), recursive = TRUE)
  }
  utils::write.table(spp_df, file = spp_file, sep = "\t", row.names = FALSE)


  return(spp_df)

}
