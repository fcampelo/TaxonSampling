#path to tabular file linking NCBI taxon IDs to sequence IDs
ids_file <- "data_files/metadata/TaxID2SeqID.txt"

# path to NCBI taxonomy files (execute "install.sh" to automatically download them)
taxonomy_path <- "data_files/taxdump/"

outlist <- get_taxonomy_data(ids_file = ids_file, taxonomy_path = taxonomy_path)



#path to multi-fasta file from where sequences should be sampled
multifasta <- "data_files/fasta/mit_vertebrata.fasta"

#path to file linking NCBI taxon IDs to sequence IDs
knownSppFile <- "data_files/metadata/TaxID2sppCounts.tsv"

#output directory
outDir <- "results/"

#path to output file
outFile <- paste0(outDir, "/output_TS.fasta")



## Parameters for TS

# where should we start looking? (root node where to start algorithm)

#mammalia
root_taxon <- 40674

#vertebrata
#root_taxon <- 7742

#whether/how to randomize (options are "yes", "no" and "after_first_round")
#randomize <- "yes"
#randomize <- "no"
randomize <- "after_first_round"

#sampling procedure (either "agnostic" or "known_species")
sampling <- "agnostic"
#sampling <- "known_species"

#number of sequences to sample
m <- 200

#IDs to be ignored during sampling procedure (either terminal or internal taxons)
ignoreIDs <- NULL

#required IDs to be present in final output file (only terminal taxa - species)
#requireIDs <- c(2026169, 8364, 57393, 241292, 61967)
requireIDs <- NULL


#from here on, users do not need to edit this file
#-----------------------------------------------------------------------------#

source("bin/TaxonSampling/TaxonSampling.R")

# These variables are highly experimental and probably
# will not to work, but future versions of TS are likely to support
# them, so we decided to keep them here for now.

#method to use (either 'diversity' or 'balance').
#"Balance" is heavily outdated, and currently does not support
#distinct sampling strategies. We haven't been using it for a
#while.
method <- "diversity"

# should we ignore the following non-leaf IDs while sampling?
ignoreNonLeafID <- NULL

#allow TS to repeat IDs if needed? ('no' is better to get a higher diversity)
replacement <- "no"

#get all nodes from NCBI taxonomy
print("Parsing NCBI Taxonomy data")
nodes <- suppressMessages(getnodes(taxondir))
print("Done.")

print("Computing sequence counts per taxon")
#count how many sequences per node
countIDs <- TS_TaxonomyData(idsFile, nodes)
print("Done.")

print("Computing species count per taxon")
#count how many known spp per node
countSpp <- TS_SpeciesData(knownSppFile, countIDs)
print("Done.")

print("Simplifying taxonomic space")
#simplify taxonomic space by
#1) removing nodes that are not children of root node
#2) removing nodes that have no sequences to be sampled in current sequence database
nodes <- Simplify_Nodes(nodes, countIDs)
print("Done.")

#list to store output IDs
outputIDs1 <- list()

#selecting user-defined IDs that must be present in output file
if (!is.null(requireIDs)) {
  total_requiredIDs <- length(requireIDs)
  if (!all(is.element(requireIDs, nodes[,1]))) {
    stop("One of the required IDs provided is not present in database")
  }
  outputIDs1 <- intersect(requireIDs, nodes[,1])
  m <- m - total_requiredIDs
}

#no required IDs anymore
requireIDs <- NULL

#creating output directory, if needed
if (dir.exists(outDir)){
} else {
  dir.create(file.path(outDir), recursive=TRUE)
}

print("Running TaxonSampling")
outputIDs2 <- TS_Algorithm(root_taxon, m, nodes, countIDs, method, randomize,
                           replacement, ignoreIDs, requireIDs,
                           ignoreNonLeafID, sampling)
print("Done.")

print("Writing output file")

#Merging user-defined IDs with the ones selectec by TS
outputIDs <- c(outputIDs1, outputIDs2)
#Writing output fasta file
WriteFasta(idsFile, multifasta, outputIDs, outFile)
print("Done.")
print("TaxonSampling analysis finished.")
