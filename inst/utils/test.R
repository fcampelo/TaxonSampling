library(dplyr)

spp_file      <- "./data_files/metadata/TaxID2sppCounts.tsv"
taxonomy_path <-  "data_files/taxdump/"
ids_file      <- "data_files/metadata/TaxID2SeqID.txt"

## Run once to save species counts to file spp_file
# spp_df <- get_taxID_spp_counts(taxonomy_path = "data_files/taxdump/",
#                                start_from_species = FALSE,
#                                spp_file = spp_file)

# t0 <- Sys.time()
taxlist <- get_counts(taxonomy_path = taxonomy_path,
                      ids_file      = ids_file,
                      spp_file      = spp_file) %>%
  run_TS(taxon            = 40674,
         m                = 200,
         seq_file         = "data_files/fasta/mit_vertebrata.fasta",
         out_file         = "results/output_TS.fasta",
         method           = "diversity",
         randomize        = "yes",
         replacement      = FALSE,
         ignoreIDs        = 10090,
         requireIDs       = 9606,
         sampling         = "known_species")

# dt <- Sys.time() - t0
# message("\nElapsed time: ", dt, " ", units(dt))

#summary(taxlist)
#plot(taxlist)
