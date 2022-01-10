spp_file      <- "data_files/metadata/TaxID2sppCounts.tsv"
taxonomy_path <- "data_files/taxdump/"
ids_file      <- "data_files/metadata/TaxID2SeqID.txt"
seq_file      <- "data_files/fasta/mit_vertebrata.fasta"

## Run once to save species counts to file spp_file
# t0 <- Sys.time()
# spp_df <- get_taxID_spp_counts(taxonomy_path = taxonomy_path,
#                                start_from_species = FALSE,
#                                out_file = NULL,
#                                what = "spp_counts")
# dt <- Sys.time() - t0
# message("\nElapsed time: ", dt, " ", units(dt))

taxlist <- get_counts(taxonomy_path = taxonomy_path,
                      ids_file      = ids_file,
                      # sff_df        = spp_df,
                      spp_file      = spp_file)
taxlist <- run_TS(taxlist,
                  taxon            = 40674,
                  m                = 200,
                  seq_file         = seq_file,
                  out_file         = NULL,
                  method           = "diversity",
                  randomize        = "yes",
                  replacement      = FALSE,
                  ignoreIDs        = NULL,
                  requireIDs       = NULL,
                  sampling         = "agnostic")

#summary(taxlist)
#plot(taxlist)
