library(dplyr)

# t0 <- Sys.time()
taxlist <- get_counts(taxonomy_path = "data_files/taxdump/",
                      ids_file      = "data_files/metadata/TaxID2SeqID.txt") %>%
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
