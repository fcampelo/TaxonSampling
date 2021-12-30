library(dplyr)

taxlist <- get_taxonomy_counts(taxonomy_path = "data_files/taxdump/",
                               ids_file      = "data_files/metadata/TaxID2SeqID.txt") %>%
  get_species_counts(spp_file = "data_files/metadata/TaxID2sppCounts.tsv") %>%
  run_TS(taxon            = 40674,                                     # root taxon: mammalia
         m                = 200,                                       # number of sequences to sample
         seq_file         = "data_files/fasta/mit_vertebrata.fasta",   # multi-fasta file from where sequences should be sampled
         out_file         = "results/output_TS.fasta",                 # output file for saving results
         method           = "diversity",                               # sampling priority ("diversity" or "balanced")
         randomize        = "after_first_round",                       # randomization strategy
         replacement      = FALSE,                                     # replacement mode
         ignoreIDs        = c(9443, 9263),                             # ignore: primates and marsupials
         requireIDs       = 10090,                                     # must have: mouse
         sampling         = "agnostic")                                # Sampling strategy


summary(taxlist)
plot(taxlist)
