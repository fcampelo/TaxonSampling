library(dplyr)

ncpus = parallel::detectCores() - 2

taxlist <- get_taxonomy_counts(taxonomy_path = "data_files/taxdump/",
                               ids_file      = "data_files/metadata/TaxID2SeqID.txt") %>%
  get_species_counts(ncpus = ncpus) %>%
  run_TS(taxon            = 40674,                                     # root taxon: mammalia
         m                = 200,                                       # number of sequences to sample
         seq_file         = "data_files/fasta/mit_vertebrata.fasta",   # multi-fasta file from where sequences should be sampled
         out_file         = "results/output_TS.fasta",                 # output file for saving results
         method           = "diversity",                               # sampling priority ("diversity" or "balanced")
         randomize        = "yes",                       # randomization strategy
         replacement      = FALSE,                                     # replacement mode
         ignoreIDs        = 10090,
         requireIDs       = c(9443, 9263, 10091, 9606),         # must have: primates, marsupials, Mus musculus, Mus musculus castaneus, H. Sapiens
         sampling         = "known_species")                                # Sampling strategy

summary(taxlist)
plot(taxlist)
