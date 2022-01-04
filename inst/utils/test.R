library(dplyr)

taxlist <- get_counts(taxonomy_path = "data_files/taxdump/",
                      ids_file      = "data_files/metadata/TaxID2SeqID.txt",
                      spp_file = "data_files/metadata/TaxID2sppCounts_2022-01-01.tsv") %>%
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
