#!/bin/sh

echo "Downloading NCBI Taxonomy"
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

# echo "Creating directory structure"
# mkdir data/validation/taxdump/

echo "Expanding file"
tar -xvzf new_taxdump.tar.gz

echo "Moving file to data/validation/taxdump/"
mv *dmp data/taxdump/

echo "Cleaning install files"
rm new_taxdump.tar.gz

echo "Done"
