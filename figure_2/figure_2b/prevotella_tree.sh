# download fasta with sequences
wget https://lpsn.dsmz.de/fasta/37185.fasta

# make an alignment
/usr/bin/mafft  --auto --reorder 37185.fasta > 37185.aln

# make the tree
fasttree -nt -gtr -pseudo 1000 37185.aln > 37185.nwk

## Tree was edited and draw using iTOL:
## itol.embl.de/

