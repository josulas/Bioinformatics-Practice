# !usr/conv/bash
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz
gunzip swissprot.gz
/bin/makeblastdb -in swissprot -dbtype prot -out swissprot_db