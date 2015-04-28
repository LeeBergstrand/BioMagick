#!/bin/bash

# Simple 1-to-1
./BioMagick.py -i Human_contigs.embl -f fasta -o out

# 4-to-1
./BioMagick.py -i apaf_phylo.xml,bcl_2_phylo.xml,phylo_distribution.xml,phylo_example.xml -f nexus -o out

# 2-to-2
./BioMagick.py -i TRBG361.embl,AAA03323.embl -f genbank,fasta -o out

# 2-to-4
./BioMagick.py -i funny.sth,simple.sth -f nexus,phylip,fasta,clustal -a ambigdna -o out

# 5-to-2 with 3 jobs/worker processes
./BioMagick.py -i funny.sth,simple.sth,cw02.aln,hedgehog.aln,opuntia.aln -f phylip,fasta -a ambigdna -o out -j 3

# Output to stdout (single-file)
./BioMagick.py -i DD231055_edited.embl -f clustal -s

# Chaining via stdout
./BioMagick.py -i simple.sth -f clustal -s | ./BioMagick.py -f nexus -a ambigdna -o out

# Shell-fu
sed 's/product=.*/product=groovy/' blank_seq.gb | ./BioMagick.py -f fasta -o out

# Just because we can
./BioMagick.py -i AAA03323.embl -f genbank -s | ./BioMagick.py -f embl -o out


