#!/usr/bin/env python

# ============================================================================================================
# Created by: Lee Bergstrand & Matt McInnes
# Description: A next generation bioinformatics file format converter and sequence feature extractor.
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
# ============================================================================================================

from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from BioID import BioID
import argparse

# -------------------------------
# Command line interface options.
# -------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', metavar='PATH', nargs='+', help='''
A list of input file paths. If not specified, input is read from stdin''')

parser.add_argument('-o', '--outdir', metavar='PATH', nargs='+', help='''
An output directory for output files. If not specified, output is piped to stdout.''')

parser.add_argument('-f', '--outfmt', metavar='FORMAT', nargs='+', help='''
A List of output file formats.''')

args = parser.parse_args()


# ----------------------------------
# Command line interface Controller.
# ----------------------------------
def cli():
	print(args.input)

if __name__ == '__main__':
	cli()