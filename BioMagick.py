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
import ntpath
import os
import argparse
import sys

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
	input_files = args.input
	out_fmt = args.outfmt
	if args.outdir:
		os.chdir(args.outdir)  # Set working directory of script to output dir.

	id_results = BioID("./formats.yml").identify(input_files)

	if out_fmt:
		direct_convert(id_results, out_fmt)

# ------------------------------------------------------------
# Generates of dictionary of sequence record objects per file.
# ------------------------------------------------------------
def generate_sequence_objects(id_results):
	seq_objects = {}

	for file_path, identity in id_results.items():
		with open(file_path, 'rU') as input_handle:
			new_seq_objects = SeqIO.parse(input_handle, identity)
			seq_objects[file_path] = new_seq_objects

	return seq_objects


# ------------------------------------------------------------
# Generates of dictionary of sequence record objects per file.
# ------------------------------------------------------------
def direct_convert(id_results, out_fmt):

	for format_type in out_fmt:
		for file_path, identity in id_results.items():

			if sys.platform == "win32":
				file_name = ntpath.basename(file_path).split('.')[0]
			else:
				file_name = os.path.basename(file_path).split('.')[0]

				file_name = file_name + "." + format_type

			print("Converting " + file_path + " to " + format_type)

			try:
				SeqIO.convert(file_path, identity.lower(), file_name, format_type)
			except ValueError as e:
				print("\nError in conversion of " + file_path + " to " + format_type + ": " + str(e))
				print("Skipping " + file_path + " ...\n")
				continue

if __name__ == '__main__':
	cli()