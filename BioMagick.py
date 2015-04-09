#!/usr/bin/env python

# ============================================================================================================
# Created by: Lee Bergstrand & Matt McInnes
# Description: A next generation bioinformatics file format converter and sequence feature extractor.
# Requirements: This script requires the Biopython module: http://biopython.org/wiki/Download
# ============================================================================================================

from __future__ import print_function
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio import AlignIO
from Bio import Phylo
from BioID import BioID
from BioID import BioIDFormat
import ntpath
import os
import argparse
import sys
import yaml


class BioMagickFormat(BioIDFormat):
	def __init__(self, name, extension, bioclass):
		super(BioMagickFormat, self).__init__(name, None, None, [])
		self.extension = extension
		self.bioclass = bioclass


# ----------------------------------
# Command line interface Controller.
# ----------------------------------
def main(args):
	input_files = args.input
	out_fmt = args.outfmt
	out_dir = args.outdir[0] if args.outdir else "."

	if out_fmt is None or out_fmt == []:
		print("Error: at laest 1 output format is needed")
		exit(1)

	if args.alphabet is not None:
		alphabet = args.alphabet[0]
		if alphabet == 'ambigdna':
			alphabet = IUPAC.IUPACAmbiguousDNA()
		elif alphabet == "unambigdna":
			alphabet = IUPAC.IUPACUnambiguousDNA()
		elif alphabet == "exdna":
			alphabet = IUPAC.ExtendedIUPACDNA()
		elif alphabet == 'ambigrna':
			alphabet = IUPAC.IUPACAmbiguousRNA()
		elif alphabet == "unambigrna":
			alphabet = IUPAC.IUPACUnambiguousRNA()
		elif alphabet == 'prot':
			alphabet = IUPAC.IUPACProtein()
		elif alphabet == "exprot":
			alphabet = IUPAC.ExtendedIUPACProtein()
		else:
			print("Error: %s is not a valid alphabet" % alphabet)
			print("Valid alphabets: ambigdna, unambigdna, exdna, ambigrna, unambigrna, prot, exprot")
			exit(1)
	else:
		alphabet = ""

	# Load and parse YAML format export settings
	with open("BioMagickFormatInfo.yml", "rU") as settings_file:
		contents = settings_file.read()

	settings = {}
	for setting in yaml.safe_load(contents):
		settings[setting["name"]] = BioMagickFormat(setting["name"], setting["extension"], setting["bioclass"])

	if input_files is None or input_files == []:
		id_results = BioID("./BioIDFormatInfo.yml").identify(input_files)
	else:
		if sys.stdin.isatty():
			print("Error: you must either specify an input file or pipe some data to stdin")
			exit(1)
		id_results = BioID("./BioIDFormatInfo.yml").identify(sys.stdin.read())

	direct_convert(settings, id_results, out_dir, out_fmt, alphabet)


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
def direct_convert(settings, id_results, out_path, out_formats, alphabet):
	for out_format in out_formats:
		for in_path, in_format in id_results.items():
			out_file = out_path
			if sys.platform == "win32":
				if out_file[-1] != "/":
					out_file += "/"

				out_file += ntpath.basename(in_path).split('.')[0]
			else:
				if out_file[-1] != "\\":
					out_file += "\\"

				out_file += os.path.basename(in_path).split('.')[0]

			out_extension = settings[out_format].extension
			out_file = out_file + "." + out_extension
			print("Converting %s file %s to %s file %s" % (in_format, in_path, out_format, out_file))

			try:
				format_setting = settings[id_results[in_path]]
				if format_setting.bioclass == "seq":
					SeqIO.convert(in_path, in_format.lower(), out_file, out_format, alphabet)
				elif format_setting.bioclass == "phylo":
					Phylo.convert(in_path, in_format.lower(), out_file, out_format)
				elif format_setting.bioclass == "align":
					AlignIO.convert(in_path, in_format.lower(), out_file, out_format)
				else:
					print("Error: invalid BioPython conversion class: %s" % format_setting.bioclass)
					exit(1)
			except ValueError as e:
				print("\nError in conversion of " + in_path + " to " + out_format + ": " + str(e))
				print("Skipping " + in_path + " ...\n")
				continue

if __name__ == '__main__':
	# -------------------------------
	# Command line interface options.
	# -------------------------------
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', metavar='INPATH', nargs='+', help='''
	A list of input file paths. If not specified, input is read from stdin''')

	parser.add_argument('-o', '--outdir', metavar='OUTPATH', nargs=1, help='''
	An output directory for output files. If not specified, output is piped to stdout.''')

	parser.add_argument('-f', '--outfmt', metavar='OUTFORMAT', nargs='+', help='''
	A List of output file formats.''')

	parser.add_argument('-a', '--alphabet', metavar='ALPHA', nargs=1, help='''
	The alphabet to use for conversion (ambigdna, unambigdna, exdna, ambigrna, unambigrna, prot, exprot).''')

	cli_args = parser.parse_args()
	cli_args.input = cli_args.input[0].split(",")
	cli_args.outfmt = cli_args.outfmt[0].split(",")
	main(cli_args)