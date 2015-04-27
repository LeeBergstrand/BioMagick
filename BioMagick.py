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
import multiprocessing
import functools


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
		sys.exit(1)

	if args.stdout is True:
		if len(input_files) > 1:
			print("Error: outputting to stdout is only possible with single-file conversions")
			sys.exit(1)
		else:
			out_dir = None  # Indicate to use stdout

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
			sys.exit(1)
	else:
		alphabet = None

	# Load and parse YAML format export settings
	with open("BioMagickFormatInfo.yml", "rU") as settings_file:
		contents = settings_file.read()

	settings = {}
	for setting in yaml.safe_load(contents):
		settings[setting["name"]] = BioMagickFormat(setting["name"], setting["extension"], setting["bioclass"])

	if input_files is None or input_files == []:
		# Convert single file from stdin
		if sys.stdin.isatty():
			print("Error: you must either specify an input file or pipe some data to stdin")
			sys.exit(1)

		id_results = BioID("./BioIDFormatInfo.yml").identify(sys.stdin.read())
		direct_convert(settings, id_results, out_dir, out_fmt, alphabet)
	elif len(input_files) == 1:
		id_results = BioID("./BioIDFormatInfo.yml").identify(input_files)
		direct_convert(settings, id_results, out_dir, out_fmt, alphabet)
	else:
		if args.JOBS:
			process_count = args.JOBS if args.JOBS >= len(input_files) else len(input_files)
		else:
			process_count = multiprocessing.cpu_count() if multiprocessing.cpu_count() > len(input_files) else len(input_files)

		pool = multiprocessing.Pool(processes=process_count)
		pool.map(functools.partial(do_conversion, format_settings=settings, output_path=out_dir, output_formats=out_fmt,
		                           input_alphabet=alphabet), input_files)


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


# -------------------------------------------------------------------------------
# A wrapper for the convert function which provides it with settings information.
# -------------------------------------------------------------------------------
def do_conversion(input_file, format_settings, output_path, output_formats, input_alphabet):
	id_results = BioID("./BioIDFormatInfo.yml").identify([input_file])
	direct_convert(format_settings, id_results, output_path, output_formats, input_alphabet)


# ---------------------------------------------------------------------------------------------
# Converts between bioinformatic formats using SeqIO's, AlignIO's and Phylo's convert function.
# ---------------------------------------------------------------------------------------------
def direct_convert(settings, id_results, out_path, out_formats, alphabet):
	if out_path is None:
		out_file = "./conv.tmp"
		in_path, in_format = list(id_results.items())[0]
		out_format = out_formats[0]

		try:
			format_setting = settings[in_format]
			if format_setting.bioclass == "seq":
				SeqIO.convert(in_path, in_format.lower(), out_file, out_format, alphabet)
			elif format_setting.bioclass == "phylo":
				Phylo.convert(in_path, in_format.lower(), out_file, out_format)
			elif format_setting.bioclass == "align":
				AlignIO.convert(in_path, in_format.lower(), out_file, out_format)
			else:
				print("Error: invalid BioPython conversion class: %s" % format_setting.bioclass)
				sys.exit(1)
		except ValueError as e:
			print("Error in conversion of " + in_path + " to " + out_format + ": " + str(e))
			sys.exit(1)

		with open(out_file, "r") as tmp_file:
			print(tmp_file.read())

		os.remove(out_file)  # Is this really necessary?
	else:
		for out_format in out_formats:
			for in_path, in_format in id_results.items():
				out_file = out_path
				if sys.platform == "win32":
					if out_file[-1] != "\\":
						out_file += "\\"

					out_file += ntpath.basename(in_path).split('.')[0]
				else:
					if out_file[-1] != "/":
						out_file += "/"

					out_file += os.path.basename(in_path).split('.')[0]

				out_extension = settings[out_format].extension
				out_file = out_file + "." + out_extension
				print("Converting %s file %s to %s file %s" % (in_format, in_path, out_format, out_file))

				try:
					format_setting = settings[in_format]
					if format_setting.bioclass == "seq":
						SeqIO.convert(in_path, in_format.lower(), out_file, out_format, alphabet)
					elif format_setting.bioclass == "phylo":
						Phylo.convert(in_path, in_format.lower(), out_file, out_format)
					elif format_setting.bioclass == "align":
						AlignIO.convert(in_path, in_format.lower(), out_file, out_format)
					else:
						print("Error: invalid BioPython conversion class: %s" % format_setting.bioclass)
						sys.exit(1)
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
	A list of input file paths. If not specified, input is read from stdin.''')

	parser.add_argument('-s', '--stdout', dest="stdout", action='store_true', help='''
	Output result of single-file conversion to stdout.''')

	parser.add_argument('-o', '--outdir', metavar='OUTPATH', nargs=1, help='''
	An output directory for output files. If not specified, output is piped to stdout.''')

	parser.add_argument('-f', '--outfmt', metavar='OUTFORMAT', nargs='+', help='''
	A List of output file formats.''')

	parser.add_argument('-a', '--alphabet', metavar='ALPHA', nargs=1, help='''
	The alphabet to use for conversion (ambigdna, unambigdna, exdna, ambigrna, unambigrna, prot, exprot).''')

	parser.add_argument('-j', '--jobs', metavar='JOBS', nargs=1, type=int, help='''
	The number of processes to use for multiple files (defaults to the number of processor cores).''')

	cli_args = parser.parse_args()
	cli_args.input = cli_args.input[0].split(",")
	cli_args.outfmt = cli_args.outfmt[0].split(",")
	main(cli_args)