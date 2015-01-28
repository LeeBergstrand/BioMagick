#!/usr/bin/env python
#
# Prototype Biological File Format Identifier
# By Lee & Matt
#

import sys
import os
import argparse

from BioID import BioID


parser = argparse.ArgumentParser(description="Prototype Biological File Format Identifier")
parser.add_argument("-f", action="store", dest="FormatsPath",
                    help="The path to the format definition file, defaults to formats.json")
parser.add_argument("inputs", metavar="InputFiles", nargs="+",
                    help="The input files whose formats are to be identified")
args = parser.parse_args()

if len(sys.argv) < 2:
	sys.exit(1)

if args.FormatsPath is None:
	# Default format definition file: formats.json
	args.FormatsPath = "formats.json"

if not os.path.exists(args.FormatsPath):
	print("Error: Failed to find format definition file " + args.FormatsPath)
	sys.exit(1)

# <shiny> Replace wildcard ("*") with list of local files </shiny>
for input_path in args.inputs:
	if input_path.endswith("*"):
		args.inputs.remove(input_path)
		for (directory_paths, directory_names, file_names) in os.walk(input_path.rstrip("/*")):
			for file_name in file_names:
				file_path = input_path.rstrip("*") + file_name
				if file_path not in args.inputs:  # Avoid dupes
					args.inputs.append(file_path)

try:
	bid = BioID(args.FormatsPath)
	formats = bid.identify(args.inputs)
	for fmt in formats:
		print(fmt + ": " + formats[fmt])
except Exception as ex:
	print("Error: " + str(ex))
	sys.exit(1)