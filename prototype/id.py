#!/usr/bin/env python
#
# Prototype Biological File Format Identifier
# By Lee & Matt
#

import sys, os, argparse, re

parser = argparse.ArgumentParser(description="Prototype Biological File Format Identifier")
parser.add_argument("-f", action="store", dest="FormatsPath", help="The path to the format definition file, defaults to formats.dat")
parser.add_argument("inputs", metavar="InputFiles", nargs="+", help="The input files whose formats are to be identified")
args = parser.parse_args()

if len(sys.argv) < 2:
	exit(1)

if args.FormatsPath == None:
	# Default format definition file: formats.dat
	args.FormatsPath = "formats.dat"

if not os.path.exists(args.FormatsPath):
	print "Error: Failed to find " + args.FormatsPath
	exit(1)

# Load definitions
fdf = open(args.FormatsPath, "r")
lines = fdf.read().splitlines()
fdf.close()

# Parse definitions
formats = []
for line in lines:
	formats.append(line.split("\t", 1)) # Format name and regex are tab-separated

# Process input files
recog = []
for filePath in args.inputs:
	with open(filePath, "r") as inFile:
		buff = inFile.read() # << Is it worth checking filesize and breaking up the reads if over a threshold? >>
		for format in formats:
			if re.findall(format[1].replace("\\n", "\n"), buff) != []:
				print filePath + ": " + format[0]
				recog.append(filePath)
				break

	# Notify if unrecognized
	if not filePath in recog:
		print filePath + ": Unrecognized"
