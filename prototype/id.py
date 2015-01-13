#!/usr/bin/env python
#
# Prototype Biological File Format Identifier
# By Lee & Matt
#

import sys, os, argparse, re, json

parser = argparse.ArgumentParser(description="Prototype Biological File Format Identifier")
parser.add_argument("-f", action="store", dest="FormatsPath", help="The path to the format definition file, defaults to formats.json")
parser.add_argument("inputs", metavar="InputFiles", nargs="+", help="The input files whose formats are to be identified")
args = parser.parse_args()

if len(sys.argv) < 2:
	exit(1)

if args.FormatsPath == None:
	# Default format definition file: formats.json
	args.FormatsPath = "formats.json"

if not os.path.exists(args.FormatsPath):
	print "Error: Failed to find " + args.FormatsPath
	exit(1)

# <Shiny> Replace wildcard ("*") with list of local files
if "*" in args.inputs:
	args.inputs.remove("*")
	for (dirPaths, dirNames, fileNames) in os.walk("."):
		for fileName in fileNames:
			if not fileName in args.inputs:
				args.inputs.append(fileName)

# Load & parse definitions
with open(args.FormatsPath, "r") as hDef:
	conts = hDef.read()
formats = json.loads(conts)["formats"]

# Process input files
recog = []
for filePath in args.inputs:
	with open(filePath, "r") as inFile:
		buff = inFile.read() # << Is it worth checking filesize and breaking up the reads if over a threshold? >>

	for format in formats:
		if format["type"] == "regex":
			if re.findall(format["regex"].replace("\\n", "\n"), buff) != []:
				print filePath + ": " + format["name"]
				recog.append(filePath)
				break
		elif format["type"] == "bytes":
			matched = True
			for bytes in format["bytes"]:
				if not bytes in buff:
					matched = False
					break

			# Move on if any byte field doesn't match
			# << Idea: Add "requirement" attribute ("required" or "optional") to byte fields >>
			if not matched:
				break

	# Notify if unrecognized
	if not filePath in recog:
		print filePath + ": Unrecognized"
