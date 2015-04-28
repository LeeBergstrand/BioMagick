#!/usr/bin/env python

# ============================================================================================================
# Created by: Lee Bergstrand & Matt McInnes
# Description: A class for auto-identification BioInformatics file formats.
# Requirements:
#   - This script requires the yaml module: pip install yaml
#   - This script requires the binaryornot module: pip intall binaryornot
# ============================================================================================================

import sys
import re
import yaml
import mmap
import codecs
from binaryornot.check import is_binary


# =========================================================================
# BioIDFormat: Defines the properties of a bioinformatic file format object
# =========================================================================
class BioIDFormat(object):
	def __init__(self, name, file_type, compression, markers):
		self.name = str(name)
		self.file_type = str(file_type)
		self.compression = str(compression)
		self.markers = list(markers)


# ===========================================================================
# BioID: Defines the properties of a bioinformatic file identification object
# ===========================================================================
class BioID(object):
	def __init__(self, path):
		self.binary_definitions = []
		self.text_definitions = []

		with open(path, "rU") as definition_file:
			contents = definition_file.read()

		# Load YAML format definitions and parse into binary and text lists
		for definition in yaml.safe_load(contents):
			if definition["type"] == "BIN":
				self.binary_definitions.append(
					BioIDFormat(definition["name"], definition["type"], definition["compression"], definition["markers"]))
			elif definition["type"] == "TEXT":
				self.text_definitions.append(
					BioIDFormat(definition["name"], definition["type"], definition["compression"], definition["markers"]))

	# ---------------------------------------------------------
	# Method used to search binary files for raw byte sequences
	# ---------------------------------------------------------
	def identify_binary(self, file_path):
		binary_input_file = open(file_path, "rb")  # Read in binary file as binary ("rb")

		if sys.platform == "win32":  # Check if os is windows
			mapped_file = mmap.mmap(binary_input_file.fileno(), 0, None, mmap.ACCESS_READ)
		else:
			mapped_file = mmap.mmap(binary_input_file.fileno(), 0, mmap.MAP_PRIVATE, mmap.PROT_READ)

		for binary_file_type in self.binary_definitions:
			pattern_match = False
			for field in binary_file_type.markers:
				# Check for byte field in binary file
				if mapped_file.find(codecs.decode(field, "unicode_escape").encode("latin1")) != -1:
					pattern_match = True
				else:
					pattern_match = False
					break  # All byte fields must be present in order for the file type check to pass
			if pattern_match:
				mapped_file.close()
				binary_input_file.close()
				return binary_file_type.name

		mapped_file.close()
		binary_input_file.close()
		return "unrecognized"

	# -------------------------------------------------------------
	# Method used to match regular expressions against "text" files
	# -------------------------------------------------------------
	def identify_text(self, input_text):
		for text_file_type in self.text_definitions:
			pattern_match = False
			for regex in text_file_type.markers:
				# Replace escaped newlines in regular expressions with actual newline characters
				if re.findall(regex.replace("\\n", "\n"), input_text, re.IGNORECASE):
					pattern_match = True
				else:
					pattern_match = False
					break
			if pattern_match:
				return text_file_type.name

		return "unrecognized"

	# ---------------------------------------------------------------------------------------------------
	# Method used for identifying the file type of a list of files or a block of text received from stdin
	# ---------------------------------------------------------------------------------------------------
	def identify(self, input_data):
		identified = {}

		if type(input_data) == list:
			# identify() was passed a list of filenames
			for file_path in input_data:
				if is_binary(file_path):
					identified[file_path] = self.identify_binary(file_path)
				else:
					text_input_file = open(file_path, "rU")  # Read in text file as text with universal newlines ("rU")
					input_text = text_input_file.read()
					text_input_file.close()
					identified[file_path] = self.identify_text(input_text)
		elif type(input_data) == str:
			# identify() was passed a string from stdin
			identified["unknown_text"] = self.identify_text(input_data)
		else:
			# indentify() was (probably) passed some useless garbage
			print("Error: identify() received unrecognized input data: %s" % str(input_data))
			sys.exit(1)

		return identified