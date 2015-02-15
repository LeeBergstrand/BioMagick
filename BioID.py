#!/usr/bin/env python

# ============================================================================================================
# Created by: Lee Bergstrand & Matt McInnes
# Description: A class for identifying bioinformatic file formats.
# Requirements: This class requires the binaryornot module: https://pypi.python.org/pypi/binaryornot/
# ============================================================================================================

import sys
import re
import yaml
import mmap
import codecs
from binaryornot.check import is_binary


# ============================================================================================================
# BioFormat: Defines the properties of a bioinformatic file format object:
class BioFormat(object):
	def __init__(self, name, file_type, compression, markers):
		self.name = str(name)
		self.file_type = str(file_type)
		self.compression = str(compression)
		self.markers = list(markers)


# ============================================================================================================
# BioFormat: Defines the properties of a bioinformatic file identification object:
class BioID(object):
	def __init__(self, path):
		self.binary_definitions = []
		self.text_definitions = []

		with open(path, "rU") as definition_file:
			contents = definition_file.read()

		# Load JSON format definitions and parse into binary and text lists
		for definition in yaml.safe_load(contents):
			if definition["type"] == "BIN":
				self.binary_definitions.append(
					BioFormat(definition["name"], definition["type"], definition["compression"], definition["markers"]))
			elif definition["type"] == "TEXT":
				self.text_definitions.append(
					BioFormat(definition["name"], definition["type"], definition["compression"], definition["markers"]))

	# ----------------------------------------------------------------------------------------------
	# Method: Used to identify a list of bioinformatic files.
	# ----------------------------------------------------------------------------------------------
	def identify(self, files):
		identified = {}

		# Check if each each file is binary or text and use the appropriate marker (regex or byte field).
		for file_path in files:
			if is_binary(file_path):
				identified[file_path] = self.identify_binary_file(file_path)
			else:
				identified[file_path] = self.identify_text_file(file_path)

		return identified

	# ---------------------------------------------------------------------------------------------------
	# Method: Used to determine if input binary file is a bioinformatic file by using raw byte sequences.
	# ---------------------------------------------------------------------------------------------------
	def identify_binary_file(self, file_path):
		binary_input_file = open(file_path, "rb")  # Read in binary file as binary ("rb")
		identity = "unrecognized"

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
				identity = binary_file_type.name
				break

		mapped_file.close()
		binary_input_file.close()
		return identity

	# ----------------------------------------------------------------------------------------------
	# Method: Used to determine if input text file is a bioinformatic file by using REGEXs.
	# ----------------------------------------------------------------------------------------------
	def identify_text_file(self, file_path):
		with open(file_path, "rU") as text_input_file:  # Read in text file as text with universal newlines ("rU")
			input_text = text_input_file.read()
			identity = self.identify_raw_text(input_text)
			return identity

	# ----------------------------------------------------------------------------------------------
	# Method: Used to determine if raw text is a bioinformatic file by using REGEXs.
	# ----------------------------------------------------------------------------------------------
	def identify_raw_text(self, input_text):
		output = "unrecognized"
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
				output = text_file_type.name
				break

		return output

