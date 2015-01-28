#!/usr/bin/env python
#
# A class for auto identifying BioInformatics file formats.
# By Lee & Matt

import sys
import re
import json
import mmap
import codecs
from binaryornot.check import is_binary


# BioFormat: Defines the properties of a bioinformatic file format object:
class BioFormat(object):
	def __init__(self, name, file_type, compression, markers):
		self.name = str(name)
		self.file_type = str(file_type)
		self.compression = str(compression)
		self.markers = list(markers)


# BioFormat: Defines the properties of a bioinformatic file identification object:
class BioID(object):
	def __init__(self, path):
		self.binary_definitions = []
		self.text_definitions = []

		with open(path, "rU") as definition_file:
			contents = definition_file.read()

		# Load JSON format definitions and parse into binary and text lists
		for definition in json.loads(contents)["formats"]:
			if definition["type"] == "BIN":
				self.binary_definitions.append(
					BioFormat(definition["name"], definition["type"], definition["compression"], definition["markers"]))
			elif definition["type"] == "TEXT":
				self.text_definitions.append(
					BioFormat(definition["name"], definition["type"], definition["compression"], definition["markers"]))

	# Method used to search binary files for raw byte sequences
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
				return binary_file_type.name

		return "unrecognized"

	# Method used to match regular expressions against "text" files
	def identify_text(self, file_path):
		text_input_file = open(file_path, "rU")  # Read in text file as text with universal newlines ("rU")
		input_text = text_input_file.read()

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

	# Method one used for identifying the file type of a list for files:
	def identify(self, files):
		identified = {}

		# [[ Preserving this (for now) for the usage of list comprehensions and the del keyword ]] #
		# binary_list = [file_format for file_format in format_list if file_format.file_type == "BIN"]
		# text_list = [file_format for file_format in format_list if file_format.file_type == "TEXT"]
		# del format_list  # Delete full file list from memory

		# Check if each each file is binary or text and use the appropriate marker (regex or byte field).
		for file_path in files:
			if is_binary(file_path):
				identified[file_path] = self.identify_binary(file_path)
			else:
				identified[file_path] = self.identify_text(file_path)

		return identified