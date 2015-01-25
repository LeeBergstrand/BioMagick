#!/usr/bin/env python
#
# A class for auto identifying BioInformatics file formats.
# By Lee & Matt

import sys
import re
import json
import mmap
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
		with open(path, "rU") as definition_file:
			contents = definition_file.read()
		self.definitions = json.loads(contents)["formats"]

	# Method one used for identifying the file type of a list for files:
	def identify(self, files):
		format_list = []
		identified = {}

		for format_type in self.definitions:
			format_list.append(
				BioFormat(format_type["name"], format_type["type"], format_type["compression"], format_type["markers"]))

		binary_list = [file_format for file_format in format_list if file_format.file_type == "BIN"]
		text_list = [file_format for file_format in format_list if file_format.file_type == "TEXT"]
		del format_list  # Delete full file list from memory

		# Check if each each file is binary or text and use the appropriate marker (regex or byte field).
		for filePath in files:
			file_type = ""
			pattern_match = bool
			if is_binary(filePath):
				binary_input_file = open(filePath, "rb")  # Read in binary file as binary ("rb")
				if sys.platform == "win32":  # Check if os is windows
					mapped_file = mmap.mmap(binary_input_file.fileno(), 0, None, mmap.ACCESS_READ)
				else:
					mapped_file = mmap.mmap(binary_input_file.fileno(), 0, mmap.MAP_PRIVATE, mmap.PROT_READ)

				for binary_file_type in binary_list:
					byte_fields = binary_file_type.markers
					for field in byte_fields:
						if mapped_file.find(field.decode("string_escape")) != -1:  # Check for byte field in binary file
							pattern_match = True
						else:
							pattern_match = False
							break  # All byte fields must be present in order for the file type check to pass
					if pattern_match:
						file_type = binary_file_type.name
						break
			else:
				text_input_file = open(filePath, "rU")  # Read in text file as text with universal newlines ("rU")
				input_text = text_input_file.read()
				for text_file_type in text_list:
					regexps = text_input_file.markers
					for regex in regexps:
						if re.findall(regex.replace("\\n", "\n"), input_text, re.IGNORECASE):
							pattern_match = True
						else:
							pattern_match = False
							break
					if pattern_match:
						file_type = text_file_type.name
						break  # All Regexps must be present in order for the file type check to pass.
			if not pattern_match:
				file_type = "unrecognized"

			if filePath not in identified:
				identified[filePath] = file_type

		return identified
