#!/usr/bin/env python
#
# A class for auto identifying BioInformatics file formats.
# By Lee & Matt

import sys
import re
import json
import mmap


class BioID:
	def __init__(self, path):
		with open(path, "rU") as definition_file:
			contents = definition_file.read()
		self.definitions = json.loads(contents)["formats"]

	def identify(self, files):
		identified = {}
		for filePath in files:
			with open(filePath, "rU") as input_file:
				if sys.platform == "win32":
					mapped_file = mmap.mmap(input_file.fileno(), 0, None, mmap.ACCESS_READ)
				else:
					mapped_file = mmap.mmap(input_file.fileno(), 0, mmap.MAP_PRIVATE, mmap.PROT_READ)

			for definition in self.definitions:
				definition_match = False
				if "regexen" in definition.keys():
					pattern_match = True
					for regex in definition["regexen"]:
						if not re.findall(regex.replace("\\n", "\n"), mapped_file, re.IGNORECASE):
							pattern_match = False
							break
					if pattern_match:
						definition_match = True
					else:
						definition_match = False
				if not definition_match and "bytes" in definition.keys():
					pattern_match = True
					for byte_field in definition["bytes"]:
						if mapped_file.find(byte_field.decode("string_escape")) == -1:
							pattern_match = False
							break
					if pattern_match:
						definition_match = True
					else:
						definition_match = False
				if definition_match:
					identified[filePath] = definition["name"]
					break

			mapped_file.close()

			if filePath not in identified:
				identified[filePath] = "unrecognized"

		return identified