#!/usr/bin/env python
#
# A class for auto identifying BioInformatics file formats.
# By Lee & Matt

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
		matched = False
		for filePath in files:
			with open(filePath, "rU") as input_file:
				infile = input_file.read()
				mapped_file = mmap.mmap(input_file.fileno(), 0, mmap.MAP_PRIVATE, mmap.PROT_READ)

			if len(infile) == 0:
				identified[filePath] = "empty"  # Empty files have no format :)
				continue

			for definition in self.definitions:
				matched = True
				if "regexen" in definition:
					for regex in definition["regexen"]:
						if not re.findall(regex.replace("\\n", "\n"), infile, re.IGNORECASE):
							matched = False
							break
			if "bytes" in definition:
				for byte_field in definition["bytes"]:
					if mapped_file.find(byte_field.decode("string_escape")) == -1:
						matched = False
						break
			if matched:
				identified[filePath] = definition["name"]
				break

			mapped_file.close()

			if file not in identified:
				identified[filePath] = "unrecognized"

		return identified