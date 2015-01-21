#!/usr/bin/env python
#
# A class for auto identifying BioInformatics file formats.
# By Lee & Matt

import re
import json
import mmap


class BioID:
    definitions = None

    def __init__(self, path):
        with open(path, "rU") as definitionfile:
            contents = definitionfile.read()
        self.definitions = json.loads(contents)["formats"]

    @classmethod
    def identify(cls, files):
        identified = {}
        for filepath in files:
            with open(filepath, "rU") as inputfile:
                filecontents = inputfile.read()
                mappedfile = mmap.mmap(inputfile.fileno(), 0, mmap.MAP_PRIVATE, mmap.PROT_READ)

            if len(filecontents) == 0:
                identified[filepath] = "empty"  # Empty files have no format :)
                continue

            for definition in cls.definitions:
                matched = True
                if "regexen" in definition:
                    for regex in definition["regexen"]:
                        if not re.findall(regex.replace("\\n", "\n"), filecontents, re.IGNORECASE):
                            matched = False
                            break
            if "bytes" in definition:
                for bytefield in definition["bytes"]:
                    if mappedfile.find(bytefield.decode("string_escape")) == -1:
                        matched = False
                        break
            if matched:
                identified[filepath] = definition["name"]
                break

            mappedfile.close()

            if file not in identified:
                identified[filepath] = "unrecognized"

        return identified