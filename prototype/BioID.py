import os, re, json, mmap

class BioID:
	def __init__(self, defpath):
		with open(defpath, "r") as deffile:
			conts = deffile.read()
		self.defs = json.loads(conts)["formats"]

	def identify(self, files):
		recog = {}
		for file in files:
			with open(file, "r") as infile:
				buff = infile.read()
				map = mmap.mmap(infile.fileno(), 0, prot=mmap.PROT_READ)

			if len(buff) == 0:
				recog[file] = "empty" # Empty files have no format :)
				continue

			for fdef in self.defs:
				matched = True
				if "regexen" in fdef:
					for regex in fdef["regexen"]:
						if re.findall(regex.replace("\\n", "\n"), buff) == []:
							matched = False
							break
				if "bytes" in fdef:
					for bytes in fdef["bytes"]:
						if map.find(bytes.decode("string_escape")) == -1:
							matched = False
							break

				if matched:
					recog[file] = fdef["name"]
					break

			map.close()
			if not file in recog:
				recog[file] = "unrecognized"

		return recog
