#!/usr/bin/env python

# ----------------------------------------------------------------------------------------------
# Created by: Lee Bergstrand
#
# Description: Contains unit test for BioID Class
# ----------------------------------------------------------------------------------------------
# ==============================================================================================

from BioID import BioID


# Nose test generator to iterate format test files defined in CSVs
class TestFormatDefinitions(object):
	def test_formats(self):
		with open("./testing/format_tests.csv", "rU") as formats_file:
			test_files = formats_file.readlines()[1:]

		for test_file in test_files:
			filename, expected_format = test_file.rstrip(",\n").split(",")
			yield self.check_format, filename, expected_format

	@staticmethod
	def check_format(test_file, expected_format):
		# Putting the test file path here saves having to specify a path for each test file in the CSV
		test_file_path = "./testing/testFiles/" + test_file
		id_results = BioID("./formats.json").identify([test_file_path])
		assert id_results[test_file_path] == expected_format