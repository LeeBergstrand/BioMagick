#!/usr/bin/env python

# ============================================================================================================
# Created by: Lee Bergstrand & Matt McInnes
# Description: A class that auto-generates unit tests for the BioID module using an input file list.
# Requirements: This class requires the BioID module.
# ============================================================================================================

from BioID import BioID


# ============================================================================================================
# Nose test generator to iterate a list of format test files as defined in the format list CSV.
class TestFormatDefinitions(object):
	def test_formats(self):
		with open("./testing/format_tests.csv", "rU") as formats_file:
			test_files = formats_file.readlines()[1:]

		for test_file in test_files:
			filename, expected_format = test_file.rstrip(",\n").split(",")
			yield self.check_format, filename, expected_format

	# ---------------------------------------------------------------------------------------------------
	# Static Method: Test if BioID has found the input file to be the correct format or is unrecognized.
	# ---------------------------------------------------------------------------------------------------
	@staticmethod
	def check_format(test_file, expected_format):
		# Putting the test file path here saves having to specify a path for each test file in the CSV
		test_file_path = "./testing/testFiles/" + test_file
		id_results = BioID("./formats.yml").identify([test_file_path])
		assert id_results[test_file_path] == expected_format