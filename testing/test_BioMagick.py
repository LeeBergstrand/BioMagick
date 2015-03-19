#!/usr/bin/env python

# ----------------------------------------------------------------------------------------------
# Created by: Lee & Matt
#
# Description: Contains unit test for BioID Class
# ----------------------------------------------------------------------------------------------
# ==============================================================================================

import os
import yaml
from subprocess import call


# Nose test generator to iterate format test files defined in CSVs
class TestConversion(object):
	def tests(self):
		with open("./testing/conversion_tests.yml", "rU") as tests_file:
			test_cases = tests_file.read()

		for test_case in yaml.safe_load(test_cases):
			inputs = test_case["inputs"]
			outputs = test_case["outputs"]
			formats = test_case["formats"]
			alphabet = test_case["alphabet"] if "alphabet" in test_case else None

			yield self.check_conversion, inputs, outputs, formats, alphabet

	@staticmethod
	def check_conversion(input_files, expected_outputs, output_formats, alphabet):
		# Set up CLI arguments
		args = "-i %s -f %s -a %s" % (",".join(input_files), ",".join(output_formats), alphabet)

		# Do conversion(s)
		ret = call("python BioMagick.py " + args)
		assert ret == 0

		# Check outputs
		files = os.listdir("./")
		for output_file in expected_outputs:
			assert output_file in files

			# Clean up each output file as it's verified
			os.remove(output_file)