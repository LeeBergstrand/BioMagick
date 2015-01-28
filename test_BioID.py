#!/usr/bin/env python

# ----------------------------------------------------------------------------------------------
# Created by: Lee Bergstrand
#
# Description: Contains unit test for BioID Class
# ----------------------------------------------------------------------------------------------
# ==============================================================================================

# Imports & Setup:
import unittest
from BioID import BioID


class FileIdTest(unittest.TestCase):
	def test_fasta(self):
		id_results = BioID('./formats.json').identify(['./testing/testFiles/NC_000932.faa'])
		test_file, first_file_type = id_results.popitem()
		self.assertEquals(first_file_type, 'FASTA')

	def test_genbank(self):
		id_results = BioID('./formats.json').identify(['./testing/testFiles/NC_000932.gb'])
		test_file, first_file_type = id_results.popitem()
		self.assertEquals(first_file_type, 'GENBANK')

suite = unittest.TestLoader().loadTestsFromTestCase(FileIdTest)
unittest.TextTestRunner(verbosity=2).run(suite)