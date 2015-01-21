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


class GenericFileIdTest(unittest.TestCase):

	def setUp(self):
		pass

	@classmethod
	def setUpClass(cls):
		cls._identifier = BioID

	def test(self):
		self.assertEquals(self._identifier.identify('./testfiles/NC_000932.faa'), 'FASTA')

if __name__ == '__main__':
	unittest.main()