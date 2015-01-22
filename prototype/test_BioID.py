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

	def test(self):
		self.assertEquals(BioID('./formats.json').identify('./testfiles/NC_000932.faa'), 'FASTA')

if __name__ == '__main__':
	unittest.main()