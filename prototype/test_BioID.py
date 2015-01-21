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
		self.identifier = BioID()

	def test(self):
		self.assertTrue(self.identifier.identify())

if __name__ == '__main__':
	unittest.main()