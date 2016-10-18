
"""test_rigidbody - Tests for `pyptools.rigidbody module."""

import unittest

from pyptools.rigidbody import RigidBody

from . import TEST_PDB

class TestRigidBody(unittest.TestCase):

    def setUp(self):
        self.rb = RigidBody(TEST_PDB)

    def test_constructor(self):
        rb = RigidBody(TEST_PDB)
        self.assertEqual(len(rb), 10)

