"""Tests for ptools.io.pdb.RobustPDBBuilder

Reference strings are taken from the PDB format documentation:
    https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
"""

import unittest

from ptools.io.pdb import RobustPDBBuilder


class TestHeader(unittest.TestCase):
    """Tests for the HEADER record creation."""

    def test_header_record_ok(self):
        rec = RobustPDBBuilder.HeaderRecord(
            classification="PHOTOSYNTHESIS",
            date="28-MAR-07",
            idcode="2UXK",
        )
        reference_string = (
            "HEADER    PHOTOSYNTHESIS                          28-MAR-07   2UXK"
        )
        self.assertEqual(rec.pdb, reference_string)


class TestHeaderDate(unittest.TestCase):
    """Tests for the HEADER record date field."""

    def setUp(self):
        self.classification = "PHOTOSYNTHESIS"
        self.idcode = "2UXK"

    def assertExceptionRaisedOnInitialization(self, datestr):
        with self.assertRaises(AssertionError):
            RobustPDBBuilder.HeaderRecord(
                classification=self.classification,
                idcode=self.idcode,
                date=datestr,
            )

    def test_too_many_tokens(self):
        self.assertExceptionRaisedOnInitialization("28-MAR-07-FOO")

    # -- Test day format ---------------------------------------------------------------
    def test_day_is_not_an_integer(self):
        self.assertExceptionRaisedOnInitialization("LU-MAR-07")

    def test_day_is_too_short(self):
        self.assertExceptionRaisedOnInitialization("2-MAR-07")

    def test_day_is_too_long(self):
        self.assertExceptionRaisedOnInitialization("288-MAR-07")

    # -- Test month format -------------------------------------------------------------
    def test_month_is_not_a_string(self):
        self.assertExceptionRaisedOnInitialization("28-03-07")

    def test_month_is_too_short(self):
        self.assertExceptionRaisedOnInitialization("28-MA-07")

    def test_month_is_too_long(self):
        self.assertExceptionRaisedOnInitialization("28-MARS-07")

    # -- Test year format --------------------------------------------------------------
    def test_year_is_not_an_integer(self):
        self.assertExceptionRaisedOnInitialization("28-MAR-UX")

    def test_year_is_too_short(self):
        self.assertExceptionRaisedOnInitialization("28-MAR-0")

    def test_year_is_too_long(self):
        self.assertExceptionRaisedOnInitialization("28-MAR-2007")
