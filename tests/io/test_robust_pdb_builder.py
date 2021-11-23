"""Tests for ptools.io.pdb.RobustPDBBuilder

Reference strings are taken from the PDB format documentation:
    https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html


Important
---------

Test detection using pytest: classes do not inherit from unittest.TestCase
because of the inheritance scheme.
"""

import unittest

import pytest

from ptools.io.pdb import RobustPDBBuilder



class TestHeaderBase:

    classification: str = "PHOTOSYNTHESIS"
    date: str = "28-MAR-07"
    idcode: str = "2UXK"

    reference_string = "HEADER    PHOTOSYNTHESIS                          28-MAR-07   2UXK"

    def assertExceptionRaisedUponInitialization(
        self, classification: str = None, date: str = None, idcode: str = None
    ):
        with pytest.raises(AssertionError):
            classification_str = self.classification if classification is None else classification
            date_str = self.date if date is None else date
            idcode_str = self.idcode if idcode is None else idcode
            RobustPDBBuilder.HeaderRecord(
                classification=classification_str,
                date=date_str,
                idcode=idcode_str,
            )


class TestHeader(TestHeaderBase):

    def test_header_ok(self):
        rec = RobustPDBBuilder.HeaderRecord(self.classification, self.date, self.idcode)
        assert rec.pdb == self.reference_string

    def test_check_format_line(self):
        RobustPDBBuilder.HeaderRecord.check_format_line(self.reference_string)

    def test_classification_string_too_long(self):
        self.assertExceptionRaisedUponInitialization(classification="X" * 80)

    def test_idcode_string_too_long(self):
        self.assertExceptionRaisedUponInitialization(idcode="X" * 5)

    def test_idcode_string_too_short(self):
        self.assertExceptionRaisedUponInitialization(idcode="X" * 3)


class TestHeaderDate(TestHeaderBase):
    """Tests for the HEADER record date field."""

    def test_too_many_tokens(self):
        self.assertExceptionRaisedUponInitialization(date="28-MAR-07-FOO")

    # -- Test day format ---------------------------------------------------------------
    def test_day_is_not_an_integer(self):
        self.assertExceptionRaisedUponInitialization(date="LU-MAR-07")

    def test_day_is_too_short(self):
        self.assertExceptionRaisedUponInitialization(date="2-MAR-07")

    def test_day_is_too_long(self):
        self.assertExceptionRaisedUponInitialization(date="288-MAR-07")

    # -- Test month format -------------------------------------------------------------
    def test_month_is_not_a_string(self):
        self.assertExceptionRaisedUponInitialization(date="28-03-07")

    def test_month_is_too_short(self):
        self.assertExceptionRaisedUponInitialization(date="28-MA-07")

    def test_month_is_too_long(self):
        self.assertExceptionRaisedUponInitialization(date="28-MARS-07")

    # -- Test year format --------------------------------------------------------------
    def test_year_is_not_an_integer(self):
        self.assertExceptionRaisedUponInitialization(date="28-MAR-UX")

    def test_year_is_too_short(self):
        self.assertExceptionRaisedUponInitialization(date="28-MAR-0")

    def test_year_is_too_long(self):
        self.assertExceptionRaisedUponInitialization(date="28-MAR-2007")


