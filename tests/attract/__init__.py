import os
from pathlib import Path


ATTRACT_DATA_TEST_DIR = Path(__file__).parent / "data"

TEST_ATTRACT_PARAMS = ATTRACT_DATA_TEST_DIR / "attract.inp"
TEST_ATTRACT_PARAMS_WITH_LIGAND = ATTRACT_DATA_TEST_DIR / "attract_ligand.inp"

TEST_RECEPTOR_RED = ATTRACT_DATA_TEST_DIR / "receptor.red"
TEST_LIGAND_RED = ATTRACT_DATA_TEST_DIR / "ligand.red"

TEST_AMINON = ATTRACT_DATA_TEST_DIR / "aminon.par"


# Example of a reduced PDB (RED) file content.
TEST_DUM_RED_CONTENT = """\
HEADER    ATTRACT1 REDUCED PDB FILE
ATOM      1 CA   CYS     1      12.025  21.956  13.016    1   0.000 0 0
ATOM      2 CSE  CYS     1      11.702  23.345  13.055    7   0.000 0 0
ATOM      3 CA   GLY     2      12.408  20.728  16.555    1   0.000 0 0
ATOM      4 CA   VAL     3      11.501  17.132  16.643    1   0.000 0 0
ATOM      5 CSE  VAL     3      10.112  16.917  16.247   29   0.000 0 0
ATOM      6 CA   PRO     4      14.215  14.528  17.062    1   0.000 0 0
ATOM      7 CSE  PRO     4      14.229  15.158  18.286   22   0.000 0 0
ATOM      8 CA   ALA     5      13.919  11.330  15.124    1   0.000 0 0
ATOM      9 CSE  ALA     5      14.424  11.148  14.581    2   0.000 0 0
"""


# Equivalent to TEST_DUM_RED_CONTENT in PDB format.
TEST_DUM_PDB_CONTENT = """\
ATOM      1  N   CYS E   1      11.377  21.513  11.770  1.00  7.18           N
ATOM      2  CA  CYS E   1      12.025  21.956  13.016  1.00  5.40           C
ATOM      3  C   CYS E   1      11.406  21.350  14.300  1.00  6.41           C
ATOM      4  O   CYS E   1      10.216  21.020  14.517  1.00  5.73           O
ATOM      5  CB  CYS E   1      12.168  23.454  12.852  1.00  3.26           C
ATOM      6  SG  CYS E   1      10.913  24.625  13.296  1.00  2.00           S
ATOM      7  N   GLY E   2      12.379  21.161  15.213  1.00  6.48           N
ATOM      8  CA  GLY E   2      12.408  20.728  16.555  1.00  5.36           C
ATOM      9  C   GLY E   2      11.698  19.535  17.075  1.00  5.75           C
"""

TEST_AMINON_CONTENT = """\
    1  2.000  1.000  0
    2  1.900  1.000  0
    3  1.950  2.000  0
    4  1.900  0.600  0
    5  1.900  0.600  0
"""
