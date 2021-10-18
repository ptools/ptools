"""testing.io - I/O testing utilities."""

from contextlib import contextmanager
from os import stat
import tempfile
from typing import Union, Tuple


class TestPDBBuilder:
    @staticmethod
    def header() -> str:
        return "HEADER    TEST PDB                                27-SEP-06   XXXX"

    @staticmethod
    def title() -> str:
        return """\
TITLE     DUMMY PDB FILE CONTAINING 10 ATOMS.
TITLE    2 THIS IS THE TITLE LINE 2.
"""

    @staticmethod
    def remark() -> str:
        return """\
REMARK   1
REMARK   1 REFERENCE 1
REMARK   1  AUTH   J.N.BREG,J.H.J.VAN  OPHEUSDEN,M.J.M.BURGERING,
REMARK   1  AUTH 2 R.BOELENS,R.KAPTEIN
REMARK   1  TITL   STRUCTURE OF ARC REPRESSOR  IN SOLUTION: EVIDENCE
REMARK   1  TITL 2 FOR A FAMILY OF B-SHEET DNA-BINDING PROTEIN
REMARK   1  REF    NATURE                        V. 346   586 1990
REMARK   1  REFN                   ISSN 0028-0836
REMARK   1  PMID   2377232
REMARK   1  DOI    10.1038/346586a0
REMARK   1 REFERENCE 2
REMARK   1  AUTH   J.N.BREG,R.BOELENS,A.V.E.GEORGE,R.KAPTEIN
REMARK   1  TITL   SEQUENCE-SPECIFIC 1H NMR  ASSIGNMENT AND SECONDARY
REMARK   1  TITL 2 STRUCTURE OF THE ARC REPRESSOR OF BACTERIOPHAGE
REMARK   1  TITL 3 P22 AS DETERMINED BY 2D 1H NMR SPECTROSCOPY
REMARK   1  REF    BIOCHEMISTRY                  V.  28  9826 1989
REMARK   1  REFN                   ISSN 0006-2960
REMARK   1  PMID   2611268

REMARK   1
REMARK   1 REFERENCE 1
REMARK   1  AUTH   J.MAREK,J.VEVODOVA,I.SMATANOVA,Y.NAGATA,
REMARK   1  AUTH 2 L.A.SVENSSON,J.NEWMAN,M.TAKAGI,J.DAMBORSKY
REMARK   1  TITL   CRYSTAL STRUCTURE OF THE  HALOALKANE DEHALOGENASE
REMARK   1  TITL 2 FROM SPHINGOMONAS PAUCIMOBILIS UT26
REMARK   1  REF    BIOCHEMISTRY                  V.  39 14082 2000
REMARK   1  REFN                   ISSN 0006-2960
REMARK   1  PMID   11087355
REMARK   1  DOI    10.1021/bi001539c
"""

    @staticmethod
    def cryst() -> str:
        return "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1"

    @staticmethod
    def orig() -> str:
        return """\
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
"""

    @staticmethod
    def scale() -> str:
        return """\
SCALE1      1.000000  0.000000  0.000000        0.00000
SCALE2      0.000000  1.000000  0.000000        0.00000
SCALE3      0.000000  0.000000  1.000000        0.00000
"""

    @staticmethod
    def atom(n_atoms: int = 10) -> str:
        if n_atoms > 10:
            raise ValueError("expecting n_atoms <= 10")

        atoms = [
            "ATOM      1  N   LYS A   1       1.000  11.000  21.000  1.00  1.40           N",
            "ATOM      2  CA  LYS A   1       2.000  12.000  22.000  1.00  0.52           C",
            "ATOM      3  C   LYS A   1       3.000  13.000  23.000  1.00  0.39           C",
            "ATOM      4  O   LYS A   1       4.000  14.000  24.000  1.00  0.33           O",
            "ATOM      5  CB  LYS A   1       5.000  15.000  25.000  1.00  1.53           C",
            "ATOM      6  CG  LYS A   1       6.000  16.000  26.000  1.00  2.38           C",
            "ATOM      7  CD  LYS A   1       7.000  17.000  27.000  1.00  3.11           C",
            "ATOM      8  CE  LYS A   1       8.000  18.000  28.000  1.00  3.58           C",
            "ATOM      9  NZ  LYS A   1       9.000  19.000  29.000  1.00  4.21           N",
            "ATOM     10  H1  LYS A   1      10.000  20.000  30.000  1.00  1.94           H",
        ]
        return "\n".join(atoms[:n_atoms])

    @staticmethod
    def atom_names():
        return ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ", "H1"]

    @staticmethod
    def model_header(model_id: int = 1):
        return f"MODEL      {model_id:4d}"

    @staticmethod
    def model_footer():
        return "ENDMDL"

    @classmethod
    def pdb(cls, has_model_header: bool = True, n_models: int = 1):
        """Returns a standard pdb.

        By default 1 model, 10 atoms.
        """
        content = [
            cls.header(),
            cls.title(),
            cls.remark(),
            cls.cryst(),
            cls.orig(),
            cls.scale(),
        ]
        if n_models > 1 and not has_model_header:
            raise ValueError("multiple model PDB must have a model header")

        if has_model_header:
            for i in range(n_models):
                content += [cls.model_header(i + 1), cls.atom(), cls.model_footer()]
        else:
            content.append(cls.atom())
        return "\n".join(content)


def random_filename() -> str:
    """Returns a random file name."""
    with tempfile.NamedTemporaryFile() as tmpfile:
        return tmpfile.name


@contextmanager
def mk_tmp_file(content: str = "", **kwargs) -> tempfile.NamedTemporaryFile:
    """Creates a temporary file."""
    try:
        tmpfile = tempfile.NamedTemporaryFile("wt", **kwargs)
        tmpfile.write(content)
        tmpfile.flush()
        tmpfile.seek(0)
        yield tmpfile
    finally:
        tmpfile.close()


def mk_empty_file(**kwargs) -> tempfile.NamedTemporaryFile:
    """Creates a temporary empty file."""
    return mk_tmp_file(content="")


def mk_pdb_no_model() -> tempfile.NamedTemporaryFile:
    """Creates a temporary file that contains 10 atoms with no 'MODEL' entry."""
    return mk_tmp_file(content=TestPDBBuilder.pdb(has_model_header=False))


def mk_pdb_10_atoms() -> tempfile.NamedTemporaryFile:
    """Creates a temporary file that contains 10 atoms with a 'MODEL' entry."""
    return mk_tmp_file(content=TestPDBBuilder.pdb())


def mk_pdb_models(n_models: int = 1) -> tempfile.NamedTemporaryFile:
    """Creates a temporary file that contains 10 atoms, no 'MODEL' entry."""
    return mk_tmp_file(content=TestPDBBuilder.pdb(n_models=n_models))
