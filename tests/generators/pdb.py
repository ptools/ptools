"""Helper classes and functions for PDB format I/O testing."""

import tempfile
from dataclasses import dataclass

from .io import generate_tmp_file


@dataclass
class PDBBuilder:
    """Creates a standard PDB file."""

    n_atoms: int = 10
    n_models: int = 1
    has_model_header: bool = True

    @staticmethod
    def header() -> str:
        return "HEADER    TEST PDB                                27-SEP-06   XXXX"

    @staticmethod
    def title() -> str:
        return """\
TITLE     DUMMY PDB FILE.
TITLE    2 THIS IS THE TITLE LINE 2.\
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
REMARK   1  DOI    10.1021/bi001539c\
"""

    @staticmethod
    def cryst() -> str:
        return "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1"

    @staticmethod
    def orig() -> str:
        return """\
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000\
"""

    @staticmethod
    def scale() -> str:
        return """\
SCALE1      1.000000  0.000000  0.000000        0.00000
SCALE2      0.000000  1.000000  0.000000        0.00000
SCALE3      0.000000  0.000000  1.000000        0.00000\
"""

    def atoms(self) -> str:
        atom_names = self.atom_names()
        atoms = [
            f"ATOM  {i+1:5d}  {atom_names[i]:3s} LYS A   1    {i+1:8.3f}{i+11:8.3f}{i+21:8.3f}  1.00  1.40           N"
            for i in range(self.n_atoms)
        ]
        return "\n".join(atoms)

    def atom_names(self):
        names = ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ", "H1"]
        return [names[i % len(names)] for i in range(self.n_atoms)]

    @staticmethod
    def model_header(model_id: int = 1):
        return f"MODEL      {model_id:4d}"

    @staticmethod
    def model_footer():
        return "ENDMDL"

    @property
    def content(self) -> str:
        """Returns a standard pdb.

        By default 1 model, 10 atoms.
        """
        content = [
            self.header(),
            self.title(),
            self.remark(),
            self.cryst(),
            self.orig(),
            self.scale(),
        ]
        if self.n_models > 1 and not self.has_model_header:
            raise ValueError("multiple model PDB must have a model header")

        if self.has_model_header:
            for i in range(self.n_models):
                content += [self.model_header(i + 1), self.atoms(), self.model_footer()]
        else:
            content.append(self.atoms())

        return "\n".join(content)


def generate_pdb_file(
    n_atoms: int = 10, n_models: int = 1, has_model_header: bool = True
) -> tempfile.NamedTemporaryFile:
    """Creates a temporary file that contains n_atoms atoms with a 'MODEL' entry."""
    return generate_tmp_file(
        content=PDBBuilder(n_atoms, n_models, has_model_header).content
    )
