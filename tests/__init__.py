import pathlib

_DATA_DIR = pathlib.Path(__file__).parent / "data"

TEST_LIGAND = _DATA_DIR / "ligand.pdb"
TEST_RECEPTOR = _DATA_DIR / "receptor.pdb"

# File containing distances between every pair of atom of TEST_RECEPTOR and TEST_LIGAND
# (calculated with VMD and a cutoff of 5 Å).
TEST_DISTANCES_RECEPTOR_LIGAND = _DATA_DIR / "dist.dat"
