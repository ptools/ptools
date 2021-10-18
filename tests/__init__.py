import os

TEST_LIGAND = os.path.join(os.path.dirname(__file__), "data", "ligand.pdb")
TEST_RECEPTOR = os.path.join(os.path.dirname(__file__), "data", "receptor.pdb")

TEST_PDB_MCOPRIGID = os.path.join(os.path.dirname(__file__), "data", "test_mcoprigid.pdb")

# File containing distances between every pair of atom of TEST_RECEPTOR
# and TEST_LIGAND
TEST_DISTANCES_RECEPTOR_LIGAND = os.path.join(
    os.path.dirname(__file__), "data", "dist.dat"
)
