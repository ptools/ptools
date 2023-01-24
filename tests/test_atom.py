"""test_atom - Tests for `ptools.atom` module."""

# Unit-test specific imports
from pytest import approx
from .generators import generate_particlecollection

from ptools.tables import atomic_masses

def test_mass_getter():
    atoms = generate_particlecollection()
    atoms.guess_masses()
    assert atoms[0].mass == atomic_masses[atoms[0].element]

def test_mass_setter():
    atoms = generate_particlecollection()
    atoms.guess_masses()
    atoms[0].mass = 42
    assert atoms[0].mass == 42

def test_equal():
    """Asserts that identical atoms from different ParticleCollection instances
    are evaluated equal when attributes are actually equal."""
    left = generate_particlecollection()
    right = left.copy()
    assert left is not right
    assert left[0] is not right[0]
    assert left[0] == right[0]
