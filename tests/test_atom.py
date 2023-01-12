"""test_atom - Tests for `ptools.atom` module."""

# Unit-test specific imports
from pytest import approx
from .testing.dummy import generate_dummy_atomcollection







def test_mass_getter():
    atoms = generate_dummy_atomcollection()
    assert atoms[0].mass == approx(12.011)

def test_mass_setter():
    atoms = generate_dummy_atomcollection()
    atoms[0].mass = 42
    assert atoms[0].mass == 42

def test_equal():
    """Asserts that identical atoms from different AtomCollection instances
    are evaluated equal when attributes are actually equal."""
    left = generate_dummy_atomcollection()
    right = generate_dummy_atomcollection()
    assert left[0] == right[0]
