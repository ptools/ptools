"""Tests for `ptools.particle.Particle`."""

from ptools.particlecollection import ParticleCollection

from .generators import generate_atoms


def test_equality():
    """Tests that two particles are equal if they have the same attributes."""
    col = ParticleCollection(generate_atoms())
    a1 = col.particles[0]
    a2 = col.particles[0]
    assert a1 == a2


def test_inequality():
    """Tests that two particles are not equal if they have different attributes."""
    col = ParticleCollection(generate_atoms())
    a1 = col.particles[0]
    a2 = col.particles[1]
    assert a1 != a2


def test_equality_non_particle():
    """Tests that a particle is equal to a non-particle object with the same attributes."""

    class Dummy:
        pass

    col = ParticleCollection(generate_atoms())
    a1 = col.particles[0]
    a2 = Dummy()
    for attr in col.atom_properties.singular_names():
        setattr(a2, attr, getattr(a1, attr))
    assert a1 == a2


def test_inequality_non_particle():
    """Tests that a particle is not equal to a non-particle object with different attributes."""

    class Dummy:
        pass

    col = ParticleCollection(generate_atoms(2))
    a1 = col.particles[0]
    a2 = Dummy()
    for attr in col.atom_properties.singular_names():
        setattr(a2, attr, "X")
    assert a1 != a2


def test_inequality_missing_attribute():
    """Tests that a particle particle equality returns False and does not fail when
    comparing two particles with different attributes (attribute not present)."""

    class Dummy:
        pass

    col = ParticleCollection(generate_atoms(2))
    a1 = col.particles[0]
    a2 = Dummy()
    assert a1 != a2
