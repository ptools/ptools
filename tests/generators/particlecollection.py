"""Generators for testing."""

from ptools.particlecollection import ParticleCollection

from .atom import generate_atoms


def generate_particlecollection(**kwargs) -> ParticleCollection:
    """Creates a dummy particle collection."""
    return ParticleCollection(generate_atoms(**kwargs))
