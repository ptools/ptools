"""heligeom.

Some more documentation coming soon.
"""

import string
import numpy as np

import ptools.rigidbody
from ptools.superpose import mat_trans_2_screw, fit_matrix


def heli_analyze(mono1, mono2):
    """Returns the screw transformation from mono1 to mono2."""
    return mat_trans_2_screw(fit_matrix(mono1, mono2))


def heli_construct(mono1, hp, N, Z=False):
    """Constructs a N-mer by repeating the screw transormation hp."""
    final = ptools.RigidBody()
    monoTest = mono1.copy()
    chain_id = 0
    origin = hp.point
    axis = hp.unit
    if Z:
        # Redefine origin and axis to Z
        origin = np.zeros(3)
        axis = np.array([0.0, 0, 1])

        # Align the screw axis on Z-axis and apply the transformation on monoTest
        monoTest.orient(hp.unit, [0.0, 0.0, 1.0])

    monoTest.set_chain(string.ascii_uppercase[chain_id % 26])
    final += monoTest
    chain_id += 1

    for _ in range(N - 1):
        monoTest.ab_rotate(origin, origin + axis, hp.angle)
        monoTest.translate(axis * hp.normtranslation)
        monoTest.set_chain(string.ascii_uppercase[chain_id % 26])
        final += monoTest
        chain_id += 1

    return final
