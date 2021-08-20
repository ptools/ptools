"""heligeom.

Some more documentation coming soon.
"""

import string
import numpy as np

from ptools.rigidbody import RigidBody
from ptools.superpose import Screw, mat_trans_2_screw, fit_matrix


def heli_analyze(mono1: RigidBody, mono2: RigidBody) -> np.ndarray:
    """Returns the screw transformation from mono1 to mono2."""
    return mat_trans_2_screw(fit_matrix(mono1, mono2))


def heli_construct(mono1: RigidBody, hp: Screw, N: int, Z: bool = False) -> RigidBody:
    """Constructs a N-mer by repeating the screw transormation hp."""
    return extend(hp, mono1, N, Z)


def extend(hp: Screw, mono1: RigidBody, N: int, Z: bool = False):
    final = RigidBody()
    mono_test = mono1.copy()
    chain_id = 0
    origin = hp.point
    axis = hp.unit
    if Z:
        # Redefine origin and axis to Z
        origin = np.zeros(3)
        axis = np.array([0.0, 0, 1])

        # Align the screw axis on Z-axis and apply the transformation on mono_test
        mono_test.orient(hp.unit, [0.0, 0.0, 1.0])

    mono_test.set_chain(string.ascii_uppercase[chain_id % 26])
    final += mono_test
    chain_id += 1

    for _ in range(N - 1):
        mono_test.ab_rotate(origin, origin + axis, hp.angle)
        mono_test.translate(axis * hp.normtranslation)
        mono_test.set_chain(string.ascii_uppercase[chain_id % 26])
        final += mono_test
        chain_id += 1

    return final


def main():  # pragma: no cover
    mono1 = ptools.RigidBody("/Users/benoist/Desktop/1BTA-new.pdb")
    mono2 = ptools.RigidBody("/Users/benoist/Desktop/1BTA-2.pdb")

    hp = heli_analyze(mono1, mono2)
    result = heli_construct(mono1, hp, 10)

    print(result.topdb())


if __name__ == "__main__":  # pragma: no cover
    main()
