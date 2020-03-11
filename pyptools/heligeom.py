"""heligeom.

Some more documentation coming soon.
"""

import pyptools
from pyptools.superpose import mat_trans_2_screw, fit_matrix

import string


def heli_analyze(mono1, mono2):
    """Returns the screw transformation from mono1 to mono2."""
    return mat_trans_2_screw(fit_matrix(mono1, mono2))


def heli_construct(mono1, hp, N, Z=False):
    """Constructs a N-mer by repeating the screw transormation hp."""
    return extend(hp, mono1, N, Z)


def extend(hp, mono1, N, Z=False):
    final = pyptools.RigidBody()
    monoTest = mono1.copy()
    i = 0
    O = hp.point
    axe = hp.unit
    if Z:
        raise NotImplementedError("not implemented yet")
    monoTest.set_chain(string.ascii_uppercase[i % 26])

    print(final.coords.shape)

    final += monoTest

    i += 1

    for j in range(N - 1):
        monoTest.ab_rotate(O, O + axe, hp.angle)
        monoTest.translate(axe * hp.normtranslation)
        monoTest.set_chain(string.ascii_uppercase[i % 26])
        final += monoTest
        i += 1

    return final





def main():
    import pyptools

    mono1 = pyptools.RigidBody("/Users/benoist/Desktop/1BTA-new.pdb")
    mono2 = pyptools.RigidBody("/Users/benoist/Desktop/1BTA-2.pdb")

    hp = heli_analyze(mono1, mono2)
    result = heli_construct(mono1, hp, 10)

    print(result.topdb())




if __name__ == '__main__':
    main()
