"""heligeom.

Some more documentation coming soon.
"""

import string
import numpy as np

from .rigidbody import RigidBody
from .superpose import Screw, mat_trans_2_screw, fit_matrix
from . import transform


def contact(receptor, ligand, cutoff=5):
    """return residues in interaction, use ptools::pairlist"""

    pl = PairList(receptor, ligand, cutoff)
    contactnat = set()  # residue list in interaction

    for id_atom_recpt, id_atom_lig in pl.contacts():

        contactnat.add((receptor[id_atom_recpt].resid,
                        ligand[id_atom_lig].resid))

    return contactnat


def fnat(receptor1, ligcrist, receptor2, ligprobe):
    """return native fraction (fnat)"""
    corig = contact(receptor1, ligcrist)
    print(len(corig))
    if len(corig) == 0:
        return 0
    cnew = contact(receptor2, ligprobe)
    print(len(cnew))
    intersect = corig & cnew
    f = float(len(intersect)) / float(len(corig))
    return f

def distAxis(mono,hp):
    """compute the distance between the axis of the screw and all the atoms of the monomer.
    Return the smallest and biggest distances
    """

    dmin, dmax = -1,-1

    for atom in mono:

        v = atom.coords - hp.point
        d = np.linalg.norm(np.cross(v, hp.unit))

        if dmin == -1:
            dmin = d
        elif d < dmin:
            dmin = d

        if dmax == -1:
            dmax = d
        elif d > dmax:
            dmax = d

    return dmin, dmax

def chain_intersect(rb1, rb2, delta_resid=0):
    """
    Return new RigidBodies containing only the common residues of the input RigidBodies.
    - Assumes RigidBodies contain a single chain. If multiple chains, split before calling
    and assemble results chain-by-chain.
    - Assumes CA atoms are present for all residues
    - Block alignment of resids is allowed through specifying a resid offset
    - Differences in the numbers of atoms in corresponding residues of the two input RigidBodies are preserved
    """
    resids1 = set([ a.resid for a in rb1.select_atom_type("CA")])
    resids2 = set([ a.resid - delta_resid for a in rb2.select_atom_type("CA")])
    resids = resids1.intersection(resids2)
    rb1new = RigidBody(atoms=[ a for a in rb1 if a.resid in resids])
    rb2new = RigidBody(atoms=[ a for a in rb2 if a.resid - delta_resid in resids])
    return rb1new, rb2new


def heli_analyze(mono1: RigidBody, mono2: RigidBody) -> Screw:
    """Returns the screw transformation from mono1 to mono2."""
    return mat_trans_2_screw(fit_matrix(mono1, mono2))


def heli_construct(mono1: RigidBody, hp: Screw, N: int, Z: bool = False) -> RigidBody:
    """Constructs a N-mer by repeating the screw transormation hp."""
    final = RigidBody()
    mono_test: RigidBody = mono1.copy()
    chain_id = 0
    origin = hp.point
    axis = hp.unit
    if Z:
        # Redefine origin and axis to Z
        origin = np.zeros(3)
        axis = np.array([0.0, 0, 1])

        # Align the screw axis on Z-axis and apply the transformation on mono_test
        transform.orient(mono_test, hp.unit, [0.0, 0.0, 1.0])

    mono_test.set_chain(string.ascii_uppercase[chain_id % 26])
    final += mono_test
    chain_id += 1

    for _ in range(N - 1):
        transform.ab_rotate(mono_test, origin, origin + axis, hp.angle, degrees=False)
        transform.translate(mono_test, axis * hp.normtranslation)
        mono_test.set_chain(string.ascii_uppercase[chain_id % 26])
        final += mono_test
        chain_id += 1

    return final
