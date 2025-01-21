"""heligeom.

Some more documentation coming soon.
"""

import string

import numpy as np

from . import transform
from .rigidbody import RigidBody
from .superpose import Screw, fit_matrix, mat_trans_2_screw


def chain_intersect(
    rb1: RigidBody, rb2: RigidBody, delta_resid: int = 0
) -> tuple[RigidBody, RigidBody]:
    """Returns new RigidBodies containing only the common residues of the input RigidBodies.

    - Assumes RigidBodies contain a single chain. If multiple chains, split before calling
    and assemble results chain-by-chain.
    - Assumes CA atoms are present for all residues
    - Block alignment of resids is allowed through specifying a resid offset
    - Differences in the numbers of atoms in corresponding residues
      of the two input RigidBodies are preserved

    Args:
        rb1: 1st RigidBody
        rb2: 2nd RigidBody
        delta_resid: the offset in residue indexes between the 2.

    Returns:
        A tuple wit the 2 new RigidBodies.
    """
    resids1 = set(rb1.select("name CA").residue_indices)
    resids2 = set(rb2.select("name CA").residue_indices - delta_resid)
    resids = list(resids1.intersection(resids2))

    atom_ids1 = np.where(np.isin(rb1.residue_indices, resids))[0]
    atom_ids2 = np.where(np.isin(rb2.residue_indices - delta_resid, resids))[0]

    return rb1[atom_ids1].copy(), rb2[atom_ids2].copy()


def heli_analyze(rb1: RigidBody, rb2: RigidBody) -> Screw:
    """Returns the screw transformation from rb1 to rb2.

    Args:
        rb1: 1st RigidBody
        rb2: 2nd RigidBody

    Returns:
        The Screw object.

    """
    return mat_trans_2_screw(fit_matrix(rb1, rb2))


def heli_construct(rb: RigidBody, hp: Screw, N: int, Z: bool = False) -> RigidBody:
    """Constructs a N-mer by repeating the screw transormation hp.

    Args:
        rb: a RigidBody
        hp: a Screw transformation
        N: the number of monomers to construct
        Z: whether or not the oligomer should be aligned to Z axis.

    Returns:
        A RigidBody of the constructed oligomer.
    """
    rb_orig = rb.copy()
    chain_id = 0
    origin = hp.point
    axis = hp.unit

    chains = [string.ascii_uppercase[chain_id % 26]] * len(rb_orig)
    rb_orig.atom_properties["chains"] = chains
    chain_id += 1

    final = rb_orig.copy()

    for _ in range(N - 1):
        transform.ab_rotate(rb_orig, origin, origin + axis, hp.angle, degrees=False)
        transform.translate(rb_orig, axis * hp.normtranslation)

        all_chain_characters = string.ascii_letters + string.digits

        # generates chains_id (based on 1 alphanumerical character) for the whole oligomer
        # 62 possibles combinaisons
        chains = [all_chain_characters[chain_id % len(all_chain_characters)]] * len(rb_orig)
        rb_orig.atom_properties["chains"] = chains
        chain_id += 1

        final += rb_orig

    if Z:
        transform.orient(final, hp.unit, [0.0, 0.0, 1.0])

    final.reset_atom_indices(start=rb.indices[0])

    return final
