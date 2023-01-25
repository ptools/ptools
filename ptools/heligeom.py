"""heligeom.

Some more documentation coming soon.
"""

import string
import numpy as np

from . import measure
from .rigidbody import RigidBody
from .superpose import Screw, mat_trans_2_screw, fit_matrix
from . import transform, linalg


def fnat(
    receptor1: RigidBody, lig1: RigidBody, receptor2: RigidBody, lig2: RigidBody
) -> float:
    """Compute the fraction of native contacts between 2 pairs of RigidBodies.

    Args:
        receptor1: 1st RigidBody of the 1st pair.
        lig1: 2nd RigidBody of the 1st pair.
        receptor2: 1st RigidBody of the 2nd pair.
        lig2: 2nd RigidBody of the 2nd pair.

    Returns:
        float: the fraction computed.
    """
    contact_cutoff = 5.0

    res_pair1 = measure.contacts_by_residue(receptor1, lig1, contact_cutoff)
    if len(res_pair1) == 0:
        return 0

    res_pair2 = measure.contacts_by_residue(receptor2, lig2, contact_cutoff)

    intersect = res_pair1 & res_pair2

    return float(len(intersect)) / float(len(res_pair1))


def dist_axis(rb: RigidBody, hp: Screw) -> tuple[float, float]:
    """compute the minimal and maximal distances between the axis of the Screw
       and all the atoms of the RigidBody rb.

    Args:
        rb: a RigidBody
        hp: a Screw transformation

    Returns:
        a tuple of float containing the minimal and maximal distances.
    """

    dmin, dmax = -1.0, -1.0

    for atom in rb:

        v = atom.coordinates - hp.point
        d = float(np.linalg.norm(np.cross(v, hp.unit)))

        if dmin == -1:
            dmin = d
        elif d < dmin:
            dmin = d

        if dmax == -1:
            dmax = d
        elif d > dmax:
            dmax = d

    return dmin, dmax


def chain_intersect(
    rb1: RigidBody, rb2: RigidBody, delta_resid: int = 0
) -> tuple[RigidBody, RigidBody]:
    """Return new RigidBodies containing only the common residues of the input RigidBodies.

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
    resids1 = set(a.resid for a in rb1.select_atom_type("CA"))
    resids2 = set(a.resid - delta_resid for a in rb2.select_atom_type("CA"))
    resids = resids1.intersection(resids2)
    rb1new = RigidBody(atoms=[a for a in rb1 if a.resid in resids])
    rb2new = RigidBody(atoms=[a for a in rb2 if a.resid - delta_resid in resids])
    return rb1new, rb2new


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
    rb_orig: RigidBody = rb.copy()
    chain_id = 0
    origin = hp.point
    axis = hp.unit
    if Z:
        # Redefine axis to Z
        axis = np.array([0.0, 0, 1])

        # Align the screw axis on Z-axis and apply the transformation on rb_orig
        transform.orient(rb_orig, hp.unit, [0.0, 0.0, 1.0])

        # Align the screw axis on Z-axis and apply the transformation on rb_orig and origin
        # "don't know why this works... but Charles does"
        t_matrix = linalg.orientation_matrix(rb_orig.coordinates, hp.unit, axis)
        transform.transform(rb_orig, t_matrix)

        origin[:] = np.inner(origin, t_matrix[:3, :3])

    chains = [string.ascii_uppercase[chain_id % 26]] * len(rb_orig)
    rb_orig.atom_properties.set("chains", chains)
    chain_id += 1

    final: RigidBody = rb_orig.copy()

    for _ in range(N - 1):
        transform.ab_rotate(rb_orig, origin, origin + axis, hp.angle, degrees=False)
        transform.translate(rb_orig, axis * hp.normtranslation)

        chains = [string.ascii_uppercase[chain_id % 26]] * len(rb_orig)
        rb_orig.atom_properties.set("chains", chains)
        chain_id += 1

        final += rb_orig

    return final
