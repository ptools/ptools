"""Ptools forcefield implementation."""

from typing import Protocol

import numpy as np
from scipy.spatial.distance import cdist


from .io.readers.attract import read_aminon
from .pairlist import PairList
from .attract import AttractRigidBody


# Name of the force fields implemented in pyattract.
ATTRACT_FORCEFIELDS = ("scorpion", "attract1", "attract2")


# Attract force field default atom type radii and amplitudes.
ATTRACT_DEFAULT_FF_PARAMS = np.array(
    [
        [2.000, 1.000],
        [1.900, 1.000],
        [1.950, 2.000],
        [1.900, 0.600],
        [1.900, 0.600],
        [1.900, 0.600],
        [1.990, 2.000],
        [1.990, 1.500],
        [1.900, 0.600],
        [1.990, 1.500],
        [1.900, 0.600],
        [1.990, 1.500],
        [1.900, 0.500],
        [2.200, 3.000],
        [2.200, 3.000],
        [2.000, 2.000],
        [1.900, 0.600],
        [1.990, 2.000],
        [1.990, 2.900],
        [2.200, 2.000],
        [2.200, 3.000],
        [1.900, 0.600],
        [1.900, 0.600],
        [1.900, 0.600],
        [1.990, 1.500],
        [2.200, 2.600],
        [2.200, 2.000],
        [2.200, 2.600],
        [1.990, 2.500],
        [2.640, 0.600],
        [2.900, 0.600],
        [2.650, 0.600],
        [2.530, 0.600],
        [2.650, 0.600],
        [2.590, 0.600],
        [2.530, 0.600],
        [2.650, 0.600],
        [2.830, 0.600],
        [2.840, 0.600],
        [2.650, 0.600],
        [2.840, 0.600],
        [2.940, 0.600],
    ],
    dtype=float,
)


class AttractForceField1:
    """The AttractForceField1."""

    receptor: AttractRigidBody
    ligand: AttractRigidBody
    cutoff: float
    paramfile: str

    def __init__(
        self,
        receptor: AttractRigidBody,
        ligand: AttractRigidBody,
        cutoff: float = 10,
        paramfile: str = "",
    ):
        self.receptor = receptor
        self.ligand = ligand
        self.cutoff = cutoff
        self.paramfile = paramfile

        self._vdw_energy = 0.0
        self._electrostatic_energy = 0.0

        _repulsive_parameters = np.empty(0)
        _attractive_parameters = np.empty(0)
        _repulsive_pairs = np.empty(0)
        _attractive_pairs = np.empty(0)

        self._initialize_parameters()

    def _initialize_parameters(self):
        if self.paramfile != "":
            params = read_aminon(self.paramfile)
        else:
            params = ATTRACT_DEFAULT_FF_PARAMS

        rad, amp = list(zip(*params))
        rad = np.array(rad)
        amp = np.array(amp)

        self._repulsive_parameters = (
            amp[:, None] * amp[:] * np.power(rad[:, None] + rad[:], 8)
        )
        self._attractive_parameters = (
            amp[:, None] * amp[:] * np.power(rad[:, None] + rad[:], 6)
        )

        # Categorie pairs.
        C = np.array(np.meshgrid(self.receptor.typeids, self.ligand.typeids)).T

        # Numpy insane trickery
        # pylint: disable=E1126
        self._attractive_pairs = self._attractive_parameters[C[..., 0], C[..., 1]]
        self._repulsive_pairs = self._repulsive_parameters[C[..., 0], C[..., 1]]


    def vdw_energy(self):
        """Returns the van der Waals energy."""
        return self._vdw_energy

    def electrostatic_energy(self):
        """Returns the electrostatic energy."""
        return self._electrostatic_energy

    def energy(self) -> float:
        """Returns the total energy between two molecules."""
        return self.non_bonded_energy()

    def update(self):
        """Calculates all energy terms (returns nothing)."""
        self.non_bonded_energy()

    def non_bonded_energy(self):
        """Non-bonded energy calculation."""
        if self.cutoff > 10:
            return self.__nb_energy_large_cutoff()
        return self.__nb_energy_small_cutoff()

    def __nb_energy_small_cutoff(self) -> float:
        """Private method for non-bonded energy calculation with small cutoffs."""

        def van_der_waals(dx, rr2):
            # pylint: disable=E1126
            a = np.array([self._attractive_pairs[i, j] for i, j in zip(*keep)])
            b = np.array([self._repulsive_pairs[i, j] for i, j in zip(*keep)])

            rr23 = np.power(rr2, 3)
            rep = b * rr2
            vlj = (rep - a) * rr23
            fb = 6.0 * vlj + 2.0 * rep * rr23
            fdb = dx * fb[:, None]

            forces = np.zeros((len(self.receptor), len(self.ligand), 3))
            forces[keep] = fdb

            self.receptor.forces += forces.sum(axis=1)
            self.ligand.forces -= forces.sum(axis=0)

            return vlj.sum()

        def electrostatics(dx, rr2):
            charge = (
                self.receptor.charges[:, None]
                * self.ligand.charges
                * (332.053986 / 20.0)
            )
            charge = charge[keep]

            et = charge * rr2
            fdb = dx * (2.0 * et)[:, None]

            forces = np.zeros((len(self.receptor), len(self.ligand), 3))
            forces[keep] = fdb

            self.receptor.forces += forces.sum(axis=1)
            self.ligand.forces -= forces.sum(axis=0)

            return et.sum()

        sq_distances = cdist(
            self.receptor.coordinates, self.ligand.coordinates, metric="sqeuclidean"
        )
        keep = np.where(sq_distances <= self.cutoff * self.cutoff)

        XA = np.take(self.receptor.coordinates, keep[0], axis=0)
        XB = np.take(self.ligand.coordinates, keep[1], axis=0)

        dx = np.asarray(XA - XB)
        rr2 = 1.0 / np.power(dx, 2).sum(axis=1)
        dx = dx + rr2[:, None]

        self._vdw_energy = van_der_waals(dx, rr2)
        self._electrostatic_energy = electrostatics(dx, rr2)
        return self._vdw_energy + self._electrostatic_energy

    def __nb_energy_large_cutoff(self) -> float:
        """Private method for non-bonded energy calculation with large cutoffs."""

        def van_der_waals(dx, rr2):
            rr23 = np.power(rr2, 3)
            rep = self._repulsive_pairs * rr2
            vlj = (rep - self._attractive_pairs) * rr23
            fb = 6.0 * vlj + 2.0 * rep * rr23
            fdb = dx * fb[:, :, None]

            self.receptor.forces += fdb.sum(axis=1)
            self.ligand.forces -= fdb.sum(axis=0)
            return vlj.sum()

        def electrostatics(dx, rr2):
            charge = (
                self.receptor.charges[:, None]
                * self.ligand.charges
                * (332.053986 / 20.0)
            )
            et = charge * rr2
            fdb = dx * (2.0 * et)[:, :, None]
            self.receptor.forces += fdb.sum(axis=1)
            self.ligand.forces -= fdb.sum(axis=0)
            return et.sum()

        dx = np.asarray(self.receptor.coordinates[:, None] - self.ligand.coordinates)
        sq_distances = cdist(
            self.receptor.coordinates, self.ligand.coordinates, metric="sqeuclidean"
        )

        exclude = np.where(sq_distances > self.cutoff * self.cutoff)

        rr2 = 1.0 / np.power(dx, 2).sum(axis=2)
        rr2[exclude] = 0.0
        dx = dx + rr2[:, :, None]

        self._vdw_energy = van_der_waals(dx, rr2)
        self._electrostatic_energy = electrostatics(dx, rr2)
        return self._vdw_energy + self._electrostatic_energy

    # pylint: disable=R0914
    # too many local variables
    def non_bonded_energy_legacy(self):
        """Old-fashioned energy calculation.

        Very slow.
        """
        raise NotImplementedError
        vdw = 0.0
        elec = 0.0

        pairlist = PairList(self.receptor, self.ligand, cutoff=self.cutoff)
        contacts = pairlist.contacts()
        norm2 = pairlist.sqdistances()

        alldx = [
            self.receptor.coordinates[ir] - self.ligand.coordinates[il]
            for ir, il in contacts
        ]

        for i, (ir, il) in enumerate(contacts):
            category_rec = self.receptor.typeids[ir]
            category_lig = self.ligand.typeids[il]

            attractive_pairs = self._attractive_parameters[category_rec][category_lig]
            repulsive_pairs = self._repulsive_parameters[category_rec][category_lig]

            dx = alldx[i]
            r2 = norm2[i]
            rr2 = 1.0 / r2
            dx = dx + rr2

            rr23 = rr2 * rr2 * rr2
            rep = repulsive_pairs * rr2
            vlj = (rep - attractive_pairs) * rr23

            vdw += vlj

            fb = 6.0 * vlj + 2.0 * rep * rr23
            fdb = fb * dx

            self.receptor.forces[ir] += fdb
            self.ligand.forces[il] -= fdb

            # Electrostatics.
            charge_rec = self.receptor.charges[ir]
            charge_lig = self.ligand.charges[il]
            charge = charge_rec * charge_lig * (332.053986 / 20.0)

            if abs(charge) > 0.0:
                et = charge * rr2
                elec += et

                fdb = (2.0 * et) * dx
                self.receptor.forces[ir] += fdb
                self.ligand.forces[il] -= fdb

        self._vdw_energy = vdw
        self._electrostatic_energy = elec

        return vdw + elec
