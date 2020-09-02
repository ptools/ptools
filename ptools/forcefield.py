
import math

import numpy as np

from .pairlist import PairList
from .pyattract.io import read_aminon


# Name of the force fields implemented in pyattract.
ATTRACT_FORCEFIELDS = ('scorpion', 'attract1', 'attract2')


# Attract force field default atom type radii and amplitudes.
ATTRACT_DEFAULT_FF_PARAMS = np.array([
    [2.000, 1.000], [1.900, 1.000], [1.950, 2.000], [1.900, 0.600],
    [1.900, 0.600], [1.900, 0.600], [1.990, 2.000], [1.990, 1.500],
    [1.900, 0.600], [1.990, 1.500], [1.900, 0.600], [1.990, 1.500],
    [1.900, 0.500], [2.200, 3.000], [2.200, 3.000], [2.000, 2.000],
    [1.900, 0.600], [1.990, 2.000], [1.990, 2.900], [2.200, 2.000],
    [2.200, 3.000], [1.900, 0.600], [1.900, 0.600], [1.900, 0.600],
    [1.990, 1.500], [2.200, 2.600], [2.200, 2.000], [2.200, 2.600],
    [1.990, 2.500], [2.640, 0.600], [2.900, 0.600], [2.650, 0.600],
    [2.530, 0.600], [2.650, 0.600], [2.590, 0.600], [2.530, 0.600],
    [2.650, 0.600], [2.830, 0.600], [2.840, 0.600], [2.650, 0.600],
    [2.840, 0.600],
    [2.940, 0.600]
], dtype=float)


class ForceField:
    """Base class for calculating the energy between two molecules.

    Methods for calculating energy terms should be overriden by children
    classes.
    """

    def __init__(self, receptor, ligand, cutoff=10):
        self.receptor = receptor
        self.ligand = ligand
        self.cutoff = cutoff
        self._vdw_energy = 0.0
        self._electrostatic_energy = 0.0

    def update(self):
        """Calculate all energy terms while returning nothing."""
        self.non_bonded_energy()

    def energy(self):
        """Return the total energy between the two molecules, which basically
        corresponds to the non-bonded energy."""
        return self.non_bonded_energy()

    def non_bonded_energy(self):
        """Calculate non-bonded energy between a receptor and a ligand."""
        return self.vdw_energy() + self.electrostatic_energy()

    def vdw_energy(self):
        """Calculate the van der Waals energy (Lennard-Jones energy)
        between a receptor and a ligand."""
        return self._vdw_energy

    def electrostatic_energy(self):
        """Calculate the van der Waals energy (Lennard-Jones energy)
        between a receptor and a ligand."""
        return self._electrostatic_energy


class AttractForceField1(ForceField):

    def __init__(self, receptor, ligand, cutoff=10, paramfile=None):
        super().__init__(receptor, ligand, cutoff)
        self._repulsive = []
        self._attractive = []
        self._init_parameters(paramfile)

    def _init_parameters(self, path=None):
        if path is not None:
            params = read_aminon(path)
        else:
            params = ATTRACT_DEFAULT_FF_PARAMS

        rad, amp = list(zip(*params))
        rad = np.array(rad)
        amp = np.array(amp)

        self._repulsive = amp[:, None] * amp[:] * np.power(rad[:, None] + rad[:], 8)
        self._attractive = amp[:, None] * amp[:] * np.power(rad[:, None] + rad[:], 6)



    def non_bonded_energy(self):
        vdw = 0.0
        elec = 0.0

        pairlist = PairList(self.receptor, self.ligand, cutoff=self.cutoff)
        contacts = pairlist.contacts()
        norm2 = pairlist.sqdistances()

        alldx = [self.receptor.coords[ir] - self.ligand.coords[il]
                 for ir, il in contacts]

        for i, (ir, il) in enumerate(contacts):
            category_rec = self.receptor.atom_categories[ir]
            category_lig = self.ligand.atom_categories[il]

            alen = self._attractive[category_rec][category_lig]
            rlen = self._repulsive[category_rec][category_lig]

            dx = alldx[i]
            r2 = norm2[i]
            rr2 = 1.0 / r2
            dx = dx + rr2

            rr23 = rr2 * rr2 * rr2
            rep = rlen * rr2
            vlj = (rep - alen) * rr23

            vdw += vlj

            fb = 6.0 * vlj + 2.0 * rep * rr23
            fdb = fb * dx

            self.receptor.atom_forces[ir] += fdb
            self.ligand.atom_forces[il] -= fdb

            # Electrostatics.
            charge_rec = self.receptor.atom_charges[ir]
            charge_lig = self.ligand.atom_charges[il]
            charge = charge_rec * charge_lig * (332.053986 / 20.0)

            if abs(charge) > 0.0:
                et = charge * rr2
                elec += et

                fdb = (2.0 * et) * dx
                self.receptor.atom_forces[ir] += fdb
                self.ligand.atom_forces[il] -= fdb

        self._vdw_energy = vdw
        self._electrostatic_energy = elec

        return vdw + elec
