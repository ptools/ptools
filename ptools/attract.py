"""Attract docking."""

import time
from dataclasses import dataclass
from typing import TypeVar

import numpy as np
from scipy.optimize import minimize
import tqdm

from . import measure, transform
from ._typing import FilePath
from .forcefield import AttractForceField1
from .io.readers.red import read_red
from .linalg import transformation_matrix
from .rigidbody import RigidBody

AttractRigidBodyType = TypeVar("AttractRigidBodyType", bound="AttractRigidBody")


@dataclass
class MinimizationParameters:

    square_cutoff: float
    maximum_iterations: int

    @property
    def cutoff(self) -> float:
        return self.square_cutoff ** 0.5


@dataclass
class MinimizationResults:
    start_energy: float
    final_energy: float
    transformation_matrix: np.ndarray
    elapsed: float


class AttractRigidBody(RigidBody):
    """AttractRigidBody is a RigidBody on which one can calculate the energy.

    It has 3 additionnal arrays compared to ParticleCollection:
        - typeids (np.ndarray(N, )):
            1 x N shaped array for atom typeids
        - charges (np.ndarray(N, )):O
            1 x N shaped array for atom charges
        - radii (np.ndarray(N, )):O
            1 x N shaped array for atom radii
        - forces (np.ndarray(N, 3)):
            N x 3 shaped array for atom forces

    Atom typeids and charges are parsed from input PDB file.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._initialize_attract_properties()

    def _initialize_attract_properties(self):
        """Initializes atom type ids, charges and forces from PDB extra field."""
        n_atoms = len(self)

        if "typeids" not in self.atom_properties:
            self.add_atom_property("typeid", "typeids", np.zeros(n_atoms, dtype=int))

        if "charges" not in self.atom_properties:
            self.add_atom_property("charge", "charges", np.zeros(n_atoms, dtype=float))

        if "radii" not in self.atom_properties:
            self.add_atom_property("radius", "radii", np.zeros(n_atoms, dtype=float))

        if "forces" not in self.atom_properties:
            self.add_atom_property("force", "forces", np.zeros((n_atoms, 3), dtype=float))

    @classmethod
    def from_red(cls: type[AttractRigidBodyType], path: FilePath) -> AttractRigidBodyType:
        rigid = cls.from_properties(read_red(path).atom_properties)
        return rigid

    @classmethod
    def from_pdb(cls: type[AttractRigidBodyType], path: FilePath) -> AttractRigidBodyType:
        raise NotImplementedError("Use AttractRigidBody.from_red instead.")

    def reset_forces(self):
        """Set all atom forces to (0, 0 0)."""
        self.forces.fill(0)

    def apply_forces(self, forces: np.ndarray):
        """Adds forces to atoms."""
        self.forces += forces  # type: ignore[attr-defined]


def _function(x: np.ndarray, ff: AttractForceField1) -> float:
    """Function to minimize.

    Args:
        x (np.array, list): 6 elements: rotation angles and translation vector
        ff (ptools.forcefield.AttractForceField1)

    Returns:
        float: energy
    """
    source_coordinates = ff.ligand.coordinates.copy()

    rotation = x[:3]
    translation = x[3:]

    transform.rotate(ff.ligand, rotation)
    transform.translate(ff.ligand, translation)

    e = ff.non_bonded_energy()
    ff.ligand.coordinates = source_coordinates
    return e


def run_attract(ligand: AttractRigidBody, receptor: AttractRigidBody, **kwargs):
    """Run the Attract docking procedure.

    Args:
        ligand (ptools.rigidbody.AttractRigidBody)
        receptor (ptools.rigidbody.AttractRigidBody)
        translations (dict)
        rotations (dict[int]->[float, float, float])
        minimlist (list[dict[str]->value])

    Example:
        >>> receptor = ptools.rigidbody.AttractRigidBody("receptor.red")
        >>> ligand = ptools.rigidbody.AttractRigidBody("ligand.red")
        >>> reference = ptools.rigidbody.AttractRigidBody("ligand.red")
        >>> nbminim, lignames, minimlist, rstk = read_attract_parameters("attract.inp")
        >>>
        >>>
        >>> options = {
        ...   translations = {0: measure.centroid(ligand)},
        ...   rotations = {0: (0, 0, 0)},
        ...   minimlist = minimlist,
        }
        >>> ptools.attract.run_attract(ligand, receptor, **options)
    """

    minimlist = kwargs.pop("minimlist", None)
    if minimlist is None:
        raise ValueError("argument 'minimlist' is required")

    _default_translation = {0: measure.centroid(ligand)}
    _default_rotation = {0: (0, 0, 0)}

    translations = kwargs.pop("translations", _default_translation)
    rotations = kwargs.pop("rotations", _default_rotation)

    docking_data = []

    for transi, transnb in enumerate(sorted(translations.keys())):
        trans = translations[transnb]

        for roti, rotnb in enumerate(sorted(rotations.keys())):
            output_data = {
                "translation": trans.tolist(),
                "rotation": rotations[rotnb],
            }

            rot = rotations[rotnb]

            transform.translate(ligand, -measure.centroid(ligand))
            transform.attract_euler_rotate(ligand, rot)
            transform.translate(ligand, trans)

            output_data["minimizations"] = []
            for i, minim in enumerate(minimlist):
                parameters = MinimizationParameters(
                    square_cutoff=minim["squarecutoff"],
                    maximum_iterations=minim["maxiter"]
                )

                results = _run_minimization(parameters, receptor, ligand)
                output_data["minimizations"].append(
                    {
                        "cutoff": parameters.cutoff,
                        "maxiter": parameters.maximum_iterations,
                        "rstk": minim["rstk"],
                        "start_energy": results.start_energy,
                        "final_energy": results.final_energy,
                        "transformation_matrix": results.transformation_matrix.tolist(),
                        "elapsed": results.elapsed,
                    }
                )

                transform.move(ligand, results.transformation_matrix)

            ff = AttractForceField1(receptor, ligand, 100.0, "aminon.par")
            output_data["final_energy"] = ff.non_bonded_energy()
            docking_data.append(output_data)



def _run_minimization(
    params: MinimizationParameters,
    receptor: AttractRigidBody,
    ligand: AttractRigidBody,
) -> MinimizationResults:
    """Run the minimization."""

    start = time.time()

    cutoff = params.cutoff
    niter = params.maximum_iterations

    ff = AttractForceField1(receptor, ligand, cutoff, "aminon.par")
    start_energy = ff.non_bonded_energy()

    x0 = np.zeros(6)
    res = minimize(_function, x0, args=(ff,), method="L-BFGS-B", options={"maxiter": niter})
    m = transformation_matrix(res.x[3:], res.x[:3])
    return MinimizationResults(
        start_energy=start_energy,
        final_energy=res.fun,
        transformation_matrix=m,
        elapsed=time.time() - start,
    )
