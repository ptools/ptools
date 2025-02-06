"""Attract docking."""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import tqdm
from scipy.optimize import minimize

from . import measure, transform
from .forcefield import AttractForceField1
from .linalg import transformation_matrix
from .rigidbody import RigidBody

vector_3d = list[float]


@dataclass
class MinimizationParameters:
    square_cutoff: float
    maximum_iterations: int
    rstk: float = 0.0

    @property
    def cutoff(self) -> float:
        return self.square_cutoff**0.5


@dataclass
class AttractDockingParameters:
    """Stores parameters for an Attract docking."""

    translations: list[vector_3d] = field(default_factory=list)
    rotations: list[vector_3d] = field(default_factory=list)
    minimizations: list[vector_3d] = field(default_factory=list)  # TODO: type is not good


@dataclass
class MinimizationResults:
    start_energy: float
    final_energy: float
    transformation_matrix: np.ndarray
    elapsed: float


def default_minimization_parameters() -> MinimizationParameters:
    """Returns default minimization parameters (cutoff=10., maxiter=100)."""
    return MinimizationParameters(square_cutoff=100.0, maximum_iterations=100)


class AttractRigidBody(RigidBody):
    """AttractRigidBody is a RigidBody on which one can calculate the energy.

    It has additionnal attributes compared to ParticleCollection.

    Attributes:
        - forcefield (str): forcefield name
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
        self.forcefield = ""

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
    def from_pdb(cls, *args, **kwargs):
        raise NotImplementedError("Use ptools.read_attract_topology instead.")

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


def run_attract(*args, **kwargs):
    """Run the Attract docking procedure."""
    return run_attract_monop(*args, **kwargs)


def run_attract_monop(
    ligand: AttractRigidBody, receptor: AttractRigidBody, parameters: AttractDockingParameters
):
    """Run the Attract docking procedure (sequential)."""

    _ligand = ligand.copy()
    _receptor = receptor.copy()

    minimlist = parameters.minimizations
    translations = parameters.translations
    rotations = parameters.rotations

    jobs = [
        (translation, rotation, minimlist) for translation in translations for rotation in rotations
    ]

    total_number_of_jobs = len(jobs) * len(minimlist)

    all_results = []

    progress = tqdm.tqdm(total=total_number_of_jobs, desc="Attract docking")
    for translation, rotation, minimlist in jobs:
        ligand = _ligand.copy()
        receptor = _receptor.copy()

        transform.translate(ligand, -measure.centroid(ligand))
        transform.attract_euler_rotate(ligand, rotation)
        transform.translate(ligand, translation)

        output_data = {
            "translation": translation,
            "rotation": rotation,
            "minimizations": [],
        }

        for minim in minimlist:
            results = _run_minimization(minim, receptor, ligand)
            new_ligand = ligand.copy()

            center = measure.centroid(new_ligand)
            transform.translate(new_ligand, -center)
            transform.transform(new_ligand, results.transformation_matrix)
            transform.translate(new_ligand, center)

            ligand = new_ligand

            output_data["minimizations"].append(
                {
                    "square_cutoff": minim.square_cutoff,
                    "maxiter": minim.maximum_iterations,
                    "rstk": minim.rstk,
                    "start_energy": results.start_energy,
                    "final_energy": results.final_energy,
                    "transformation_matrix": results.transformation_matrix.tolist(),
                    "elapsed": results.elapsed,
                }
            )
            progress.update()

        ff = AttractForceField1(receptor, ligand, 100.0)
        output_data["final_energy"] = ff.non_bonded_energy()

        all_results.append(output_data)

    return all_results


def run_attract_parallel(
    ligand: AttractRigidBody, receptor: AttractRigidBody, parameters: AttractDockingParameters
):
    """Run the Attract docking procedure (parallel)."""
    ligand = ligand.copy()
    receptor = receptor.copy()

    minimlist = parameters.minimizations
    translations = parameters.translations
    rotations = parameters.rotations


    start = time.perf_counter()

    jobs = [
        (_single_attract_job, (ligand.copy(), receptor.copy(), minimlist, translation, rotation))
        for translation in translations
        for rotation in rotations
    ]

    print("Created jobs in", time.perf_counter() - start)
    exit()




    import multiprocessing
    with multiprocessing.Pool(2) as pool:
        results = pool.starmap(_single_attract_job, jobs)
    


def _single_attract_job(
    ligand: AttractRigidBody,
    receptor: AttractRigidBody,
    minimlist: MinimizationParameters,
    translation: vector_3d,
    rotation: vector_3d,
):
    transform.translate(ligand, -measure.centroid(ligand))
    transform.attract_euler_rotate(ligand, rotation)
    transform.translate(ligand, translation)

    output_data = {
        "translation": translation,
        "rotation": rotation,
        "minimizations": [],
    }

    for minim in minimlist:
        results = _run_minimization(minim, receptor, ligand)
        new_ligand = ligand.copy()

        center = measure.centroid(new_ligand)
        transform.translate(new_ligand, -center)
        transform.transform(new_ligand, results.transformation_matrix)
        transform.translate(new_ligand, center)

        ligand = new_ligand

        output_data["minimizations"].append(
            {
                "square_cutoff": minim.square_cutoff,
                "maxiter": minim.maximum_iterations,
                "rstk": minim.rstk,
                "start_energy": results.start_energy,
                "final_energy": results.final_energy,
                "transformation_matrix": results.transformation_matrix.tolist(),
                "elapsed": results.elapsed,
            }
        )

        ff = AttractForceField1(receptor, ligand, 100.0)
        output_data["final_energy"] = ff.non_bonded_energy()
    return output_data



def _run_minimization(
    params: MinimizationParameters,
    receptor: AttractRigidBody,
    ligand: AttractRigidBody,
) -> MinimizationResults:
    """Run the minimization."""

    start = time.perf_counter()

    cutoff = params.cutoff
    niter = params.maximum_iterations

    ff = AttractForceField1(receptor, ligand, cutoff)
    start_energy = ff.non_bonded_energy()

    x0 = np.zeros(6)
    res = minimize(_function, x0, args=(ff,), method="L-BFGS-B", options={"maxiter": niter})
    m = transformation_matrix(res.x[3:], res.x[:3])
    results = MinimizationResults(
        start_energy=start_energy,
        final_energy=res.fun,
        transformation_matrix=m,
        elapsed=time.perf_counter() - start,
    )
    results.x = res.x
    return results
