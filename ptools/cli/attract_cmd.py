"""PTools attract command."""

import json
import time
from pathlib import Path

import click
from loguru import logger

import ptools
from ptools.attract import (
    AttractDockingParameters,
    AttractRigidBody,
    default_minimization_parameters,
)
from ptools.io import read_attract_docking_parameters, read_attract_topology

from .header import print_header

ExistingFile = click.Path(exists=True, dir_okay=False)

__COMMAND__ = "attract"


def assert_forcefield_match(receptor: AttractRigidBody, ligand: AttractRigidBody):
    """Asserts that the force fields of the receptor and ligand match."""
    if receptor.forcefield != ligand.forcefield:
        raise ValueError(
            f"Receptor and ligand force fields do not match: "
            f"'{receptor.forcefield}' != '{ligand.forcefield}'"
        )


def assert_forcefield_is_attract1(topology: ptools.AttractRigidBody):
    """Asserts that the force field of the topology is 'ATTRACT1'."""
    if topology.forcefield != "ATTRACT1":
        raise NotImplementedError(f"Force field {topology.forcefield!r} not implemented yet")


def read_docking_parameters_file(file_path: Path) -> AttractDockingParameters:
    """Reads docking parameters from a file."""
    logger.info(f"Reading parameters file: {file_path}")
    parameters = read_attract_docking_parameters(file_path)
    logger.info(f"  - {len(parameters.translations)} translations")
    logger.info(f"  - {len(parameters.rotations)} rotations per translation")
    logger.info(f"  - {len(parameters.minimizations)} series of minimizations")
    N = len(parameters.translations) * len(parameters.rotations) * len(parameters.minimizations)
    logger.info(f"Total number of minimizations: {N:,}")
    return parameters


def format_seconds(seconds: float) -> str:
    """Format seconds to a human-readable string."""
    if seconds < 60:
        return f"{seconds:.2f} s"
    minutes, seconds = divmod(seconds, 60)
    if minutes < 60:
        return f"{minutes:.0f} m {seconds:.2f} s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours:.0f} h {minutes:.0f} m {seconds:.2f} s"


@click.command()
@click.option(
    "-r",
    "--receptor",
    "receptor_path",
    required=True,
    type=ExistingFile,
    help="path to receptor file (attract pdb)",
)
@click.option(
    "-l",
    "--ligand",
    "ligand_path",
    required=True,
    type=ExistingFile,
    help="path to ligand file (attract pdb)",
)
@click.option(
    "--ref",
    "reference_path",
    type=ExistingFile,
    help="path to reference ligand file (pdb)",
)
@click.option(
    "-c",
    "--conf",
    "configuration_path",
    type=ExistingFile,
    help="path to attract configuration file (json)",
)
def attract(receptor_path, ligand_path, reference_path, configuration_path):
    """Runs attract."""
    print_header(__COMMAND__)
    time_start = time.perf_counter()

    # Load receptor and ligand.
    receptor = read_attract_topology(receptor_path)
    ligand = read_attract_topology(ligand_path)
    logger.info(f"Read receptor (fixed): {receptor_path} with {len(receptor)} particules")
    logger.info(f"Read ligand (mobile): {ligand_path} with {len(ligand)} particules")

    # Sanity checks for forcefield: should be 'ATTRACT1' and match between receptor and ligand.
    assert_forcefield_match(receptor, ligand)
    assert_forcefield_is_attract1(receptor)
    logger.info(f"Detected forcefield: {receptor.forcefield!r}")

    # Reads reference topology if provided.
    reference_topology = None
    if reference_path:
        reference_topology = ptools.read_pdb(reference_path)
        logger.info(
            f"Read reference file: {reference_path} with {len(reference_topology)} particules"
        )

    # If a configuration file is provided, read it.
    # Otherwise, minimize from the starting configuration.
    parameters = None
    if configuration_path:
        parameters = read_docking_parameters_file(configuration_path)
    else:
        logger.info("Minimize from starting configuration")
        translations = [ptools.measure.centroid(ligand)]
        rotations = [(0, 0, 0)]
        minimlist = [default_minimization_parameters()]
        parameters = AttractDockingParameters(translations, rotations, minimlist)

    # Run attract.
    # results = ptools.attract.run_attract(ligand, receptor, parameters)
    results = ptools.attract.run_attract_parallel(ligand, receptor, parameters)

    # Converts results to JSON-serializable format.
    logger.info("Writing results to 'results.json'")
    with open("results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    # print end and elapsed time
    elapsed = time.perf_counter() - time_start
    logger.info(f"Elapsed time: {format_seconds(elapsed)}")
