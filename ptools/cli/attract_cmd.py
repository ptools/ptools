"""PTools attract command."""

import argparse
import datetime
from pathlib import Path

import ptools
from ptools.attract import AttractRigidBody, AttractDockingParameters, default_minimization_parameters
from ptools.io import read_attract_docking_parameters, read_attract_topology

from .header import print_header

__COMMAND__ = "attract"


def create_subparser(parent):
    """Creates command-line parser."""
    parser = parent.add_parser(__COMMAND__, help=__doc__)
    parser.set_defaults(func=run)
    parser.add_argument(
        "-r",
        "--receptor",
        dest="receptor_path",
        required=True,
        type=Path,
        help="name of the receptor file",
    )
    parser.add_argument(
        "-l",
        "--ligand",
        dest="ligand_path",
        required=True,
        type=Path,
        help="name of the ligand file",
    )
    parser.add_argument("--ref", dest="reference_path", type=Path, help="reference ligand for rmsd")
    parser.add_argument(
        "-c",
        "--conf",
        default=None,
        type=Path,
        help="attract configuration file (json)",
    )
    # parser.add_argument(
    #     "-p",
    #     "--param",
    #     dest="parameters_path",
    #     type=Path,
    #     help="attract force field parameter file " "(default=default force field file)",
    # )
    # parser.add_argument(
    #     "--ngroups",
    #     action="store",
    #     type=int,
    #     default=1,
    #     help="Desired number of divisions of translations file",
    # )
    # parser.add_argument(
    #     "--ngroup",
    #     action="store",
    #     type=int,
    #     default=1,
    #     help="Which translation group (1 <= ngroup <= ngroups) " "to run (requires --ngroups)",
    # )
    # parser.add_argument(
    #     "--translation",
    #     type=int,
    #     dest="transnb",
    #     default=None,
    #     help="minimize for the provided translation number",
    # )
    # parser.add_argument(
    #     "--rotation",
    #     type=int,
    #     dest="rotnb",
    #     default=None,
    #     help="minimize for the given rotation number",
    # )


def parse_args(args: argparse.Namespace):
    """Parses command-line arguments."""
    assert args.receptor_path.exists(), f"Receptor file not found: {args.receptor_path}"
    assert args.ligand_path.exists(), f"Ligand file not found: {args.ligand_path}"
    assert args.conf.exists(), f"Configuration file not found: {args.conf}"

    if args.reference_path:
        assert args.reference_path.exists(), f"Reference file not found: {args.reference_path}"

    if args.parameters_path:
        assert args.parameters_path.exists(), f"Parameter file not found: {args.parameters_path}"


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
    log(f"Reading parameters file: {file_path}")
    parameters = read_attract_docking_parameters(file_path)
    log(f"  - {len(parameters.translations)} translations")
    log(f"  - {len(parameters.rotations)} rotations per translation")
    log(f"  - {len(parameters.minimizations)} series of minimizations")
    N = len(parameters.translations) * len(parameters.rotations) * len(parameters.minimizations)
    log(f"Total number of minimizations: {N:,}")
    return parameters


def log(*args, **kwargs):
    """Logs a message."""
    print(*args, **kwargs)


def run(args):
    """Runs attract."""
    print_header(__COMMAND__)
    time_start = datetime.datetime.now()
    log("Start time:", time_start)

    # Load receptor and ligand.
    receptor = read_attract_topology(args.receptor_path)
    ligand = read_attract_topology(args.ligand_path)
    log(f"Read receptor (fixed): {args.receptor_path} with {len(receptor)} particules")
    log(f"Read ligand (mobile): {args.ligand_path} with {len(ligand)} particules")

    # Sanity checks for forcefield: should be 'ATTRACT1' and match between receptor and ligand.
    assert_forcefield_match(receptor, ligand)
    assert_forcefield_is_attract1(receptor)
    log(f"Detected forcefield: {receptor.forcefield!r}")

    # Reads reference topology if provided.
    reference_topology = None
    if args.reference_path:
        reference_topology = ptools.read_pdb(args.reference_path)
        log(f"Read reference file: {args.reference_path} with {len(reference_topology)} particules")

    # If a configuration file is provided, read it.
    # Otherwise, minimize from the starting configuration.
    parameters = None
    if args.conf:
        parameters = read_docking_parameters_file(args.conf)
    else:
        log("Minimize from starting configuration")
        translations = {0: ptools.measure.centroid(ligand)}
        rotations = {0: (0, 0, 0)}
        minimlist = [default_minimization_parameters()]
        parameters = AttractDockingParameters(translations, rotations, minimlist)

    # Run attract.
    _ = ptools.attract.run_attract(ligand, receptor, parameters)

    # output compressed ligand and receptor:
    # if not single and print_files:
    #     print(docking.compress_file(args.receptor_path))
    #     print(docking.compress_file(args.ligand_path))
    #     print(docking.compress_file(ff_specs['ff_file']))
    #     print(docking.compress_file('translation.dat'))
    #     print(docking.compress_file('rotation.dat'))
    #     print(docking.compress_file('attract.inp'))

    # print end and elapsed time
    time_end = datetime.datetime.now()
    # print "Finished at: ",now.strftime("%A %B %d %Y, %H:%M")
    log("End time:", time_end)
    log("Elapsed time:", time_end - time_start)
