"""PTools reduce command."""

import argparse
import datetime
import logging
from pathlib import Path

from .header import print_header

from .. import reduce
from ..io import assert_file_exists, write_reduced_pdb


__COMMAND__ = "reduce"


logging.basicConfig(format="%(name)s:%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)


def create_subparser(parent):
    """Creates command-line parser."""
    parser = parent.add_parser(__COMMAND__, help=__doc__)
    parser.set_defaults(func=run)

    parser.add_argument(
        "topology", help="input topology file in atomistic resolution.", type=Path
    )

    parser.add_argument(
        "--reduction-parameters",
        help="path to reduction parameters file.",
        type=Path,
    )

    parser.add_argument(
        "--name-conversion-rules",
        help="path to atom/residue name conversion rules.",
        type=Path,
        default=reduce.DEFAULT_ATOM_RENAME_RULES_PATH,
    )

    parser.add_argument(
        "--ignore-errors",
        help="ignore errors of the specified type(s).",
        action="append",
        nargs="?",
        default=[],
        choices=reduce.exceptions.all_exceptions_names() + ["all"],
    )

    parser.add_argument(
        "--ff",
        help="force field to use for reduction.",
        choices=[name.lower() for name in reduce.FORCEFIELDS.keys()],
        default="attract1",
    )

    parser.add_argument(
        "-o",
        "--output",
        help="output file name.",
        type=Path,
        default="reduced.pdb",
    )

    group = parser.add_argument_group("scorpion force field specific options")
    group.add_argument(
        "--optimize-charges",
        help="optimize charges of the reduced model (scorpion force field only).",
        action="store_true",
    )

    # group.add_argument(
    #     "--cgopt",
    #     help="alias for ``--optimize-charges``.",
    #     dest="optimize_charges",
    #     action="store_true",
    # )


def parse_args(args: argparse.Namespace):
    """Parses command-line arguments."""
    if not args.reduction_parameters:
        args.reduction_parameters = reduce.FORCEFIELDS[args.ff]

    if args.optimize_charges and args.ff != "scorpion":
        raise ValueError(
            "Charge optimization is only supported for the scorpion force field."
        )

    assert_file_exists(args.topology)
    assert_file_exists(args.reduction_parameters)
    assert_file_exists(args.name_conversion_rules)


def run(args: argparse.Namespace):
    """Runs reduce."""
    print_header(__COMMAND__)
    parse_args(args)

    logger.info("Start time: %s", str(datetime.datetime.now()))

    if "all" in args.ignore_errors:
        args.ignore_errors = reduce.exceptions.all_exceptions_names()
    ignore_exceptions = reduce.exceptions.exceptions_from_names(args.ignore_errors)

    reducer = reduce.Reducer(
        args.topology, args.reduction_parameters, args.name_conversion_rules
    )
    reducer.reduce(ignore_exceptions)

    logger.info(
        "Reduced atomistic model from %d to %d",
        reducer.number_of_atoms(),
        reducer.number_of_beads(),
    )

    if args.ff == "scorpion":
        # Sets charges of the first and last "CA" bead to 1 and -1, respectively
        ca_beads = [bead for bead in reducer.beads if bead.type == "CA"]
        ca_beads[0].charge = 1.0
        ca_beads[-1].charge = -1.0

        if args.optimize_charges:
            raise NotImplementedError("Charge optimization is not yet implemented.")
        #     reducer.optimize_charges()

    logger.info("Writing reduced model to %s", args.output)
    write_reduced_pdb(reducer.beads, args.output)

    logger.info("End time: %s", str(datetime.datetime.now()))
