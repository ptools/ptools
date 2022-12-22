"""PTools reduce command."""

import argparse
import datetime
import logging
from pathlib import Path

from .header import print_header
from .._typing import PathLike

from .. import reduce
from ..io import assert_file_exists


logging.basicConfig(format="%(name)s:%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)


def create_subparser(parent):
    """Creates command-line parser."""
    parser = parent.add_parser("reduce", help=__doc__)
    parser.set_defaults(func=run)
    parser.add_argument(
        "topology", help="input topology file in atomistic resolution.", type=Path
    )
    parser.add_argument(
        "--reduction-parameters",
        help="path to reduction parameters file.",
        type=Path,
        default=reduce.ATTRACT1_PROTEIN_REDUCTION_PARAMETERS_PATH,
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


def run(args: argparse.Namespace):
    """Runs reduce."""
    print_header("reduce")
    # logger.info(f"Start time: %s", str(datetime.datetime.now()))


    assert_file_exists(args.topology)
    assert_file_exists(args.reduction_parameters)
    assert_file_exists(args.name_conversion_rules)

    if args.ignore_errors == ["all"]:
        args.ignore_errors = reduce.exceptions.all_exceptions_names()
    ignore_exceptions = reduce.exceptions.exceptions_from_names(args.ignore_errors)

    reducer = reduce.Reducer(args.topology, args.reduction_parameters, args.name_conversion_rules)
    reducer.reduce(ignore_exceptions=ignore_exceptions)


    logger.info("Reduced atomistic model from %d to %d", reducer.number_of_atoms(), reducer.number_of_beads())

    # logger.info(f"End time: %s", str(datetime.datetime.now()))
