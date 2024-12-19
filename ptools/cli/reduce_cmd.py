"""PTools reduce command."""

import datetime
import logging

import click

from .. import reduce as Reduce
from ..io import write_reduced_pdb
from .header import print_header

__COMMAND__ = "reduce"

ExistingFile = click.Path(exists=True, dir_okay=False)

logging.basicConfig(format="%(name)s:%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)


@click.command()
@click.argument(
    "topology_path",
    type=ExistingFile,
)
@click.option(
    "--reduction-parameters",
    "reduction_parameters_path",
    help="path to reduction parameters file.",
    type=ExistingFile,
)
@click.option(
    "--name-conversion-rules",
    "name_conversion_rules_path",
    help="path to atom/residue name conversion rules.",
    type=ExistingFile,
    default=Reduce.DEFAULT_ATOM_RENAME_RULES_PATH,
)
@click.option(
    "-i",
    "--ignore-errors",
    help="ignore errors of the specified type(s).",
    multiple=True,
    default=[],
    type=click.Choice(Reduce.exceptions.all_exceptions_names() + ["all"]),
)
@click.option(
    "--ff",
    help="force field to use for reduction.",
    type=click.Choice([name.lower() for name in Reduce.FORCEFIELDS.keys()]),
    default="attract1",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    help="output file name.",
    default="reduced.red",
)
@click.option(
    "--optimize-charges",
    help="optimize charges of the reduced model (scorpion force field only).",
    is_flag=True,
)
def reduce(
    topology_path,
    reduction_parameters_path,
    name_conversion_rules_path,
    ignore_errors,
    ff,
    output_path,
    optimize_charges,
):
    """Runs reduce.

    Args:
        topology_path: input topology file in atomistic resolution.
    """
    print_header(__COMMAND__)
    logger.info("Start time: %s", str(datetime.datetime.now()))

    if "all" in ignore_errors:
        ignore_errors = Reduce.exceptions.all_exceptions_names()
    ignore_exceptions = Reduce.exceptions.exceptions_from_names(ignore_errors)

    if not reduction_parameters_path:
        reduction_parameters_path = Reduce.FORCEFIELDS[ff]

    reducer = Reduce.Reducer(topology_path, reduction_parameters_path, name_conversion_rules_path)
    reducer.reduce(ignore_exceptions)

    logger.info(
        "Reduced atomistic model from %d to %d",
        reducer.number_of_atoms(),
        reducer.number_of_beads(),
    )

    if ff == "scorpion":
        # Sets charges of the first and last "CA" bead to 1 and -1, respectively
        ca_beads = [bead for bead in reducer.beads if bead.type == "CA"]
        ca_beads[0].charge = 1.0
        ca_beads[-1].charge = -1.0

        if optimize_charges:
            _optimize_charges(reducer)

    logger.info("Writing reduced model to %s", output_path)
    write_reduced_pdb(reducer.reduced_model, output_path)

    logger.info("End time: %s", str(datetime.datetime.now()))


def _optimize_charges(reducer: Reduce.Reducer):
    """Optimizes charges of the reduced model."""
    raise NotImplementedError("Charge optimization is not implemented yet.")
