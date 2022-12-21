"""PTools reduce command."""

import datetime
import logging

from .header import print_header

logging.basicConfig(format="%(name)s:%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

def create_subparser(parent):
    """Creates command-line parser."""
    parser = parent.add_parser("reduce", help=__doc__)
    parser.set_defaults(func=run)


def run(args):
    """Runs reduce."""
    print_header("reduce")
    # logger.info(f"Start time: %s", str(datetime.datetime.now()))















    # logger.info(f"End time: %s", str(datetime.datetime.now()))
