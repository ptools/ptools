"""PTools command-line interface.

For more specific help, provide the name of a positional argument, e.g.,

   ptools --help
"""

import click

from .attract_cmd import attract
from .logger import configure_logger, logger
from .reduce_cmd import reduce


@click.group(help=__doc__)
def main():
    logger.enable("ptools")
    configure_logger(debug=False)


main.add_command(attract)
main.add_command(reduce)


if __name__ == "__main__":
    main()
