"""PTools command-line interface.

For more specific help, provide the name of a positional argument, e.g.,

   ptools --help
"""

import click

from .attract_cmd import attract
from .reduce_cmd import reduce


@click.group(help=__doc__)
def main():
    pass


main.add_command(attract)
main.add_command(reduce)


if __name__ == "__main__":
    main()
