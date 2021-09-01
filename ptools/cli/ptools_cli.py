"""PTools command-line interface.

For more specific help, provide the name of a positional argument, e.g.,

   ptools --help
"""


import argparse

from . import attract_cmd


def parse_command_line(args=None):
    """Main ptools command-line parsing."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version")
    subparsers = parser.add_subparsers()
    attract_cmd.create_subparser(subparsers)
    return parser.parse_args(args)


def main(args=None):
    """Main ptools command."""
    args = parse_command_line(args)
    args.func(args)


if __name__ == "__main__":
    main()
