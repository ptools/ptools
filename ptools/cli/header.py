"""Defines the ``print_header`` function.

``print_header`` prints a header to the console.

It is intented to be used in the ``run`` function of each command-line.
"""

import ptools


COMMAND_TEMPLATE = "PTools {command_name}"


def print_header(command_name: str):
    """Prints a header to the console.

    Args:
        command_name: Name of the command.
    """
    width = 70
    command = COMMAND_TEMPLATE.format(command_name=command_name)
    version_str = f"based on the PTools library v{ptools.__version__}"

    print("*" * width)
    for line in ("", command, "", "Python edition", "", version_str, ""):
        print(f"* {line.center(width - 4)} *")
    print("*" * width)
    print()
