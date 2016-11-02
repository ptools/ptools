
"""PyAttract script."""


import argparse

# from pyptools.rigidbody import AttractRigidBody


def parse_command_line(args=None):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-r', '--receptor', required=True,
                        help="path to receptor file")
    parser.add_argument('-l', '--ligand', required=True,
                        help="path to ligand file")
    parser.add_argument('-s', '--single', action="store_true",
                        help="single minimization mode")
    parser.add_argument('--ref',
                        help="reference ligand for rmsd")
    return parser.parse_args(args)


def main(args=None):
    return

#     args = parse_command_line(args)
#     # trjname = 'minimization.trj'  # save minimization variables to trjname

#     # TODO: check file exits:
#     #   - attract.inp
#     #   - aminon.par

#     # Load receptor and ligand.
#     receptor = AttractRigidBody(args.receptor)
#     ligand = AttractRigidBody(args.ligand)

#     # Single calculation mode.
#     translations = [[1, ligand.get_center()]]
#     rotations = [(0, 0, 0)]

#     # Attract core algorithm.
#     for i, trans in enumerate(translations):
#         transnb = i + 1
#         print("Translation {}/{}".format(transnb, len(translations)))
#         for j, rotnb in enumerate(rotations):
#             rotnb = j + 1
#             print("Rotation {}/{}".format(rotnb, len(rotations)))


if __name__ == '__main__':
    main()
