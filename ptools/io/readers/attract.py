from __future__ import annotations

import json

from ..._typing import FilePath
from ...attract import AttractDockingParameters, AttractRigidBody
from .red import read_forcefield_from_red, read_red


def read_topology(path: FilePath) -> AttractRigidBody:
    """Reads a topology file."""
    rigid = AttractRigidBody.from_properties(read_red(path).atom_properties)
    rigid.forcefield = read_forcefield_from_red(path)
    return rigid


def read_docking_parameters(json_file: FilePath) -> AttractDockingParameters:
    """Reads parameters from an attract parameter file."""
    with open(json_file, encoding="utf-8") as file:
        data = json.load(file)
        translations = data["translations"]
        rotations = data["rotations"]
        minimizations = data["minimizations"]
        return AttractDockingParameters(translations, rotations, minimizations)
