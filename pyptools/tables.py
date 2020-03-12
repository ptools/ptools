"""Tables.

Used to guess atom types, masses, etc.
"""


# Masses for elements in atomic units (u).
# Taken from MDAnalysis (which took it from CHARMM and Gromacs atommass.dat)

masses = {
    "Ac": 227.028,
    "Al": 26.981539,
    "Am": 243,
    "Sb": 121.757,
    "Ar": 39.948,
    "As": 74.92159,
    "At": 210,
    "Ba": 137.327,
    "Bk": 247,
    "Be": 9.012182,
    "Bi": 208.98037,
    "Bh": 262,
    "B" : 10.811,
    "BR": 79.90400,
    "Cd": 112.411,
    "CA": 40.08000,
    "Cf": 251,
    "C" : 12.01100,
    "Ce": 140.11600,
    "CS": 132.90000,
    "CL": 35.45000,
    "Cr": 51.9961,
    "Co": 58.9332,
    "CU": 63.54600,
    "Cm": 247,
    "Db": 262,
    "Dy": 162.5,
    "Es": 252,
    "Er": 167.26,
    "Eu": 151.965,
    "Fm": 257,
    "F" : 18.99800,
    "Fr": 223,
    "Gd": 157.25,
    "Ga": 69.723,
    "Ge": 72.61,
    "Au": 196.96654,
    "Hf": 178.49,
    "Hs": 265,
    "HE": 4.00260,
    "Ho": 164.93032,
    "H" : 1.00800,
    "In": 114.82,
    "I" : 126.90450,
    "Ir": 192.22,
    "FE": 55.84700,
    "Kr": 83.8,
    "La": 138.9055,
    "Lr": 262,
    "Pb": 207.2,
    "Li": 6.941,
    "Lu": 174.967,
    "MG": 24.30500,
    "Mn": 54.93805,
    "Mt": 266,
    "Md": 258,
    "Hg": 200.59,
    "Mo": 95.94,
    "N" : 14.00700,
    "NA": 22.98977,
    "Nd": 144.24,
    "NE": 20.17970,
    "Np": 237.048,
    "Ni": 58.6934,
    "Nb": 92.90638,
    "No": 259,
    "Os": 190.2,
    "O" : 15.99900,
    "Pd": 106.42,
    "P" : 30.97400,
    "Pt": 195.08,
    "Pu": 244,
    "Po": 209,
    "K" : 39.10200,
    "Pr": 140.90765,
    "Pm": 145,
    "Pa": 231.0359,
    "Ra": 226.025,
    "Rn": 222,
    "Re": 186.207,
    "Rh": 102.9055,
    "RB": 85.46780,
    "Ru": 101.07,
    "Rf": 261,
    "Sm": 150.36,
    "Sc": 44.95591,
    "Sg": 263,
    "Se": 78.96,
    "Si": 28.0855,
    "Ag": 107.8682,
    "Na": 22.989768,
    "Sr": 87.62,
    "S" : 32.06000,
    "Ta": 180.9479,
    "Tc": 98,
    "Te": 127.6,
    "Tb": 158.92534,
    "Tl": 204.3833,
    "Th": 232.0381,
    "Tm": 168.93421,
    "Sn": 118.71,
    "Ti": 47.88,
    "W" : 183.85,
    "U" : 238.0289,
    "V" : 50.9415,
    "Xe": 131.29,
    "Yb": 173.04,
    "Y" : 88.90585,
    "ZN": 65.37000,
    "Zr": 91.224,
}
