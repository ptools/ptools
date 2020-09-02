
import time

import numpy as np
from scipy.optimize import minimize

from .forcefield import AttractForceField1
from .spatial import transformation_matrix


def _function(x, ff):
    """Function to minimize.

    Args:
        x (np.array, list): 6 elements: rotation angles and translation vector
        ff (ptools.forcefield.ForceField)

    Returns:
        float: energy
    """
    X = ff.ligand.coords.copy()

    rotation = x[:3]
    translation = x[3:]

    ff.ligand.rotate(rotation)
    ff.ligand.translate(translation)

    e = ff.non_bonded_energy()
    ff.ligand.coords = X
    return e


def run_attract(ligand, receptor, translations, rotations, minimlist):
    """Run the Attract docking procedure.

    Args:
        ligand (ptools.rigidbody.RigidBody)
        receptor (ptools.rigidbody.RigidBody)
        translations (dict)
        rotations (dict[int]->[float, float, float])
        minimlist (list[dict[str]->value])
    """

    for transi, transnb in enumerate(sorted(translations.keys())):
        trans = translations[transnb]
        print(f"@@ Translation #{transnb} {transi}/{len(translations)}")
        for roti, rotnb in enumerate(sorted(rotations.keys())):
            print(f"@@ Rotation #{rotnb} {roti + 1}/{len(rotations)}")
            rot = rotations[rotnb]

            ligand.translate(-ligand.center())
            ligand.attract_euler_rotate(rot[0], rot[1], rot[2])
            ligand.translate(trans)

            for i, minim in enumerate(minimlist):
                start = time.time()
                cutoff = minim["squarecutoff"] ** 0.5
                niter = minim["maxiter"]

                ff = AttractForceField1(receptor, ligand, cutoff, "aminon.par")

                print(f"- Minimization {i + 1}/{len(minimlist)}:")
                print(f"  - cutoff: {cutoff:.2f} A")
                print(f"  - maxiter: {niter}")
                print(f"  - start energy: {ff.non_bonded_energy():.2f}", flush=True)


                x0 = np.zeros(6)
                res = minimize(_function, x0, args=(ff, ), method="L-BFGS-B",
                            options={"maxiter": niter})


                print("  - results:")
                print(f"    - energy: {res.fun:6.2f}")
                print(f"    - transformation matrix:")
                m = transformation_matrix(res.x[3:], res.x[:3])
                print(m)

                # Moving ligand accordingly.
                ligand.transform(m)

                print(f"    - elapsed: {time.time() - start:.1f} seconds")
                print("=" * 40, flush=True)

            ff = AttractForceField1(receptor, ligand, 100.0, "aminon.par")
            print(f"  - Final energy: {ff.non_bonded_energy(): 6.2f}")

