# Scientific libraries.
import numpy as np

# PTools imports.
from ptools.linalg import transform

# More test-specific imports.
from ..testing import assert_array_almost_equal


def test_rotate_by():
    # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
    coords = np.array(list([i + 1.0, i + 11.0, i + 21.0] for i in range(10)))

    # Rotate by 12° along X-axis.
    transform.rotate_by(coords, [12, 0, 0])

    # Reference coordinates calculated with VMD:
    #     set m [transaxis x 12]
    #     set sel [atomselect 0 "all"]
    #     $sel move $m
    expected = [
        [1.0000000, 6.393478, 22.828129],
        [2.0000000, 7.163715, 24.014187],
        [3.0000000, 7.933951, 25.200247],
        [4.0000000, 8.704186, 26.386307],
        [5.0000000, 9.474422, 27.572367],
        [6.0000000, 10.244658, 28.758427],
        [7.0000000, 11.014894, 29.944485],
        [8.0000000, 11.785130, 31.130545],
        [9.0000000, 12.555366, 32.316605],
        [10.0000000, 13.325602, 33.502663],
    ]

    assert_array_almost_equal(coords, expected, decimal=5)


def test_transform_rotation():
    """Test transformation with a rotation matrix.

    In this test, the last row and columns are identity.
    """
    # Matrix calculated with VMD that correspond to a rotation
    # in the X-axis by 10°.
    matrix = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.984807753012208, -0.17364817766693033, 0.0],
            [0.0, 0.17364817766693033, 0.984807753012208, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )

    # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
    coords = np.array(list([i + 1.0, i + 11.0, i + 21.0] for i in range(10)))

    ref = coords.copy()
    transform.rotate_by(ref, [10, 0, 0])

    transform.transform(coords, matrix)
    assert_array_almost_equal(coords, ref)


def test_transform():
    """Test transformation with a random matrix."""
    # Input transformation matrix (just a random matrix).
    matrix = np.array(
        [
            [0.72898284, 0.43879558, 0.31011219, 0.32696626],
            [0.74291592, 0.00244003, 0.20881428, 0.91014385],
            [0.71045798, 0.85828462, 0.71570462, 0.96008097],
            [0.43315630, 0.30949433, 0.93486660, 0.48177513],
        ]
    )

    # Coordinates are [(1, 11, 21), (2, 12, 22), ..., (10, 20, 30)]
    coords = np.array(list([i + 1.0, i + 11.0, i + 21.0] for i in range(10)))

    # Reference calculated with VMD
    # with input topology = tests/data/test_10atoms.pdb.
    expected = np.array(
        [
            [0.5175050497055054, 0.253219336271286, 1.0914303064346313],
            [0.5412969589233398, 0.2738751769065857, 1.1091270446777344],
            [0.5621657967567444, 0.291993111371994, 1.1246496438980103],
            [0.5806189775466919, 0.3080138564109802, 1.1383754014968872],
            [0.5970528721809387, 0.3222815990447998, 1.1505991220474243],
            [0.611781895160675, 0.33506909012794495, 1.1615549325942993],
            [0.6250582337379456, 0.3465954065322876, 1.171429991722107],
            [0.6370866894721985, 0.35703831911087036, 1.1803770065307617],
            [0.6480352282524109, 0.3665436804294586, 1.1885206699371338],
            [0.6580431461334229, 0.37523239850997925, 1.1959648132324219],
        ]
    )

    transform.transform(coords, matrix)
    assert_array_almost_equal(coords, expected)
