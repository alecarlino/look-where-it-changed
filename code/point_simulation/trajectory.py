# Title: Camera motion on the locus
# Author: Alessandro Carlino
# Date: 2026/09/19

'''
Given two poses on the camera locus, it builds the path the camera follows
from one to the other, measures it, and finds the poses of the locus it
passes by on the way.

The arm is abstract for now. The camera stays on the locus sphere, always
looking at its centre, and moves by interpolating azimuth and elevation
rather than along a great circle: it resembles a base turning and a shoulder
lifting, and it never leaves the band of elevations the locus was built in,
which a great circle between two low poses on opposite sides could.

The poses passed by are what makes a path worth more than its end points:
the camera can slow down and take them on the way, far cheaper than
travelling to them on purpose.

See notes.md for the motion model and its limits.
'''

# Import libraries
import numpy as np
from scipy.spatial import KDTree

# --- Constants -------------------------------------------------------------

STEP = np.radians(1.0) # Angular spacing of the samples along a path
PASSING = 0.5          # Distance counted as passing by, in lattice spacings

# --- Path ------------------------------------------------------------------

def leg(start: np.ndarray, end: np.ndarray) -> np.ndarray:
    """
    Camera positions from one pose to another, on the locus sphere.
    Args:
        start: Centre of the first pose, shape (3,).
        end: Centre of the last pose, shape (3,), at the same radius.
    Returns:
        Array (n_samples, 3), from start to end inclusive.
    """

    radius = float(np.linalg.norm(start))
    azimuth = np.arctan2([start[1], end[1]], [start[0], end[0]])
    elevation = np.arcsin(np.clip([start[2] / radius, end[2] / radius],
                                  -1.0, 1.0))

    # Turn the short way round, across the seam at plus or minus pi
    turn = (azimuth[1] - azimuth[0] + np.pi) % (2.0 * np.pi) - np.pi
    lift = elevation[1] - elevation[0]

    count = max(2, int(np.ceil(max(abs(turn), abs(lift)) / STEP)) + 1)
    share = np.linspace(0.0, 1.0, count)
    a = azimuth[0] + share * turn
    e = elevation[0] + share * lift

    return radius * np.stack(
        [np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)], axis=1
    )

def length(path: np.ndarray) -> float:
    """
    Length of a path, as the sum of its segments.
    """

    return float(np.linalg.norm(np.diff(path, axis=0), axis=1).sum())

def passed_poses(path: np.ndarray, cameras) -> np.ndarray:
    """
    Poses of the locus the path passes by, in the order it meets them.
    Args:
        path: Output of leg.
        cameras: A CameraGrid, the locus.
    Returns:
        Indices into the locus, end points included, each once.
    """

    # Passing by means coming within half a lattice spacing of a pose
    unit = cameras.centre / np.linalg.norm(cameras.centre, axis=1,
                                           keepdims=True)
    tree = KDTree(unit)
    spacing = np.median(tree.query(unit, k=2)[0][:, 1])

    near = tree.query_ball_point(
        path / np.linalg.norm(path, axis=1, keepdims=True),
        PASSING * spacing,
    )

    order = []
    for sample in near:
        for index in sorted(sample):
            if index not in order:
                order.append(index)

    return np.array(order, dtype=int)
