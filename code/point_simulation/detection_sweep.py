# Title: Detection sweep
# Author: Alessandro Carlino
# Date: 2026/09/19

'''
Given the camera locus and a scene before and after a change, it simulates
the cheap detection pass that precedes planning: a few poses spread as far
apart as possible, looking at the scene as it is now, and the part of the
change they actually observe. What they observe is the change the planner
knows about; what they miss is the change exploration still has to find.

The poses are picked by farthest point sampling over the existing locus
rather than from a fresh small lattice, so that the sweep is a subset of the
candidates and the planner can tell which views have already been taken. The
first pose is the highest one, which makes the choice deterministic without
a seed.

What the sweep observes is decided by the detection proxy in knowledge.py,
the same one that later updates the known change with every view the
planner takes: an addition is known once a view sees it, a removal once a
view sees past it.

See notes.md for the reasoning and the limits.
'''

# Import libraries
import numpy as np

from knowledge import Knowledge, no_knowledge, observe
from synthetic_point_cloud import ScenePair
from visibility_filter import CloudGeometry

# --- Detection sweep -------------------------------------------------------

def detection_sweep(scene: ScenePair, old_geometry: CloudGeometry,
                    new_geometry: CloudGeometry, cameras,
                    n_sweep: int) -> tuple[np.ndarray, Knowledge]:
    """
    Pick a few spread out poses and learn the change they observe.
    Args:
        scene: The scene before and after the change.
        old_geometry: Output of describe_cloud for the old scene.
        new_geometry: Output of describe_cloud for the new scene.
        cameras: A CameraGrid, the locus the sweep is drawn from.
        n_sweep: Number of poses in the detection pass.
    Returns:
        The poses in the order taken, and the Knowledge after the sweep,
        starting from none.
    """

    # Inputs check
    if not 1 <= n_sweep <= len(cameras.centre):
        raise ValueError(
            f"n_sweep must lie in [1, {len(cameras.centre)}], got {n_sweep}"
        )

    views = _farthest_points(cameras.centre, n_sweep)

    return views, observe(no_knowledge(scene, cameras), scene, old_geometry,
                          new_geometry, cameras, views)

# --- Helper functions ------------------------------------------------------

def _farthest_points(centre: np.ndarray, count: int) -> np.ndarray:
    """Greedy farthest point sampling, starting from the highest centre."""

    chosen = [int(np.argmax(centre[:, 2]))]
    distance = np.linalg.norm(centre - centre[chosen[0]], axis=1)

    # Each new pose is the one furthest from all the poses already taken
    for _ in range(count - 1):
        chosen.append(int(np.argmax(distance)))
        distance = np.minimum(
            distance, np.linalg.norm(centre - centre[chosen[-1]], axis=1)
        )

    return np.array(chosen)
