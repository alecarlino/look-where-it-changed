# Title: Change region definition
# Author: Alessandro Carlino
# Date: 2026/09/18

'''
Given a point cloud, it marks at random the points belonging to the changed
region of the scene. The status is assigned after the generation and
independently from it.
A region is a geodesic ball on the surface: a seed point is drawn at random
and every point within reach of it, measured along a nearest neighbour graph
rather than through empty space, is marked as changed.
'''

# Import libraries
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import KDTree

# --- Constants -------------------------------------------------------------

NEIGHBOURS = 10  # Neighbours per point in the graph the geodesic runs on
SEPARATION = 2.0 # Minimum seed distance, in units of the region reach
ATTEMPTS = 200   # Seed candidates drawn before giving up on a region


# --- Define change points --------------------------------------------------

def define_change_points(
    points: np.ndarray,
    n_regions: int = 1,
    extent: float = 0.25,
    seed: int | None = None,
) -> np.ndarray:
    """
    Mark the points that belong to the changed region of the scene.
    Args:
        points: Point cloud of shape (n_points, 3).
        n_regions: Number of separate changed regions.
        extent: Geodesic reach of a region, as a fraction of cloud radius.
        seed: Seed for reproducibility, independent from the scene one.
    Returns:
        Boolean array (n_points,), True where the point has changed.
    Raises:
        RuntimeError: If the regions cannot be placed far enough apart.
    """

    # Inputs check
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"points must have shape (n, 3), got {points.shape}")
    if n_regions < 0:
        raise ValueError(f"n_regions must be non negative, got {n_regions}")
    if extent <= 0:
        raise ValueError(f"extent must be positive, got {extent}")

    changed = np.zeros(len(points), dtype=bool)
    if n_regions == 0:
        return changed

    rng = np.random.default_rng(seed)
    graph = _neighbour_graph(points)

    # The reach is a fraction of the radius of the cloud, so that the size of
    # a region does not depend on the scale the scene was generated at
    radius = np.linalg.norm(points - points.mean(axis=0), axis=1).max()
    reach = extent * radius

    placed = []
    for _region in range(n_regions):
        origin = _draw_seed(rng, len(points), placed, reach)
        distance = dijkstra(graph, indices=origin, limit=SEPARATION * reach)

        placed.append(distance)
        changed |= distance <= reach

    return changed

# --- Helper functions ------------------------------------------------------

def _neighbour_graph(points: np.ndarray) -> csr_matrix:
    """
    Symmetric nearest neighbour graph, weighted by Euclidean distance.
    """

    distance, index = KDTree(points).query(points, k=NEIGHBOURS + 1)

    # Column zero is the point itself, at distance zero, and is dropped
    weight = np.asarray(distance)[:, 1:].ravel()
    column = np.asarray(index)[:, 1:].ravel()
    row = np.repeat(np.arange(len(points)), NEIGHBOURS)

    graph = csr_matrix((weight, (row, column)),
                       shape=(len(points), len(points)))

    # Symmetrise: being a neighbour of is not a symmetric relation
    return graph.maximum(graph.T)

def _draw_seed(rng: np.random.Generator, n_points: int, placed: list,
               reach: float) -> int:
    """
    Draw a seed far enough from the regions already placed.
    """

    for _attempt in range(ATTEMPTS):
        candidate = int(rng.integers(n_points))

        # Reject a candidate the previous regions already reach, so that the
        # regions stay distinct instead of merging into fewer, larger ones
        if all(d[candidate] > SEPARATION * reach for d in placed):
            return candidate

    raise RuntimeError(
        f"could not place {len(placed) + 1} regions at least "
        f"{SEPARATION:.3g} reaches apart after {ATTEMPTS} attempts; "
        f"reduce n_regions or extent"
    )