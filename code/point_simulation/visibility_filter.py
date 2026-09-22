# Title: Point cloud visibility filter
# Author: Alessandro Carlino
# Date: 2026/09/18

'''
Given a point cloud and a camera, it tells which points the camera actually
sees and which ones a nearer surface hides. It needs the points alone: no
mesh, no reconstruction, no normals given from outside, so the same code runs
on a synthetic scene and on a real scan.

The method is a depth buffer with footprints. Points are projected through
the calibration and each one covers a disk instead of a single pixel, because
a cloud is a set of samples and without a footprint the background leaks
through the gaps between them. A point is seen when nothing closer than its
own depth, by more than the thickness of its own surface patch, covers its
pixel.

Every length is derived from the local sampling spacing of the cloud, so
there is no scene dependent parameter to tune: the footprint of a sample is
how far it sits from its neighbours, and the depth tolerance is how thick its
patch looks from the camera. The only knob left is how many neighbours the
spacing is estimated over, which is a statistical choice and not a scale.

See notes.md for the derivation and the failure modes on real scans.
'''

# Import libraries
from typing import NamedTuple

import numba
import numpy as np
from scipy.spatial import KDTree

# --- Constants -------------------------------------------------------------

NEIGHBOURS = 12 # Neighbours the spacing and the normal are estimated over
NEAR = 1e-6     # Depth below which a point is behind the camera
FOOTPRINT = 0.7 # Footprint radius, in local spacings
GRAZING = 0.1   # Floor on the cosine of incidence, caps the patch thickness
MARGIN = 1.0    # Tolerance in patch thicknesses, absorbs the pixel rounding

class CloudGeometry(NamedTuple):
    """Local shape of the cloud, estimated once and reused by every camera."""

    spacing: np.ndarray # (n_points,), distance to the neighbours
    normal: np.ndarray  # (n_points, 3), unoriented unit surface normal

# --- Describe the cloud ----------------------------------------------------

def describe_cloud(points: np.ndarray,
                   neighbours: int = NEIGHBOURS) -> CloudGeometry:
    """Estimate the sampling spacing and the surface normal of every point.

    Args:
        points: Point cloud of shape (n_points, 3).
        neighbours: Neighbours used for both estimates.

    Returns:
        A CloudGeometry, to be passed to visible_points for every camera.
    """

    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"points must have shape (n, 3), got {points.shape}")
    if not 3 <= neighbours < len(points):
        raise ValueError(
            f"neighbours must lie in [3, {len(points)}), got {neighbours}"
        )

    distance, index = KDTree(points).query(points, k=neighbours + 1)
    neighbourhood = points[np.asarray(index)[:, 1:]]

    # Median neighbour distance
    spacing = np.median(np.asarray(distance)[:, 1:], axis=1)

    # Normal as the least spread direction of the neighbourhood, the standard
    # local plane fit: eigenvector of the smallest eigenvalue
    centred = neighbourhood - neighbourhood.mean(axis=1, keepdims=True)
    covariance = np.einsum("nki,nkj->nij", centred, centred)
    normal = np.linalg.eigh(covariance)[1][:, :, 0]

    return CloudGeometry(spacing, normal)

# --- Visibility ------------------------------------------------------------

def visible_points(points: np.ndarray, geometry: CloudGeometry,
                   centre: np.ndarray, rotation: np.ndarray,
                   calibration: np.ndarray) -> np.ndarray:
    """
    Mark the points a camera sees, hiding those a nearer surface covers.
    Args:
        points: Point cloud of shape (n_points, 3).
        geometry: Output of describe_cloud for the same cloud.
        centre: Camera centre in world coordinates, shape (3,).
        rotation: World to camera rotation, shape (3, 3).
        calibration: Intrinsic matrix, shape (3, 3).
    Returns:
        Boolean array (n_points,), True where the camera sees the point.
    """

    view = _to_image(points, geometry, centre, rotation, calibration)
    depth, row, column, inside, thickness = view
    buffer = _render(view, geometry, calibration)

    seen = np.zeros(len(points), dtype=bool)
    seen[inside] = depth[inside] <= (
        buffer[row[inside], column[inside]] + MARGIN * thickness[inside]
    )

    return seen

def vacated_points(query: np.ndarray, query_geometry: CloudGeometry,
                   points: np.ndarray, geometry: CloudGeometry,
                   centre: np.ndarray, rotation: np.ndarray,
                   calibration: np.ndarray) -> np.ndarray:
    """
    Mark the query points a camera sees past, so no longer there.
    A removal is detected by looking through the place a surface used to be:
    when the surface the camera observes at the pixel of a query point lies
    farther than the point, by more than the thickness of its patch, the
    point is empty space now.
    Args:
        query: Points of the old scene to check, shape (n_query, 3).
        query_geometry: Output of describe_cloud for the query points.
        points: Point cloud the camera actually observes, shape (n_points, 3).
        geometry: Output of describe_cloud for the observed cloud.
        centre: Camera centre in world coordinates, shape (3,).
        rotation: World to camera rotation, shape (3, 3).
        calibration: Intrinsic matrix, shape (3, 3).
    Returns:
        Boolean array (n_query,), True where the camera sees past the point.
    """

    observed = _to_image(points, geometry, centre, rotation, calibration)
    buffer = _render(observed, geometry, calibration)

    depth, row, column, inside, thickness = _to_image(
        query, query_geometry, centre, rotation, calibration
    )

    # An empty pixel holds infinity, so seeing nothing at all counts as past
    vacated = np.zeros(len(query), dtype=bool)
    vacated[inside] = buffer[row[inside], column[inside]] > (
        depth[inside] + MARGIN * thickness[inside]
    )

    return vacated

# --- Helper functions ------------------------------------------------------

def _to_image(points: np.ndarray, geometry: CloudGeometry,
              centre: np.ndarray, rotation: np.ndarray,
              calibration: np.ndarray) -> tuple:
    """
    Depth, pixel and patch thickness of every point seen from a camera.
    """

    width = int(round(2 * calibration[0, 2]))
    height = int(round(2 * calibration[1, 2]))
    focal = np.array([calibration[0, 0], calibration[1, 1]])
    principal = np.array([calibration[0, 2], calibration[1, 2]])

    local = (points - centre) @ rotation.T
    depth = local[:, 2]

    # Behind the camera, or on its plane, is never seen
    front = depth > NEAR
    pixel = np.zeros((len(points), 2))
    pixel[front] = local[front, :2] * focal / depth[front, None] + principal

    column = np.floor(pixel[:, 0]).astype(int)
    row = np.floor(pixel[:, 1]).astype(int)
    inside = front & (0 <= column) & (column < width)
    inside &= (0 <= row) & (row < height)

    # A tilted patch spans a range of depths, and a point must not be hidden
    # by its own surface. The floor on the cosine stops a patch seen edge on
    # from claiming an unbounded thickness
    towards = points - centre
    towards /= np.linalg.norm(towards, axis=1, keepdims=True)
    cosine = np.abs(np.einsum("ni,ni->n", towards, geometry.normal))
    thickness = geometry.spacing / np.maximum(cosine, GRAZING)

    return depth, row, column, inside, thickness

def _render(view: tuple, geometry: CloudGeometry,
            calibration: np.ndarray) -> np.ndarray:
    """
    Depth buffer of a cloud already projected by _to_image.
    """

    depth, row, column, inside, _ = view
    width = int(round(2 * calibration[0, 2]))
    height = int(round(2 * calibration[1, 2]))
    focal = np.array([calibration[0, 0], calibration[1, 1]])

    # A sample stands for a patch as wide as its distance to the neighbours,
    # so its footprint is that width seen from here
    radius = np.zeros(len(depth))
    radius[inside] = (
        FOOTPRINT * focal.mean() * geometry.spacing[inside] / depth[inside]
    )

    return _depth_buffer(column, row, depth, radius, inside, width, height)

def _depth_buffer(column: np.ndarray, row: np.ndarray, depth: np.ndarray,
                  radius: np.ndarray, inside: np.ndarray,
                  width: int, height: int) -> np.ndarray:
    """
    Nearest depth per pixel, every point covering a disk of its radius.
    """

    buffer = np.full(height * width, np.inf)
    _splat(buffer, column, row, depth, radius, np.flatnonzero(inside),
           width, height)

    return buffer.reshape(height, width)

@numba.njit(cache=True)
def _splat(buffer: np.ndarray, column: np.ndarray, row: np.ndarray,
           depth: np.ndarray, radius: np.ndarray, index: np.ndarray,
           width: int, height: int) -> None:
    """
    Write every point into the flat buffer over the disk of its footprint.
    A plain loop over the points and the pixels of their own disk, compiled:
    a sparse cloud needs footprints of a dozen pixels to seal its surface, and
    building index arrays for millions of writes was most of the cost.
    """

    for point in index:
        limit = radius[point] * radius[point]
        extent = int(np.ceil(radius[point]))

        for shift_row in range(-extent, extent + 1):
            target_row = row[point] + shift_row
            if target_row < 0 or target_row >= height:
                continue

            for shift_column in range(-extent, extent + 1):
                if shift_row * shift_row + shift_column * shift_column > limit:
                    continue
                target_column = column[point] + shift_column
                if target_column < 0 or target_column >= width:
                    continue

                pixel = target_row * width + target_column
                if depth[point] < buffer[pixel]:
                    buffer[pixel] = depth[point]
