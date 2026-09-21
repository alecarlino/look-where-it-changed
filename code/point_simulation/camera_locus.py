# Title: Camera locus generator
# Author: Alessandro Carlino
# Date: 2026/09/18

'''
Given a number of cameras and the radius of the locus, it creates a Fibonacci
lattice containing the camera centers of the cameras and orients every one of
them towards yhe center of the sphere.
The lattice can be restricted to a band of elevations, to keep the locus inside
the feasible workspace.
'''

# Import libraries
from typing import NamedTuple

import numpy as np

# --- Constants -------------------------------------------------------------

GOLDEN = np.pi * (3.0 - np.sqrt(5.0)) # Golden angle, spacing of the lattice
UP = np.array([0.0, 0.0, 1.0])        # World up, used to fix the camera roll
ASIDE = np.array([1.0, 0.0, 0.0])     # Fallback up for a camera over the pole


class CameraGrid(NamedTuple):
    """Centres, world to camera rotations, and the shared calibration."""

    centre: np.ndarray   # (n_cameras, 3), in world coordinates
    rotation: np.ndarray # (n_cameras, 3, 3), rows right, down, forward
    calibration: np.ndarray # (3, 3)


# --- Generate the camera locus ---------------------------------------------

def camera_grid(
    n_cameras: int,
    radius: float,
    calibration: np.ndarray,
    min_elevation: float = -90.0,
    max_elevation: float = 90.0,
) -> CameraGrid:
    """
    Place cameras over a sphere, all looking at the origin.
    Args:
        n_cameras: Number of cameras on the locus.
        radius: Radius of the locus, independent of the scene radius.
        calibration: Intrinsic matrix (3, 3), shared by every camera.
        min_elevation: Lowest elevation of the band over the equator.
        max_elevation: Highest elevation of the band.
    Returns:
        A CameraGrid. Rotations take world points to camera coordinates in the
        OpenCV convention, so that x_cam = rotation @ (x_world - centre).
    """

    # Inputs check
    if n_cameras < 1:
        raise ValueError(f"n_cameras must be at least 1, got {n_cameras}")
    if radius <= 0:
        raise ValueError(f"radius must be positive, got {radius}")
    calibration = np.asarray(calibration, dtype=float)
    if calibration.shape != (3, 3):
        raise ValueError(
            f"calibration must be (3, 3), got {calibration.shape}"
        )
    if not -90.0 <= min_elevation < max_elevation <= 90.0:
        raise ValueError(
            f"need -90 <= min_elevation < max_elevation <= 90, got "
            f"{min_elevation} and {max_elevation}"
        )

    # Fibonacci lattice: equal steps in height, golden angle in azimuth
    low = np.sin(np.radians(min_elevation))
    high = np.sin(np.radians(max_elevation))

    index = np.arange(1, n_cameras + 1)
    height = low + (high - low) * (index - 0.5) / n_cameras
    ring = np.sqrt(np.maximum(1.0 - height ** 2, 0.0))
    azimuth = index * GOLDEN

    direction = np.stack(
        [ring * np.cos(azimuth), ring * np.sin(azimuth), height], axis=1
    )
    centre = radius * direction

    # Looking at the origin from the locus means looking inwards
    forward = -direction

    # Right stays horizontal to fix the roll
    right = np.cross(forward, UP)
    degenerate = np.linalg.norm(right, axis=1) < 1e-6
    right[degenerate] = np.cross(forward[degenerate], ASIDE)
    right /= np.linalg.norm(right, axis=1, keepdims=True)

    down = np.cross(forward, right)
    rotation = np.stack([right, down, forward], axis=1)

    return CameraGrid(centre, rotation, calibration)
