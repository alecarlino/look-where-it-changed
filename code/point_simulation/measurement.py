# Title: RGB measurement model
# Author: Alessandro Carlino
# Date: 2026/09/19

'''
Given points and a camera, it returns the Fisher information that one image
carries about the position of each point.

The sensor is an RGB camera, so a view measures where a point falls in the
image and not how far it is: the information lies across the viewing ray and
is zero along it. Each block is rank two, a single view cannot place a point
in depth, and only two views looking at it from different directions make
its position well defined. That is triangulation, and it is what makes the
planner look for baseline rather than for close, frontal views.

See notes.md for the derivation.
'''

# Import libraries
import numpy as np

# --- Constants -------------------------------------------------------------

PIXEL_NOISE = 1.0 # Standard deviation of an image measurement, in pixels

# --- Information -----------------------------------------------------------

def projection_information(points: np.ndarray, centre: np.ndarray,
                           rotation: np.ndarray,
                           calibration: np.ndarray) -> np.ndarray:
    """
    Fisher information of one image about the position of each point.
    Args:
        points: Points seen by the camera, shape (n_points, 3).
        centre: Camera centre in world coordinates, shape (3,).
        rotation: World to camera rotation, shape (3, 3).
        calibration: Intrinsic matrix, shape (3, 3).
    Returns:
        Array (n_points, 3, 3) of symmetric rank two blocks.
    """

    local = (points - centre) @ rotation.T
    x, y, z = local[:, 0], local[:, 1], local[:, 2]
    fx, fy = calibration[0, 0], calibration[1, 1]

    # Derivative of the pixel with respect to the point in camera frame
    jacobian = np.zeros((len(points), 2, 3))
    jacobian[:, 0, 0] = fx / z
    jacobian[:, 0, 2] = -fx * x / z ** 2
    jacobian[:, 1, 1] = fy / z
    jacobian[:, 1, 2] = -fy * y / z ** 2

    # Back to world frame, where the camera frame is a rotation away
    jacobian = jacobian @ rotation

    return np.einsum("nak,nal->nkl", jacobian, jacobian) / PIXEL_NOISE ** 2
