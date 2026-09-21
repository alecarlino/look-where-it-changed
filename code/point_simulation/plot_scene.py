# Title: Scene and camera locus plotting
# Author: Alessandro Carlino
# Date: 2026/09/18

'''
Draws a synthetic scene, its changed region and the camera locus, as a three
dimensional view plus the three orthogonal projections. The projections are
what actually show whether the shape has concavities, which perspective
hides.

Cameras are drawn as frustums rather than as points with an arrow, so that
position, orientation and field of view are all visible at once. The corners
come from the calibration, assuming the principal point sits at the centre of
the image, so the image is twice the principal point in each direction.

The figure is returned and not shown, so the caller decides between an
interactive window and viz.save_fig.
'''

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection

# --- Constants -------------------------------------------------------------

UNCHANGED = "#4c72b0" # Colour of the unchanged points
CHANGED = "#d62728"   # Colour of the changed region
CAMERA = "#2ca02c"    # Colour of the camera frustums
FRUSTUM = 0.12        # Frustum depth, as a fraction of the locus radius


# --- Plot the scene --------------------------------------------------------

def plot_scene(points: np.ndarray, changed: np.ndarray,
               cameras=None) -> plt.Figure:
    """Draw the scene, the changed region and the camera locus.

    Args:
        points: Point cloud of shape (n_points, 3).
        changed: Boolean mask (n_points,), True where the point has changed.
        cameras: A CameraGrid, or None to draw the scene alone.

    Returns:
        The figure, for the caller to show or save.
    """

    limit = np.abs(points).max()
    if cameras is not None:
        limit = max(limit, np.abs(cameras.centre).max())

    figure = plt.figure(figsize=(12, 10), facecolor="white")

    # Three dimensional view, with the frustums
    axis = figure.add_subplot(2, 2, 1, projection="3d")
    axis.scatter(*points[~changed].T, s=0.5, c=UNCHANGED, alpha=0.30)
    axis.scatter(*points[changed].T, s=2.0, c=CHANGED, alpha=0.90)

    if cameras is not None:
        axis.add_collection3d(Line3DCollection(
            _frustums(cameras), colors=CAMERA, linewidths=0.4, alpha=0.45
        ))

    axis.set_xlim(-limit, limit)
    axis.set_ylim(-limit, limit)
    axis.set_zlim(-limit, limit)
    axis.set_box_aspect((1, 1, 1))
    axis.set_xticklabels([])
    axis.set_yticklabels([])
    axis.set_zticklabels([])
    axis.set_title("3D view")

    # Orthogonal projections, where the concavities and the elevation band of
    # the locus are visible. Cameras are centres only, frustums would hide
    # the very structure the projections are there to show.
    projections = [((0, 1), "x", "y"), ((0, 2), "x", "z"), ((1, 2), "y", "z")]

    for position, ((first, second), name_a, name_b) in enumerate(
            projections, 2):
        axis = figure.add_subplot(2, 2, position)

        if cameras is not None:
            axis.scatter(cameras.centre[:, first], cameras.centre[:, second],
                         s=3, c=CAMERA, alpha=0.7, edgecolors="none")

        axis.scatter(points[~changed][:, first], points[~changed][:, second],
                     s=0.4, c=UNCHANGED, alpha=0.25, edgecolors="none")
        axis.scatter(points[changed][:, first], points[changed][:, second],
                     s=1.5, c=CHANGED, alpha=0.8, edgecolors="none")

        axis.set_xlim(-limit, limit)
        axis.set_ylim(-limit, limit)
        axis.set_aspect("equal")
        axis.set_xlabel(name_a)
        axis.set_ylabel(name_b)
        axis.set_title(f"{name_a}{name_b} projection")
        axis.grid(alpha=0.25)

    figure.tight_layout()

    return figure


# --- Helper functions ------------------------------------------------------

def _frustums(cameras) -> np.ndarray:
    """Line segments of every camera frustum, in world coordinates."""

    K = cameras.calibration
    focal = np.array([K[0, 0], K[1, 1]])
    principal = np.array([K[0, 2], K[1, 2]])
    depth = FRUSTUM * np.linalg.norm(cameras.centre, axis=1).max()

    # Image corners, the image being twice the principal point
    corner = np.array([[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]])
    corner = corner * principal

    # Back project the corners to the frustum depth
    local = np.ones((4, 3)) * depth
    local[:, :2] = (corner - principal) / focal * depth

    # Camera to world is the transpose of the world to camera rotation
    tip = cameras.centre[:, None, :]
    far = tip + np.einsum("nji,kj->nki", cameras.rotation, local)

    # Four edges from the centre to the corners, and the four of the rectangle
    edges = np.concatenate([
        np.stack([np.broadcast_to(tip, far.shape), far], axis=2),
        np.stack([far, np.roll(far, -1, axis=1)], axis=2),
    ], axis=1)

    return edges.reshape(-1, 2, 3)
