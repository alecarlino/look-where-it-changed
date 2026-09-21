# Title: Known change and merged cloud
# Author: Alessandro Carlino
# Date: 2026/09/19

'''
Keeps track of what the system knows about the change, and turns it into the
single cloud the planner works on.

Detection is not simulated but replaced by a proxy: a detector that makes no
mistake on what the acquired views observe, and knows nothing about the
rest. A view reveals the added points it sees and confirms the removed points
it sees past, and it also measures the added points it sees, so the Fisher
information on each of them builds up with every view that looks at it. The
same proxy serves the detection sweep and every view the planner acquires,
so the known change and the information on it grow with the views taken.

The merged cloud is the old cloud, which is what the outdated model knows,
plus the added points observed so far, each with a status. A removed point
not yet confirmed is still believed to be there, and an added part not yet
seen is not there at all: the planner sees the scene through what it knows,
and is wrong exactly where the change has not been observed.

See notes.md for the reasoning and the limits.
'''

# Import libraries
from typing import NamedTuple

import numpy as np

from measurement import projection_information
from synthetic_point_cloud import ScenePair
from visibility_filter import CloudGeometry, vacated_points, visible_points

# --- Constants -------------------------------------------------------------

UNCHANGED = 0 # Old point believed to be still there
REMOVED = 1   # Old point a view has seen past
ADDED = 2     # New point a view has seen

class Knowledge(NamedTuple):
    """
    Views taken so far and the change they have observed.
    """

    acquired: np.ndarray    # (n_cameras,), True where the pose has been taken
    added: np.ndarray       # (n_new,), added points observed so far
    removed: np.ndarray     # (n_old,), removed points confirmed so far
    information: np.ndarray # (n_new, 3, 3), Fisher information on new points

class MergedCloud(NamedTuple):
    """
    The cloud the planner works on, with the status of every point.
    """

    points: np.ndarray # (n_points, 3), the old cloud then the added points
    status: np.ndarray # (n_points,), UNCHANGED, REMOVED or ADDED

# --- Known change ----------------------------------------------------------

def no_knowledge(scene: ScenePair, cameras) -> Knowledge:
    """
    The state before any view: nothing taken, no change known.
    """

    return Knowledge(np.zeros(len(cameras.centre), dtype=bool),
                     np.zeros(len(scene.new), dtype=bool),
                     np.zeros(len(scene.old), dtype=bool),
                     np.zeros((len(scene.new), 3, 3)))

def observe(knowledge: Knowledge, scene: ScenePair,
            old_geometry: CloudGeometry, new_geometry: CloudGeometry,
            cameras, views: np.ndarray) -> Knowledge:
    """
    Update the known change with what some more views observe.
    Args:
        knowledge: What is known before these views.
        scene: The scene before and after the change, the ground truth.
        old_geometry: Output of describe_cloud for the old scene.
        new_geometry: Output of describe_cloud for the new scene.
        cameras: A CameraGrid.
        views: Indices into the locus of the poses just acquired.
    Returns:
        The updated Knowledge, which only ever grows.
    """

    acquired = knowledge.acquired.copy()
    added = knowledge.added.copy()
    removed = knowledge.removed.copy()
    information = knowledge.information.copy()

    # A view looks at the scene as it is now; the old one is only what the
    # outdated model remembers, checked against what the view sees
    for index in views:
        pose = (cameras.centre[index], cameras.rotation[index],
                cameras.calibration)
        seen = scene.added & visible_points(scene.new, new_geometry, *pose)
        removed |= scene.removed & vacated_points(
            scene.old, old_geometry, scene.new, new_geometry, *pose
        )

        # Every view that sees a point measures it again, retakes included
        information[seen] += projection_information(scene.new[seen], *pose)
        added |= seen
        acquired[index] = True

    return Knowledge(acquired, added, removed, information)

# --- Merged cloud ----------------------------------------------------------

def merge(scene: ScenePair, knowledge: Knowledge) -> MergedCloud:
    """
    The old cloud plus the added points observed so far, with their status.
    """

    status = np.where(knowledge.removed, REMOVED, UNCHANGED)
    points = np.vstack([scene.old, scene.new[knowledge.added]])
    status = np.concatenate([
        status, np.full(knowledge.added.sum(), ADDED)
    ]).astype(np.int8)

    return MergedCloud(points, status)

def occluding(merged: MergedCloud) -> np.ndarray:
    """
    The points the planner believes are solid: all but confirmed removals.
    """

    return merged.status != REMOVED
