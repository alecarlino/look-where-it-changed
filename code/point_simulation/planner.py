# Title: View planner
# Author: Alessandro Carlino
# Date: 2026/09/19

'''
Given what is known about the change and where the camera stands, it plans
the next leg of the path: where to stop, and which poses to photograph on
the way.

The information objective is the D optimal one, summed over the known added
points,

    J(S) = sum_j log det( A_j + sum_{i in S} F_ij )

where A_j is the information already gathered on point j plus a weak prior,
and F_ij is what an image from pose i would add. A leg takes time, the
travel along it plus a fixed time for every image, since the camera slows
down or stops to take it, and it takes images. Its cost is the share of each
budget it uses up, each measured against what is left of it, so that the
scarcer budget weighs more: with time to spare the planner maximises the
gain per image and travel is almost free, with images to spare it maximises
the gain per second. A leg is scored by its efficiency, gain over cost, and
the planner picks the most efficient leg that fits both budgets.

The poses passed by on a leg are cheap, they cost only the time of their
image and no detour. One is photographed only if it raises the efficiency of
the leg, so a pose whose gain is not worth the stop is let go by. Poses
already taken are never chosen as a stop, since going back on purpose to a
view already held is a waste, but they may be photographed again when
passed by: a retake is an independent measurement and still adds some
information.

Log det is submodular, but with a cost that depends on the order of the
stops the greedy on efficiency no longer carries the 1 - 1/e guarantee of
the cardinality constrained case; it is a heuristic here.

The planner sees the scene only through the merged cloud: what an image
would see is predicted on the occluders it believes in, never on the true
scene. It only refines change that is already known; a leg that happens to
see unknown change reveals it to the next one.

See notes.md for the reasoning and the limits.
'''

# Import libraries
from typing import NamedTuple

import numpy as np

from knowledge import ADDED, Knowledge, MergedCloud, occluding
from measurement import projection_information
from trajectory import leg, length, passed_poses
from visibility_filter import describe_cloud, visible_points

# --- Constants -------------------------------------------------------------

PRIOR = 1.0 # Prior standard deviation of a point position, in cloud radii

class Leg(NamedTuple):
    """
    One leg of the path: where it stops and what it photographs.
    """

    stop: int             # Index into the locus of the pose it ends at
    captures: np.ndarray  # Indices of the poses photographed, in path order
    path: np.ndarray      # (n_samples, 3), camera positions along the leg
    time: float           # Travel plus capture time

# --- Plan a leg ------------------------------------------------------------

def plan_leg(merged: MergedCloud, knowledge: Knowledge, cameras,
             position: int, time_left: float, images_left: int,
             speed: float, image_time: float) -> Leg | None:
    """
    The most efficient next leg that fits the budgets, or None if none does.
    Args:
        merged: The merged cloud, what the planner knows of the scene.
        knowledge: The views taken and the information they gathered.
        cameras: A CameraGrid, the locus the camera moves on.
        position: Index into the locus of the pose the camera is at.
        time_left: Time budget still available.
        images_left: Image budget still available.
        speed: Speed of the camera along the locus.
        image_time: Time spent slowing down or stopping for an image.
    Returns:
        The chosen Leg.
    """

    if speed <= 0 or image_time < 0:
        raise ValueError("speed must be positive and image_time not negative")
    if images_left < 1:
        return None

    current = knowledge.information[knowledge.added] + _prior(merged)
    effect = _predict(merged, cameras)

    best, best_rate = None, 0.0
    for stop in np.flatnonzero(~knowledge.acquired):
        if stop == position:
            continue

        path = leg(cameras.centre[position], cameras.centre[stop])
        travel = length(path) / speed
        if travel + image_time > time_left:
            continue

        # The stop is always photographed, the poses on the way only when
        # their image makes the leg more efficient
        trial = current.copy()
        gain = _take(trial, effect[stop])
        captures = [stop]
        passed = passed_poses(path, cameras)
        on_way = [p for p in passed if p not in (position, stop)]

        def cost(images: int) -> float:
            spent = travel + images * image_time
            return spent / time_left + images / images_left

        while on_way and len(captures) < images_left:
            if travel + (len(captures) + 1) * image_time > time_left:
                break
            extra = [_gain(trial[effect[p][0]], effect[p][1]) for p in on_way]
            pick = int(np.argmax(extra))
            now, then = cost(len(captures)), cost(len(captures) + 1)
            if (gain + extra[pick]) / then <= gain / now:
                break
            gain += _take(trial, effect[on_way[pick]])
            captures.append(on_way.pop(pick))

        rate = gain / cost(len(captures))
        if rate > best_rate:
            order = [p for p in passed if p in captures]
            best_rate = rate
            best = Leg(int(stop), np.array(order, dtype=int), path,
                       travel + len(captures) * image_time)

    return best

# --- Evaluate --------------------------------------------------------------

def worst_deviation(information: np.ndarray, radius: float) -> np.ndarray:
    """
    Standard deviation of each point along its least known direction.
    Args:
        information: Fisher information blocks, shape (n_points, 3, 3).
        radius: Radius of the cloud, which sets the prior.
    Returns:
        Array (n_points,), equal to the prior where nothing was measured.
    """

    prior = np.eye(3) / (PRIOR * radius) ** 2
    return 1.0 / np.sqrt(np.linalg.eigvalsh(information + prior)[:, 0])

# --- Helper functions ------------------------------------------------------

def _prior(merged: MergedCloud) -> np.ndarray:
    """
    Weak prior information, so that a point seen from one side is invertible.
    """

    radius = float(np.linalg.norm(merged.points, axis=1).max())
    return np.eye(3) / (PRIOR * radius) ** 2

def _predict(merged: MergedCloud, cameras) -> list:
    """
    For every pose, the known added points it would see and what it would add.
    The prediction runs on the occluders the planner believes in, never on the
    true scene. The known added points come last in the merged cloud, in the
    same order as their information in the knowledge.
    """

    solid = occluding(merged)
    believed = merged.points[solid]
    geometry = describe_cloud(believed)
    target = merged.status[solid] == ADDED

    effect = []
    for centre, rotation in zip(cameras.centre, cameras.rotation):
        pose = (centre, rotation, cameras.calibration)
        visible = visible_points(believed, geometry, *pose)[target]
        seen = np.flatnonzero(visible)
        effect.append((seen, projection_information(
            believed[target][seen], *pose
        )))

    return effect

def _take(current: np.ndarray, effect: tuple) -> float:
    """
    Add the information of one image in place and return what it gained.
    """

    seen, blocks = effect
    gain = _gain(current[seen], blocks)
    current[seen] += blocks
    return gain

def _gain(current: np.ndarray, blocks: np.ndarray) -> float:
    """
    Increase of the log det objective if these blocks were added.
    """

    if len(blocks) == 0:
        return 0.0

    before = np.linalg.slogdet(current)[1]
    after = np.linalg.slogdet(current + blocks)[1]
    return float((after - before).sum())
