# Title: Scene generation and plotting
# Author: Alessandro Carlino
# Date: 2026/09/18

'''
Generates a synthetic scene before and after objects are added, removed or
moved, builds the camera locus around it, runs the detection pass on the new
scene, then lets the planner move the camera leg after leg, each leg planned
on the merged cloud of what is known so far, until the time or the image
budget runs out, and shows every stage in rerun as one interactive
timeline.

Run it from the point_simulation folder, with the virtual environment that
carries the rerun sdk:

    ../../.venv/bin/python run_scene.py
'''

# Import libraries
import numpy as np

from camera_locus import camera_grid
from detection_sweep import detection_sweep
from knowledge import ADDED, REMOVED, merge, observe
from planner import plan_leg, worst_deviation
from scene_viewer import view_scene
from synthetic_point_cloud import generate_scene_pair
from visibility_filter import describe_cloud

# --- Scene parameters ------------------------------------------------------

n_points = 2000   # Points on the background surface
radius = 1.00     # Maximum radius of the scene
n_modes = 60      # Number of modes of the field
decay = 1.60      # Spectral decay (complexity)
fill = 0.30       # Fraction of the scene volume inside the background
seed = 0          # Randomization seed of the background

# --- Change parameters -----------------------------------------------------

n_added = 1          # Objects present only after the change
n_removed = 1        # Objects present only before the change
n_moved = 1          # Objects present in both, at a different pose
object_radius = 0.25 # Radius of an object, as a fraction of the scene one
change_seed = 0      # Randomization seed of the objects and their poses

# --- Camera parameters -----------------------------------------------------

n_cameras = 120     # Cameras over the locus
locus_radius = 3.00 # Radius of the locus, independent of the scene one
focal = 800.0       # Focal length in pixels
resolution = 500    # Square image side in pixels
min_elevation = -10 # Lowest camera elevation over the equator, in degrees
max_elevation = 85  # Highest camera elevation, in degrees
n_sweep = 4         # Poses in the detection pass, spread over the locus

# --- Mission parameters ----------------------------------------------------

speed = 1.0         # Camera speed along the locus, in scene units per second
image_time = 2.0    # Seconds spent slowing down or stopping for an image
time_budget = 40.0  # Seconds available after the sweep
image_budget = 12   # Images available after the sweep

calibration = np.array([
    [focal, 0.0, resolution / 2],
    [0.0, focal, resolution / 2],
    [0.0, 0.0, 1.0],
])

# --- Build the scene -------------------------------------------------------

scene = generate_scene_pair(n_points, radius, n_modes=n_modes, decay=decay,
                            fill=fill, seed=seed, n_added=n_added,
                            n_removed=n_removed, n_moved=n_moved,
                            object_radius=object_radius,
                            change_seed=change_seed)
cameras = camera_grid(n_cameras, locus_radius, calibration,
                      min_elevation=min_elevation,
                      max_elevation=max_elevation)
old_geometry = describe_cloud(scene.old)
new_geometry = describe_cloud(scene.new)
views, knowledge = detection_sweep(scene, old_geometry, new_geometry,
                                   cameras, n_sweep)
rounds = [(views, knowledge, None)]

# One leg at a time: plan it on what is known, take its images, let them
# update what is known, and plan the next one from where the camera stopped
position = int(views[-1])
time_left, images_left = time_budget, image_budget
while True:
    step = plan_leg(merge(scene, knowledge), knowledge, cameras, position,
                    time_left, images_left, speed, image_time)
    if step is None:
        break

    knowledge = observe(knowledge, scene, old_geometry, new_geometry,
                        cameras, step.captures)
    rounds.append((step.captures, knowledge, step.path))

    time_left -= step.time
    images_left -= len(step.captures)
    position = step.stop

# --- Print a summary -------------------------------------------------------

elevation = np.degrees(np.arcsin(cameras.centre[:, 2] / locus_radius))

print(f"Scene: {n_points} background points, n_modes={n_modes}, "
      f"decay={decay}, fill={fill}, seed={seed}")
print(f"  change:      {n_added} added, {n_removed} removed, {n_moved} moved")
print(f"  old scene:   {len(scene.old)} points, {scene.removed.sum()} removed")
print(f"  new scene:   {len(scene.new)} points, {scene.added.sum()} added")
print(f"  cameras:     {len(cameras.centre)} at radius {locus_radius}, "
      f"elevation {elevation.min():.1f} to {elevation.max():.1f} degrees")
radius = float(np.linalg.norm(scene.new, axis=1).max())
for number, (taken, known, _) in enumerate(rounds):
    deviation = worst_deviation(known.information[scene.added], radius)
    label = "sweep" if number == 0 else f"leg {number}"
    print(f"  {label:<11}  images {taken.tolist()}: "
          f"{known.added.sum()} of {scene.added.sum()} added known, "
          f"{known.removed.sum()} of {scene.removed.sum()} removed, "
          f"median worst deviation {np.median(deviation) / radius:.4f} R")

print(f"  spent:       {time_budget - time_left:.1f} of {time_budget:.0f} s, "
      f"{image_budget - images_left} of {image_budget} images")

merged = merge(scene, knowledge)
print(f"  merged:      {len(merged.points)} points, "
      f"{(merged.status == REMOVED).sum()} removed and "
      f"{(merged.status == ADDED).sum()} added known")
print(f"  still wrong: {(scene.removed & ~knowledge.removed).sum()} removed "
      f"points believed present, {(scene.added & ~knowledge.added).sum()} "
      f"added points unknown")

# --- Visualization ---------------------------------------------------------

view_scene(scene, cameras, rounds)
