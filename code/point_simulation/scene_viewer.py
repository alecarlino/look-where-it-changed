# Title: Interactive scene viewer
# Author: Alessandro Carlino
# Date: 2026/09/18

'''
Logs the scene before and after the change, the camera locus, and every
leg of the camera path with the images it took and the change they made
known, to rerun, which shows them in a single interactive window.

The stages are steps of one timeline, in the order things happen: the old
scene, the new one, the change between them, the camera locus, then one step
per round: the detection sweep first, then every leg of the path, drawn as
the line the camera followed. At each round the images just taken are added
and the known change is logged again, so scrubbing the timeline shows the
camera moving and the knowledge growing. Nothing is ever
cleared, every layer stays visible at the end of the timeline, which is
where the viewer opens, and any of them can be switched off by hand from the
side panel. Switching the old scene off leaves the world as it is now.

The two scenes share their background points, which are drawn in the same
colour so that the layers coincide where nothing changed. Added points are
drawn over the new scene and removed ones over the old; the change the views
have made known is drawn over both, so what keeps its first colour is change
nobody has observed yet.

Every candidate pose is logged, and the sweep and the planned poses, both
subsets of it, are logged over it in their own colours, at the same frustum
size so the highlight sits on the pose it refers to rather than framing it.
Two frustums of exactly the same geometry leave the renderer to pick one,
and it does not reliably pick the later, so a highlight is drawn with a
thicker line and a barely longer plane distance: enough to win the tie, not
enough to see.

Cameras are logged as a transform plus a pinhole, so the frustums are drawn
by the viewer from the calibration itself instead of being built by hand as
matplotlib requires. Rerun wants the world from camera transform, which is
the transpose of the world to camera rotation the locus produces, and its
pinhole assumes the right down forward convention, the same OpenCV one used
throughout.

A blueprint is sent with every recording. Rerun remembers the layout it built
for an application the first time it saw one, so an entity added or renamed
later can stay missing from the side panel until the stored layout is reset
by hand. Sending one makes the tree depend on this run only.

See notes.md for the conventions.
'''

# Import libraries
import numpy as np
import rerun as rr
import rerun.blueprint as rrb

# --- Constants -------------------------------------------------------------

CLOUD = [130, 150, 175]     # Colour of the cloud as a whole
ADDED = [220, 50, 50]       # Colour of the added surface
REMOVED = [150, 70, 200]    # Colour of the removed surface
CANDIDATE = [120, 120, 120] # Colour of a candidate camera
SWEEP = [60, 120, 230]      # Colour of a detection sweep pose
PLANNED = [60, 180, 90]     # Colour of a pose chosen by the planner
KNOWN_ADDED = [250, 190, 40]   # Colour of an addition already seen
KNOWN_REMOVED = [40, 200, 220] # Colour of a removal already seen past
FRUSTUM = 0.10              # Frustum depth, as a fraction of the locus radius
TIEBREAK = 1.02             # Extra depth per highlight level, wins overlaps
STROKE = -2.5               # Line width of the highlight, negative for pixels
DOT = 0.005                 # Point radius, as a fraction of the cloud radius
TIMELINE = "stage"          # Timeline the steps live on


# --- Show the scene --------------------------------------------------------

def view_scene(scene, cameras, rounds: list,
               save_path: str | None = None) -> None:
    """
    Log the scene to rerun as steps of a single timeline.
    Args:
        scene: A ScenePair, the scene before and after the change.
        cameras: A CameraGrid.
        rounds: Triples of the images taken in a round, the Knowledge after
            it, and the camera path of the leg, None for the detection sweep,
            which comes first.
        save_path: Write a recording there instead of opening the viewer.
            A recording can be reopened later or passed to someone else.
    """

    rr.init("point_simulation", spawn=save_path is None)
    if save_path is not None:
        rr.save(save_path)

    # One view holding everything, overriding whatever layout the viewer
    # stored for this application when the entities last looked different
    rr.send_blueprint(rrb.Blueprint(
        rrb.Spatial3DView(origin="/", contents="/**", name="scene"),
        rrb.TimePanel(state="expanded"),
    ))

    # Z up, so the elevation band of the locus reads the way it is defined
    rr.log("/", rr.ViewCoordinates.RIGHT_HAND_Z_UP, static=True)

    # A drawing size, not a measurement: the points are a surface, so a
    # small fraction of the cloud radius reads well at any scale
    dot = DOT * float(np.linalg.norm(scene.new, axis=1).max())

    # Steps one and two, the scene as the model knows it and as it is now
    rr.set_time(TIMELINE, sequence=0)
    rr.log("scene/old", rr.Points3D(scene.old, colors=CLOUD, radii=dot))
    rr.set_time(TIMELINE, sequence=1)
    rr.log("scene/new", rr.Points3D(scene.new, colors=CLOUD, radii=dot))

    # Step three, the change, wider so it draws over the scenes
    rr.set_time(TIMELINE, sequence=2)
    rr.log("scene/added", rr.Points3D(
        scene.new[scene.added], colors=ADDED, radii=2.0 * dot
    ))
    rr.log("scene/removed", rr.Points3D(
        scene.old[scene.removed], colors=REMOVED, radii=2.0 * dot
    ))

    # Step four, the whole locus
    rr.set_time(TIMELINE, sequence=3)
    _log_cameras("cameras/candidate", cameras,
                 np.ones(len(cameras.centre), bool), CANDIDATE)

    # One step per round: the views just taken, over the candidates at the
    # same size so the highlight covers the pose it refers to, and the known
    # change after them, over the ground truth and wider, so the colour of
    # the change still showing is what nobody has observed yet
    for number, (views, knowledge, path) in enumerate(rounds):
        rr.set_time(TIMELINE, sequence=4 + number)

        taken = np.zeros(len(cameras.centre), dtype=bool)
        taken[views] = True
        if number == 0:
            _log_cameras("cameras/sweep", cameras, taken, SWEEP, highlight=1)
        else:
            _log_cameras("cameras/planned", cameras, taken, PLANNED,
                         highlight=1)
            rr.log(f"trajectory/leg_{number:02d}", rr.LineStrips3D(
                [path], colors=PLANNED, radii=STROKE
            ))

        rr.log("known/added", rr.Points3D(
            scene.new[knowledge.added], colors=KNOWN_ADDED, radii=3.0 * dot
        ))
        rr.log("known/removed", rr.Points3D(
            scene.old[knowledge.removed], colors=KNOWN_REMOVED,
            radii=3.0 * dot
        ))

# --- Helper functions ------------------------------------------------------

def _log_cameras(path: str, cameras, keep: np.ndarray, colour: list,
                 highlight: int = 0) -> None:
    """Log the kept cameras as a transform and a pinhole each."""

    width = int(round(2 * cameras.calibration[0, 2]))
    height = int(round(2 * cameras.calibration[1, 2]))

    # The default image plane sits one world unit away, which at this scale
    # buries the scene under a cage of overlapping frustums
    depth = FRUSTUM * float(np.linalg.norm(cameras.centre, axis=1).max())
    depth *= TIEBREAK ** highlight

    for index in np.flatnonzero(keep):
        entity = f"{path}/{index:04d}"

        # Rerun wants the transform of the entity in its parent space, that is
        # world from camera, the transpose of the world to camera rotation
        rr.log(entity, rr.Transform3D(
            translation=cameras.centre[index],
            mat3x3=cameras.rotation[index].T,
        ))

        rr.log(entity, rr.Pinhole(
            image_from_camera=cameras.calibration,
            resolution=[width, height],
            camera_xyz=rr.ViewCoordinates.RDF,
            image_plane_distance=depth,
            color=colour,
            line_width=STROKE if highlight else None,
        ))
