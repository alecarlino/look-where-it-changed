# Title: Synthetic point cloud generator
# Author: Alessandro Carlino
# Date: 2026/09/15

'''
Given the desired number of points and the maximum radius of the scene, it
generates a random tridimensional point cloud lying on the surface of a
procedurally generated object. Points lie on a surface to reproduce a depth
sensor measurements.
The object is the excursion set of a Gaussian process with a Matérn
covariance and a mean falling towards the border, the region where a sample
of it exceeds a level, cut to a ball. The falling mean keeps the object
inside the ball, which then only closes the rare part that reaches it.
The sample is drawn with random Fourier features, a sum of cosines whose
frequencies come from the spectral density of the kernel. The length scale
sets the size of the features and the smoothness how rough they are; the
level follows from the expected fraction of the ball to be inside, since
the process is Gaussian with unit variance. Every piece is a published
construction: random Fourier features (Rahimi and Recht, 2007), the Matérn
kernel, its spectral density and a non zero mean function (Rasmussen and
Williams, 2006), excursion sets (Adler and Taylor, 2007) and intersection by
minimum (Ricci, 1973).
The scene is generated two times, before and after a change with objects
added, removed or moved.
Every part is a field shifted, solid where positive. As a consequence a scene
is the maximum of its parts, so an object is added or removed by adding or
dropping  a term. The movement is given by move is a removal at one pose plus
an addition at another.
The change is measured as the geometric difference of the two scenes, so that
an object resting on a surface removes the surface it covers, and taking it
away reveals it.
'''

# Import libraries
from typing import NamedTuple

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm

# --- Constants -------------------------------------------------------------

MEAN_DROP = 4.0     # Fall of the mean from centre to border, in deviations
SHELL = 0.02        # Half thickness of the sampling shell, over the radius
STEPS = 4           # Newton steps used to land on the isosurface
LANDED = 1e-9       # Distance from the surface counted as landed
MAX_BATCHES = 200   # Candidate batches drawn before giving up
OBJECT_FILL = 0.45  # Fraction of its own ball an object fills
OBJECT_PROBE = 256  # Points drawn first to measure the area of an object
MIN_OBJECT = 50     # Fewest points an object is sampled with
UPWARD = 0.3        # Least vertical normal of a surface an object rests on
CONTACT = (1.0, 0.8, 0.6, 0.4, 0.2) # Offsets tried to rest an object, in radii
BURIED = 0.02       # Fraction of an object sunk into the surface for contact
ATTEMPTS = 200      # Poses tried before giving up on placing an object
TOLERANCE = 0.5     # Distance counted as change, in point spacings

class ScenePair(NamedTuple):
    """
    The scene as the outdated model knows it and as it is now.
    """

    old: np.ndarray     # (n_old, 3), the scene before the change
    new: np.ndarray     # (n_new, 3), the scene after the change
    added: np.ndarray   # (n_new,), new points away from the old surface
    removed: np.ndarray # (n_old,), old points away from the new surface

# --- Generate a scene pair -------------------------------------------------

def generate_scene_pair(
    n_points: int,
    radius: float,
    n_modes: int = 256,
    length: float = 0.3,
    smoothness: float = 2.5,
    fill: float = 0.30,
    seed: int | None = None,
    n_added: int = 1,
    n_removed: int = 0,
    n_moved: int = 0,
    object_radius: float = 0.25,
    change_seed: int | None = None,
) -> ScenePair:
    """
    Generate a scene before and after objects are added, removed or moved.
    Args:
        n_points: Points on the background surface.
        radius: Radius of the sphere the scene lives in.
        n_modes: Random Fourier features summed to build each field.
        length: Length scale of the features, over the radius.
        smoothness: Matérn smoothness, lower is rougher.
        fill: Fraction of the scene volume inside the background.
        seed: Seed of the background, the same in both scenes.
        n_added: Objects present only after the change.
        n_removed: Objects present only before the change.
        n_moved: Objects present in both, at a different pose.
        object_radius: Radius of an object, as a fraction of the scene one.
        change_seed: Seed of the objects and their poses.
    Returns:
        A ScenePair, the change being the geometric difference of the two.
    Raises:
        RuntimeError: If an object finds no free surface to rest on.
    """

    # Inputs check
    _check_scene(n_points, radius, n_modes, length, smoothness, fill)
    if min(n_added, n_removed, n_moved) < 0:
        raise ValueError("object counts must be non negative")
    if not 0 < object_radius < 0.5:
        raise ValueError(
            f"object_radius must lie in (0, 0.5), got {object_radius}"
        )

    # Background, identical in both scenes and sampled once for both
    rng = np.random.default_rng(seed)
    background = _component(rng, n_modes, length, smoothness, radius, fill)
    base, area = _sample_surface(rng, background, n_points, radius)
    spacing = np.sqrt(area / n_points)

    # Objects share the density of the background, so every point stands for
    # the same patch of surface wherever it lies
    change = np.random.default_rng(change_seed)
    size = object_radius * radius
    shapes = [
        _object_shape(change, n_modes, length, smoothness, size,
                      n_points / area)
        for _ in range(n_added + n_removed + n_moved)
    ]

    # A move is a removal at one pose plus an addition at another
    roles = (
        [(k, True, False) for k in range(n_removed)]
        + [(n_removed + k, False, True) for k in range(n_added)]
        + [(n_removed + n_added + k, flag, not flag)
           for k in range(n_moved) for flag in (True, False)]
    )

    gradient = background(base)[1]
    outward = -gradient / np.linalg.norm(gradient, axis=1, keepdims=True)

    old_parts, new_parts, centres = [], [], []
    for shape, in_old, in_new in roles:
        placed = _place_object(change, background, base, outward,
                               shapes[shape], size, radius, centres)
        centres.append(placed[2])
        if in_old:
            old_parts.append(placed[:2])
        if in_new:
            new_parts.append(placed[:2])

    old = _assemble(background, base, old_parts)
    new = _assemble(background, base, new_parts)

    # The change is geometric: a point is changed when the other scene has no
    # surface within tolerance of it, so a short move flags only the part of
    # the object that does not overlap its old pose
    old_field = _union([background] + [part[0] for part in old_parts])
    new_field = _union([background] + [part[0] for part in new_parts])
    added = _distance(old_field, new) > TOLERANCE * spacing
    removed = _distance(new_field, old) > TOLERANCE * spacing

    return ScenePair(old, new, added, removed)

# --- Helper functions ------------------------------------------------------

def _check_scene(n_points: int, radius: float, n_modes: int,
                 length: float, smoothness: float,
                 fill: float) -> None:
    """
    Validate the parameters shared by both generators.
    """

    if n_points <= 0:
        raise ValueError(f"n_points must be positive, got {n_points}")
    if radius <= 0:
        raise ValueError(f"radius must be positive, got {radius}")
    if n_modes < 1:
        raise ValueError(f"n_modes must be at least 1, got {n_modes}")
    if length <= 0:
        raise ValueError(f"length must be positive, got {length}")
    if smoothness <= 0:
        raise ValueError(f"smoothness must be positive, got {smoothness}")
    if not 0 < fill < 1:
        raise ValueError(f"fill must lie in (0, 1), got {fill}")

def _component(rng: np.random.Generator, n_modes: int, length: float,
               smoothness: float, radius: float, fill: float):
    """
    An excursion set cut to a ball, as a field solid where it is positive.
    """

    field = _random_field(rng, n_modes, length, smoothness, radius)

    # The mean falls as the square of the distance from the centre, so the
    # process sits above the level mostly inside and the object stays there
    def mean(r: float) -> float:
        return -MEAN_DROP * r ** 2

    # A point at fraction r of the radius is inside with probability
    # 1 - Phi(level - mean), and the level is the one whose average over the
    # ball, weighted by the volume 3 r^2 dr, is the fraction fill
    def inside(level: float) -> float:
        return 3.0 * quad(
            lambda r: (1.0 - norm.cdf(level - mean(r))) * r ** 2, 0.0, 1.0
        )[0]

    level = brentq(lambda level: inside(level) - fill, -10.0, 10.0)

    def component(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        value, gradient = field(x)
        distance = np.linalg.norm(x, axis=1)
        value = value + mean(distance / radius) - level
        gradient = gradient - 2.0 * MEAN_DROP * x / radius ** 2

        # Intersection with the ball is the minimum of the two, so the solid
        # is closed by the ball wherever the excursion set reaches it
        wall = radius - distance
        cut = wall < value
        outward = x / np.maximum(distance, 1e-12)[:, None]

        return (np.where(cut, wall, value),
                np.where(cut[:, None], -outward, gradient))

    return component

def _sample_surface(rng: np.random.Generator, field, count: int,
                    radius: float) -> tuple[np.ndarray, float]:
    """
    Points on the zero level of a field, and the area it measures on the way.
    """

    points = np.empty((count, 3))
    filled = drawn = near = 0

    # Candidates come from a ball one shell wider, so that the cap on the
    # sphere gets a full shell around it and not only its inner half
    outer = radius * (1.0 + SHELL)

    # Rejection sampling in batches, the acceptance rate not being known
    for _batch in range(MAX_BATCHES):
        missing = count - filled
        x = _uniform_in_ball(rng, max(8 * missing, 4096), outer)
        value, gradient = field(x)

        # Constant thickness shell, so the surface density stays uniform
        slope = np.linalg.norm(gradient, axis=1)
        thickness = SHELL * radius * np.maximum(slope, 1e-12)
        shell = np.abs(value) <= thickness
        drawn += len(x)
        near += int(shell.sum())
        x = _project(x[shell], field)

        # Near the crease where the excursion set meets the cap, Newton can
        # bounce between the two and never land, so those points are dropped
        value, gradient = field(x)
        landed = np.abs(value) <= LANDED * radius * np.maximum(
            np.linalg.norm(gradient, axis=1), 1e-12
        )

        # Points on the cap sit on the sphere itself, up to rounding
        x = x[landed & (np.linalg.norm(x, axis=1) <= radius * (1.0 + 1e-9))]

        take = min(len(x), missing)
        points[filled:filled + take] = x[:take]
        filled += take

        if filled == count:
            # The shell holds a volume twice its half thickness times the area
            volume = 4.0 / 3.0 * np.pi * outer ** 3
            area = near / drawn * volume / (2.0 * SHELL * radius)
            return points, area

    raise RuntimeError(
        f"could only place {filled} of {count} points on the isosurface"
    )

def _object_shape(rng: np.random.Generator, n_modes: int, length: float,
                  smoothness: float, size: float, density: float):
    """
    A closed random object in its own ball, sampled at the given density.
    """

    field = _component(rng, n_modes, length, smoothness, size,
                       OBJECT_FILL)
    points, area = _sample_surface(rng, field, OBJECT_PROBE, size)

    # The first batch measures the area, the second tops up to the density
    count = max(MIN_OBJECT, int(round(density * area)))
    if count > OBJECT_PROBE:
        extra = _sample_surface(rng, field, count - OBJECT_PROBE, size)[0]
        points = np.vstack([points, extra])

    return field, points[:count]

def _place_object(rng: np.random.Generator, background, base: np.ndarray,
                  outward: np.ndarray, shape, size: float, radius: float,
                  centres: list) -> tuple:
    """
    Rest an object on the background, clear of the objects already placed.
    """

    field, local = shape

    # Surfaces facing up, where an object would rest, or any if there is none
    anchors = np.flatnonzero(outward[:, 2] > UPWARD)
    if len(anchors) == 0:
        anchors = np.arange(len(base))

    for _attempt in range(ATTEMPTS):
        anchor = anchors[rng.integers(len(anchors))]
        rotation = np.linalg.qr(rng.standard_normal((3, 3)))[0]
        rotation *= np.sign(np.linalg.det(rotation))

        # Slide it down along the normal until part of it sinks in the surface
        for offset in CONTACT:
            # The background fills the scene sphere up to its border, so an
            # object resting on top may stick out of it by up to its radius
            centre = base[anchor] + offset * size * outward[anchor]
            if np.linalg.norm(centre) > radius:
                continue
            if any(np.linalg.norm(centre - c) < 2.0 * size for c in centres):
                continue

            points = local @ rotation.T + centre
            if (background(points)[0] > 0.0).mean() >= BURIED:
                return _placed(field, rotation, centre), points, centre

    raise RuntimeError(
        f"could not rest an object of radius {size:.3g} after {ATTEMPTS} "
        f"attempts; reduce the number of objects or object_radius"
    )

def _placed(field, rotation: np.ndarray, centre: np.ndarray):
    """
    A field moved to a pose, evaluated in the coordinates of the object.
    """

    def moved(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        value, gradient = field((x - centre) @ rotation)
        return value, gradient @ rotation.T

    return moved

def _union(fields: list):
    """
    Solid union of fields, which is their maximum.
    """

    def union(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        results = [field(x) for field in fields]
        values = np.stack([result[0] for result in results])
        best = np.argmax(values, axis=0)
        column = np.arange(len(x))
        gradient = np.stack([result[1] for result in results])[best, column]
        return values[best, column], gradient

    return union

def _assemble(background, base: np.ndarray, parts: list) -> np.ndarray:
    """
    Surface points of a scene: those of each part outside every other part.
    """

    fields = [background] + [part[0] for part in parts]
    clouds = [base] + [part[1] for part in parts]

    kept = []
    for index, cloud in enumerate(clouds):
        outside = np.ones(len(cloud), dtype=bool)
        for other, field in enumerate(fields):
            if other != index:
                outside &= field(cloud)[0] < 0.0
        kept.append(cloud[outside])

    return np.vstack(kept)

def _distance(field, x: np.ndarray) -> np.ndarray:
    """
    First order distance from points to the zero level of a field.
    """

    value, gradient = field(x)
    return np.abs(value) / np.maximum(np.linalg.norm(gradient, axis=1), 1e-12)

def _random_field(rng: np.random.Generator, n_modes: int, length: float,
                  smoothness: float, radius: float):
    """
    A sample of a unit variance Gaussian process with a Matérn covariance.
    Random Fourier features: a sum of cosines with uniform random phases and
    frequencies drawn from the spectral density of the kernel, which for a
    Matérn kernel is a multivariate t with twice the smoothness as degrees of
    freedom. Returns a function giving the value and the gradient of the
    sample at a set of points, with the features drawn once.
    """

    # Multivariate t: a Gaussian divided by an independent chi square root
    scale = np.sqrt(2.0 * smoothness / rng.chisquare(2.0 * smoothness,
                                                      n_modes))
    waves = (rng.standard_normal((n_modes, 3)) * scale[:, None]
             / (length * radius))
    phase = rng.random(n_modes) * 2.0 * np.pi
    weight = np.sqrt(2.0 / n_modes)

    def field(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        projection = x @ waves.T + phase
        return (weight * np.cos(projection).sum(axis=1),
                -weight * np.sin(projection) @ waves)

    return field

def _project(x: np.ndarray, field) -> np.ndarray:
    """
    Newton iterations along the gradient, to land on the isosurface.
    """

    # Newton on the scalar equation, stepping along the level set normal
    for _ in range(STEPS):
        value, gradient = field(x)
        squared = np.einsum("ij,ij->i", gradient, gradient)
        step = value / np.maximum(squared, 1e-12)
        x = x - step[:, None] * gradient

    return x

def _uniform_in_ball(rng: np.random.Generator, count: int,
                     radius: float) -> np.ndarray:
    """
    Points uniform in the volume of a ball centred on the origin.
    """

    directions = rng.standard_normal((count, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)

    # Radius from inverting the cumulative (r / radius) ^ 3
    return directions * (radius * rng.random(count) ** (1.0 / 3.0))[:, None]