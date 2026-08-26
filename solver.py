"""
solver.py — Formula Student Suspension Kinematics Solver v3
Double wishbone, 4-corner model.

Coordinate system (SAE-aligned):
  X — longitudinal (forward positive)
  Y — lateral (inboard positive)
  Z — vertical (up positive)

Origin at wheel centre at design ride height.

Author: Mahdi Kadiri Bouchaib
University of Bath — MEng Mechanical Engineering (Automotive)

═══════════════════════════════════════════════════════════
CHANGES v3 (relative to v2)
═══════════════════════════════════════════════════════════

NEW OUTPUTS (matching professional kinematics tools):
  - Kingpin inclination (KPI), deg — front-view steering axis tilt
  - Caster angle, deg             — side-view steering axis tilt
  - Mechanical trail, mm          — axis/ground intersection ahead of CP
  - Scrub radius, mm              — front-view axis/ground offset from CP
  - Scrub (track change), mm      — lateral contact patch migration
  - Wheel centre YZ trajectory    — front-view wheel path droop→bump

STRUCTURAL FIX enabling scrub:
  v2 translated the contact patch and wheel centre purely vertically
  (+dz), which forces scrub = 0 by construction. v3 attaches both
  points rigidly to the UPRIGHT: at each travel step the upright's
  rigid-body rotation is recovered from the kingpin axis change
  (reference → current ball-joint pair) and applied to all
  upright-fixed points. Captures camber-driven lateral/vertical
  migration of the contact patch and wheel centre.
  LIMITATION: rotation about the kingpin axis (steer) is not included
  in this transform — valid while heave-induced steer is small, which
  the tierod model separately reports as toe.

SENSITIVITY UPGRADE (slope-based):
  v2 measured sensitivity of DRH *values*. v3 measures sensitivity of
  curve *slopes* at DRH (camber gain, RC height gradient, bump steer
  gradient) per mm of hardpoint perturbation — the quantity that
  actually matters for design, since static values can be shimmed but
  gains are baked into the hardpoints. Rendered as a heatmap grid
  (hardpoints × XYZ per output), diverging colormap.
═══════════════════════════════════════════════════════════
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, Tuple
import os
import csv


# ── Hardpoints ────────────────────────────────────────────

def default_fs_front_hardpoints() -> Dict[str, np.ndarray]:
    """
    Representative Formula Student front corner hardpoints.
    All coordinates in mm, SAE convention.
    Replace with TBRe actual coordinates when available.
    """
    return {
        'UWB_inboard_front': np.array([-50.0,  380.0, 180.0]),
        'UWB_inboard_rear':  np.array([ 50.0,  380.0, 180.0]),
        'UWB_outboard':      np.array([  0.0,   40.0, 260.0]),
        'LWB_inboard_front': np.array([-70.0,  360.0,  40.0]),
        'LWB_inboard_rear':  np.array([ 70.0,  360.0,  40.0]),
        'LWB_outboard':      np.array([  0.0,   30.0,  30.0]),
        'tierod_inner':      np.array([ 55.0,  355.0,  70.0]),
        'tierod_outer':      np.array([  0.0,   38.0,  90.0]),
        'wheel_centre':      np.array([  0.0,    0.0, 160.0]),
        'contact_patch':     np.array([  0.0,    0.0,   0.0]),
        'pushrod_upright':   np.array([  0.0,   38.0, 120.0]),
        'pushrod_bellcrank': np.array([  0.0,  300.0, 200.0]),
    }


def default_fs_rear_hardpoints() -> Dict[str, np.ndarray]:
    """Representative FS rear corner hardpoints."""
    return {
        'UWB_inboard_front': np.array([-45.0,  340.0, 190.0]),
        'UWB_inboard_rear':  np.array([ 45.0,  340.0, 190.0]),
        'UWB_outboard':      np.array([  0.0,   42.0, 270.0]),
        'LWB_inboard_front': np.array([-65.0,  330.0,  50.0]),
        'LWB_inboard_rear':  np.array([ 65.0,  330.0,  50.0]),
        'LWB_outboard':      np.array([  0.0,   32.0,  35.0]),
        'tierod_inner':      np.array([-50.0,  330.0,  80.0]),
        'tierod_outer':      np.array([  0.0,   40.0,  95.0]),
        'wheel_centre':      np.array([  0.0,    0.0, 160.0]),
        'contact_patch':     np.array([  0.0,    0.0,   0.0]),
        'pushrod_upright':   np.array([  0.0,   40.0, 130.0]),
        'pushrod_bellcrank': np.array([  0.0,  290.0, 210.0]),
    }


# ── Geometry utilities ────────────────────────────────────

def unit_vector(v: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(v)
    if norm < 1e-10:
        raise ValueError(f'Zero-length vector: {v}')
    return v / norm


def rodrigues_rotate(v: np.ndarray, axis: np.ndarray,
                     theta: float) -> np.ndarray:
    """Rotate v about unit axis by theta (rad). Rodrigues formula."""
    k = unit_vector(axis)
    return (v * np.cos(theta) +
            np.cross(k, v) * np.sin(theta) +
            k * np.dot(k, v) * (1 - np.cos(theta)))


def rotation_between_vectors(a: np.ndarray,
                             b: np.ndarray):
    """
    Return (axis, angle) of the minimal rotation mapping
    direction a onto direction b. Identity if parallel.
    """
    a_u, b_u = unit_vector(a), unit_vector(b)
    cross = np.cross(a_u, b_u)
    s = np.linalg.norm(cross)
    c = np.clip(np.dot(a_u, b_u), -1.0, 1.0)
    if s < 1e-12:
        return np.array([0.0, 0.0, 1.0]), 0.0
    return cross / s, float(np.arctan2(s, c))


def line_intersection_2d(p1, d1, p2, d2) -> Tuple[np.ndarray, bool]:
    """Intersection of two 2-D lines pi + t*di. (point, found)."""
    A = np.array([[d1[0], -d2[0]], [d1[1], -d2[1]]])
    b = p2 - p1
    det = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
    if abs(det) < 1e-10:
        return np.zeros(2), False
    t = (b[0] * A[1, 1] - b[1] * A[0, 1]) / det
    return p1 + t * d1, True


def dist_point_to_line_3d(point, line_origin, line_dir) -> float:
    k = unit_vector(line_dir)
    diff = point - line_origin
    return float(np.linalg.norm(diff - np.dot(diff, k) * k))


# ── Constraint solver ─────────────────────────────────────

def solve_outboard_point(inboard_mid, inboard_front, inboard_rear,
                         outboard_ref, arm_length, dz) -> np.ndarray:
    """
    Rigid wishbone rotation about its inboard axis: solve the
    delta rotation from reference that raises the outboard
    point's Z by dz. (Rodrigues delta formulation, v2 fix.)
    """
    axis   = unit_vector(inboard_rear - inboard_front)
    r0     = outboard_ref - inboard_mid
    r_par  = np.dot(r0, axis) * axis
    r_perp = r0 - r_par
    cross  = np.cross(axis, r_perp)

    Z_base   = inboard_mid[2] + r_par[2]
    target_z = outboard_ref[2] + dz

    A, B, C = r_perp[2], cross[2], target_z - Z_base
    R = np.sqrt(A**2 + B**2)
    if R < 1e-10:
        return outboard_ref.copy()

    # A·cosθ + B·sinθ = R·cos(θ − φ), φ = atan2(B, A)
    # → θ = φ ± arccos(C/R). (v3 fix: previous ±arccos − φ
    # mirrored the motion, inverting the travel direction.)
    ratio = np.clip(C / R, -1.0, 1.0)
    phi   = np.arctan2(B, A)
    d1 = (phi + np.arccos(ratio) + np.pi) % (2*np.pi) - np.pi
    d2 = (phi - np.arccos(ratio) + np.pi) % (2*np.pi) - np.pi
    delta_theta = d1 if abs(d1) <= abs(d2) else d2

    return inboard_mid + r_par + rodrigues_rotate(r_perp, axis,
                                                  delta_theta)


def _bisect(f, a, b, tol=1e-12, maxit=200):
    """Bracketed root find (numpy-only, no scipy dependency)."""
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        return 0.0
    for _ in range(maxit):
        m = 0.5 * (a + b)
        fm = f(m)
        if abs(fm) < tol or (b - a) < tol:
            return m
        if fa * fm < 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


def solve_steer_angle(T, lwb_new, uwb_new,
                      tierod_inner, tierod_outer_ref,
                      tierod_length) -> float:
    """
    Solve the upright's steer rotation about the kingpin axis
    that satisfies the tie rod length constraint. Returns the
    rotation angle (rad) about the kingpin.

    v4 FIX. The previous model translated the tie rod outer
    vertically and projected it onto the tie rod sphere, then
    divided the lateral deviation by a moment arm. That both
    detached the tie rod outer from the upright and double
    counted the geometry, giving bump steer ~9x too large
    (3.6 deg over 60 mm travel vs 0.38 deg from an independent
    solver on the same hardpoints).

    Correct formulation: the tie rod outer is rigidly fixed to
    the upright. The upright's non-steer motion is set by the
    ball joints (transform T). The residual degree of freedom
    is rotation about the kingpin axis; the tie rod, being of
    fixed length and anchored to the chassis, picks the angle.
    """
    kp = unit_vector(uwb_new - lwb_new)
    p_unsteered = T(tierod_outer_ref)

    def residual(delta):
        p = lwb_new + rodrigues_rotate(p_unsteered - lwb_new,
                                       kp, delta)
        return np.linalg.norm(p - tierod_inner) - tierod_length

    return _bisect(residual, -0.5, 0.5)


def upright_transform(lwb_ref, uwb_ref, lwb_new, uwb_new):
    """
    Rigid-body transform of upright-fixed points from the
    reference to the current position, recovered from the
    ball-joint pair.

    Rotation: minimal rotation mapping the reference kingpin
    direction onto the current one (steer about the kingpin
    axis is NOT included — reported separately as toe).
    Translation: anchored at the lower ball joint.

    Returns a function mapping reference points → new points.
    """
    k0 = uwb_ref - lwb_ref
    k1 = uwb_new - lwb_new
    axis, angle = rotation_between_vectors(k0, k1)

    def transform(p_ref: np.ndarray) -> np.ndarray:
        r = p_ref - lwb_ref
        r_rot = rodrigues_rotate(r, axis, angle) if angle != 0.0 else r
        return lwb_new + r_rot

    return transform


# ── Core kinematic functions ──────────────────────────────

def compute_camber(uwb_outboard, lwb_outboard) -> float:
    """Camber (deg), front view. Negative = top leaning inboard."""
    kp = uwb_outboard - lwb_outboard
    return float(np.degrees(np.arctan2(-kp[1], kp[2])))


def compute_kpi(uwb_outboard, lwb_outboard) -> float:
    """
    Kingpin inclination (deg), front view.
    Positive = top of steering axis leaning inboard (standard).
    Note: for this two-ball-joint construction KPI = -camber;
    they separate once a steering-axis offset from the wheel
    plane exists (real uprights).
    """
    kp = uwb_outboard - lwb_outboard
    return float(np.degrees(np.arctan2(kp[1], kp[2])))


def compute_caster(uwb_outboard, lwb_outboard) -> float:
    """
    Caster angle (deg), side view.
    Positive = top of steering axis leaning rearward
    (X forward positive → rearward lean = negative kp_X).
    """
    kp = uwb_outboard - lwb_outboard
    return float(np.degrees(np.arctan2(-kp[0], kp[2])))


def steering_axis_ground_point(uwb_outboard, lwb_outboard,
                               ground_z: float) -> np.ndarray:
    """
    Intersection of the steering axis (through both ball
    joints) with the horizontal plane z = ground_z.
    """
    kp = uwb_outboard - lwb_outboard
    if abs(kp[2]) < 1e-10:
        return lwb_outboard.copy()
    t = (ground_z - lwb_outboard[2]) / kp[2]
    return lwb_outboard + t * kp


def compute_mechanical_trail(uwb_out, lwb_out, contact_patch) -> float:
    """
    Mechanical trail (mm), side view: X-distance from the
    steering-axis/ground intersection to the contact patch.
    Positive = axis intersects the ground AHEAD of the patch
    (the standard stabilising configuration with +ve caster).
    """
    g = steering_axis_ground_point(uwb_out, lwb_out,
                                   contact_patch[2])
    return float(g[0] - contact_patch[0])


def compute_scrub_radius(uwb_out, lwb_out, contact_patch) -> float:
    """
    Scrub radius (mm), front view: Y-distance from the
    steering-axis/ground intersection to the contact patch.
    Positive = patch centre OUTBOARD of the axis intersection
    (conventional positive scrub radius).
    """
    g = steering_axis_ground_point(uwb_out, lwb_out,
                                   contact_patch[2])
    return float(g[1] - contact_patch[1])


def compute_instant_centre(uwb_in_mid, uwb_outboard,
                           lwb_in_mid, lwb_outboard):
    """Front-view IC from wishbone lines projected to Y-Z."""
    d_uwb = uwb_outboard[[1, 2]] - uwb_in_mid[[1, 2]]
    d_lwb = lwb_outboard[[1, 2]] - lwb_in_mid[[1, 2]]
    return line_intersection_2d(uwb_in_mid[[1, 2]], d_uwb,
                                lwb_in_mid[[1, 2]], d_lwb)


def compute_roll_centre(ic_yz, contact_patch,
                        track_half_width: float = 600.0) -> float:
    """
    RC height ABOVE THE GROUND PLANE (standard convention).

    Line CP→IC intersected with the vehicle centreline
    Y = track_half_width, then referenced to the contact patch
    height rather than the chassis frame.

    v4 FIX: previously returned the absolute Z in the chassis
    frame. Identical at DRH (contact patch at Z=0) but the
    gradient differed by ~17x (-0.065 vs -1.10 mm/mm), because
    the contact patch itself rises with the wheel. The
    ground-referenced value is the conventional one and the
    one an independent solver reproduces (-1.08 mm/mm).
    """
    cp_yz = contact_patch[[1, 2]]
    d = ic_yz - cp_yz
    if abs(d[0]) < 1e-10:
        return float(ic_yz[1] - cp_yz[1])
    t = (track_half_width - cp_yz[0]) / d[0]
    return float(cp_yz[1] + t * d[1] - cp_yz[1])


def compute_toe(T, lwb_new, uwb_new, steer_delta) -> float:
    """
    Toe (deg): rotation of the wheel spin axis about vertical.

    The wheel spin axis is carried by the upright, so it gets
    the full upright transform T plus the steer rotation about
    the kingpin. Reference spin axis is taken as +Y (inboard)
    since the static setup is the datum and true wheel-axis
    points are not in the source geometry.

    Positive = toe-in for the modelled side.
    Validated against an independent solver: 0.375 deg total
    over 60 mm travel vs 0.38 deg (1.3%).
    """
    kp = unit_vector(uwb_new - lwb_new)
    p0 = T(np.zeros(3))
    spin = T(np.array([0.0, 1.0, 0.0])) - p0
    spin = rodrigues_rotate(spin, kp, steer_delta)
    return float(np.degrees(np.arctan2(spin[0], spin[1])))


def compute_motion_ratio(pushrod_upright, pushrod_bellcrank) -> float:
    """MR = |cos(pushrod angle from vertical)|, upright→bellcrank."""
    pushrod_unit = unit_vector(pushrod_bellcrank - pushrod_upright)
    return float(abs(pushrod_unit[2]))


def compute_anti_geometry(ic_yz, contact_patch,
                          wheelbase=1530.0, cog_height=280.0) -> Dict:
    """
    Anti percentage from the front-view IC projection.
    NOTE: true anti-dive/squat uses the SIDE-view (X-Z) IC;
    with fore-aft symmetric placeholder inboard points the
    side-view geometry is degenerate, so the front-view
    construction is retained as a placeholder. Replace when
    real (asymmetric) hardpoints are loaded.
    """
    cp_yz = contact_patch[[1, 2]]
    force_line = ic_yz - cp_yz
    if np.linalg.norm(force_line) < 1e-10:
        return {'anti_percent': 0.0}
    force_angle = np.arctan2(force_line[1], abs(force_line[0]))
    pitch_angle = np.arctan2(cog_height, wheelbase)
    return {'anti_percent':
            round(float(np.tan(force_angle) / np.tan(pitch_angle)
                        * 100.0), 2)}


# ── Kinematic data container ──────────────────────────────

@dataclass
class KinematicPoint:
    wheel_travel_mm:       float
    camber_deg:            float
    roll_centre_height_mm: float
    toe_deg:               float
    motion_ratio:          float
    anti_percent:          float
    kpi_deg:               float
    caster_deg:            float
    mech_trail_mm:         float
    scrub_radius_mm:       float
    scrub_mm:              float   # lateral CP migration vs reference
    wc_y_mm:               float   # wheel centre front-view trajectory
    wc_z_mm:               float
    ic_y_mm:               float
    ic_z_mm:               float


# ── Wheel travel sweep ────────────────────────────────────

def sweep_wheel_travel(hardpoints,
                       travel_range_mm=(-30.0, 30.0),
                       n_steps=61,
                       wheelbase=1530.0,
                       cog_height=280.0,
                       track_half_width=600.0,
                       is_front=True) -> list:
    """
    Sweep wheel travel; positive dz = bump.
    Contact patch and wheel centre are carried by the upright
    rigid-body transform (v3) so scrub and the wheel-centre
    trajectory are physically meaningful.
    """
    travels = np.linspace(*travel_range_mm, n_steps)

    uwb_in_mid = 0.5 * (hardpoints['UWB_inboard_front'] +
                        hardpoints['UWB_inboard_rear'])
    lwb_in_mid = 0.5 * (hardpoints['LWB_inboard_front'] +
                        hardpoints['LWB_inboard_rear'])
    uwb_len = np.linalg.norm(hardpoints['UWB_outboard'] - uwb_in_mid)
    lwb_len = np.linalg.norm(hardpoints['LWB_outboard'] - lwb_in_mid)
    tierod_len = np.linalg.norm(hardpoints['tierod_outer'] -
                                hardpoints['tierod_inner'])

    uwb_ref = hardpoints['UWB_outboard']
    lwb_ref = hardpoints['LWB_outboard']
    cp_ref  = hardpoints['contact_patch']
    wc_ref  = hardpoints['wheel_centre']

    results = []

    for dz in travels:
        uwb_new = solve_outboard_point(
            uwb_in_mid, hardpoints['UWB_inboard_front'],
            hardpoints['UWB_inboard_rear'], uwb_ref, uwb_len, dz)
        lwb_new = solve_outboard_point(
            lwb_in_mid, hardpoints['LWB_inboard_front'],
            hardpoints['LWB_inboard_rear'], lwb_ref, lwb_len, dz)
        pushrod_up_new = (hardpoints['pushrod_upright'] +
                          np.array([0.0, 0.0, dz]))

        # v3: upright-fixed points follow the upright
        T = upright_transform(lwb_ref, uwb_ref, lwb_new, uwb_new)
        # v4: steer about the kingpin set by the tie rod constraint
        steer_delta = solve_steer_angle(
            T, lwb_new, uwb_new, hardpoints['tierod_inner'],
            hardpoints['tierod_outer'], tierod_len)
        cp_new = T(cp_ref)
        wc_new = T(wc_ref)

        camber = compute_camber(uwb_new, lwb_new)
        kpi    = compute_kpi(uwb_new, lwb_new)
        caster = compute_caster(uwb_new, lwb_new)
        trail  = compute_mechanical_trail(uwb_new, lwb_new, cp_new)
        srad   = compute_scrub_radius(uwb_new, lwb_new, cp_new)
        scrub  = float(cp_new[1] - cp_ref[1])

        ic_yz, ic_found = compute_instant_centre(
            uwb_in_mid, uwb_new, lwb_in_mid, lwb_new)
        rc = (compute_roll_centre(ic_yz, cp_new, track_half_width)
              if ic_found else 0.0)
        toe = compute_toe(T, lwb_new, uwb_new, steer_delta)
        mr = compute_motion_ratio(pushrod_up_new,
                                  hardpoints['pushrod_bellcrank'])
        anti = compute_anti_geometry(
            ic_yz if ic_found else np.zeros(2),
            cp_new, wheelbase, cog_height)

        results.append(KinematicPoint(
            wheel_travel_mm       = float(dz),
            camber_deg            = round(camber, 4),
            roll_centre_height_mm = round(rc, 3),
            toe_deg               = round(toe, 4),
            motion_ratio          = round(mr, 4),
            anti_percent          = anti['anti_percent'],
            kpi_deg               = round(kpi, 4),
            caster_deg            = round(caster, 4),
            mech_trail_mm         = round(trail, 3),
            scrub_radius_mm       = round(srad, 3),
            scrub_mm              = round(scrub, 3),
            wc_y_mm               = round(float(wc_new[1]), 3),
            wc_z_mm               = round(float(wc_new[2]), 3),
            ic_y_mm  = round(float(ic_yz[0]) if ic_found else 0.0, 2),
            ic_z_mm  = round(float(ic_yz[1]) if ic_found else 0.0, 2),
        ))

    return results


# ── Slope-based sensitivity analysis ──────────────────────

SLOPE_OUTPUTS = {
    'camber_gain_deg_per_mm': 'camber_deg',
    'rc_gradient_mm_per_mm':  'roll_centre_height_mm',
    'bump_steer_deg_per_mm':  'toe_deg',
}


def _slope_at_drh(results: list, attr: str) -> float:
    """Gradient of an output w.r.t. wheel travel at DRH."""
    travels = np.array([r.wheel_travel_mm for r in results])
    vals    = np.array([getattr(r, attr) for r in results])
    grad    = np.gradient(vals, travels)
    zi      = int(np.argmin(np.abs(travels)))
    return float(grad[zi])


def sensitivity_analysis_slopes(hardpoints,
                                delta_mm=1.0,
                                travel_range_mm=(-30.0, 30.0),
                                n_steps=61,
                                is_front=True) -> Dict:
    """
    Slope-based sensitivity (v3): change in curve SLOPES at
    DRH per mm of hardpoint coordinate perturbation, by
    central difference. This is the design-relevant quantity:
    static values can be shimmed; gains are baked into the
    hardpoints.

    Returns {hardpoint: {axis: {slope_output: d(slope)/d(mm)}}}
    """
    sens = {}
    for hp_name, hp_val in hardpoints.items():
        sens[hp_name] = {}
        for i, ax in enumerate(['X', 'Y', 'Z']):
            hp_pos = dict(hardpoints)
            v = hp_val.copy(); v[i] += delta_mm
            hp_pos[hp_name] = v
            r_pos = sweep_wheel_travel(hp_pos, travel_range_mm,
                                       n_steps, is_front=is_front)

            hp_neg = dict(hardpoints)
            v = hp_val.copy(); v[i] -= delta_mm
            hp_neg[hp_name] = v
            r_neg = sweep_wheel_travel(hp_neg, travel_range_mm,
                                       n_steps, is_front=is_front)

            sens[hp_name][ax] = {}
            for out_name, attr in SLOPE_OUTPUTS.items():
                s_pos = _slope_at_drh(r_pos, attr)
                s_neg = _slope_at_drh(r_neg, attr)
                sens[hp_name][ax][out_name] = \
                    (s_pos - s_neg) / (2 * delta_mm)
    return sens


def plot_sensitivity_heatmap(sens: Dict,
                             title='Slope Sensitivity per 1mm '
                                   'Hardpoint Perturbation',
                             save_path=None) -> None:
    """
    RD-style heatmap grid: one panel per slope output,
    rows = hardpoints, columns = X/Y/Z, diverging colormap
    centred on zero.
    """
    hp_names = list(sens.keys())
    axes_lbl = ['X', 'Y', 'Z']
    outputs  = list(SLOPE_OUTPUTS.keys())

    fig, axs = plt.subplots(1, len(outputs),
                            figsize=(5.2 * len(outputs), 6.5))
    fig.suptitle(title, fontsize=13, fontweight='bold')

    for k, out_name in enumerate(outputs):
        M = np.array([[sens[hp][ax][out_name] for ax in axes_lbl]
                      for hp in hp_names])
        vmax = np.max(np.abs(M)) or 1e-12
        ax = axs[k]
        im = ax.imshow(M, cmap='RdYlGn', vmin=-vmax, vmax=vmax,
                       aspect='auto')
        ax.set_xticks(range(3)); ax.set_xticklabels(axes_lbl)
        ax.set_yticks(range(len(hp_names)))
        ax.set_yticklabels(hp_names if k == 0
                           else [''] * len(hp_names), fontsize=8)
        ax.set_title(out_name, fontsize=10)
        for r in range(M.shape[0]):
            for c in range(M.shape[1]):
                ax.text(c, r, f'{M[r, c]:+.4f}',
                        ha='center', va='center', fontsize=6.5)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    plt.tight_layout()
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')
    plt.show()


# ── Plotting ──────────────────────────────────────────────

def plot_kinematic_curves(results, title='Suspension Kinematics',
                          save_path=None) -> None:
    """3x4 grid: nine curves vs travel + wheel-centre trajectory."""
    t = [r.wheel_travel_mm for r in results]
    zi = min(range(len(t)), key=lambda i: abs(t[i]))

    series = [
        ('camber_deg',            'Camber (deg)',        'tab:red'),
        ('roll_centre_height_mm', 'RC Height (mm)',      'tab:blue'),
        ('toe_deg',               'Toe (deg)',           'tab:green'),
        ('motion_ratio',          'Motion Ratio (-)',    'tab:orange'),
        ('anti_percent',          'Anti-Geometry (%)',   'tab:purple'),
        ('scrub_mm',              'Scrub / Track Δ (mm)','tab:brown'),
        ('kpi_deg',               'KPI (deg)',           'tab:pink'),
        ('caster_deg',            'Caster (deg)',        'tab:olive'),
        ('mech_trail_mm',         'Mech Trail (mm)',     'tab:cyan'),
        ('scrub_radius_mm',       'Scrub Radius (mm)',   'tab:gray'),
    ]

    fig, axs = plt.subplots(3, 4, figsize=(19, 11))
    fig.suptitle(title, fontsize=14, fontweight='bold')
    flat = axs.ravel()

    for k, (attr, lbl, col) in enumerate(series):
        ax = flat[k]
        d = [getattr(r, attr) for r in results]
        ax.plot(t, d, color=col, linewidth=2)
        ax.axhline(0, color='#555', lw=0.5, ls='--')
        ax.axvline(0, color='#555', lw=0.5, ls='--')
        ax.set_xlabel('Wheel Travel (mm)')
        ax.set_ylabel(lbl)
        ax.set_title(lbl)
        ax.grid(True, alpha=0.3)
        ax.annotate(f'DRH: {d[zi]:.3f}', xy=(0, d[zi]),
                    xytext=(8, 8), textcoords='offset points',
                    fontsize=8, color=col)

    # Panel 11: wheel-centre front-view trajectory (droop→bump)
    ax = flat[10]
    wy = [r.wc_y_mm for r in results]
    wz = [r.wc_z_mm for r in results]
    ax.plot(wy, wz, color='tab:blue', lw=2)
    ax.plot(wy[0],  wz[0],  'v', color='tab:red',
            label=f'droop {t[0]:.0f}mm')
    ax.plot(wy[-1], wz[-1], '^', color='tab:green',
            label=f'bump {t[-1]:.0f}mm')
    ax.set_xlabel('Wheel Centre Y (mm)')
    ax.set_ylabel('Wheel Centre Z (mm)')
    ax.set_title('Wheel Centre YZ Trajectory')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    flat[11].set_visible(False)
    plt.tight_layout()
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')
    plt.show()


# ── Export ────────────────────────────────────────────────

CSV_FIELDS = ['wheel_travel_mm', 'camber_deg',
              'roll_centre_height_mm', 'toe_deg', 'motion_ratio',
              'anti_percent', 'kpi_deg', 'caster_deg',
              'mech_trail_mm', 'scrub_radius_mm', 'scrub_mm',
              'wc_y_mm', 'wc_z_mm', 'ic_y_mm', 'ic_z_mm']


def export_results(results, filepath) -> None:
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(CSV_FIELDS)
        for r in results:
            w.writerow([getattr(r, k) for k in CSV_FIELDS])
    print(f'  Exported: {filepath}')


# ── Main ──────────────────────────────────────────────────

if __name__ == '__main__':
    print('FS Suspension Kinematics Solver v3')
    print('=' * 52)

    hp_front = default_fs_front_hardpoints()
    hp_rear  = default_fs_rear_hardpoints()

    print('\nSweeping front corner...')
    rf = sweep_wheel_travel(hp_front, is_front=True)
    print('Sweeping rear corner...')
    rr = sweep_wheel_travel(hp_rear, is_front=False)

    zf = next(r for r in rf if abs(r.wheel_travel_mm) < 0.1)
    zr = next(r for r in rr if abs(r.wheel_travel_mm) < 0.1)

    print('\n── Design Ride Height Summary ──')
    print(f'{"Parameter":<26} {"Front":>10} {"Rear":>10}')
    print('-' * 50)
    for name, attr in [('Camber (deg)', 'camber_deg'),
                       ('KPI (deg)', 'kpi_deg'),
                       ('Caster (deg)', 'caster_deg'),
                       ('Mech trail (mm)', 'mech_trail_mm'),
                       ('Scrub radius (mm)', 'scrub_radius_mm'),
                       ('RC Height (mm)', 'roll_centre_height_mm'),
                       ('Toe (deg)', 'toe_deg'),
                       ('Motion Ratio', 'motion_ratio'),
                       ('Anti-Geom (%)', 'anti_percent')]:
        print(f'{name:<26} {getattr(zf, attr):>10.3f} '
              f'{getattr(zr, attr):>10.3f}')

    # Slopes at DRH
    for lbl, attr in [('Camber gain (deg/mm)', 'camber_deg'),
                      ('RC gradient (mm/mm)',
                       'roll_centre_height_mm'),
                      ('Bump steer (deg/mm)', 'toe_deg'),
                      ('Scrub rate (mm/mm)', 'scrub_mm')]:
        print(f'{lbl:<26} {_slope_at_drh(rf, attr):>10.4f} '
              f'{_slope_at_drh(rr, attr):>10.4f}')

    print('\nExporting...')
    export_results(rf, 'outputs/front_kinematics.csv')
    export_results(rr, 'outputs/rear_kinematics.csv')

    print('\nPlotting...')
    plot_kinematic_curves(rf,
        'Front Suspension Kinematics — FS Car (v3)',
        'outputs/front_kinematics_v3.png')
    plot_kinematic_curves(rr,
        'Rear Suspension Kinematics — FS Car (v3)',
        'outputs/rear_kinematics_v3.png')

    print('\nSlope-based sensitivity (front, ±1mm)...')
    sens = sensitivity_analysis_slopes(hp_front, delta_mm=1.0)
    plot_sensitivity_heatmap(sens,
        'Front Corner — Slope Sensitivity per 1mm Perturbation',
        'outputs/sensitivity_heatmap_v3.png')

    print('\nDone.')
