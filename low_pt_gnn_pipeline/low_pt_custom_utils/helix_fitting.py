"""
Helix fitting utilities for segment matching.

Fits helical trajectories to sets of detector hits. In a uniform magnetic field B
along z, charged particles follow helices: circular in x-y, linear z vs arc-length.

Circle fitting uses the Kasa algebraic method (least squares on ALL hits).
Pitch estimation fits z vs cumulative arc length along the fitted circle.

Reference: Kasa, I. "A circle fitting procedure and its error analysis."
           IEEE Trans. Instrum. Meas., 25(1):8-14, 1976.
"""

from dataclasses import dataclass
from typing import Tuple, Optional

import numpy as np


@dataclass
class HelixParams:
    """Fitted helix parameters for a segment."""
    xc: float               # Circle center x (mm)
    yc: float               # Circle center y (mm)
    R: float                # Circle radius (mm)
    pitch: Optional[float]  # dz/ds (mm/mm), None if not computed
    z0: Optional[float]     # z-intercept at s=0, None if not computed
    pT: float               # Transverse momentum (GeV/c), 0.3 * B * R / 1000
    phi_center: float       # atan2(yc, xc) — azimuthal angle of circle center
    fit_quality: str        # "good" (3+ hits), "poor" (2 hits), "none" (1 hit)
    residual_rms: float     # RMS of circle fit residuals (mm)
    nhits: int              # Number of hits used in fit


def kasa_circle_fit(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float, np.ndarray]:
    """
    Algebraic circle fit using the Kasa method.

    Minimizes sum of (x_i^2 + y_i^2 - 2*xc*x_i - 2*yc*y_i - c)^2
    by solving the linear system:
        [2x, 2y, 1] * [xc, yc, c]^T = x^2 + y^2

    Uses ALL points via least squares — more points = more robust fit.

    Args:
        x, y: Arrays of hit coordinates (mm). Must have len >= 3.

    Returns:
        xc, yc: Circle center coordinates (mm)
        R: Circle radius (mm)
        residuals: Per-point residuals |sqrt((x_i-xc)^2 + (y_i-yc)^2) - R| (mm)

    Raises:
        ValueError: If fewer than 3 points provided or fit is degenerate.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = len(x)

    if n < 3:
        raise ValueError(f"Need at least 3 points for circle fit, got {n}")

    # Build linear system: A * [xc, yc, c]^T = b
    A = np.column_stack([2.0 * x, 2.0 * y, np.ones(n)])
    b = x**2 + y**2

    # Solve the normal equations (A^T A) x = A^T b — a 3x3 solve, ~5x cheaper than the
    # SVD-based lstsq and identical for this well-conditioned over-determined fit.
    try:
        result = np.linalg.solve(A.T @ A, A.T @ b)
    except np.linalg.LinAlgError:
        raise ValueError("Degenerate circle fit (points may be collinear)")

    xc, yc, c = result
    R_squared = c + xc**2 + yc**2

    if R_squared <= 0:
        raise ValueError(f"Invalid circle fit: R^2 = {R_squared:.4f} <= 0")

    R = np.sqrt(R_squared)

    # Compute residuals: |distance_to_center - R|
    distances = np.sqrt((x - xc)**2 + (y - yc)**2)
    residuals = np.abs(distances - R)

    return xc, yc, R, residuals


def fit_circle_with_outlier_rejection(
    x: np.ndarray,
    y: np.ndarray,
    max_iterations: int = 3,
    sigma_cut: float = 5.0,
    min_residual_mm: float = 5.0,
) -> Tuple[float, float, float, np.ndarray]:
    """
    Iterative Kasa circle fit with outlier rejection.

    After each fit, points with residual > sigma_cut * median_residual
    AND residual > min_residual_mm are rejected and the fit is repeated.

    Args:
        x, y: Arrays of hit coordinates (mm).
        max_iterations: Maximum number of fit-reject cycles.
        sigma_cut: Multiplier for median residual threshold.
        min_residual_mm: Absolute minimum residual to be considered an outlier.

    Returns:
        xc, yc: Circle center coordinates (mm)
        R: Circle radius (mm)
        inlier_mask: Boolean array, True for inliers.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    mask = np.ones(len(x), dtype=bool)

    for _ in range(max_iterations):
        if mask.sum() < 3:
            break

        xc, yc, R, residuals_masked = kasa_circle_fit(x[mask], y[mask])

        # Compute residuals for ALL points (including previously masked)
        distances = np.sqrt((x - xc)**2 + (y - yc)**2)
        residuals_all = np.abs(distances - R)

        # Identify outliers among currently included points
        median_res = np.median(residuals_masked)
        threshold = max(sigma_cut * median_res, min_residual_mm)

        new_mask = residuals_all <= threshold
        # Don't unmask previously rejected points
        new_mask = new_mask & mask

        if np.array_equal(new_mask, mask):
            break  # No change, converged

        # Ensure we still have enough points
        if new_mask.sum() < 3:
            break

        mask = new_mask

    # Final fit on inliers
    xc, yc, R, _ = kasa_circle_fit(x[mask], y[mask])

    return xc, yc, R, mask


def compute_arc_lengths(
    x: np.ndarray,
    y: np.ndarray,
    xc: float,
    yc: float,
    R: float,
) -> np.ndarray:
    """
    Compute cumulative arc lengths along a fitted circle.

    Points must be in trajectory order (e.g., sorted by radial position outward).
    Arc length between consecutive points is R * |delta_theta|.

    Args:
        x, y: Hit coordinates in trajectory order (mm).
        xc, yc: Circle center (mm).
        R: Circle radius (mm).

    Returns:
        arc_lengths: Cumulative arc length from first point (mm).
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    # Angles relative to circle center
    theta = np.arctan2(y - yc, x - xc)

    # Compute angular differences and unwrap
    dtheta = np.diff(theta)
    # Wrap to [-pi, pi]
    dtheta = (dtheta + np.pi) % (2 * np.pi) - np.pi

    # Cumulative arc length
    arc_lengths = np.zeros(len(x))
    arc_lengths[1:] = np.cumsum(np.abs(dtheta) * R)

    return arc_lengths


def fit_pitch(arc_lengths: np.ndarray, z: np.ndarray) -> Tuple[float, float, np.ndarray]:
    """
    Linear fit of z vs arc-length: z = z0 + pitch * s.

    Args:
        arc_lengths: Cumulative arc lengths (mm).
        z: z-coordinates of hits (mm).

    Returns:
        pitch: dz/ds (dimensionless, mm/mm).
        z0: z-intercept at s=0 (mm).
        residuals: Per-point residuals |z_i - (z0 + pitch * s_i)| (mm).
    """
    arc_lengths = np.asarray(arc_lengths, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)

    if len(arc_lengths) < 2:
        return 0.0, z[0] if len(z) > 0 else 0.0, np.array([0.0])

    # Closed-form least-squares line z = pitch*s + z0 (identical to lstsq, no SVD).
    s_mean = arc_lengths.mean(); z_mean = z.mean()
    ds = arc_lengths - s_mean
    denom = float(ds @ ds)
    if denom < 1e-12:
        pitch, z0 = 0.0, z_mean
    else:
        pitch = float(ds @ (z - z_mean)) / denom
        z0 = z_mean - pitch * s_mean

    z_predicted = z0 + pitch * arc_lengths
    residuals = np.abs(z - z_predicted)

    return pitch, z0, residuals


def fit_helix_to_segment(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    B_field: float = 2.0,
    outlier_rejection: bool = False,
) -> HelixParams:
    """
    Full helix fit to a segment's hits.

    For 3+ hits: Kasa circle fit in x-y, then pitch fit in z vs arc-length.
    For 2 hits: No circle fit possible, returns HelixParams with fit_quality="poor".
    For 1 hit: Returns HelixParams with fit_quality="none".

    Args:
        x, y, z: Hit coordinates in trajectory order (mm).
        B_field: Magnetic field strength (Tesla). Default 2.0 T.
        outlier_rejection: Whether to use iterative outlier rejection.

    Returns:
        HelixParams with all fitted parameters.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    n = len(x)

    if n == 0:
        return HelixParams(
            xc=0.0, yc=0.0, R=0.0, pitch=None, z0=None, pT=0.0,
            phi_center=0.0, fit_quality="none", residual_rms=0.0, nhits=0,
        )

    if n == 1:
        return HelixParams(
            xc=0.0, yc=0.0, R=0.0, pitch=None, z0=z[0], pT=0.0,
            phi_center=0.0, fit_quality="none", residual_rms=0.0, nhits=1,
        )

    if n == 2:
        # Can't fit a circle with 2 points. Store endpoint info only.
        return HelixParams(
            xc=0.0, yc=0.0, R=0.0, pitch=None, z0=None, pT=0.0,
            phi_center=0.0, fit_quality="poor", residual_rms=0.0, nhits=2,
        )

    # 3+ hits: full circle fit
    try:
        if outlier_rejection:
            xc, yc, R, inlier_mask = fit_circle_with_outlier_rejection(x, y)
            # Compute residuals on inliers
            distances = np.sqrt((x[inlier_mask] - xc)**2 + (y[inlier_mask] - yc)**2)
            circle_residuals = np.abs(distances - R)
        else:
            xc, yc, R, circle_residuals = kasa_circle_fit(x, y)
            inlier_mask = np.ones(n, dtype=bool)
    except ValueError:
        # Degenerate fit (e.g., collinear points)
        return HelixParams(
            xc=0.0, yc=0.0, R=0.0, pitch=None, z0=None, pT=0.0,
            phi_center=0.0, fit_quality="poor", residual_rms=0.0, nhits=n,
        )

    residual_rms = np.sqrt(np.mean(circle_residuals**2))

    # Very large radius likely means near-straight track segment
    # (R > 10m is effectively infinite for our detector at ~1m scale)
    if R > 10000:
        R = 10000  # Cap to avoid numerical issues in pT calculation

    pT = 0.3 * B_field * R / 1000.0  # R in mm → R/1000 in m → pT in GeV/c
    phi_center = np.arctan2(yc, xc)

    # Pitch fit: z vs arc-length (using inlier hits only)
    x_fit = x[inlier_mask]
    y_fit = y[inlier_mask]
    z_fit = z[inlier_mask]

    arc_lengths = compute_arc_lengths(x_fit, y_fit, xc, yc, R)

    # Only fit pitch if we have sufficient arc length span
    arc_span = arc_lengths[-1] - arc_lengths[0] if len(arc_lengths) > 1 else 0.0
    if arc_span > 1e-6:
        pitch, z0, _ = fit_pitch(arc_lengths, z_fit)
    else:
        pitch = None
        z0 = np.mean(z_fit)

    return HelixParams(
        xc=xc, yc=yc, R=R, pitch=pitch, z0=z0, pT=pT,
        phi_center=phi_center, fit_quality="good",
        residual_rms=residual_rms, nhits=n,
    )
