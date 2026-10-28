"""Experiment-agnostic monitor gamma fitting, validation, and JSON I/O."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import socket
import sys
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.optimize import minimize_scalar


SCHEMA_VERSION = "portable-monitor-calibration-v1"
DEFAULT_LEVEL_COUNT = 17
PLAUSIBLE_GAMMA_RANGE = (1.0, 3.5)
MINIMUM_R_SQUARED = 0.99
MAXIMUM_VERIFICATION_ERROR = 0.05


def calibration_levels(level_count: int = DEFAULT_LEVEL_COUNT) -> np.ndarray:
    """Return evenly spaced display-code values from black through white."""

    level_count = int(level_count)
    if level_count < 3:
        raise ValueError("a gamma fit requires at least three luminance levels")
    return np.linspace(0.0, 1.0, level_count, dtype=np.float64)


def luminance_model(
    levels: Sequence[float],
    minimum_luminance: float,
    maximum_luminance: float,
    gamma: float,
) -> np.ndarray:
    """Evaluate Lmin + (Lmax - Lmin) * code**gamma."""

    code = np.clip(np.asarray(levels, dtype=np.float64), 0.0, 1.0)
    return float(minimum_luminance) + (
        float(maximum_luminance) - float(minimum_luminance)
    ) * code ** float(gamma)


def fit_gamma(
    levels: Sequence[float], luminances_cd_m2: Sequence[float]
) -> dict[str, float | str]:
    """Fit a power-law gamma with directly measured black and white endpoints."""

    code = np.asarray(levels, dtype=np.float64)
    luminance = np.asarray(luminances_cd_m2, dtype=np.float64)
    if code.shape != luminance.shape:
        raise ValueError("levels and luminance readings must have the same shape")
    if code.ndim != 1 or code.size < 3:
        raise ValueError("a gamma fit requires at least three one-dimensional readings")
    if not np.isfinite(code).all() or not np.isfinite(luminance).all():
        raise ValueError("levels and luminance readings must be finite")
    if np.any(luminance < 0):
        raise ValueError("luminance readings cannot be negative")
    order = np.argsort(code)
    code = code[order]
    luminance = luminance[order]
    if code[0] != 0.0 or code[-1] != 1.0:
        raise ValueError("levels must include the black (0) and white (1) endpoints")
    minimum, maximum = float(luminance[0]), float(luminance[-1])
    if maximum <= minimum:
        raise ValueError("white luminance must exceed black luminance")

    def squared_error(candidate_gamma: float) -> float:
        predicted = luminance_model(code, minimum, maximum, candidate_gamma)
        return float(np.sum(np.square(luminance - predicted)))

    result = minimize_scalar(squared_error, bounds=(0.1, 10.0), method="bounded")
    if not result.success:
        raise RuntimeError(f"gamma fitting failed: {result.message}")
    gamma = float(result.x)
    predicted = luminance_model(code, minimum, maximum, gamma)
    residuals = luminance - predicted
    residual_sum = float(np.sum(np.square(residuals)))
    total_sum = float(np.sum(np.square(luminance - luminance.mean())))
    return {
        "model": "power",
        "gamma": gamma,
        "min_luminance_cd_m2": minimum,
        "max_luminance_cd_m2": maximum,
        "r_squared": 1.0 - residual_sum / total_sum if total_sum > 0 else float("nan"),
        "rmse_cd_m2": float(np.sqrt(residual_sum / code.size)),
        "max_abs_residual_cd_m2": float(np.max(np.abs(residuals))),
    }


def linearity_error(
    levels: Sequence[float], luminances_cd_m2: Sequence[float]
) -> dict[str, float]:
    """Measure corrected-display departures from a linear luminance ramp."""

    code = np.asarray(levels, dtype=np.float64)
    luminance = np.asarray(luminances_cd_m2, dtype=np.float64)
    if code.shape != luminance.shape or code.ndim != 1:
        raise ValueError("verification levels and readings must be matching vectors")
    if not np.isfinite(code).all() or not np.isfinite(luminance).all():
        raise ValueError("verification levels and readings must be finite")
    minimum, maximum = float(luminance.min()), float(luminance.max())
    if maximum <= minimum:
        raise ValueError("verification readings must span a luminance range")
    ideal = minimum + (maximum - minimum) * code
    residuals = luminance - ideal
    luminance_range = maximum - minimum
    return {
        "max_abs_error_fraction": float(np.max(np.abs(residuals)) / luminance_range),
        "rms_error_fraction": float(
            np.sqrt(np.mean(np.square(residuals))) / luminance_range
        ),
    }


def software_provenance() -> dict[str, Any]:
    """Record the host and numerical runtime used for calibration."""

    record: dict[str, Any] = {
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
    }
    try:
        import psychopy

        record["psychopy_version"] = psychopy.__version__
    except Exception:
        record["psychopy_version"] = None
    return record


def build_calibration_record(
    display_name: str,
    levels: Sequence[float],
    luminances_cd_m2: Sequence[float],
    *,
    display_profile: Mapping[str, Any],
    screen_index: int,
    photometer_model: str = "",
    operator: str = "",
    warmup_minutes: float = 0.0,
    corrected_luminances_cd_m2: Sequence[float] | None = None,
    measured_refresh_hz: float | None = None,
) -> dict[str, Any]:
    """Build a portable record accepted by strict PsychoPy experiments."""

    fit = fit_gamma(levels, luminances_cd_m2)
    created = datetime.now(timezone.utc).isoformat()
    record: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "display_name": str(display_name),
        "gamma": float(fit["gamma"]),
        "created": created,
        "created_at": created,
        "correction_mode": "gamma",
        "luminance_cd_m2": {
            "min": float(fit["min_luminance_cd_m2"]),
            "max": float(fit["max_luminance_cd_m2"]),
        },
        "fit": {
            key: fit[key]
            for key in (
                "model",
                "r_squared",
                "rmse_cd_m2",
                "max_abs_residual_cd_m2",
            )
        },
        "measurements": {
            "levels": np.asarray(levels, dtype=np.float64).tolist(),
            "luminance_cd_m2": np.asarray(
                luminances_cd_m2, dtype=np.float64
            ).tolist(),
            "luminance_post_cd_m2": (
                None
                if corrected_luminances_cd_m2 is None
                else np.asarray(
                    corrected_luminances_cd_m2, dtype=np.float64
                ).tolist()
            ),
        },
        "display_profile": dict(display_profile),
        "provenance": {
            **software_provenance(),
            "screen_index": int(screen_index),
            "photometer_model": str(photometer_model),
            "operator": str(operator),
            "warmup_minutes": float(warmup_minutes),
            "measured_refresh_hz": (
                None if measured_refresh_hz is None else float(measured_refresh_hz)
            ),
        },
    }
    if corrected_luminances_cd_m2 is not None:
        record["verification"] = linearity_error(
            levels, corrected_luminances_cd_m2
        )
    record["quality_control"] = quality_control(record)
    return record


def quality_control(
    record: Mapping[str, Any], *, require_verification: bool = True
) -> dict[str, Any]:
    """Return explicit pass/fail gates without discarding diagnostic records."""

    reasons: list[str] = []
    try:
        gamma = float(record["gamma"])
    except (KeyError, TypeError, ValueError):
        gamma = float("nan")
        reasons.append("missing or invalid gamma")
    if not PLAUSIBLE_GAMMA_RANGE[0] <= gamma <= PLAUSIBLE_GAMMA_RANGE[1]:
        reasons.append(
            f"gamma must be within {PLAUSIBLE_GAMMA_RANGE[0]}-{PLAUSIBLE_GAMMA_RANGE[1]}"
        )
    try:
        r_squared = float(record["fit"]["r_squared"])
    except (KeyError, TypeError, ValueError):
        r_squared = float("nan")
        reasons.append("missing or invalid fit R-squared")
    if not r_squared >= MINIMUM_R_SQUARED:
        reasons.append(f"fit R-squared must be at least {MINIMUM_R_SQUARED}")
    verification = record.get("verification")
    if not isinstance(verification, Mapping):
        if require_verification:
            reasons.append("a corrected verification pass is required")
    else:
        for key in ("max_abs_error_fraction", "rms_error_fraction"):
            try:
                error = float(verification[key])
            except (KeyError, TypeError, ValueError):
                reasons.append(f"verification is missing {key}")
                continue
            if error > MAXIMUM_VERIFICATION_ERROR:
                reasons.append(
                    f"{key} exceeds {MAXIMUM_VERIFICATION_ERROR:.0%} of luminance range"
                )
    return {
        "passed": not reasons,
        "require_verification": bool(require_verification),
        "minimum_r_squared": MINIMUM_R_SQUARED,
        "plausible_gamma_range": list(PLAUSIBLE_GAMMA_RANGE),
        "maximum_verification_error_fraction": MAXIMUM_VERIFICATION_ERROR,
        "failure_reasons": reasons,
    }


def save_calibration(path: str | Path, record: Mapping[str, Any]) -> Path:
    """Save calibration JSON atomically."""

    destination = Path(path).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(f".{destination.name}.partial-{os.getpid()}")
    partial.write_text(json.dumps(record, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(partial, destination)
    return destination


def load_calibration(
    path: str | Path,
    *,
    expected_display_name: str | None = None,
    require_verification: bool = True,
) -> dict[str, Any]:
    """Load a calibration and reject records that fail the frozen quality gates."""

    source = Path(path).resolve()
    record = json.loads(source.read_text(encoding="utf-8"))
    if expected_display_name is not None and record.get("display_name") != expected_display_name:
        raise ValueError(
            f"calibration is for {record.get('display_name')!r}, not {expected_display_name!r}"
        )
    qc = quality_control(record, require_verification=require_verification)
    if not qc["passed"]:
        raise ValueError("calibration failed quality control: " + "; ".join(qc["failure_reasons"]))
    record["quality_control"] = qc
    record["source_path"] = str(source)
    return record
