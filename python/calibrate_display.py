#!/usr/bin/env python3
"""Calibrate a presentation display with PsychoPy and a hand-held photometer."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

import numpy as np

from calibration_core import (
    DEFAULT_LEVEL_COUNT,
    build_calibration_record,
    calibration_levels,
    fit_gamma,
    luminance_model,
    save_calibration,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_PROFILE = "ASUS B"
ENTRY_KEYS = {str(number): str(number) for number in range(10)}
ENTRY_KEYS.update({"period": ".", "num_decimal": ".", "minus": "-"})
ENTRY_KEYS.update({f"num_{number}": str(number) for number in range(10)})


def load_profiles(path: str | Path) -> dict[str, dict[str, Any]]:
    profiles = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(profiles, dict) or DEFAULT_PROFILE not in profiles:
        raise ValueError(f"display profiles must contain {DEFAULT_PROFILE!r}")
    required = {
        "width_px",
        "height_px",
        "width_cm",
        "viewing_distance_cm",
        "expected_refresh_hz",
        "fullscreen",
        "win_type",
    }
    for name, profile in profiles.items():
        missing = required - set(profile)
        if missing:
            raise ValueError(f"display profile {name!r} is missing {sorted(missing)}")
    return profiles


def activate_window_backend(win_type: str) -> str | None:
    """Activate PsychoPy's external window-backend plugin when required."""

    if str(win_type).lower() != "glfw":
        return None
    from psychopy import plugins

    plugin_name = "psychopy-glfw"
    plugins.scanPlugins()
    if not plugins.loadPlugin(plugin_name):
        raise RuntimeError(
            "The ASUS display profile requires PsychoPy's GLFW plugin, but "
            "psychopy-glfw could not be activated. In the active environment run "
            "`python -m pip install psychopy-glfw`, then retry."
        )
    return plugin_name


def apply_hardware_gamma(window: Any, gamma: float) -> None:
    """Apply a hardware gamma ramp or report an actionable Linux display error."""

    try:
        window.gamma = float(gamma)
    except OSError as error:
        raise RuntimeError(
            "PsychoPy cannot access this display session's hardware gamma ramp. "
            "On Linux, check `echo $XDG_SESSION_TYPE`; Wayland and some compositor "
            "configurations do not expose XF86VidMode gamma ramps. Use an X11/Xorg "
            "session for hardware-linearized acquisition, or use a separately "
            "validated software-correction workflow."
        ) from error


def simulated_readings(
    levels: Sequence[float], gamma: float, *, seed: int = 0
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    clean = luminance_model(levels, 0.3, 250.0, gamma)
    return clean + rng.normal(0.0, 0.0005 * 250.0, size=len(levels))


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", default=str(ROOT / "display_profiles.json"))
    parser.add_argument("--profile", default=DEFAULT_PROFILE)
    parser.add_argument(
        "--window-backend",
        choices=("pyglet", "glfw"),
        default="pyglet",
        help="PsychoPy window backend; pyglet avoids known Linux EGL/GLFW failures",
    )
    parser.add_argument("--screen-index", type=int, default=0)
    parser.add_argument("--levels", type=int, default=DEFAULT_LEVEL_COUNT)
    parser.add_argument("--photometer-model", default="")
    parser.add_argument("--operator", default="")
    parser.add_argument("--output-dir", default=str(ROOT / "calibrations"))
    parser.add_argument("--windowed", action="store_true")
    parser.add_argument("--no-verification", action="store_true")
    parser.add_argument("--no-gui", action="store_true")
    parser.add_argument(
        "--simulate-photometer",
        type=float,
        metavar="GAMMA",
        help="Rehearse the complete workflow with deterministic simulated readings",
    )
    return parser.parse_args(argv)


def _selected(value: Any) -> Any:
    return value[0] if isinstance(value, (list, tuple)) else value


def collect_settings(args: argparse.Namespace, profiles: Mapping[str, Any]) -> dict[str, Any] | None:
    from psychopy import gui

    names = list(profiles)
    if args.profile in names:
        names = [args.profile, *(name for name in names if name != args.profile)]
    values = {
        "Display Profile": names,
        "Window Backend": [
            args.window_backend,
            *(name for name in ("pyglet", "glfw") if name != args.window_backend),
        ],
        "Screen Index": str(args.screen_index),
        "Levels": str(args.levels),
        "Photometer Model": args.photometer_model,
        "Operator": args.operator,
        "Verification Pass": ["Yes", "No"] if not args.no_verification else ["No", "Yes"],
        "Windowed Rehearsal": ["Yes", "No"] if args.windowed else ["No", "Yes"],
        "Simulated Gamma (blank = photometer)": (
            "" if args.simulate_photometer is None else str(args.simulate_photometer)
        ),
        "Output Directory": args.output_dir,
    }
    dialog = gui.DlgFromDict(
        dictionary=values,
        title="Standalone Monitor Calibration",
        order=list(values),
        sortKeys=False,
    )
    if not dialog.OK:
        return None
    simulated = str(values["Simulated Gamma (blank = photometer)"]).strip()
    return {
        "profile": str(_selected(values["Display Profile"])),
        "window_backend": str(_selected(values["Window Backend"])),
        "screen_index": int(values["Screen Index"]),
        "levels": int(values["Levels"]),
        "photometer_model": str(values["Photometer Model"]),
        "operator": str(values["Operator"]),
        "verification": _selected(values["Verification Pass"]) == "Yes",
        "windowed": _selected(values["Windowed Rehearsal"]) == "Yes",
        "simulate_photometer": None if not simulated else float(simulated),
        "output_dir": str(Path(str(values["Output Directory"])).expanduser()),
    }


def command_line_settings(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "profile": args.profile,
        "window_backend": args.window_backend,
        "screen_index": args.screen_index,
        "levels": args.levels,
        "photometer_model": args.photometer_model,
        "operator": args.operator,
        "verification": not args.no_verification,
        "windowed": args.windowed,
        "simulate_photometer": args.simulate_photometer,
        "output_dir": str(Path(args.output_dir).expanduser()),
    }


def validate_settings(settings: Mapping[str, Any], profiles: Mapping[str, Any]) -> None:
    if settings["profile"] not in profiles:
        raise ValueError(f"unknown display profile: {settings['profile']!r}")
    if settings["window_backend"] not in {"pyglet", "glfw"}:
        raise ValueError(f"unknown window backend: {settings['window_backend']!r}")
    if int(settings["screen_index"]) < 0:
        raise ValueError("screen index must be zero or greater")
    if int(settings["levels"]) < 3:
        raise ValueError("at least three luminance levels are required")
    if settings["simulate_photometer"] is not None and float(settings["simulate_photometer"]) <= 0:
        raise ValueError("simulated gamma must be positive")
    if not str(settings["output_dir"]).strip():
        raise ValueError("output directory cannot be empty")


def save_figure(
    path: Path,
    levels: np.ndarray,
    luminances: np.ndarray,
    fit: Mapping[str, Any],
    corrected: np.ndarray | None,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panel_count = 3 if corrected is not None else 2
    figure, axes = plt.subplots(1, panel_count, figsize=(5.0 * panel_count, 4.2))
    dense = np.linspace(0.0, 1.0, 300)
    predicted = luminance_model(
        dense,
        fit["min_luminance_cd_m2"],
        fit["max_luminance_cd_m2"],
        fit["gamma"],
    )
    axes[0].plot(levels, luminances, "o", label="measured")
    axes[0].plot(dense, predicted, label=f"fit gamma={fit['gamma']:.3f}")
    axes[0].set(xlabel="Display code", ylabel="Luminance (cd/m2)", title=f"R2={fit['r_squared']:.5f}")
    axes[0].legend()
    fitted_at_levels = luminance_model(
        levels,
        fit["min_luminance_cd_m2"],
        fit["max_luminance_cd_m2"],
        fit["gamma"],
    )
    axes[1].axhline(0, color="0.6", linewidth=0.8)
    axes[1].plot(levels, luminances - fitted_at_levels, "o-")
    axes[1].set(xlabel="Display code", ylabel="Residual (cd/m2)", title=f"RMSE={fit['rmse_cd_m2']:.3f}")
    if corrected is not None:
        ideal = corrected.min() + (corrected.max() - corrected.min()) * levels
        axes[2].plot(levels, corrected, "o-", label="corrected")
        axes[2].plot(levels, ideal, "--", color="0.6", label="ideal")
        axes[2].set(xlabel="Display code", ylabel="Luminance (cd/m2)", title="Verification")
        axes[2].legend()
    figure.tight_layout()
    figure.savefig(path)
    plt.close(figure)


def _safe_display_filename(display_name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._ -]+", "_", display_name).strip() or "display"


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    profiles = load_profiles(args.profiles)
    settings = command_line_settings(args) if args.no_gui else collect_settings(args, profiles)
    if settings is None:
        return 0
    validate_settings(settings, profiles)
    profile = dict(profiles[settings["profile"]])
    profile["win_type"] = settings["window_backend"]
    levels = calibration_levels(settings["levels"])

    activate_window_backend(str(profile["win_type"]))
    from psychopy import core, event, visual

    width, height = int(profile["width_px"]), int(profile["height_px"])
    full_screen = bool(profile["fullscreen"] and not settings["windowed"])
    window_size = (width, height) if full_screen else (min(width, 1200), min(height, 800))
    window = None
    original_ramp = None
    luminances = None
    corrected = None
    measured_refresh = None
    try:
        window = visual.Window(
            size=window_size,
            fullscr=full_screen,
            screen=int(settings["screen_index"]),
            color=0,
            colorSpace="rgb",
            units="pix",
            winType=str(profile["win_type"]),
            gamma=None,
            waitBlanking=True,
            checkTiming=True,
            autoLog=False,
        )
        try:
            original_ramp = np.asarray(window.backend.getGammaRamp(), dtype=float).copy()
        except Exception:
            original_ramp = None
        window.mouseVisible = False
        measured_refresh = window.getActualFrameRate(
            nIdentical=30, nMaxFrames=180, nWarmUpFrames=60, threshold=1
        )
        uniform_field = visual.Rect(
            window,
            width=window_size[0] + 4,
            height=window_size[1] + 4,
            lineWidth=0,
            units="pix",
            autoLog=False,
        )
        prompt = visual.TextStim(
            window,
            color="white",
            height=20,
            units="pix",
            pos=(0, -window_size[1] / 2 + 90),
            wrapWidth=0.85 * window_size[0],
            autoLog=False,
        )

        def show_uniform_field(level: float) -> None:
            uniform_field.fillColor = [2 * float(level) - 1] * 3
            uniform_field.draw()
            window.flip()

        def show_prompt(message: str) -> None:
            prompt.text = message
            prompt.draw()
            window.flip()

        def read_luminance(level: float, index: int, label: str) -> float | None:
            show_uniform_field(level)
            core.wait(0.25)
            event.clearEvents(eventType="keyboard")
            keys = event.waitKeys(keyList=["space", "escape"])
            if "escape" in keys:
                return None
            typed = ""
            while True:
                show_prompt(
                    f"{label}  {index + 1}/{len(levels)}  code={level:.4f}\n\n"
                    f"Photometer reading (cd/m2): {typed}_\n\n"
                    "RETURN accepts, BACKSPACE edits, ESCAPE aborts",
                )
                for key in event.waitKeys():
                    if key == "escape":
                        return None
                    if key == "backspace":
                        typed = typed[:-1]
                    elif key in ("return", "num_enter"):
                        try:
                            return float(typed)
                        except ValueError:
                            typed = ""
                    elif key in ENTRY_KEYS:
                        typed += ENTRY_KEYS[key]

        def measure(label: str, simulation_seed: int) -> np.ndarray | None:
            if settings["simulate_photometer"] is not None:
                return simulated_readings(
                    levels, float(settings["simulate_photometer"]), seed=simulation_seed
                )
            readings = []
            for index, level in enumerate(levels):
                reading = read_luminance(float(level), index, label)
                if reading is None:
                    return None
                readings.append(reading)
            return np.asarray(readings, dtype=np.float64)

        luminances = measure("Uncorrected", 0)
        if luminances is None:
            return 1
        fit = fit_gamma(levels, luminances)
        print(
            f"gamma={fit['gamma']:.4f}  R2={fit['r_squared']:.5f}  "
            f"RMSE={fit['rmse_cd_m2']:.4f} cd/m2"
        )
        if settings["verification"]:
            if settings["simulate_photometer"] is not None:
                # A corrected simulated display is linear apart from measurement noise.
                rng = np.random.default_rng(1)
                corrected = 0.3 + (250.0 - 0.3) * levels
                corrected += rng.normal(0, 0.0005 * 250.0, size=levels.size)
            else:
                apply_hardware_gamma(window, float(fit["gamma"]))
                corrected = measure("Corrected verification", 1)
            if corrected is None:
                print("Verification aborted; no calibration was saved.")
                return 1
    finally:
        if window is not None:
            try:
                if original_ramp is not None:
                    window.gammaRamp = original_ramp
            finally:
                window.close()

    assert luminances is not None
    is_simulation = settings["simulate_photometer"] is not None
    record_display_name = (
        f"SIMULATION {settings['profile']}" if is_simulation else settings["profile"]
    )
    record = build_calibration_record(
        record_display_name,
        levels,
        luminances,
        display_profile=profile,
        screen_index=int(settings["screen_index"]),
        photometer_model=settings["photometer_model"],
        operator=settings["operator"],
        warmup_minutes=0.0,
        corrected_luminances_cd_m2=corrected,
        measured_refresh_hz=measured_refresh,
    )
    record["simulation"] = {
        "enabled": is_simulation,
        "target_display_name": settings["profile"],
        "requested_gamma": (
            None
            if settings["simulate_photometer"] is None
            else float(settings["simulate_photometer"])
        ),
        "scientific_use_allowed": False if is_simulation else True,
    }
    record["measurement_field"] = {
        "mode": "full_display_uniform",
        "measurement_text_visible": False,
        "entry_screen_separate": True,
    }
    output = Path(settings["output_dir"]).resolve()
    if is_simulation:
        output = output / "simulations"
    output.mkdir(parents=True, exist_ok=True)
    safe_name = _safe_display_filename(record_display_name)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    timestamped_json = save_calibration(output / f"{safe_name}_{timestamp}.json", record)
    latest_json = save_calibration(output / f"{safe_name}.json", record)
    fit = fit_gamma(levels, luminances)
    timestamped_figure = output / f"{safe_name}_{timestamp}_residuals.pdf"
    save_figure(timestamped_figure, levels, luminances, fit, corrected)
    latest_figure = output / f"{safe_name}_residuals.pdf"
    latest_figure.write_bytes(timestamped_figure.read_bytes())
    print(f"Saved timestamped calibration: {timestamped_json}")
    print(f"Saved current calibration:     {latest_json}")
    qc_label = "Rehearsal QC passed" if is_simulation else "Quality control passed"
    print(f"{qc_label + ':':31s}{record['quality_control']['passed']}")
    if not record["quality_control"]["passed"]:
        for reason in record["quality_control"]["failure_reasons"]:
            print(f"  - {reason}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
