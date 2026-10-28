from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import calibrate_display
from calibrate_display import activate_window_backend, apply_hardware_gamma, parse_args
import pytest

from calibration_core import (
    build_calibration_record,
    calibration_levels,
    fit_gamma,
    linearity_error,
    load_calibration,
    luminance_model,
    save_calibration,
)


def test_non_plugin_backend_needs_no_activation():
    assert activate_window_backend("pyglet") is None


def test_pyglet_is_the_portable_default_backend():
    arguments = parse_args([])
    assert arguments.window_backend == "pyglet"
    assert not hasattr(arguments, "patch_fraction")


def test_warmup_cannot_be_requested_from_the_command_line():
    with pytest.raises(SystemExit) as error:
        parse_args(["--warmup-minutes", "30"])
    assert error.value.code == 2


def test_asus_e_gui_and_measurement_start_without_warmup(monkeypatch, tmp_path):
    profiles = calibrate_display.load_profiles(calibrate_display.ROOT / "display_profiles.json")
    profile = profiles["ASUS E"]
    assert (profile["width_px"], profile["height_px"]) == (1920, 1080)
    assert profile["width_cm"] == 58.8
    assert profile["viewing_distance_cm"] == 73.0
    assert profile["expected_refresh_hz"] == 60.0
    assert all(p["width_cm"] == 58.8 for n, p in profiles.items() if n.startswith("ASUS "))
    levels = calibration_levels(3)
    raw = luminance_model(levels, 0.3, 250.0, 2.2)
    corrected = 0.3 + 249.7 * levels
    key_events = []
    for reading in np.concatenate((raw, corrected)):
        key_events.append(["space"])
        key_events.append([
            "period" if key == "." else key for key in f"{reading:.8f}"
        ] + ["return"])
    events = iter(key_events)
    colors = []
    waits = []
    closed = []
    original_ramp = np.linspace(0.0, 1.0, 256)[None, :].repeat(3, axis=0)
    window = SimpleNamespace(
        backend=SimpleNamespace(getGammaRamp=lambda: original_ramp),
        getActualFrameRate=lambda **kwargs: 60.0,
        flip=lambda: None,
        close=lambda: closed.append(True),
    )
    field = SimpleNamespace()
    field.draw = lambda: colors.append(tuple(field.fillColor))

    def dialog(dictionary, **kwargs):
        assert "Warm-up Minutes" not in dictionary
        assert dictionary["Display Profile"][0] == "ASUS E"
        return SimpleNamespace(OK=True)

    def unexpected_warmup(*args, **kwargs):
        raise AssertionError("Unexpected warm-up timer or polling")

    monkeypatch.setitem(sys.modules, "psychopy", SimpleNamespace(
        __version__="test",
        gui=SimpleNamespace(DlgFromDict=dialog),
        core=SimpleNamespace(wait=lambda delay: waits.append(delay), getTime=unexpected_warmup),
        event=SimpleNamespace(
            clearEvents=lambda **kwargs: None,
            waitKeys=lambda **kwargs: next(events),
            getKeys=unexpected_warmup,
        ),
        visual=SimpleNamespace(
            Window=lambda **kwargs: window,
            Rect=lambda *args, **kwargs: field,
            TextStim=lambda *args, **kwargs: SimpleNamespace(draw=lambda: None),
        ),
    ))
    monkeypatch.setattr(calibrate_display, "save_figure", lambda path, *args: Path(path).write_bytes(b"mock figure"))
    assert calibrate_display.main([
        "--profile", "ASUS E", "--levels", "3", "--output-dir", str(tmp_path),
    ]) == 0
    record = load_calibration(tmp_path / "ASUS E.json", expected_display_name="ASUS E")
    assert record["display_profile"] == profile
    assert record["provenance"]["warmup_minutes"] == 0.0
    assert record["provenance"]["measured_refresh_hz"] == 60.0
    assert colors[0] == (-1.0, -1.0, -1.0)
    assert waits == [0.25] * 6
    assert closed == [True]
    np.testing.assert_allclose(window.gammaRamp, original_ramp)


def test_hardware_gamma_failure_has_actionable_message():
    class WindowWithoutGammaRamp:
        @property
        def gamma(self):
            return None

        @gamma.setter
        def gamma(self, value):
            raise OSError("XF86VidModeGetGammaRamp failed")

    try:
        apply_hardware_gamma(WindowWithoutGammaRamp(), 2.2)
    except RuntimeError as error:
        assert "X11/Xorg" in str(error)
    else:
        raise AssertionError("hardware gamma failure was not reported")


def test_glfw_plugin_is_explicitly_activated(monkeypatch):
    requested = []
    scans = []
    plugins = SimpleNamespace(
        scanPlugins=lambda: scans.append(True) or 1,
        loadPlugin=lambda name: requested.append(name) or True,
    )
    monkeypatch.setitem(sys.modules, "psychopy", SimpleNamespace(plugins=plugins))
    assert activate_window_backend("glfw") == "psychopy-glfw"
    assert scans == [True]
    assert requested == ["psychopy-glfw"]


def test_gamma_fit_recovers_known_power():
    levels = calibration_levels(17)
    luminances = luminance_model(levels, 0.3, 250.0, 2.2)
    fit = fit_gamma(levels, luminances)
    assert abs(fit["gamma"] - 2.2) < 1e-4
    assert fit["r_squared"] > 0.999999


def test_linearity_error_is_zero_for_linear_corrected_values():
    levels = calibration_levels(17)
    error = linearity_error(levels, 0.3 + 249.7 * levels)
    assert error["max_abs_error_fraction"] < 1e-12
    assert error["rms_error_fraction"] < 1e-12


def test_saved_verified_record_passes_strict_loader(tmp_path):
    levels = calibration_levels(17)
    raw = luminance_model(levels, 0.3, 250.0, 2.2)
    corrected = 0.3 + 249.7 * levels
    record = build_calibration_record(
        "ASUS B",
        levels,
        raw,
        display_profile={
            "width_px": 1920,
            "height_px": 1080,
            "width_cm": 58.0,
            "viewing_distance_cm": 45.0,
            "expected_refresh_hz": 60.0,
            "fullscreen": True,
            "win_type": "glfw",
        },
        screen_index=0,
        corrected_luminances_cd_m2=corrected,
    )
    path = save_calibration(tmp_path / "ASUS B.json", record)
    loaded = load_calibration(path, expected_display_name="ASUS B")
    assert loaded["quality_control"]["passed"]
    assert loaded["display_name"] == "ASUS B"
    assert json.loads(path.read_text())["schema_version"] == "portable-monitor-calibration-v1"


def test_unverified_record_is_saved_for_diagnosis_but_rejected(tmp_path):
    levels = calibration_levels(17)
    raw = luminance_model(levels, 0.3, 250.0, 2.2)
    record = build_calibration_record(
        "ASUS B",
        levels,
        raw,
        display_profile={},
        screen_index=0,
    )
    assert not record["quality_control"]["passed"]
    path = save_calibration(tmp_path / "unverified.json", record)
    try:
        load_calibration(path)
    except ValueError as error:
        assert "verification" in str(error)
    else:
        raise AssertionError("strict loader accepted an unverified calibration")
