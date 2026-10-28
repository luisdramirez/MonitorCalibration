# Standalone monitor calibration

This folder measures and verifies the luminance transfer function of a visual
stimulus display. It is intentionally independent of any experiment repository:
there are no imports from `sd_continuity`, `naturalness_tuning`, or the pSFT
toolbox.

## Defaults inherited from the lab workflow

- Display: **ASUS B**
- Resolution: **1920 x 1080**
- Physical width: **58.8 cm**
- Viewing distance: **45 cm**
- Refresh rate: **60 Hz**
- PsychoPy screen index: **0**
- PsychoPy window backend: **pyglet**
- Luminance samples: **17**, evenly spaced from code 0 through 1
- Measurement field: **uniform across the full display**
- Scripted display warm-up: **none**
- Corrected verification pass: **required by default**
- Quality gates: gamma 1.0-3.5, fit R2 >= .99, verification maximum and RMS
  errors each <= 5% of luminance range

ASUS C, ASUS D, ASUS E, Macbook, and S32D850 precedents are also included in
`display_profiles.json`. Edit that JSON to add a new lab display; do not embed
experiment-specific settings in the calibration code.

ASUS E uses **1920 x 1080 at 60 Hz**, **58.8 cm** active width and **73 cm**
viewing distance. Set the mirrored displays to `low-res-60` before calibration.
The profile does not change the operating system's resolution or refresh rate.

## Install

Python 3.10 is recommended. In a dedicated environment:

```bash
python -m pip install -r requirements.txt
```

The default `pyglet` backend is built into PsychoPy and avoids Linux EGL
failures observed with GLFW. GLFW remains available through
`--window-backend glfw`; PsychoPy 2026 then requires the separately installed
`psychopy-glfw` plugin, which the script scans for and activates explicitly.

The calibration software is portable, but every generated calibration belongs
to the exact monitor, GPU, cable, operating-system display path, monitor menu
settings, resolution, and refresh rate on which it was measured.

## Rehearse without a photometer

```bash
python calibrate_display.py --simulate-photometer 2.2 --windowed
```

The startup GUI retains all normal defaults. Simulation exercises display
creation, fitting, synthetic verification, figures, and JSON output without
altering the computer's real hardware gamma ramp.

Simulations are written under `calibrations/simulations/`, use a display name
such as `SIMULATION ASUS B`, and are explicitly marked as unsuitable for
scientific use. They never replace the acquisition file `calibrations/ASUS B.json`.

On Linux, real calibration with hardware gamma correction requires a display
session that exposes the gamma lookup table. If PsychoPy reports an
`XF86VidModeGetGammaRamp` failure, check `echo $XDG_SESSION_TYPE`; an X11/Xorg
session is generally required. Do not suppress the error for scientific
acquisition, because doing so would run without the fitted correction.

## Calibrate a display

```bash
python calibrate_display.py
```

The GUI collects the display profile, screen index, number of levels,
photometer model, operator, verification choice, and
output directory.

To start with ASUS E selected, run from the repository root:

```bash
python python/calibrate_display.py --profile "ASUS E"
```

Before measuring:

1. Disable HDR, adaptive brightness, Night Shift/blue-light filtering,
   automatic brightness limiting, dynamic contrast, and eco modes.
2. Use the exact resolution, refresh rate, GPU, cable, and monitor settings that
   will be used for experiments.
3. Darken the room, shield the sensor from ambient light, and place the
   photometer perpendicular to the center of the display.

The script has no warm-up timer or delay setting. After display setup and the
short refresh-rate measurement, it shows a centered instruction screen. The
instructions explain photometer placement, reading the uniform fields, number
entry and the selected verification procedure. Press Space to begin or Escape
to cancel without saving. Saved
provenance records zero minutes of scripted warm-up; it does not measure how
long the physical monitor was powered on.

GUI-based simulated rehearsals show their own instructions. With both
`--no-gui` and `--simulate-photometer`, instructions are skipped so the rehearsal
can run unattended.

Each level appears immediately as an unmarked, full-display uniform field.
Once the photometer reading is stable, press Space; the uniform field is then
replaced by a separate data-entry screen where the typed value is visible. Type
the reading in cd/m2 and press Return to record it and advance directly to the
next uniform level. Backspace edits and Escape aborts without saving. During
verification, the fitted gamma correction is applied and the levels are
measured again.

## Outputs

The default `calibrations/` directory receives:

- `ASUS B.json`: current stable calibration path
- `ASUS B_<UTC timestamp>.json`: immutable measurement record
- Matching residual/verification PDFs

The JSON includes display geometry, raw and corrected measurements, fit,
verification errors, software/host provenance, and explicit quality-control
results. It is directly consumable by the naturalness calibration task when the
display is named `ASUS B`.

Validate a file independently with:

```bash
python verify_calibration.py "calibrations/ASUS B.json" \
  --display-name "ASUS B"
```

Do not use a file whose quality control fails. Keep historical timestamped
records; do not silently replace one display's calibration with another's.

## Google Drive or GitHub

The complete source folder can be transferred through Google Drive, email, or
GitHub. Install dependencies locally on each presentation computer. Running the
code from a synced folder is acceptable, but copy the chosen calibration JSON
to a stable local experiment path before acquisition so syncing cannot alter it
during a session.
