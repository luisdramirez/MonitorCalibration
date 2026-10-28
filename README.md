# Monitor calibration

This folder keeps the lab's two monitor-calibration implementations separate:

- `python/` contains the portable, experiment-agnostic PsychoPy workflow. This
  is the recommended implementation for new calibration measurements.
- `matlab/` contains the legacy MATLAB calibration and timing scripts.

See each subfolder's README for requirements and usage. Calibration files made
by the Python workflow are written to `python/calibrations/` by default.
