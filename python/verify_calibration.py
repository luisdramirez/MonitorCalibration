#!/usr/bin/env python3
"""Validate a portable monitor-calibration JSON against frozen quality gates."""

from __future__ import annotations

import argparse
import json

from calibration_core import load_calibration


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("calibration")
    parser.add_argument("--display-name")
    parser.add_argument("--allow-unverified", action="store_true")
    args = parser.parse_args()
    record = load_calibration(
        args.calibration,
        expected_display_name=args.display_name,
        require_verification=not args.allow_unverified,
    )
    print(
        json.dumps(
            {
                "display_name": record["display_name"],
                "gamma": record["gamma"],
                "fit": record["fit"],
                "verification": record.get("verification"),
                "quality_control": record["quality_control"],
                "source_path": record["source_path"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
