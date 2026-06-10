"""
Explicit detector / region-map selection for the CSV->PyG conversion.

The convert config (`convert_csv_to_pyg_sets.yaml`) holds one region map per detector
geometry under `region_labels_by_detector` (keys: 'generic' = GenericDetector/FATRAS,
'odd' = OpenDataDetector/Geant4). The caller picks one explicitly via `--detector`; there
is intentionally NO auto-detection. `resolve_region_labels` writes the chosen map into
`config['region_labels']`, which is what `add_region_labels` consumes.
"""

import sys

VALID_DETECTORS = ("generic", "odd")


def resolve_region_labels(config, detector=None):
    """Set ``config['region_labels']`` to the map for the chosen detector.

    detector: 'generic' or 'odd' (case-insensitive). If None, prompt interactively when a
    TTY is attached, otherwise raise (so batch jobs fail loudly instead of guessing).
    Backward-compat: if the config already has a plain ``region_labels`` and no
    ``region_labels_by_detector``, it is left untouched.
    """
    maps = config.get("region_labels_by_detector")
    if maps is None:
        if "region_labels" in config:
            return config  # legacy single-map config; nothing to select
        raise KeyError(
            "config has neither 'region_labels_by_detector' nor 'region_labels'"
        )

    if detector is None:
        if sys.stdin.isatty():
            detector = input(f"Select detector {VALID_DETECTORS}: ").strip()
        else:
            raise SystemExit(
                f"ERROR: --detector is required (one of {VALID_DETECTORS}) "
                "in non-interactive mode."
            )

    key = str(detector).strip().lower()
    if key not in maps:
        raise SystemExit(
            f"ERROR: --detector must be one of {list(maps)}, got '{detector}'."
        )

    config["region_labels"] = maps[key]
    print(f"[region_labels] using '{key}' detector map "
          f"({len(maps[key])} regions)")
    return config
