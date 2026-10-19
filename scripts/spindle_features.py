"""Measure spindles of tracked nuclei before division and store them as features.

Usage:
    python scripts/spindle_features.py TRACKS TUBULIN DNA OUT.geff
        [--config FILE] [--spacing 1,0.26,0.26] [--max-frame T]
        [--timepoints-before 6] [--correction ignore|none]
        [--threshold-exclude-um 1.0] [--pole-profile-radius java|um]
        [--half-size-um 20] [--z-extra-um 4] [--workers N]
        [--membrane LABELS] [--qc-dir DIR] [--qc-debug | --no-qc]
    python scripts/spindle_features.py --write-default-config FILE

TRACKS is a geff store (corrected tracks; node id = label in the segmentation).
TUBULIN (microtubules, e.g. 488) and DNA (e.g. 561) are zarr arrays, tiffs or
folders of one tiff per timepoint, with the same shape as the segmentation.

Parameters: the defaults (Spindle3D v0.8.0 for [spindle3d]), then the values of
--config FILE (a TOML file with sections [spindle3d] and [pipeline], see
src/motile_tracker/spindle/default_config.toml), then the explicit options
below. Unknown keys in the file are an error. The voxel size is --spacing, else
[pipeline] spacing, else the one stored with the tracks.

OUT.geff gets a copy of the tracks with the spindle_* node features; next to it
<stem>_config.toml (the complete resolved config: pass it back with --config to
reproduce the run; spindle_view.py reads it), <stem>_spindle.csv (one row per
candidate), <stem>_threshold_diagnostic.csv and, unless --no-qc, QC PNGs in
--qc-dir (default <out dir>/qc_<stem>, one per candidate + overview.png).
See src/motile_tracker/spindle/README.md.
"""

import argparse
import csv
import functools
import time
from pathlib import Path

from funtracks.import_export import import_from_geff

from motile_tracker.import_export.geff_io import write_geff_over
from motile_tracker.membrane.io import (
    load_membrane_labels,
    read_label_file,
    save_membrane_labels,
)
from motile_tracker.spindle import add_spindle_features, measure_tracks
from motile_tracker.spindle.config import (
    DEFAULT_CONFIG_HEADER,
    SpindleConfig,
    load_config,
    write_config,
)
from motile_tracker.spindle.core import resolve_spacing

print = functools.partial(print, flush=True)  # noqa: A001

# CLI option -> (config section, key); only options given explicitly override
OVERRIDES = {
    "spacing": ("pipeline", "spacing"),
    "max_frame": ("pipeline", "max_frame"),
    "timepoints_before": ("pipeline", "timepoints_before"),
    "correction": ("pipeline", "correction"),
    "threshold_exclude_um": ("pipeline", "threshold_exclude_um"),
    "half_size_um": ("pipeline", "half_size_um"),
    "z_extra_um": ("pipeline", "z_extra_um"),
    "pole_profile_radius": ("spindle3d", "pole_profile_radius"),
}


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        epilog="Options left out take their value from --config, else the default.",
    )
    p.add_argument("tracks", type=Path, nargs="?")
    p.add_argument("tubulin", type=Path, nargs="?")
    p.add_argument("dna", type=Path, nargs="?")
    p.add_argument("out_geff", type=Path, nargs="?")
    p.add_argument("--config", type=Path, default=None, help="config TOML file")
    p.add_argument(
        "--write-default-config",
        type=Path,
        default=None,
        metavar="FILE",
        help="write the default config to FILE and exit",
    )
    p.add_argument(
        "--spacing", default=None, help="voxel size z,y,x in um ([pipeline] spacing)"
    )
    p.add_argument(
        "--max-frame", type=int, default=None, help="only divisions up to this frame"
    )
    p.add_argument(
        "--timepoints-before", type=int, default=None, help="frames per division"
    )
    p.add_argument("--correction", choices=["ignore", "none"], default=None)
    p.add_argument("--threshold-exclude-um", type=float, default=None)
    p.add_argument("--pole-profile-radius", choices=["java", "um"], default=None)
    p.add_argument("--half-size-um", type=float, default=None)
    p.add_argument("--z-extra-um", type=float, default=None)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument(
        "--membrane",
        type=Path,
        default=None,
        help="membrane cell labels (zarr/tiff/folder) for the spindle_cell_* "
        "features; default: the input geff's membrane_labels, if any",
    )
    p.add_argument(
        "--qc",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="QC images: one spindle-aligned PNG per candidate (axis left to "
        "right; DNA, tubulin, DNA/spindle/cell outlines) + overview.png "
        "(--no-qc: none, no matplotlib needed)",
    )
    p.add_argument(
        "--qc-debug",
        action="store_true",
        help="also the threshold panels (MIPs of the raw crop) per candidate",
    )
    p.add_argument("--qc-dir", type=Path, default=None)
    args = p.parse_args(argv)

    if args.write_default_config is not None:
        write_config(SpindleConfig(), args.write_default_config, DEFAULT_CONFIG_HEADER)
        print(f"wrote {args.write_default_config}")
        return
    if None in (args.tracks, args.tubulin, args.dna, args.out_geff):
        p.error("TRACKS, TUBULIN, DNA and OUT.geff are required")

    in_geff, out_geff = args.tracks.resolve(), args.out_geff.resolve()
    if in_geff == out_geff:
        raise SystemExit("OUT.geff must differ from the input")
    overrides: dict = {"pipeline": {}, "spindle3d": {}}
    for opt, (section, key) in OVERRIDES.items():
        value = getattr(args, opt)
        if value is None:
            continue
        if opt == "spacing":
            value = [float(v) for v in value.split(",")]
        overrides[section][key] = value
    config = load_config(args.config, overrides)
    tracks = import_from_geff(in_geff)
    config.pipeline = resolve_spacing(config.pipeline, tracks)
    if all(s == 1.0 for s in config.pipeline.spacing):
        # The spindle values are in um (and um3), so a voxel size of 1 taken
        # over from a geff that stores none would make them pixel values.
        raise SystemExit(
            "voxel size 1,1,1 (the tracks store no physical scale): pass "
            "--spacing z,y,x in um (or [pipeline] spacing), e.g. --spacing 1,0.26,0.26"
        )
    params, settings = config.pipeline, config.spindle3d
    run = out_geff.stem
    out_dir = out_geff.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_dir = (args.qc_dir or out_dir / f"qc_{run}") if args.qc else None
    if qc_dir is not None:
        from motile_tracker.spindle.qc import contact_sheet, draw_qc, draw_rotated

        qc_dir.mkdir(parents=True, exist_ok=True)
        for old in qc_dir.glob("*.png"):
            old.unlink()
    write_config(
        config,
        out_dir / f"{run}_config.toml",
        f"Resolved config of the spindle run that wrote {out_geff.name}\n"
        f"tracks: {in_geff}\ntubulin: {args.tubulin}\ndna: {args.dna}\n"
        "Reproduce: scripts/spindle_features.py TRACKS TUBULIN DNA OUT.geff "
        f"--config {run}_config.toml",
    )

    tubulin = read_label_file(args.tubulin)
    dna = read_label_file(args.dna)
    in_labels = load_membrane_labels(in_geff)
    if args.membrane is not None:
        membrane = read_label_file(args.membrane)
        print(f"cell features from {args.membrane}")
    elif in_labels is not None:
        membrane = in_labels
        print(f"cell features from the membrane labels in {in_geff.name}")
    else:
        membrane = None
        print("no membrane labels: spindle_cell_* features are left NaN")
    if membrane is not None and tuple(membrane.shape) != tuple(
        tracks.segmentation.shape
    ):
        raise SystemExit(
            f"membrane labels {tuple(membrane.shape)} do not match the "
            f"segmentation {tuple(tracks.segmentation.shape)}"
        )
    t0 = time.perf_counter()
    rows, diag, thumbs = [], [], []
    for row, case, res in measure_tracks(
        tracks,
        tubulin,
        dna,
        params,
        settings,
        workers=args.workers,
        membrane=membrane,
    ):
        rows.append(row)
        if res.get("diag"):
            diag.append({"node": row["node"], "t": row["t"], **res["diag"]})
        if res.get("traceback"):
            (out_dir / "logs").mkdir(exist_ok=True)
            (out_dir / "logs" / f"{run}_node{row['node']}_t{row['t']}.txt").write_text(
                res["traceback"], encoding="utf-8"
            )
        if qc_dir is not None and case:
            tag = row["spindle_status"].split(":")[0]
            name = f"node{row['node']}_t{row['t']}_{tag}"
            drawn = draw_rotated(
                case, res, row, params.spacing, qc_dir / f"{name}.png", thumbs
            )
            if args.qc_debug or not drawn:  # no poles: the raw-crop panels
                draw_qc(
                    case,
                    res,
                    row,
                    params.spacing,
                    qc_dir / f"{name}_{'debug' if drawn else 'crop'}.png",
                    None if drawn else thumbs,
                )
        print(
            f"  node {row['node']} t{row['t']} -{row['spindle_frames_before_division']}: "
            f"{row['spindle_status']} | {row.get('spindle_stage', '')} | "
            f"L={row['spindle_length_um']:.1f} W={row['spindle_width_avg_um']:.1f} "
            f"V={row['spindle_volume_um3']:.0f} vsN={row['spindle_axis_dna_offset_deg']:.1f} "
            f"nC={row['spindle_n_centrosomes']}"
        )

    rows.sort(key=lambda r: (r["t"], r["node"]))
    if qc_dir is not None:
        order = {r["node"]: i for i, r in enumerate(rows)}
        thumbs.sort(key=lambda th: order.get(int(th[0].split()[0]), 0))
        contact_sheet(thumbs, qc_dir / "overview.png")
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(out_dir / f"{run}_spindle.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    if diag:
        with open(
            out_dir / f"{run}_threshold_diagnostic.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as fh:
            w = csv.DictWriter(fh, fieldnames=list(diag[0]))
            w.writeheader()
            w.writerows(sorted(diag, key=lambda r: (r["t"], r["node"])))
    add_spindle_features(tracks, rows)
    write_geff_over(tracks, out_geff)
    if in_labels is not None:
        # keep the membrane cells, so spindle and membrane steps can run in
        # either order
        print("copied membrane labels:", save_membrane_labels(out_geff, in_labels))
    n_ok = sum(r["spindle_status"] == "ok" for r in rows)
    n_meta = sum(
        r["spindle_status"] == "ok" and r.get("spindle_is_metaphase") == 1 for r in rows
    )
    print(
        f"[{run}] {n_ok} ok ({n_meta} metaphase), {len(rows) - n_ok} failed of {len(rows)} "
        f"candidates; {time.perf_counter() - t0:.0f}s; geff {out_geff}"
    )


if __name__ == "__main__":
    main()
