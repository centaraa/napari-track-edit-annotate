"""QC images of the spindle measurement (optional; needs matplotlib).

matplotlib is imported lazily: it is only a transitive dependency (via napari).
"""

from __future__ import annotations

import numpy as np

from .core import ISO as _ISO


def _extent(shape2, sc2):
    (ny, nx), (sy, sx) = shape2, sc2
    return [-sx / 2, (nx - 0.5) * sx, (ny - 0.5) * sy, -sy / 2]


def _contour(ax, mask2d, sc, color, lw=1.0):
    if mask2d is None or not mask2d.any():
        return
    ax.contour(
        np.arange(mask2d.shape[1]) * sc[1],
        np.arange(mask2d.shape[0]) * sc[0],
        mask2d.astype(float),
        levels=[0.5],
        colors=color,
        linewidths=lw,
    )


def draw_qc(case: dict, res: dict, row: dict, spacing, png_path, thumbs=None) -> None:
    """One PNG per node: MIPs along z and y of 488, 561 and the masks/poles.

    thumbs: if a list, a thumbnail of the overlay panel is appended for
    `contact_sheet`.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sp = np.asarray(spacing, float)
    qc = res.get("qc", {}) if res else {}
    ISO = qc.get("iso_voxel_um", _ISO)  # noqa: N806 - isotropic analysis voxel (um)
    fig, axes = plt.subplots(2, 4, figsize=(17, 9))
    chans = [
        ("488 (bg-subtracted, used as is)", case["a_b"], "gray"),
        ("488 + threshold exclusion (orange) / shell (yellow)", case["tub"], "gray"),
        ("561 DNA (other nuclei removed)", case["dna"], "magma"),
        (
            (
                "overlay: DNA mask (magenta), spindle mask (cyan), threshold exclusion "
                "(orange),\nrefined poles (lime), centrosomes (red x)"
            ),
            case["tub"],
            "gray",
        ),
    ]
    views = [(0, "MIP along z (XY)", (1, 2)), (1, "MIP along y (XZ)", (0, 2))]
    for r, (pax, lab, keep) in enumerate(views):
        for c, (title, arr, cmap) in enumerate(chans):
            ax = axes[r, c]
            mip = arr.max(axis=pax)
            ext = _extent(mip.shape, sp[list(keep)])
            vmax = max(float(np.percentile(mip, 99.7)), 1e-3)
            ax.imshow(
                mip, cmap=cmap, vmin=0, vmax=vmax, extent=ext, interpolation="nearest"
            )
            ax.set_title(f"{title}\n{lab}", fontsize=8)
            ax.set_xlabel("x (um)")
            ax.set_ylabel("y (um)" if r == 0 else "z (um)")
            if c == 1:
                for key, col in (
                    ("threshold_exclusion", "orange"),
                    ("threshold_shell", "yellow"),
                ):
                    if key in qc:
                        _contour(ax, qc[key].max(axis=pax), (ISO, ISO), col, 0.8)
            if c == 3:
                if "dna_mask" in qc:
                    _contour(ax, qc["dna_mask"].max(axis=pax), (ISO, ISO), "magenta")
                elif case.get("stage_mask") is not None:
                    _contour(
                        ax, case["stage_mask"].max(axis=pax), sp[list(keep)], "magenta"
                    )
                if "spindle_mask" in qc:
                    _contour(ax, qc["spindle_mask"].max(axis=pax), (ISO, ISO), "cyan")
                if (
                    "threshold_exclusion" in qc
                    and case["meta"]["correction"] == "ignore"
                ):
                    _contour(
                        ax,
                        qc["threshold_exclusion"].max(axis=pax),
                        (ISO, ISO),
                        "orange",
                    )
                if "poles" in qc:
                    P = qc["poles"] * ISO
                    ax.plot(
                        P[:, keep[1]],
                        P[:, keep[0]],
                        "o-",
                        ms=7,
                        mfc="none",
                        mec="lime",
                        color="lime",
                        lw=0.8,
                        mew=2,
                    )
                if "centrosomes" in qc:
                    C = qc["centrosomes"] * ISO
                    ax.plot(C[:, keep[1]], C[:, keep[0]], "x", ms=9, color="red", mew=2)
            ax.set_xlim(ext[:2])
            ax.set_ylim(ext[2:])
    st = f"{row['spindle_status']} | stage: {row.get('spindle_stage', '')}"
    fig.suptitle(
        f"node {row['node']}  frame {row['t']}  frames before division "
        f"{row['spindle_frames_before_division']}  (division node {row['division_node']})\n"
        f"length {row['spindle_length_um']:.1f} um  width {row['spindle_width_avg_um']:.1f} um  "
        f"volume {row['spindle_volume_um3']:.0f} um3  vs plate normal "
        f"{row['spindle_axis_dna_offset_deg']:.1f} deg  centrosomes {row['spindle_n_centrosomes']}  "
        f"plate flatness {row['spindle_plate_flatness']:.2f}\nstatus: {st}",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(png_path, dpi=70)
    if thumbs is not None:
        fig.canvas.draw()
        bb = axes[0, 3].get_window_extent()
        buf = np.asarray(fig.canvas.buffer_rgba())[..., :3]
        H = buf.shape[0]
        label = (
            f"{row['node']} t{row['t']} -{row['spindle_frames_before_division']} "
            f"{str(row.get('spindle_stage', '')).upper()}\n{row['spindle_status'][:40]}"
        )
        thumbs.append(
            (
                label,
                buf[int(H - bb.y1) : int(H - bb.y0), int(bb.x0) : int(bb.x1)].copy(),
            )
        )
    plt.close(fig)


def contact_sheet(thumbs, path, ncol: int = 8) -> None:
    """All thumbnails of `draw_qc` on one page."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not thumbs:
        return
    nrow = int(np.ceil(len(thumbs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol, 2.9 * nrow))
    axes = np.atleast_2d(axes)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, (lab, img) in zip(axes.ravel(), thumbs, strict=False):
        ax.imshow(img)
        ax.set_title(lab, fontsize=6)
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


# ---------------------------------------------------------------------------
# spindle-aligned view (as Java Spindle3D's output image)
# ---------------------------------------------------------------------------
def spindle_frame(poles_um) -> tuple[np.ndarray, np.ndarray]:
    """Centre and orthonormal basis (rows ex, ey, ez; z, y, x components) of the
    spindle-aligned frame: ex along the axis pole A -> pole B, ez as close to
    the stack z axis as possible (so the side view shows the coverslip tilt)."""
    a, b = np.asarray(poles_um, float)[:2]
    ex = (b - a) / max(float(np.linalg.norm(b - a)), 1e-9)
    ref = np.array([1.0, 0.0, 0.0])
    if abs(ex @ ref) > 0.98:  # axis along z: use y as "up"
        ref = np.array([0.0, 1.0, 0.0])
    ez = ref - (ref @ ex) * ex
    ez /= np.linalg.norm(ez)
    ey = np.cross(ez, ex)
    return (a + b) / 2, np.stack([ex, ey, ez])


def sample_aligned(vol, voxel_um, center, basis, u, v, w, order=1) -> np.ndarray:
    """`vol` (voxel size `voxel_um`, z,y,x) sampled at center + u*ex + v*ey +
    w*ez (um, crop frame; u, v, w same shape); 0 outside."""
    from scipy import ndimage as ndi

    pts = (
        center[:, None]
        + basis[0][:, None] * np.ravel(u)[None]
        + basis[1][:, None] * np.ravel(v)[None]
        + basis[2][:, None] * np.ravel(w)[None]
    )
    idx = pts / np.asarray(voxel_um, float)[:, None]
    out = ndi.map_coordinates(
        np.asarray(vol, np.float32), idx, order=order, mode="constant", cval=0.0
    )
    return out.reshape(np.shape(u))


def rotated_panels(case, res, spacing, half_um=18.0, slab_um=1.5, step_um=0.25):
    """Slab maximum projections of the spindle-aligned cube, per view.

    Returns {view: {layer: 2D array}} with views "top" (rows y' bottom-up,
    cols x'), "side" (rows z' top-down, cols x') and "end" (rows z' top-down,
    cols y', slab across the plate), and the projected poles / centrosomes as
    aligned (x', y', z') coordinates; None without poles.
    """
    qc = res.get("qc", {}) if res else {}
    if "poles" not in qc or len(qc["poles"]) < 2:
        return None
    sp = np.asarray(spacing, float)
    iso = qc.get("iso_voxel_um", _ISO)
    poles = np.asarray(qc["poles"], float) * iso
    center, basis = spindle_frame(poles)
    g = np.arange(-half_um, half_um + step_um / 2, step_um)
    s = np.arange(-slab_um, slab_um + step_um / 2, step_um)
    gu, gv, gs = np.meshgrid(g, g, s, indexing="xy")  # rows: 2nd axis ascending
    grids = {
        "top": (gu, gv, gs),
        "side": (gu, gs, gv[::-1]),
        "end": (gs, gu, gv[::-1]),
    }
    layers = {
        "tub": (case.get("tub"), sp, 1),
        "dna": (case.get("dna"), sp, 1),
        "dna_mask": (qc.get("dna_mask"), (iso,) * 3, 0),
        "spindle_mask": (qc.get("spindle_mask"), (iso,) * 3, 0),
        "cell": (case.get("cell_crop"), sp, 0),
    }
    panels = {}
    for name, (u, v, w) in grids.items():
        panels[name] = {
            key: sample_aligned(vol, vs, center, basis, u, v, w, order).max(axis=2)
            for key, (vol, vs, order) in layers.items()
            if vol is not None
        }

    def project(points):
        return (np.asarray(points, float).reshape(-1, 3) - center) @ basis.T

    marks = {"poles": project(poles)}
    if "centrosomes" in qc and len(qc["centrosomes"]):
        marks["centrosomes"] = project(np.asarray(qc["centrosomes"], float) * iso)
    return panels, marks


def draw_rotated(
    case: dict,
    res: dict,
    row: dict,
    spacing,
    png_path,
    thumbs=None,
    half_um: float = 18.0,
    slab_um: float = 1.5,
) -> bool:
    """One PNG per node with the spindle axis left to right (like the output
    image of Java Spindle3D: a 36 um cube aligned on the spindle).

    Three slab maximum projections (+-`slab_um` um) around the spindle centre:
    top view (x' = axis, y'), side view (x' = axis, z' up = towards the stack
    z) and end-on through the plate (y', z'). DNA magenta, tubulin green;
    outlines: DNA mask yellow, spindle mask cyan, the cell white (with
    membrane labels); poles lime, centrosomes red. All in um. Returns False
    (nothing drawn) when the spindle has no poles.
    """
    out = rotated_panels(case, res, spacing, half_um, slab_um)
    if out is None:
        return False
    panels, marks = out
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def norm(a, ref):
        return np.clip(a / max(float(np.percentile(ref, 99.7)), 1e-6), 0, 1)

    tub_ref = np.concatenate([p["tub"].ravel() for p in panels.values()])
    dna_ref = np.concatenate([p["dna"].ravel() for p in panels.values()])
    h = half_um
    views = (
        ("top", "top view: x' = spindle axis", (0, 1), "y' (um)", [-h, h, -h, h]),
        (
            "side",
            "side view: x' = spindle axis, z' up",
            (0, 2),
            "z' (um)",
            [-h, h, -h, h],
        ),
        ("end", "end-on through the plate", (1, 2), "z' (um)", [-h, h, -h, h]),
    )
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.8))
    for ax, (name, title, cols, ylab, ext) in zip(axes, views, strict=True):
        p = panels[name]
        rgb = np.zeros((*p["tub"].shape, 3), np.float32)
        d, t = norm(p["dna"], dna_ref), norm(p["tub"], tub_ref)
        rgb[..., 0], rgb[..., 1], rgb[..., 2] = d, t, d
        origin = "lower" if name == "top" else "upper"
        ax.imshow(rgb, extent=ext, origin=origin, interpolation="nearest")
        n = p["tub"].shape[0]
        yy = np.linspace(-h, h, n) if origin == "lower" else np.linspace(h, -h, n)
        xx = np.linspace(-h, h, p["tub"].shape[1])
        for key, col, lw in (
            ("cell", "white", 1.2),
            ("dna_mask", "yellow", 1.0),
            ("spindle_mask", "cyan", 1.0),
        ):
            if key in p and p[key].any():
                ax.contour(xx, yy, p[key], levels=[0.5], colors=col, linewidths=lw)
        pp = marks["poles"]
        ax.plot(
            pp[:, cols[0]], pp[:, cols[1]], "o", ms=7, mfc="none", mec="lime", mew=2
        )
        if "centrosomes" in marks:
            cc = marks["centrosomes"]
            ax.plot(cc[:, cols[0]], cc[:, cols[1]], "x", ms=9, color="red", mew=2)
        ax.set_xlim(-h, h)
        ax.set_ylim(-h, h)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("y' (um)" if name == "end" else "x' (um)")
        ax.set_ylabel(ylab)
    cell = ""
    if np.isfinite(row.get("spindle_cell_diameter_um", np.nan)):
        cell = (
            f"\ncell {row['spindle_cell_volume_um3']:.0f} um3 "
            f"(diameter {row['spindle_cell_diameter_um']:.1f} um)  "
            f"length / cell diameter {row['spindle_length_to_cell_diameter']:.2f}  "
            f"centre offset {row['spindle_center_offset_rel']:.2f} x radius  "
            f"axis vs cell long axis {row['spindle_axis_vs_cell_long_axis_deg']:.0f} deg"
        )
    fig.suptitle(
        f"node {row['node']}  t{row['t']}  -{row['spindle_frames_before_division']}  "
        f"{str(row.get('spindle_stage', '')).upper()}  ({row['spindle_status']})\n"
        f"length {row['spindle_length_um']:.1f} um  width "
        f"{row['spindle_width_avg_um']:.1f} um  tilt {row['spindle_angle_deg']:.0f} deg  "
        f"vs plate normal {row['spindle_axis_dna_offset_deg']:.0f} deg  "
        f"centrosomes {row['spindle_n_centrosomes']}{cell}\n"
        "DNA magenta, tubulin green; outlines DNA yellow, spindle cyan, cell white; "
        "poles lime, centrosomes red",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(png_path, dpi=70)
    if thumbs is not None:
        fig.canvas.draw()
        bb = axes[1].get_window_extent()
        buf = np.asarray(fig.canvas.buffer_rgba())[..., :3]
        hh = buf.shape[0]
        label = (
            f"{row['node']} t{row['t']} -{row['spindle_frames_before_division']} "
            f"{str(row.get('spindle_stage', '')).upper()}\n{row['spindle_status'][:40]}"
        )
        crop = buf[int(hh - bb.y1) : int(hh - bb.y0), int(bb.x0) : int(bb.x1)]
        thumbs.append((label, crop.copy()))
    plt.close(fig)
    return True
