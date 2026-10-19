import numpy as np
import pytest

from motile_tracker.spindle.core import cell_features, cell_stats
from motile_tracker.spindle.qc import draw_rotated, rotated_panels, spindle_frame

SP = (1.0, 0.5, 0.5)


def _ellipsoid(shape, center, semi_um, spacing):
    zz, yy, xx = np.indices(shape).astype(float)
    pts = [
        (a - c) * s / r
        for a, c, s, r in zip((zz, yy, xx), center, spacing, semi_um, strict=True)
    ]
    return sum(p * p for p in pts) <= 1.0


def test_cell_stats_volume_diameter_axis():
    mem = np.zeros((40, 80, 80), np.uint32)
    mem[_ellipsoid(mem.shape, (20, 40, 40), (8.0, 8.0, 16.0), SP)] = 7  # long in x
    cell = cell_stats(mem, 7, SP)
    expected = 4 / 3 * np.pi * 8 * 8 * 16
    assert cell["volume_um3"] == pytest.approx(expected, rel=0.05)
    assert cell["diameter_um"] == pytest.approx(
        (6 * expected / np.pi) ** (1 / 3), rel=0.03
    )
    np.testing.assert_allclose(cell["centroid_um"], [20, 20, 20], atol=0.3)
    assert abs(cell["long_axis"][2]) > 0.99  # x
    assert cell["elongation"] == pytest.approx(2.0, rel=0.1)
    assert cell_stats(mem, 99, SP) is None


def test_cell_features_relate_spindle_to_cell():
    cell = {
        "cell_id": 7,
        "volume_um3": 1000.0,
        "diameter_um": 20.0,
        "centroid_um": np.array([20.0, 20.0, 20.0]),
        "long_axis": np.array([0.0, 0.0, 1.0]),
        "elongation": 2.0,
    }
    row = {
        "spindle_status": "ok",
        "spindle_length_um": 10.0,
        "spindle_center_z_um": 20.0,
        "spindle_center_y_um": 23.0,
        "spindle_center_x_um": 24.0,
        "spindle_axis_z": 0.0,
        "spindle_axis_y": 1.0,
        "spindle_axis_x": 0.0,
    }
    f = cell_features(cell, row)
    assert f["spindle_cell_id"] == 7 and f["spindle_cell_volume_um3"] == 1000.0
    assert f["spindle_length_to_cell_diameter"] == pytest.approx(0.5)
    assert f["spindle_center_to_cell_center_um"] == pytest.approx(5.0)
    assert f["spindle_center_offset_rel"] == pytest.approx(0.5)
    assert f["spindle_axis_vs_cell_long_axis_deg"] == pytest.approx(90.0)
    # failed spindle: cell size only; no cell: all missing
    g = cell_features(cell, {"spindle_status": "failed: x"})
    assert g["spindle_cell_volume_um3"] == 1000.0
    assert np.isnan(g["spindle_length_to_cell_diameter"])
    h = cell_features(None, row)
    assert h["spindle_cell_id"] == -1 and np.isnan(h["spindle_cell_volume_um3"])


def test_spindle_frame_is_orthonormal_with_axis_first():
    a, b = np.array([5.0, 10.0, 10.0]), np.array([9.0, 16.0, 22.0])
    center, basis = spindle_frame([a, b])
    np.testing.assert_allclose(center, (a + b) / 2)
    np.testing.assert_allclose(basis @ basis.T, np.eye(3), atol=1e-9)
    np.testing.assert_allclose(basis[0], (b - a) / np.linalg.norm(b - a))
    assert basis[2][0] > 0.5  # "up" leans towards the stack z


def _oblique_case():
    """A rod-shaped spindle along an oblique axis, as a case/result pair."""
    iso = 0.25
    shape = (96, 160, 160)  # 24 x 40 x 40 um on the iso grid
    a = np.array([10.0, 14.0, 12.0])
    b = np.array([14.0, 24.0, 26.0])
    zz, yy, xx = np.indices(shape).astype(float) * iso
    p = np.stack([zz, yy, xx], -1)
    u = (b - a) / np.linalg.norm(b - a)
    t = np.clip((p - a) @ u, 0, np.linalg.norm(b - a))
    d = np.linalg.norm(p - (a + t[..., None] * u), axis=-1)
    rod = d < 1.5
    center = (a + b) / 2
    plate = (np.abs((p - center) @ u) < 0.6) & (np.linalg.norm(p - center, axis=-1) < 4)
    raw_shape = (24, 80, 80)  # spacing 1, 0.5, 0.5
    rz, ry, rx = np.indices(raw_shape).astype(float) * np.array(SP)[:, None, None, None]
    rp = np.stack([rz, ry, rx], -1)
    rt = np.clip((rp - a) @ u, 0, np.linalg.norm(b - a))
    tub = (np.linalg.norm(rp - (a + rt[..., None] * u), axis=-1) < 1.5).astype(
        np.float32
    )
    dna = (
        (np.abs((rp - center) @ u) < 0.6) & (np.linalg.norm(rp - center, axis=-1) < 4)
    ).astype(np.float32)
    case = {"tub": tub, "dna": dna}
    res = {
        "qc": {
            "iso_voxel_um": iso,
            "poles": np.stack([a, b]) / iso,
            "spindle_mask": rod,
            "dna_mask": plate,
        }
    }
    return case, res, np.linalg.norm(b - a)


def test_rotated_panels_put_the_axis_left_to_right():
    case, res, length = _oblique_case()
    panels, marks = rotated_panels(case, res, SP, half_um=12, slab_um=1.0, step_um=0.25)
    top = panels["top"]["spindle_mask"]
    rows, cols = np.nonzero(top)
    # the rod is horizontal: long along the columns (x'), thin along rows
    assert np.ptp(cols) * 0.25 == pytest.approx(length + 3.0, abs=1.0)  # + rounded ends
    assert np.ptp(rows) * 0.25 < 4.0
    # poles project onto the x' axis at -L/2 and +L/2
    np.testing.assert_allclose(
        marks["poles"][:, 0], [-length / 2, length / 2], atol=1e-6
    )
    np.testing.assert_allclose(marks["poles"][:, 1:], 0, atol=1e-6)
    # the plate is seen face-on in the end view, edge-on in the top view
    assert panels["end"]["dna_mask"].sum() > 3 * panels["top"]["dna_mask"].sum()
    assert rotated_panels(case, {"qc": {}}, SP) is None


def test_draw_rotated_writes_png(tmp_path):
    pytest.importorskip("matplotlib")
    case, res, _ = _oblique_case()
    row = {
        "node": 3,
        "t": 5,
        "spindle_frames_before_division": 1,
        "spindle_stage": "metaphase",
        "spindle_status": "ok",
        "spindle_length_um": 10.0,
        "spindle_width_avg_um": 3.0,
        "spindle_angle_deg": 15.0,
        "spindle_axis_dna_offset_deg": 2.0,
        "spindle_n_centrosomes": 0,
    }
    thumbs = []
    assert draw_rotated(case, res, row, SP, tmp_path / "q.png", thumbs, half_um=12)
    assert (tmp_path / "q.png").stat().st_size > 0 and len(thumbs) == 1
    assert not draw_rotated(case, {"qc": {}}, row, SP, tmp_path / "n.png")
