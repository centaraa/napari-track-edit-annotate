

# Napari Track Edit

[![tests](https://github.com/live-image-tracking-tools/napari-track-edit/workflows/tests/badge.svg)](https://github.com/live-image-tracking-tools/napari-track-edit/actions)
[![codecov](https://codecov.io/gh/live-image-tracking-tools/napari-track-edit/branch/main/graph/badge.svg)](https://codecov.io/gh/live-image-tracking-tools/napari-track-edit)
[![PyPI version](https://badge.fury.io/py/napari-track-edit.svg)](https://badge.fury.io/py/napari-track-edit)
[![License: BSD 3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

<img src="docs/source/images/logo_transparent.png" align="right" width="180" alt="Napari Track Edit logo" />

A napari plugin for interactive visualization, navigation, and editing of object tracking results. You can open and edit existing tracking data, manually create new tracking results, or run automatic tracking with [motile](https://github.com/funkelab/motile).
The full documentation of the plugin can be found [here](https://liveimagetrackingtools.org/napari-track-edit/).

Motile is a library that makes it easy to solve tracking problems using optimization
by framing the task as an Integer Linear Program (ILP).
See the motile [documentation](https://funkelab.github.io/motile)
for more details on the concepts and method.

----------------------------------

## This fork (centaraa/motile_tracker_annotate)

This repository is a fork of
[napari-track-edit](https://github.com/live-image-tracking-tools/napari-track-edit),
kept in sync with its `main`. It adds:

- editable per-node annotation columns in the tracks table (check boxes,
  numbers, comments of up to 150 characters), saved with the geff;
- measurement columns computed from the masks, refreshed when a mask is edited;
- membrane labels stored in the geff, shown as an extra labels layer (see
  [Membrane labels and other analysis results](#membrane-labels-and-other-analysis-results));
- int32 labels, so node ids above 255 can be selected and painted.

`pip install napari-track-edit` installs the upstream package without these.
To run the fork, clone this repository and start it with
`uv run napari-track-edit` (see [`DEVELOPER.md`](DEVELOPER.md)).

## Installation

Users can download and install an executable application from the github release, or
install from `pypi` in the environment of their choice (e.g. `venv`, `conda`) with the command
`pip install napari-track-edit`.
Currently, napari-track-edit requires python >=3.11.

### Recommended extras

For better performance, you can install optional extras:

- **numba**: Speeds up candidate graph construction significantly.
  ```bash
  pip install napari-track-edit[numba]
  ```

- **gurobi**: Uses the Gurobi solver instead of the default open-source solver. Gurobi is
  much faster but requires a license (free for academics).
  ```bash
  pip install napari-track-edit[gurobi]
  ```

You can install multiple extras at once: `pip install napari-track-edit[numba,gurobi]`.

### Gurobi license version mismatch

If you have a Gurobi license and encounter an error about license version mismatch,
you may need to install a specific version of `gurobipy` that matches your license.
Use one of the version-specific extras:

```bash
pip install napari-track-edit[gurobi12]  # For Gurobi 12.x licenses
pip install napari-track-edit[gurobi13]  # For Gurobi 13.x licenses
```

### For developers

Developers can clone the GitHub repository and then use `uv` to install and run the code.
See the developer guide in [`DEVELOPER.md`](DEVELOPER.md) for more information.

## Usage
Start napari and open all widgets via `Plugins` > `Napari Track Edit` > `Open all widgets`.

From here you can:
- load your own Labels or Points data to [track with motile](https://liveimagetrackingtools.org/napari-track-edit/tracking.html#tracking-with-motile),
- load an image to [manually track objects](https://liveimagetrackingtools.org/napari-track-edit/tracking.html#tracking-from-scratch),
- or [load existing tracking data](https://liveimagetrackingtools.org/napari-track-edit/saving_loading.html#loading-tracks) to explore or edit.

If you would like to see an example first, go to `Plugins` > `Napari Track Edit` > `Widget - Getting started`, and click on one of the two examples, `Hela cells (2D)` or `Mouse embryo (3D)`, to download and view them:

- Fluo-N2DL-HeLa is a 2D dataset of images and segmentations of HeLa cells from the [`cell tracking challenge`](https://celltrackingchallenge.net/2d-datasets/).
- Mouse Embryo Membrane is a 3D dataset of images and segmentations of a membrane labeled developing early mouse embryo (4-26 cells)
from [`Fabrèges et al (2024)`](https://www.science.org/doi/10.1126/science.adh1145) available [`here`](https://zenodo.org/records/13903500).

For details, please read the [documentation](https://liveimagetrackingtools.org/napari-track-edit/).
If you are new to napari-track-edit, you can follow this [tutorial](https://liveimagetrackingtools.org/napari-track-edit/napari-track-edit_tutorial.html) to learn the basics.

https://github.com/user-attachments/assets/cd23271d-bbe6-40c2-80cb-8404136a564a

## Membrane labels and other analysis results

The analysis of the tracked embryos runs headless in
[DestiNuc](https://github.com/centaraa/DestiNuc) (`pixi run -e analysis ...`):
`merge-membranes` turns membrane segmentations into one cell per tracked
nucleus, `spindle-features` measures mitotic spindles. They read the corrected
tracks (geff) and write a new geff; napari-track-edit only shows and saves what
they wrote. Load the output via "Tracks (geff)".

What napari-track-edit reads from such a geff (the contract with DestiNuc):

- **Node features**: any node property is loaded as a feature, shown in the
  table and written to CSV exports, e.g. `membrane_id`, `membrane_volume`
  (voxels), `membrane_qc`, `membrane_boundary_signal`, and the `spindle_*`
  features (µm).
- **Membrane cells**: a uint32 zarr array `membrane_labels` inside the geff
  store, next to `nodes`/`edges`, with the shape of the nuclei segmentation
  (t, z, y, x); a cell's value is the node id of its nucleus. It is registered
  in the geff metadata under `related_objects` (`type: labels`,
  `path: membrane_labels`, `label_prop: membrane_id`). napari-track-edit shows it
  as a read-only `<name>_membrane` overlay coloured like the nuclei, and keeps it
  when the tracks are saved again (`napari_track_edit.membrane.io`).

## Issues

If you encounter any problems, please
[file an issue](https://github.com/live-image-tracking-tools/napari-track-edit/issues)
along with a detailed description.
