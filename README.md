# Cross-domain evaluation of pretrained 3D tree detectors on airborne LiDAR

Evaluation code and aggregated results for the manuscript *Cross-Domain Evaluation of
Pretrained 3D Tree Detectors: Model Configuration and Evaluation Rules in Airborne LiDAR*,
by Manuel Josue Malla Campoverde, Walter Andrés Osorio-Tinitana and Eduardo Alejandro Tusa
Jumbo (Universidad Técnica de Machala).

Three public pretrained networks for tree instance segmentation (SegmentAnyTree,
ForestFormer3D, TreeLite3D) are applied without retraining to airborne laser scanning (ALS)
of 22 alpine plots in France and Italy, and compared with two classical methods tuned on each
plot (an adaptive 3-D mean shift, AMS3D, and a canopy height model watershed). The reference
is a field inventory of 1437 stems. Detections are matched one to one with stem positions
within 3.5 m inside the inventoried region, intervals come from resampling plots, and the
recall that a random detector would reach with the same number of detections is reported
as a chance floor.

## What is and is not here

This repository holds **code and aggregated results only**.

The ALS point clouds and the field inventories belong to INRAE (formerly IRSTEA) and the
Fondazione Edmund Mach and are **not redistributed here**, nor is anything from which tree
positions could be reconstructed: no point clouds, no per-tree inventory records and no
detection coordinates. Request the data from the original teams; every script reads it
from the path given in `TREELEARN_DATA`.

The adaptive mean shift is third-party work and is not redistributed either: the scripts
that score it read its stored output from `ams_port/out/` under `TREELEARN_DATA`.

## Layout

```
analysis/   the evaluation: one protocol for all five methods
modal/      inference runners for the three networks and the watershed baseline
viewer/     exports of the per-plot reference used by the analysis, and a local viewer
figures/    scripts that write the manuscript tables and the plot maps
results/    aggregated and per-plot numbers, no coordinates
```

## Where each result of the manuscript comes from

| Result | Scripts | Results file |
|---|---|---|
| Detection table: recall, precision, F1 and excess over chance, with intervals | `analysis/extra_analyses.py`, `figures/manuscript_tables.py` | `extra_analyses.json` |
| Paired differences against SegmentAnyTree, plain and stratified by site | `analysis/extra_analyses.py`, `analysis/paired_comparison.py` | `extra_analyses.json`, `paired_comparison.json` |
| Crown strata and species, with the chance floor of each stratum | `analysis/main_evaluation.py`, `analysis/extra_analyses.py`, `analysis/species_recall.py`, `analysis/species_table.py` | `main_evaluation.json`, `extra_analyses.json`, `species_recall.json` |
| TreeLite3D: published and retuned clustering, leave one plot out and leave one site out | `analysis/treelite_results.py`, `analysis/extra_analyses.py`, `modal/run_treelite.py` | `treelite_evaluation.json`, `extra_analyses.json` |
| AMS3D chosen by its own protocol, by an oracle, and on the other plots of the same site or of both sites | `analysis/review_checks.py`, `analysis/ams3d_grid.py` | `review_checks.json`, `ams3d_grid.json` |
| Watershed grid | `modal/watershed_sweep.py` | `watershed_sweep.json` |
| Matching tolerance, optimal assignment and scored region | `analysis/main_evaluation.py`, `analysis/review_checks2.py`, `modal/single_protocol_metrics.py` | `main_evaluation.json`, `consistent_results.json` |
| Chance floor by tolerance | `analysis/chance_floor.py` | `chance_floor.json` |
| Crown Jaccard index | `analysis/crown_jaccard.py`, `modal/jaccard.py` | `crown_jaccard.json`, `jaccard.json` |
| Computing cost | `analysis/measure_runtimes.py`, `analysis/runtime_table.py` | `measure_runtimes.json`, `runtime_summary.json` |
| Per-plot appendix table | `analysis/extra_analyses.py`, `figures/manuscript_tables.py` | `extra_analyses.json` |

`analysis/measure_runtimes.py` times the mean shift itself, so it needs the AMS3D
implementation, which is not here; its output is in `results/measure_runtimes.json`.

## Reproducing the numbers

Set up the environment once:

```
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Two variables locate what is not in the repository: `TREELEARN_DATA` points at the
directory that holds `IRSTEA_dataset/`, `trento_dataset/` and `ams_port/out/`, and
`TREELEARN_ROOT` at the one that holds each model's output folder (`sat_output_all/`,
`ff3d_output/`, `treelite_output/`). Both default to the clone and its parent.

First rebuild the per-plot reference under `viewer/data/`. It is not committed, because
it carries tree positions, and everything downstream reads it:

```
python3 viewer/export_chm.py                  # canopy model, region and per-plot statistics
python3 viewer/export_trees.py                # inventory, detections, pairs and metrics
python3 viewer/patch_exports.py               # coverage, and the detections outside the region
python3 viewer/export_species.py              # species of every inventory tree
```

Then, in this order:

```
python3 analysis/gather_detections.py         # all five methods' detections in one file
python3 analysis/evaluate_ff3d.py             # ForestFormer3D, plot by plot
python3 analysis/treelite_results.py          # the TreeLite3D sweep and its chosen setting
python3 viewer/add_viewer_methods.py          # adds those two to the per-plot reference
python3 analysis/main_evaluation.py           # main counts, matching, region, tolerance, strata
python3 analysis/species_recall.py            # recall by species
python3 analysis/crown_jaccard.py             # crown overlap index
python3 analysis/paired_comparison.py         # paired differences against SegmentAnyTree
python3 analysis/chance_floor.py              # what a uniform detector recovers, per method
python3 analysis/ams3d_grid.py                # the mean shift's 25 configurations per plot
python3 analysis/review_checks.py             # AMS3D choices, equivalence, optimal matching
python3 analysis/extra_analyses.py            # intervals, chance per plot and stratum, per-plot table
python3 figures/manuscript_tables.py          # LaTeX rows of the detection and per-plot tables
```

Each script writes its result next to `TREELEARN_ROOT`, and reads the ones before it from
there, falling back to the copy published in `results/`. So you can run any single step
against the published inputs without rerunning the whole chain. No number in the
manuscript tables is typed by hand.

## The three networks

Run from the authors' released pipelines and weights:

- SegmentAnyTree: the authors' Docker image, with a fix to the export of results to LAS
  files that does not change the predictions (`modal/run_sat.py`).
- ForestFormer3D: the authors' inference repository and the weights deposited on Zenodo
  (`modal/download_ff3d_weights.py`, `modal/run_ff3d.py`).
- TreeLite3D: the authors' repository on Pointcept 1.5.0 (`modal/run_treelite.py`). Its
  clustering radius and minimum proposal size were swept; the sweep is in the same file.

SegmentAnyTree and ForestFormer3D use their published defaults.

## License

Code released under the MIT License (`LICENSE`). The data are not covered by it.
