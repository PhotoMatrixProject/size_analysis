# size_analysis

Analysis code, data, and results accompanying **Shifting Dimensions: The Transformation of Reproduction Sizes in Art Magazines (1900–1950)** by **Viktorie Vítů and Petr Tureček**.

The project studies changes in reproduction sizes using relative image area, absolute image width, and relative image width.

## Technical documentation

[Technical details of the statistical analysis](Technical%20details.pdf) describes the temporal models, identification and tracking of size regimes, and comparisons between magazines using Jensen–Shannon divergence. It also includes the WAIC model comparison table.

## Repository structure

| Location | Contents |
| --- | --- |
| `pipeline_area/` | Relative image area analysis, input data, Stan source, figures, and summaries |
| `pipeline_widthabs/` | Absolute image width analysis, input data, and results |
| `pipeline_width_rel/` | Relative image width analysis, input data, and results |
| `Technical details.pdf` | Supplementary statistical methods |
| `combine images.R` | Code for assembling figure panels |
| `Figure*.png` | Combined figures |
| `image_files_inventory.csv`, `image_paths_for_R.txt` | Figure inventory and paths |
| `size.Rproj` | Root RStudio project |
| `IMPORT_EXCLUSIONS.csv` | Manifest of artifacts omitted from the initial import |

Each pipeline includes clustering/standardization and temporal-model scripts. Output folders contain figures and summary tables, including temporal size tracks, standardization/divergence summaries, and model-comparison results.

## Analyses

- **Temporal models:** Bayesian models with no time effect, a linear trend, a spline trend, or step changes. Relative area and relative width are modelled on the logit scale; absolute width is modelled on its original scale.
- **Size regimes:** Density modes are identified in moving time windows and linked across windows to track recurring size regimes.
- **Standardization and divergence:** Binned size distributions are compared within and between magazines, including comparisons with the previous year.

See the technical document for definitions, parameter settings, and interpretation.

## Working with the code

Use R/RStudio. Packages referenced by the scripts include `cmdstanr`, `posterior`, `loo`, `nnet`, `splines`, `randomcoloR`, `rethinking`, and `stringi`. Figure assembly uses `magick`. Bayesian temporal models also require a working CmdStan installation and compiler toolchain.

1. Open the relevant pipeline folder or its RStudio project and inspect the script you intend to run.
2. Check the working directory and input filenames. Each pipeline contains `size_data_2026_10_03.csv`, but some scripts retain references to older data filenames or local paths.
3. For clustering and standardization, use the pipeline's `*clustering_standardization.R` script.
4. For temporal models, inspect the corresponding `time_model_pipeline*.R` script and ensure its referenced `.stan` files are available relative to the working directory. Some Stan sources and script copies are in `time_results/`.
5. Inspect the existing output folders for figures and tables. Use `combine images.R` to assemble figure panels after checking its paths.

The scripts are preserved as supplied, including alternate copies. This repository is an analysis snapshot; running all scripts from the root without adjusting paths is not supported. Analyses have not been rerun as part of the upload.

## Included and excluded files

The initial import from `pipelines_v3` on 2026-10-14 includes 366 original files: R and Stan source, CSV inputs, figures, result tables, project files, and inventories. All 366 copies were verified against the originals using SHA-256 hashes.

Compiled Stan executables, RStudio session/history files, and RStudio export ZIP archives were omitted. Their paths are recorded in `IMPORT_EXCLUSIONS.csv`; `.gitignore` helps prevent these artifacts from being committed later. The `.stan` source files are included so models can be compiled locally through `cmdstanr`.
