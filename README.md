# size_analysis

R and Stan pipelines for image area, absolute width, and relative width analysis.

Imported from `pipelines_v3` on 2026-10-14, preserving the original folder structure and file contents.

- `pipeline_area/`: area analysis, data, Stan source, figures, and summaries.
- `pipeline_widthabs/`: absolute width analysis and results.
- `pipeline_width_rel/`: relative width analysis and results.
- Root files include combined figures, `combine images.R`, image inventories, and the RStudio project.

The import includes CSV inputs and existing output figures/tables. Compiled Stan executables, RStudio session state/history, and RStudio export ZIPs are excluded. `IMPORT_EXCLUSIONS.csv` records all omitted files; `.gitignore` prevents those artifacts from being added later. Stan executables can be rebuilt from the `.stan` source using CmdStan through `cmdstanr`.

Scripts are preserved as supplied. Run them from their expected working directories and check input paths first; some scripts refer to local paths or older dataset filenames. This import has been checked for file-copy integrity, but the analyses have not been rerun.
