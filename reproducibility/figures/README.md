# Figure reproduction

Run `python reproducibility/figures/build_figures.py` and `python reproducibility/figures/build_prisma_figure.py` from the source root with Matplotlib, NumPy and Pillow installed. Set `PROVFOLD_FIGURE_OUTPUT` to a separate destination. Figures are generated from the bundled analytical tables; no generative image model is used.

| Figure | Content |
|---|---|
| 1 | Input, operations and output workflow |
| 2 | HCC counting policy and support provenance |
| 3 | Simulation estimators and support trade-offs |
| 4 | Separate public-resource arms and concordance |
| S1 | AH review flow |
| S2 | AH Shannon forest plot |
| S3 | AH threshold and group-omission analysis |

The builders write PDF, SVG, PNG and TIFF. Multi-panel figures compare related analytical quantities, and do not treat different counting units as sequential attrition.
