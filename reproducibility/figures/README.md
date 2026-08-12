# Figure reproduction

The four article figures are generated from the registered input tables in `inputs/`. The script performs no manual or generative-image post-processing. Install Matplotlib, NumPy and Pillow, then run:

```bash
python reproducibility/figures/build_figures.py
```

Set `PROVFOLD_FIGURE_OUTPUT` to direct the PDF, SVG, PNG, TIFF and provenance-manifest outputs to another directory.
