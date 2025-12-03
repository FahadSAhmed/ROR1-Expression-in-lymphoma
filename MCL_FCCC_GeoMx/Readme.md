# ROR1 Expression in Mantle Cell Lymphoma (GeoMx DSP, CD20 Compartment)

This repository documents a Python/Jupyter workflow for visualizing **ROR1 RNA expression** from **CD20⁺ regions of interest (ROIs)** in a **mantle cell lymphoma (MCL)** cohort profiled on the **NanoString GeoMx DSP** platform.

The notebook reads a normalized RNA Excel file exported from GeoMx, filters for a given diagnosis and compartment, extracts the ROR1 row, and produces a **sorted bar plot** of ROR1 expression across ROIs.

---

## Overview

Given a GeoMx-normalized RNA Excel file with at least:

- Sheet: `SegmentProperties`
- Sheet: `TargetCountMatrix`

the code:

1. Filters segments by:
   - `ROIComments` == desired diagnosis (e.g., `"MCL"`)
   - `SegmentLabel` == desired compartment (e.g., `"CD20"`)
2. Extracts the **ROR1** row from `TargetCountMatrix` using `TargetName`.
3. Selects the columns corresponding to the selected segments (`SegmentDisplayName`).
4. Sorts ROR1 expression **in descending order**.
5. Generates a **red bar chart** with one bar per segment and saves it as a PNG file.

This is intended for quick QC/visualization of ROR1 expression in a specific compartment (e.g., CD20⁺ B-cell mask) across MCL ROIs.

---

## Input Requirements

### Excel file structure

**Sheet: `SegmentProperties`**

Must include at least:

- `ROIComments` – diagnosis label (e.g., `MCL`, `IF`, `HL`).
- `SegmentLabel` – compartment/mask (e.g., `CD20`, `CD3`, `CD30`).
- `SegmentDisplayName` – unique segment identifier; must match column names in `TargetCountMatrix`.

**Sheet: `TargetCountMatrix`**

Must include at least:

- `TargetName` – target (gene/protein) identifier; must include `ROR1`.
- One column **per segment**, with names matching `SegmentDisplayName` from `SegmentProperties`.

The notebook raises clear errors if required columns or the ROR1 row are missing.

---

## Dependencies

- Python ≥ 3.8
- `pandas`
- `numpy`
- `matplotlib`

Install:

```bash
pip install pandas numpy matplotlib
```

---

## Configuration

At the top of the notebook/script, define:

```python
FILE_PATH   = "/path/to/Normalized RNA HK.xlsx"
DIAGNOSIS   = "MCL"
COMPARTMENT = "CD20"   # e.g. CD20, CD3, CD30
```

Edit:

- `FILE_PATH` → path to your GeoMx-normalized RNA Excel file.
- `DIAGNOSIS` → label in `ROIComments`.
- `COMPARTMENT` → label in `SegmentLabel`.

---

## Usage

### 1. In Jupyter / Google Colab

1. Open `ROR1_MCL_GeoMx_FCCC.ipynb`.
2. Update `FILE_PATH`, `DIAGNOSIS`, and `COMPARTMENT` as needed.
3. Run all cells.

Outputs:

- PNG bar plot:
  - `ROR1_<DIAGNOSIS>_<COMPARTMENT>_bar_sorted_desc_red.png`
- Printed `pandas.Series` of ROR1 values (indexed by segment), sorted high → low.

### 2. As a standalone Python function (optional)

You can wrap the logic as:

```python
from ror1_plot import plot_ror1_by_diag_compartment

out_png = f"ROR1_{DIAGNOSIS}_{COMPARTMENT}_bar_sorted_desc_red.png"
ror1_series = plot_ror1_by_diag_compartment(
    FILE_PATH,
    DIAGNOSIS,
    COMPARTMENT,
    out_png=out_png,
    show_plot=False,
)
print("Saved plot to:", out_png)
print(ror1_series)
```

Run:

```bash
python ror1_plot.py
```

---

## Output Description

- **Figure (PNG)**  
  - Bar chart of normalized ROR1 counts per segment.
  - Bars sorted in **descending** order.
  - X-axis: segment names (derived from `SegmentDisplayName`).
  - Y-axis: normalized ROR1 expression.

- **Tabular / console output**  
  - Sorted ROR1 expression values for all included segments, as a `pandas.Series`.

---

## Adapting to Other Targets or Cohorts

- To analyze a different diagnosis: change `DIAGNOSIS`.
- To analyze a different compartment: change `COMPARTMENT`.
- To analyze another marker instead of ROR1, modify the filter on `TargetName`, e.g.:

```python
gene_symbol = "ROR1"
row = tcm[tcm["TargetName"].astype(str).str.upper() == gene_symbol.upper()]
```

Replace `"ROR1"` with your target of interest.

---

## Notes

- The workflow assumes de-identified data with no PHI.
- This notebook is intended for exploratory analysis and visualization; it is **not** a clinical diagnostic tool.
 
