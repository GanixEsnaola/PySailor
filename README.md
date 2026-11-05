# PySailor

Python implementation of the Sailor diagram, ported from the R package SailoR 1.2. PySailor is developed by Jon Sáenz and Ganix Esnaola. The method is described in (with its 2021 corrigendum):

Sáenz, J., Carreno-Madinabeitia, S., Esnaola, G., González-Rojí, S. J., Ibarra-Berastegi, G., and Ulazia, A.: The Sailor diagram – A new diagram for the verification of two-dimensional vector data from multiple models, *Geosci. Model Dev.*, 13, 3221–3240, [https://doi.org/10.5194/gmd-13-3221-2020](https://doi.org/10.5194/gmd-13-3221-2020), 2020.

## Setup

Two independent options; both create an environment named `pysailor` with Python 3.12.

### Option 1: uv (recommended)

Uses the exact versions pinned in `requirements.txt`. Prerequisite: [uv](https://docs.astral.sh/uv/), a Python package and environment manager. If `uv --version` fails, install it (no root needed; it goes into `~/.local/bin`):

```sh
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Alternatively: `pip install --user uv`.

Then create the virtual environment (uv downloads Python 3.12 if needed) and install the dependencies:

```sh
uv venv --seed --python 3.12 pysailor
uv pip install --python pysailor/bin/python -r requirements.txt
source pysailor/bin/activate
```

To leave the environment: `deactivate`. To use it again later: `source pysailor/bin/activate`.

### Option 2: conda

Needs a recent conda (23.10 or later, with the libmamba solver); older versions can run out of memory while solving. conda-forge does not always have the exact versions pinned in `requirements.txt`, so the packages are given without versions (tested with pandas 3.0.6 and the other versions as pinned):

```sh
conda create -n pysailor -c conda-forge python=3.12 numpy scipy pandas matplotlib pytest
conda activate pysailor
```

To leave the environment: `conda deactivate`. To use it again later: `conda activate pysailor`.

## Usage

```python
import pandas as pd
from SailorIndices import sailorIndices, sailorTable
from SailorPlot import sailorPlot

ref = pd.read_csv("data/Synthetic_ref.csv")   # columns mod, U, V
mod = pd.read_csv("data/Synthetic_mod.csv")   # one or more models, paired row by row with ref
print(sailorTable(sailorIndices(ref, mod), roundDigits=4))
fig, ax = sailorPlot(ref, mod, xlab="Ux (m/s)", ylab="Uy (m/s)", rmseUnits=" m/s")
fig.savefig("sailor.png")
```

- `Sailor.py`: `UVError(U, V)`, all the diagnostics for one model: means, EOFs, σ, total variance, eccentricity, θu, θv, θvu, congruence g11, bias magnitude, RMSE, vector R² and the canonical correlations it is built from. U and V must be paired N×2 arrays (row t is the same time in both).
- `SailorIndices.py`: `sailorIndices` (one `UVError` per model) and `sailorTable` (summary table). With `pooled=True`, ensemble members stacked under one model name are compared with the reference repeated once per member.
- `SailorPlot.py`: `sailorPlot`, the Sailor diagram with matplotlib. Returns `(fig, ax)`; accepts a matplotlib `style` and an existing `ax`. With `refOnce=True` the reference ellipse is drawn once, at the reference mean, instead of once per model.

## Examples and tests

- `examples/`: one script per dataset, reproducing the paper's figures and tables (`python examples/Synthetic.py`). Their tables and figures are in `examples/output/`.
- `tests/`: `pytest` (run from the repository root). Checks the diagnostics' properties, and that every example reproduces `examples/output/`.

## Contents

- `data/`: the SailoR 1.2 datasets as CSV.

## License

`Sailor.py` is licensed under the MIT License; see [LICENSE-MIT](LICENSE-MIT).

`SailorIndices.py`, `SailorPlot.py`, the scripts in `examples/` and the datasets in `data/` are derived from the R package SailoR 1.2 (CRAN, 2020): the code is ported to Python and modified (2026), and the datasets are redistributed as CSV. As SailoR, they are licensed under the GNU General Public License version 3; see [LICENSE](LICENSE).

SailoR 1.2: Copyright (C) 2020 J. Sáenz, S. Carreno-Madinabeitia, S. J. González-Rojí, G. Esnaola, G. Ibarra-Berastegi, A. Ulazia.
