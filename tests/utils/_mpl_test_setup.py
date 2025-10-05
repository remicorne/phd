import matplotlib

# Use a headless, stable rasterizer everywhere
matplotlib.use("agg")

import matplotlib as mpl

# Pin runtime config params to avoid cross-OS drift (fonts/AA/DPI/etc.)
mpl.rcParams.update(
    {
        "figure.figsize": (8, 5),
        "figure.dpi": 120,
        "savefig.dpi": 120,
        "savefig.facecolor": "white",
        "savefig.transparent": False,
        "savefig.bbox": None,
        "font.family": "DejaVu Sans",
        "text.antialiased": False,
        "path.simplify": False,
    }
)

# Stabilise seaborn defaults
try:
    import seaborn as sns

    sns.set_theme(
        rc={
            "figure.figsize": (8, 5),
            "figure.dpi": 120,
        }
    )
except Exception:
    pass

# Avoid randomness
try:
    import numpy as np

    np.random.seed(0)
except Exception:
    pass
