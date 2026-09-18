"""
Export_spatialdata_object_summary.py

Renders the SpatialData object's repr as a supplementary methods figure
(not tied to a single manuscript Figure N).
"""

import os
import matplotlib.pyplot as plt
import spatialdata as sd

# Set to your own cluster storage path
CLUSTER_BASE_DIR = os.environ.get("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")

supp_plot_outdir = "../../../../Manuscript/Figures/Supp/Methods"
os.makedirs(supp_plot_outdir, exist_ok=True)


def save_text_figure(text, out_path, fontsize=11, dpi=300):
    """Render a block of text (e.g. an object's repr) as a monospace-text figure."""
    lines = text.splitlines()
    fig_width = 0.11 * max(len(line) for line in lines) + 0.4
    fig_height = 0.22 * len(lines) + 0.4

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.axis("off")
    ax.text(0, 1, text, family="monospace", fontsize=fontsize,
            va="top", ha="left", transform=ax.transAxes)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# SpatialData object summary ------------------------------------------------

sdata = sd.SpatialData().read(
    f"{CLUSTER_BASE_DIR}/Results/Upstream/SpatialData/Final_SpatialData_obj_April_22_2025"
)

save_text_figure(repr(sdata), os.path.join(supp_plot_outdir, "SpatialData_summary.png"))
