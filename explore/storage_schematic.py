"""Where HSA keeps its data, and what each store is trusted for."""
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BLUE, ORANGE, AQUA, YELLOW, MAG, VIOLET, RED = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6", "#e34948")
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
mpl.rcParams.update({"font.size": 9, "figure.facecolor": "white"})

fig, ax = plt.subplots(figsize=(13.2, 8.2))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, title, lines, col, fill, note=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.4",
                                fc=fill, ec=col, lw=1.8, zorder=2))
    ax.text(x + w / 2, y + h - 4.0, title, ha="center", va="center", fontsize=9.8,
            fontweight="bold", color=INK, zorder=3)
    for i, t in enumerate(lines):
        ax.text(x + w / 2, y + h - 9.0 - i * 3.9, t, ha="center", va="center",
                fontsize=8.2, color=INK2, zorder=3)
    if note:
        ax.text(x + w / 2, y + 2.6, note, ha="center", va="center", fontsize=7.6,
                color=col, zorder=3, style="italic")

def arrow(x1, y1, x2, y2, label="", col=INK2, ls="-", dx=0, dy=1.8, fs=7.9):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14,
                                 lw=1.6, color=col, linestyle=ls, zorder=1,
                                 shrinkA=1, shrinkB=1))
    if label:
        ax.text((x1 + x2) / 2 + dx, (y1 + y2) / 2 + dy, label, ha="center", va="center",
                fontsize=fs, color=col, zorder=4,
                bbox=dict(fc="white", ec="none", pad=1.5))

ax.text(0, 97, "Where HSA stores its data", fontsize=15.5, fontweight="bold", color=INK)
ax.text(0, 92.6, "each arrow is a step the agents run unattended; the box it writes to is what "
                 "that step is trusted for", fontsize=9.2, color=INK2)

# ---------------------------------------------------------------- row 1: the pipeline
Y1, H1 = 68, 19
box(0, Y1, 18, H1, "17 public sources",
    ["GEO · Zenodo · HuBMAP", "figshare · 10x · CELLxGENE"],
    VIOLET, "#f6f5fd", note="scouted weekly")
box(25, Y1, 21, H1, "data/   1.4 TB",
    ["raw deposits, as downloaded", "376 GB of it is .rds we cannot read"],
    ORANGE, "#fdf2ed", note="disposable input")
box(53, Y1, 22, H1, "standardised/   309 GB",
    ["7,333 .h5ad  —  THE ATLAS", "one gzip file per sample"],
    AQUA, "#edfaf5", note="the product")
box(82, Y1, 18, H1, "quarantine/   1.7 GB",
    ["56 files pulled out", "50 mouse · 6 synthetic"],
    RED, "#fdeeee", note="kept, not deleted")

arrow(18, Y1 + 9.5, 25, Y1 + 9.5, "download")
arrow(46, Y1 + 9.5, 53, Y1 + 9.5, "materialise")
arrow(75, Y1 + 9.5, 82, Y1 + 9.5, "audit")
ax.text(49.5, Y1 - 3.4, "raw is deleted only after the .h5ad is written AND re-opened",
        fontsize=8.0, color=ORANGE, ha="center", style="italic")

# ---------------------------------------------------------------- row 2: the record
Y2, H2 = 38, 19
box(0, Y2, 28, H2, "hsa_catalogue.sqlite   44 MB",
    ["sample_meta · datasets · usage", "$31 of LLM annotation"],
    BLUE, "#eef5fd", note="the authoritative metadata")
box(36, Y2, 28, H2, "runs/   99 MB",
    ["sample_manifest.csv (74 col)", "ledger · audit · sex_refresh"],
    BLUE, "#eef5fd", note="every step's receipt")
box(72, Y2, 28, H2, "exports/   26 MB",
    ["hsa_master_sheet.csv", "dashboard JSON"],
    YELLOW, "#fdf7e9", note="local only — never pushed")

arrow(9, Y1, 9, Y2 + H2, "curate (LLM)", BLUE, dx=-7.5, dy=0)
arrow(60, Y1, 52, Y2 + H2, "QC + audit", BLUE, dx=5.5, dy=0)
arrow(28, Y2 + 9.5, 36, Y2 + 9.5, "join", BLUE)
arrow(64, Y2 + 9.5, 72, Y2 + 9.5, "export", YELLOW)

# ---------------------------------------------------------------- row 3: derived + public
Y3, H3 = 12, 16
box(0, Y3, 28, H3, "visium_hd_cells/   44 GB",
    ["bin2cell segmentations"], MAG, "#fdf1f6", note="derived, reproducible")
box(36, Y3, 28, H3, "git repo   (code only)",
    ["hsa/*.py · explore/*.ipynb"], AQUA, "#edfaf5", note="no data, no keys, no sheet")
box(72, Y3, 28, H3, "GitHub Pages   docs/",
    ["public dashboard"], AQUA, "#edfaf5", note="headline counts, no sample rows")

arrow(56, Y1, 30, Y2 + 4, "", MAG, ls="--")
arrow(30, Y2 + 2, 16, Y3 + H3, "segment (Visium HD only)", MAG, ls="--", dx=-3, dy=2.5)
arrow(64, Y3 + 8, 72, Y3 + 8, "publish", AQUA)

# everything above the git box is gitignored
ax.add_patch(FancyBboxPatch((33.2, Y3 - 1.2), 33.6, H3 + 50, boxstyle="round,pad=0,rounding_size=1",
                            fc="none", ec=MUTED, lw=0, zorder=0))
ax.text(50, 5.2, "ALL of data/, standardised/, runs/, exports/, quarantine/, visium_hd_cells/ "
                 "and the .sqlite are gitignored — the repository carries code only",
        ha="center", fontsize=8.4, color=INK2,
        bbox=dict(fc="#faf9f6", ec=MUTED, lw=1, pad=5, boxstyle="round,pad=0.6"))

ax.text(50, 0.6, "if standardised/ were lost it is re-downloadable from the catalogue's URLs;   "
                 "if the .sqlite were lost the curated metadata is NOT recoverable",
        ha="center", fontsize=8.2, color=RED)

fig.savefig("explore/storage_schematic.png", dpi=170, bbox_inches="tight")
print("wrote explore/storage_schematic.png")
