"""HSA mark: 19 capture spots in a hexagon (a Visium-like array in miniature) - indigo tissue, a pink
region on one edge, two empty spots. Writes logo.svg; svg(inline=True) uses the page's CSS colour tokens."""
from .config import ROOT

# axial hex coordinates (q, r) within radius 2, with a region label
EOSIN = {(2, -2), (2, -1), (1, -2)}
EMPTY = {(-2, 2), (-1, 2)}
SPACING, SPOT, PAD = 22.0, 9.0, 12.0


def spots():
    out = []
    for q in range(-2, 3):
        for r in range(-2, 3):
            if abs(q + r) > 2:
                continue
            x = SPACING * (q + r / 2)
            y = SPACING * (r * 3 ** 0.5 / 2)
            kind = "eosin" if (q, r) in EOSIN else "bg" if (q, r) in EMPTY else "accent"
            out.append((x, y, kind))
    return out


def svg(inline=False):
    fill = {"accent": "var(--accent)", "eosin": "var(--eosin)", "bg": "var(--line)"} if inline else \
           {"accent": "#3a3f9a", "eosin": "#c24d7a", "bg": "#d9d9e6"}
    half = 2 * SPACING + SPOT + PAD
    body = "".join(f'<circle cx="{x + half:.1f}" cy="{y + half:.1f}" r="{SPOT}" fill="{fill[k]}"/>'
                   for x, y, k in spots())
    size = 2 * half
    attrs = 'class="logo" role="img" aria-label="HSA logo"' if inline else 'xmlns="http://www.w3.org/2000/svg"'
    return f'<svg {attrs} viewBox="0 0 {size:.0f} {size:.0f}" width="{size:.0f}" height="{size:.0f}">{body}</svg>'


if __name__ == "__main__":
    (ROOT / "logo.svg").write_text(svg())
    print(len(spots()), "spots ->", ROOT / "logo.svg")
