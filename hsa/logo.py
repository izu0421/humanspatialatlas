"""Generate the HSA mark: a hexagonal spot lattice (Visium-like capture area) forming a tissue section,
spots coloured by 'region' in the dashboard's H&E palette. Writes logo.svg (standalone) and returns an inline
version that uses CSS variables so it follows the page theme."""
import math

from .config import ROOT

R, PITCH, SPOT = 46.0, 7.0, 2.45          # canvas radius, spot spacing, spot radius


def _tissue(x, y):
    """Signed 'inside' value for a lobed tissue section centred at origin."""
    a = math.atan2(y, x)
    r = 27 + 6 * math.sin(2 * a + 0.9) + 3 * math.cos(3 * a - 0.4)
    return r - math.hypot(x, y)


def spots():
    out = []
    rows = int(2 * R / (PITCH * math.sqrt(3) / 2)) + 2
    for j in range(-rows // 2, rows // 2 + 1):
        y = j * PITCH * math.sqrt(3) / 2
        off = PITCH / 2 if j % 2 else 0
        for i in range(-12, 13):
            x = i * PITCH + off
            if math.hypot(x, y) > R - 4:
                continue
            inside = _tissue(x, y)
            if inside < -1:
                kind = "bg"
            else:
                # a crescent "region" along the upper-left boundary + a small focal cluster
                focal = math.hypot(x - 9, y - 6) < 6.5
                kind = "eosin" if (inside < 8 and x < -4) or focal else "accent"
            out.append((round(x + R, 2), round(y + R, 2), kind))
    return out


def svg(inline=False):
    fill = {"accent": "var(--accent)", "eosin": "var(--eosin)", "bg": "var(--line)"} if inline else \
           {"accent": "#3a3f9a", "eosin": "#c24d7a", "bg": "#d9d9e6"}
    body = "".join(f'<circle cx="{x}" cy="{y}" r="{SPOT if k != "bg" else SPOT * 0.62}" fill="{fill[k]}"/>'
                   for x, y, k in spots())
    size = 2 * R
    attrs = 'class="logo" role="img" aria-label="HSA logo"' if inline else 'xmlns="http://www.w3.org/2000/svg"'
    return f'<svg {attrs} viewBox="0 0 {size:.0f} {size:.0f}" width="{size:.0f}" height="{size:.0f}">{body}</svg>'


if __name__ == "__main__":
    (ROOT / "logo.svg").write_text(svg())
    print(len(spots()), "spots ->", ROOT / "logo.svg")
