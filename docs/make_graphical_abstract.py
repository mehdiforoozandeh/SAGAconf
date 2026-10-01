"""Draw the README graphical abstract in light and dark variants.

    python docs/make_graphical_abstract.py

Writes docs/graphical_abstract_light.svg and docs/graphical_abstract_dark.svg. The picture
is schematic: two replicate annotations -> SAGAconf -> confident annotation + r-values.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

THEMES = {
    "light": dict(card="#F6F8FA", ink="#1F2328", muted="#59636E", line="#D1D9E0",
                  accent="#0969DA", quie="#D5DAE0", faded="#FFFFFF"),
    "dark": dict(card="#161B22", ink="#E6EDF3", muted="#9198A1", line="#3D444D",
                 accent="#4493F8", quie="#3D444D", faded="#0D1117"),
}
STATE = {"Prom": "#D1495B", "Enha": "#EDAE49", "Tran": "#00798C", "Facu": "#7B5EA7"}

BASE = [(40, "Prom"), (30, "Enha"), (60, "Tran"), (25, "Quie"), (45, "Enha"),
        (20, "Facu"), (50, "Tran"), (50, "Quie")]
VERIF = [(38, "Prom"), (34, "Enha"), (55, "Tran"), (30, "Facu"), (40, "Enha"),
         (25, "Prom"), (48, "Tran"), (50, "Quie")]
IRREPRODUCIBLE = [(130, 155), (200, 220)]  # base-track spans the replicates disagree on
BOUNDARIES = [40, 70, 175, 270]
TRACK_W, TRACK_H = 320, 22


def fill(state, t):
    return t["quie"] if state == "Quie" else STATE[state]


def track(x, y, segs, t, hide=()):
    out, pos = [], 0
    for w, s in segs:
        hidden = any(a <= pos and pos + w <= b for a, b in hide)
        if hidden:
            out.append('<rect x="%g" y="%g" width="%g" height="%g" fill="url(#hatch)" '
                       'stroke="%s" stroke-dasharray="3 2"/>' % (x + pos, y, w, TRACK_H, t["muted"]))
        else:
            out.append('<rect x="%g" y="%g" width="%g" height="%g" fill="%s"/>'
                       % (x + pos, y, w, TRACK_H, fill(s, t)))
        pos += w
    out.append('<rect x="%g" y="%g" width="%g" height="%g" rx="3" fill="none" stroke="%s"/>'
               % (x, y, TRACK_W, TRACK_H, t["line"]))
    return "\n".join(out)


def r_value(p):
    if any(a - 3 <= p <= b + 3 for a, b in IRREPRODUCIBLE):
        return 0.28
    if any(abs(p - b) <= 6 for b in BOUNDARIES):
        return 0.84
    return 0.93 + 0.05 * ((p * 7) % 3) / 3


def svg(t):
    W, H = 1200, 440
    o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" '
         'font-family="-apple-system, BlinkMacSystemFont, \'Segoe UI\', Helvetica, Arial, sans-serif" '
         'role="img" aria-label="SAGAconf compares two replicate chromatin state annotations, '
         'assigns each genomic bin an r-value, and keeps the reproducible calls.">' % (W, H, W, H),
         '<defs><pattern id="hatch" width="6" height="6" patternUnits="userSpaceOnUse" '
         'patternTransform="rotate(45)"><rect width="6" height="6" fill="%s"/>'
         '<line x1="0" y1="0" x2="0" y2="6" stroke="%s" stroke-width="1.5"/></pattern>'
         '<marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" '
         'orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="%s"/></marker></defs>'
         % (t["faded"], t["line"], t["muted"]),
         '<rect x="1" y="1" width="%d" height="%d" rx="16" fill="%s" stroke="%s"/>'
         % (W - 2, H - 2, t["card"], t["line"]),
         '<text x="40" y="62" font-size="34" font-weight="700" fill="%s">SAGAconf</text>' % t["ink"],
         '<text x="40" y="90" font-size="16" fill="%s">Which chromatin state calls can you trust? '
         'Reproducibility across replicates, turned into a calibrated score per genomic bin.</text>'
         % t["muted"]]

    def heading(x, y, n, text):
        o.append('<circle cx="%g" cy="%g" r="11" fill="%s"/>' % (x + 11, y - 5, t["accent"]))
        o.append('<text x="%g" y="%g" font-size="13" font-weight="700" fill="#FFFFFF" '
                 'text-anchor="middle">%d</text>' % (x + 11, y - 0.5, n))
        o.append('<text x="%g" y="%g" font-size="16" font-weight="600" fill="%s">%s</text>'
                 % (x + 30, y, t["ink"], text))

    # Panel 1: two replicate annotations.
    x1, y0 = 40, 150
    heading(x1, y0, 1, "Two replicate annotations")
    o.append('<text x="%g" y="%g" font-size="12" fill="%s">base replicate</text>' % (x1, y0 + 40, t["muted"]))
    o.append(track(x1, y0 + 48, BASE, t))
    o.append('<text x="%g" y="%g" font-size="12" fill="%s">verification replicate</text>' % (x1, y0 + 104, t["muted"]))
    o.append(track(x1, y0 + 112, VERIF, t))
    o.append('<text x="%g" y="%g" font-size="12" fill="%s">from any SAGA method '
             '(ChromHMM, Segway, …): posteriors, K states</text>' % (x1, y0 + 166, t["muted"]))
    lx = x1
    for name in ["Prom", "Enha", "Tran", "Facu", "Quie"]:
        o.append('<rect x="%g" y="%g" width="12" height="12" rx="2" fill="%s"/>' % (lx, y0 + 188, fill(name, t)))
        o.append('<text x="%g" y="%g" font-size="12" fill="%s">%s</text>' % (lx + 17, y0 + 198, t["muted"], name))
        lx += 64

    # Arrow 1 -> 2.
    o.append('<line x1="378" y1="265" x2="438" y2="265" stroke="%s" stroke-width="2" '
             'marker-end="url(#arrow)"/>' % t["muted"])

    # Panel 2: SAGAconf.
    x2 = 455
    o.append('<rect x="%g" y="120" width="300" height="280" rx="12" fill="none" stroke="%s" '
             'stroke-width="2"/>' % (x2, t["accent"]))
    heading(x2 + 20, y0, 2, "SAGAconf")
    # Mini IoU overlap heatmap: strong diagonal = matched states.
    hx, hy, c = x2 + 24, y0 + 30, 16
    iou = [[.9, .05, 0, 0, .05], [.1, .8, .05, 0, 0], [0, .05, .85, .1, 0],
           [0, .1, .45, .4, 0], [.05, 0, 0, 0, .95]]
    for i, row in enumerate(iou):
        for j, v in enumerate(row):
            o.append('<rect x="%g" y="%g" width="%g" height="%g" fill="%s" fill-opacity="%.2f" '
                     'stroke="%s" stroke-width="0.5"/>' % (hx + j * c, hy + i * c, c, c, t["accent"],
                                                           0.08 + 0.92 * v, t["line"]))
    steps = [("Match states", "IoU overlap between replicates"),
             ("Score every bin", "r-value: reproducibility score"),
             ("Keep confident calls", "r-value ≥ α (default 0.8)")]
    for k, (a, b) in enumerate(steps):
        ty = hy + 8 + k * 58
        o.append('<text x="%g" y="%g" font-size="14" font-weight="600" fill="%s">%s</text>'
                 % (hx + 5 * c + 16, ty + 6, t["ink"], a))
        o.append('<text x="%g" y="%g" font-size="12" fill="%s">%s</text>'
                 % (hx + 5 * c + 16, ty + 24, t["muted"], b))
    o.append('<text x="%g" y="%g" font-size="12" fill="%s">+ granularity, misalignment, and '
             'calibration reports</text>' % (hx, 376, t["muted"]))

    # Arrow 2 -> 3.
    o.append('<line x1="772" y1="265" x2="832" y2="265" stroke="%s" stroke-width="2" '
             'marker-end="url(#arrow)"/>' % t["muted"])

    # Panel 3: confident annotation + r-values.
    x3 = 845
    heading(x3, y0, 3, "Confident annotation")
    o.append('<text x="%g" y="%g" font-size="12" fill="%s">r-value per bin</text>' % (x3, y0 + 40, t["muted"]))
    by, bh = y0 + 48, 70
    thr_y = by + bh * (1 - 0.8)
    for k in range(40):
        p = k * 8 + 4
        r = r_value(p)
        col = t["accent"] if r >= 0.8 else t["muted"]
        op = 1.0 if r >= 0.8 else 0.45
        o.append('<rect x="%g" y="%g" width="6" height="%g" rx="1" fill="%s" fill-opacity="%.2f"/>'
                 % (x3 + k * 8 + 1, by + bh * (1 - r), bh * r, col, op))
    o.append('<line x1="%g" y1="%g" x2="%g" y2="%g" stroke="%s" stroke-dasharray="5 4"/>'
             % (x3, thr_y, x3 + TRACK_W, thr_y, t["ink"]))
    o.append('<text x="%g" y="%g" font-size="12" fill="%s" text-anchor="end">α</text>'
             % (x3 - 6, thr_y + 4, t["ink"]))
    o.append('<text x="%g" y="%g" font-size="12" fill="%s">base annotation, reproducible calls only</text>'
             % (x3, y0 + 146, t["muted"]))
    o.append(track(x3, y0 + 154, BASE, t, hide=IRREPRODUCIBLE))
    o.append('<rect x="%g" y="%g" width="12" height="12" rx="2" fill="url(#hatch)" stroke="%s" '
             'stroke-dasharray="3 2"/>' % (x3, y0 + 188, t["muted"]))
    o.append('<text x="%g" y="%g" font-size="12" fill="%s">irreproducible, filtered out</text>'
             % (x3 + 17, y0 + 198, t["muted"]))
    o.append("</svg>")
    return "\n".join(o)


if __name__ == "__main__":
    for name, t in THEMES.items():
        with open(os.path.join(HERE, "graphical_abstract_%s.svg" % name), "w") as fh:
            fh.write(svg(t) + "\n")
