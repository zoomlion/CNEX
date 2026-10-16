#!/usr/bin/env python3
"""Panels b and c of the CNEX overview figure, side by side in one file
(183 mm x 82.65 mm, 1240 x 560 units).

b (left)  - one read position, looked up: the 13-mer at x = 4 gives three
            (element, locus) assignments; the arrow points at the three
            corresponding dots in c (circled).
c (right) - the complex window as a dot plot: element 9898 has 12 matches at
            offset 213 and a 25-match staircase at offsets 291/293/295; one
            of the 36 links fails the 2-nt tolerance; 55 other candidates
            scatter as light dots.

Real data: docs/figures/validate_example.json (extract_validate_example.py).
"""

import json
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = json.load(open(os.path.join(HERE, "validate_example.json")))

INK = "#1f2937"
MUTED = "#6b7280"
FAINT = "#9ca3af"
AXIS = "#4b5563"
RULE = "#cbd5e1"
CAPSULE = "#f3f4f6"
GRAY_DOT = "#b9c0c9"
GRID = "#eef1f4"
FRAME = "#d7dbe0"
ORANGE = "#dd6b20"
ORANGE_L = "#f6ad55"
SANS = "Arial, Liberation Sans, Helvetica, sans-serif"
MONO = "DejaVu Sans Mono, monospace"

WIN_ELE = DATA["emitted"]["element"]
pts = sorted(DATA["chain_points"])
others = [(x, l) for x, e, l in DATA["matches"] if e != WIN_ELE]
n_cand = len(set(e for _, e, _ in DATA["matches"]))
lookup = DATA["lookup_example"]
align = DATA["alignment"]
ax0 = align["x_range"][0]
read_letters = align["read_letters"]
passing = DATA["rc" if DATA["emitted"]["strand"] == -1 else "fwd"]["passing"]

offsets = defaultdict(list)
for x, l in pts:
    offsets[l - x].append(x)
keys = sorted(offsets)
o1, n1 = keys[0], len(offsets[keys[0]])
o2a, o2b = keys[1], keys[-1]
n2 = sum(len(offsets[k]) for k in keys[1:])

fails = [(p, q) for p, q in zip(pts, pts[1:])
         if abs((q[0] - p[0]) - (q[1] - p[1])) > 2]
assert len(fails) == 1
fp, fq = fails[0]

# ── panel b (left column) ──────────────────────────────────────────────────
L0, ADV, FS_LET = 62.6, 9.34, 15.5
Y_BASE = 148
N_LET = 25
MER_X = lookup["x"]

letters = "".join(
    f'<text x="{L0 + i * ADV:.1f}" y="{Y_BASE}">{read_letters[i]}</text>'
    for i in range(N_LET))
box_x = L0 + (MER_X - ax0) * ADV - 2
box_w = 13 * ADV + 4
mer_cx = box_x + box_w / 2

assign_txt = ""
for i, (e, l) in enumerate(lookup["assignments"]):
    if i:
        assign_txt += f'<tspan fill="{MUTED}"> · </tspan>'
    col = ORANGE if e == WIN_ELE else MUTED
    assign_txt += f'<tspan fill="{col}">({e},{l})</tspan>'

# ── panel c (right column) ─────────────────────────────────────────────────
CX0, CX1 = 655.0, 1220.0
CY0, CY1 = 470.0, 90.0
X_MAX, Y_MAX = 150, 900


def px(x):
    return CX0 + x * (CX1 - CX0) / X_MAX


def py(l):
    return CY0 - l * (CY0 - CY1) / Y_MAX


grid = []
for gl in range(100, Y_MAX, 100):
    grid.append(f'<line x1="{CX0}" y1="{py(gl):.1f}" x2="{CX1}" '
                f'y2="{py(gl):.1f}" stroke="{GRID}" stroke-width="0.6"/>')
for gx in range(25, X_MAX, 25):
    grid.append(f'<line x1="{px(gx):.1f}" y1="{CY1}" x2="{px(gx):.1f}" '
                f'y2="{CY0}" stroke="{GRID}" stroke-width="0.6"/>')
grid = "\n    ".join(grid)

xticks = []
for t in range(0, X_MAX + 1, 25):
    h = 6 if t % 50 == 0 else 3.5
    xticks.append(f'<line x1="{px(t):.1f}" y1="{CY0}" x2="{px(t):.1f}" '
                  f'y2="{CY0 + h}" stroke="{AXIS}" stroke-width="0.9"/>')
    if t % 50 == 0:
        xticks.append(f'<text x="{px(t):.1f}" y="{CY0 + 22}" font-size="10" '
                      f'fill="{FAINT}" text-anchor="middle">{t}</text>')
yticks = []
for t in range(0, Y_MAX + 1, 200):
    yticks.append(f'<line x1="{CX0 - 6}" y1="{py(t):.1f}" x2="{CX0}" '
                  f'y2="{py(t):.1f}" stroke="{AXIS}" stroke-width="0.9"/>'
                  f'<text x="{CX0 - 10}" y="{py(t) + 3.5:.1f}" font-size="10" '
                  f'fill="{FAINT}" text-anchor="end">{t}</text>')
xticks = "\n    ".join(xticks)
yticks = "\n    ".join(yticks)

guide = "\n    ".join([
    f'<line x1="{px(0):.1f}" y1="{py(o1):.1f}" x2="{px(24):.1f}" '
    f'y2="{py(o1 + 24):.1f}" stroke="{ORANGE_L}" stroke-width="1" '
    f'stroke-dasharray="4 3"/>',
    f'<line x1="{px(105):.1f}" y1="{py(105 + o2a):.1f}" '
    f'x2="{px(145):.1f}" y2="{py(145 + o2a):.1f}" stroke="{ORANGE_L}" '
    f'stroke-width="1" stroke-dasharray="4 3"/>',
    f'<line x1="{px(MER_X):.1f}" y1="{py(265):.1f}" '
    f'x2="{px(MER_X):.1f}" y2="{py(185):.1f}" stroke="{ORANGE}" '
    f'stroke-width="1" stroke-dasharray="3 3" stroke-opacity="0.7"/>',
])
fail_line = (f'<line x1="{px(fp[0]):.1f}" y1="{py(fp[1]):.1f}" '
             f'x2="{px(fq[0]):.1f}" y2="{py(fq[1]):.1f}" stroke="{RULE}" '
             f'stroke-width="1" stroke-dasharray="3 3"/>')

other_dots = "\n      ".join(
    f'<circle cx="{px(x):.1f}" cy="{py(l):.1f}" r="1.5" fill="{GRAY_DOT}"/>'
    for x, l in others)

chain_segs = []
for (x1, l1), (x2, l2) in zip(pts, pts[1:]):
    if abs((x2 - x1) - (l2 - l1)) <= 2:
        chain_segs.append(
            f'<line x1="{px(x1):.1f}" y1="{py(l1):.1f}" x2="{px(x2):.1f}" '
            f'y2="{py(l2):.1f}" stroke="{ORANGE}" stroke-width="1.4" '
            f'stroke-opacity="0.75"/>')
chain_segs = "\n    ".join(chain_segs)

win_dots = "\n      ".join(
    f'<circle cx="{px(x):.1f}" cy="{py(l):.1f}" r="2.4" fill="{ORANGE}"/>'
    for x, l in pts)

SVG = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1240 560" width="183mm" height="82.65mm">
  <rect width="1240" height="560" fill="#ffffff"/>
  <g font-family="{SANS}" fill="{INK}">

    <!-- ════ panel b (left) ════ -->
    <text x="22" y="43" font-size="24" font-weight="bold">b</text>
    <text x="52" y="42" font-size="16.5" font-weight="bold">One position → three candidates</text>

    <text x="30" y="112" font-size="11.5" font-weight="bold" fill="{MUTED}">read (RC)</text>
    <rect x="52" y="125" width="252" height="33" rx="16.5" fill="{CAPSULE}"/>
    <rect x="{box_x:.1f}" y="127" width="{box_w:.1f}" height="30" rx="2"
          fill="#fee8d8" stroke="{ORANGE_L}" stroke-width="1.2"/>
    <g font-family="{MONO}" font-size="{FS_LET}" fill="{INK}">
    {letters}
    </g>
    <text x="{L0 + N_LET * ADV + 2:.1f}" y="{Y_BASE}" font-size="14" fill="{FAINT}">…</text>

    <path d="M{box_x:.1f},163 V169 M{box_x + box_w:.1f},163 V169 M{box_x:.1f},166 H{box_x + box_w:.1f}"
          fill="none" stroke="{MUTED}" stroke-width="1"/>
    <text x="{mer_cx:.1f}" y="184" font-size="10.5" fill="{ORANGE}"
          text-anchor="middle">13-mer at x = {MER_X}</text>

    <line x1="{mer_cx:.1f}" y1="192" x2="{mer_cx:.1f}" y2="340"
          stroke="{MUTED}" stroke-width="1.2"/>
    <polygon points="{mer_cx - 4:.1f},340 {mer_cx + 4:.1f},340 {mer_cx:.1f},348" fill="{MUTED}"/>

    <line x1="30" y1="348" x2="534" y2="348" stroke="{RULE}" stroke-width="1"/>
    <line x1="30" y1="396" x2="534" y2="396" stroke="{RULE}" stroke-width="1"/>
    <text x="534" y="340" font-size="10" fill="{FAINT}" text-anchor="end">13-mer table</text>
    <text x="98" y="377" font-size="15" font-family="{MONO}">{lookup["mer"]}</text>
    <text x="228" y="376" font-size="14" fill="{MUTED}">→</text>
    <text x="252" y="377" font-size="12.5" font-family="{MONO}">{assign_txt}</text>

    <line x1="540" y1="375" x2="644" y2="375" stroke="{MUTED}" stroke-width="1.2"/>
    <polygon points="644,371 644,379 652,375" fill="{MUTED}"/>

    <!-- ════ panel c (right) ════ -->
    <text x="602" y="43" font-size="24" font-weight="bold">c</text>
    <text x="632" y="42" font-size="16.5" font-weight="bold">Two offset groups, one chain</text>

    <g>
    {grid}
    </g>
    <line x1="{CX0}" y1="{CY0}" x2="{CX1}" y2="{CY0}" stroke="{AXIS}" stroke-width="0.9"/>
    <line x1="{CX0}" y1="{CY0}" x2="{CX0}" y2="{CY1}" stroke="{AXIS}" stroke-width="0.9"/>
    <line x1="{CX0}" y1="{CY1}" x2="{CX1}" y2="{CY1}" stroke="{FRAME}" stroke-width="0.8"/>
    <line x1="{CX1}" y1="{CY1}" x2="{CX1}" y2="{CY0}" stroke="{FRAME}" stroke-width="0.8"/>
    {xticks}
    {yticks}
    <text x="{(CX0 + CX1) / 2:.1f}" y="{CY0 + 42}" font-size="11.5" fill="{MUTED}"
          text-anchor="middle">read position</text>
    <text x="612" y="{(CY0 + CY1) / 2:.1f}" font-size="11.5" fill="{MUTED}"
          text-anchor="middle" transform="rotate(-90 612 {(CY0 + CY1) / 2:.1f})">element position</text>

    <g>
    {guide}
    {fail_line}
    </g>
    <g>
    {chain_segs}
    </g>
    <g>
      {other_dots}
    </g>
    <g>
      {win_dots}
    </g>

    <circle cx="663" cy="106" r="2.6" fill="{ORANGE}"/>
    <text x="673" y="109.5" font-size="10.5" fill="{MUTED}">element {WIN_ELE}</text>
    <circle cx="663" cy="123" r="1.8" fill="{GRAY_DOT}"/>
    <text x="673" y="126.5" font-size="10.5" fill="{MUTED}">{n_cand - 1} other candidates</text>

    <text x="670" y="414" font-size="10.5" fill="{ORANGE}" text-anchor="middle">x = {MER_X}</text>
    <text x="697" y="358" font-size="10.5" fill="{FAINT}">offset {o1} ({n1} matches)</text>
    <text x="1095" y="264" font-size="10.5" fill="{FAINT}" text-anchor="middle">offsets {o2a}→{o2b} ({n2} matches)</text>
    <text x="898" y="350" font-size="10.5" fill="{FAINT}" text-anchor="middle">{len(fails)} failed link</text>
  </g>
</svg>
'''

out = os.path.join(HERE, "validate_principle.svg")
with open(out, "w") as fh:
    fh.write(SVG)
print("wrote", out, f"| b: x={MER_X} assigns={len(lookup['assignments'])}"
      f" | c: pts={len(pts)} offsets={dict(offsets)} fail={fails} passing={passing}")
