#!/usr/bin/env python3
"""Generate Figure 2: Toroidal Electron SVG with FEA mesh, swapped fields, cross-section, photon path.
Run this script to regenerate the SVG file."""
import math
import os

# === PARAMETERS ===
R = 160        # major radius
r_minor = 55   # minor radius
CX, CY = 450, 340   # center in SVG coords
tilt_deg = 25
tilt = math.radians(tilt_deg)
SIN_T = math.sin(tilt)
COS_T = math.cos(tilt)

N_TOR = 24    # toroidal divisions
N_POL = 12    # poloidal divisions

def torus_3d(theta, phi):
    x3 = (R + r_minor * math.cos(phi)) * math.cos(theta)
    y3 = (R + r_minor * math.cos(phi)) * math.sin(theta)
    z3 = r_minor * math.sin(phi)
    return x3, y3, z3

def project(x3, y3, z3):
    sx = CX + x3
    sy = CY + y3 * SIN_T - z3 * COS_T
    return sx, sy

def torus_2d(theta, phi):
    return project(*torus_3d(theta, phi))

def is_front(theta, phi):
    nx = math.cos(phi) * math.cos(theta)
    ny = math.cos(phi) * math.sin(theta)
    nz = math.sin(phi)
    dot = -ny * SIN_T + nz * COS_T
    return dot > 0

def F(v):
    return f"{v:.1f}"

lines = []
def L(s):
    lines.append(s)

L('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 800" width="1000" height="800">')
L('  <defs>')
L('    <marker id="arrowB" markerWidth="10" markerHeight="7" refX="9" refY="3.5" orient="auto" markerUnits="strokeWidth">')
L('      <polygon points="0 0, 10 3.5, 0 7" fill="#2a9d8f"/>')
L('    </marker>')
L('    <marker id="arrowE" markerWidth="10" markerHeight="7" refX="9" refY="3.5" orient="auto" markerUnits="strokeWidth">')
L('      <polygon points="0 0, 10 3.5, 0 7" fill="#e76f51"/>')
L('    </marker>')
L('    <marker id="arrowDim" markerWidth="8" markerHeight="6" refX="7" refY="3" orient="auto" markerUnits="strokeWidth">')
L('      <polygon points="0 0, 8 3, 0 6" fill="#e9c46a"/>')
L('    </marker>')
L('    <marker id="arrowDimRev" markerWidth="8" markerHeight="6" refX="1" refY="3" orient="auto" markerUnits="strokeWidth">')
L('      <polygon points="8 0, 0 3, 8 6" fill="#e9c46a"/>')
L('    </marker>')
L('    <filter id="glowB" x="-20%" y="-20%" width="140%" height="140%">')
L('      <feGaussianBlur stdDeviation="2" result="blur"/>')
L('      <feMerge><feMergeNode in="blur"/><feMergeNode in="SourceGraphic"/></feMerge>')
L('    </filter>')
L('    <filter id="glowE" x="-20%" y="-20%" width="140%" height="140%">')
L('      <feGaussianBlur stdDeviation="2" result="blur"/>')
L('      <feMerge><feMergeNode in="blur"/><feMergeNode in="SourceGraphic"/></feMerge>')
L('    </filter>')
L('    <filter id="glowGold" x="-20%" y="-20%" width="140%" height="140%">')
L('      <feGaussianBlur stdDeviation="3" result="blur"/>')
L('      <feMerge><feMergeNode in="blur"/><feMergeNode in="SourceGraphic"/></feMerge>')
L('    </filter>')
L('  </defs>')

# Background
L('  <rect width="1000" height="800" fill="#181a20"/>')

# Title
L('  <text x="500" y="38" font-family="Georgia, serif" font-size="20" fill="#e0e0e0" text-anchor="middle" font-weight="bold">Toroidal Electron Model</text>')
L('  <text x="500" y="60" font-family="Georgia, serif" font-size="13" fill="#999" text-anchor="middle" font-style="italic">Electromagnetic energy confined in a torus with linked E and B fields</text>')

# ============================================================
# TORUS FILLED SURFACE (front-facing cells)
# ============================================================
L('  <!-- Torus filled surface (front-facing cells) -->')
L('  <g id="torus-fill" opacity="0.4">')

for i in range(N_TOR):
    theta0 = 2 * math.pi * i / N_TOR
    theta1 = 2 * math.pi * (i + 1) / N_TOR
    for j in range(N_POL):
        phi0 = 2 * math.pi * j / N_POL
        phi1 = 2 * math.pi * (j + 1) / N_POL
        tc = (theta0 + theta1) / 2
        pc = (phi0 + phi1) / 2
        if not is_front(tc, pc):
            continue
        p00 = torus_2d(theta0, phi0)
        p01 = torus_2d(theta0, phi1)
        p10 = torus_2d(theta1, phi0)
        p11 = torus_2d(theta1, phi1)
        d = f"M {F(p00[0])},{F(p00[1])} L {F(p01[0])},{F(p01[1])} L {F(p11[0])},{F(p11[1])} L {F(p10[0])},{F(p10[1])} Z"
        L(f'    <path d="{d}" fill="#2a1e2e" stroke="none"/>')

L('  </g>')

# ============================================================
# FEA MESH — back lines (dashed, low opacity)
# ============================================================
L('  <!-- FEA Mesh - back lines -->')
L('  <g id="mesh-back">')

def draw_curve_segments(pts_with_front, front_sw, front_op, back_sw, back_op):
    """Given list of (sx, sy, is_front), draw front and back segments."""
    seg = []
    cur_front = None
    for sx, sy, f in pts_with_front:
        if cur_front is None:
            cur_front = f
            seg = [(sx, sy)]
        elif f == cur_front:
            seg.append((sx, sy))
        else:
            if len(seg) > 1:
                d = "M " + " L ".join(f"{F(x)},{F(y)}" for x, y in seg)
                if cur_front:
                    L(f'    <path d="{d}" fill="none" stroke="#a855f7" stroke-width="{front_sw}" opacity="{front_op}"/>')
                else:
                    L(f'    <path d="{d}" fill="none" stroke="#a855f7" stroke-width="{back_sw}" opacity="{back_op}" stroke-dasharray="3,3"/>')
            cur_front = f
            seg = [(sx, sy)]
    if len(seg) > 1:
        d = "M " + " L ".join(f"{F(x)},{F(y)}" for x, y in seg)
        if cur_front:
            L(f'    <path d="{d}" fill="none" stroke="#a855f7" stroke-width="{front_sw}" opacity="{front_op}"/>')
        else:
            L(f'    <path d="{d}" fill="none" stroke="#a855f7" stroke-width="{back_sw}" opacity="{back_op}" stroke-dasharray="3,3"/>')

# Toroidal lines (constant phi)
for j in range(N_POL):
    phi = 2 * math.pi * j / N_POL
    pts = []
    for i in range(N_TOR * 4 + 1):
        theta = 2 * math.pi * i / (N_TOR * 4)
        sx, sy = torus_2d(theta, phi)
        front = is_front(theta, phi)
        pts.append((sx, sy, front))
    draw_curve_segments(pts, "0.7", "0.45", "0.5", "0.18")

# Poloidal lines (constant theta)
for i in range(N_TOR):
    theta = 2 * math.pi * i / N_TOR
    pts = []
    for j in range(N_POL * 4 + 1):
        phi = 2 * math.pi * j / (N_POL * 4)
        sx, sy = torus_2d(theta, phi)
        front = is_front(theta, phi)
        pts.append((sx, sy, front))
    draw_curve_segments(pts, "0.7", "0.45", "0.5", "0.18")

L('  </g>')

# ============================================================
# B-FIELD (POLOIDAL) — small loops around tube cross-section
# ============================================================
L('  <!-- B-field: Poloidal loops -->')
L('  <g id="B-field" filter="url(#glowB)">')

b_thetas = [0, math.pi/3, 2*math.pi/3, math.pi, 4*math.pi/3, 5*math.pi/3]
b_loop_frac = 0.72

for bt in b_thetas:
    cx3 = R * math.cos(bt)
    cy3 = R * math.sin(bt)
    center_depth = cy3 * SIN_T
    is_back = center_depth > 40

    loop_pts = []
    n_pts = 48
    lr = r_minor * b_loop_frac
    for k in range(n_pts + 1):
        phi = 2 * math.pi * k / n_pts
        x3 = (R + lr * math.cos(phi)) * math.cos(bt)
        y3 = (R + lr * math.cos(phi)) * math.sin(bt)
        z3 = lr * math.sin(phi)
        sx, sy = project(x3, y3, z3)
        loop_pts.append((sx, sy))

    if is_back:
        opacity = "0.2"
        sw = "1.5"
        dash = ' stroke-dasharray="4,3"'
    else:
        opacity = "0.85"
        sw = "2.2"
        dash = ""

    d = "M " + " L ".join(f"{F(x)},{F(y)}" for x, y in loop_pts)
    L(f'    <path d="{d}" fill="none" stroke="#2a9d8f" stroke-width="{sw}" opacity="{opacity}"{dash}/>')

    if not is_back:
        ai = n_pts // 8
        ax, ay = loop_pts[ai]
        ax2, ay2 = loop_pts[ai + 2]
        angle = math.atan2(ay2 - ay, ax2 - ax)
        adx = 7 * math.cos(angle)
        ady = 7 * math.sin(angle)
        px = -ady * 0.5
        py = adx * 0.5
        L(f'    <polygon points="{F(ax+adx)},{F(ay+ady)} {F(ax+px)},{F(ay+py)} {F(ax-px)},{F(ay-py)}" fill="#2a9d8f" opacity="0.9"/>')

L('  </g>')

# B-field label
bx, by = torus_2d(0, math.pi/4)
L(f'  <text x="{F(bx+55)}" y="{F(by-15)}" font-family="Georgia, serif" font-size="15" fill="#2a9d8f" font-weight="bold" font-style="italic">B</text>')
L(f'  <text x="{F(bx+55)}" y="{F(by+2)}" font-family="Georgia, serif" font-size="12" fill="#2a9d8f">(poloidal)</text>')
L(f'  <line x1="{F(bx+53)}" y1="{F(by-8)}" x2="{F(bx+15)}" y2="{F(by-5)}" stroke="#2a9d8f" stroke-width="0.8" opacity="0.5" stroke-dasharray="3,2"/>')

# ============================================================
# E-FIELD (TOROIDAL) — flowing along the tube
# ============================================================
L('  <!-- E-field: Toroidal arrows -->')
L('  <g id="E-field" filter="url(#glowE)">')

e_phis = [0, math.pi/2, math.pi, 3*math.pi/2]
n_seg = 8

for ep in e_phis:
    seg_len = 2 * math.pi / n_seg
    for s in range(n_seg):
        ts = s * seg_len + seg_len * 0.1
        te = (s + 1) * seg_len - seg_len * 0.1
        seg_pts = []
        any_visible = False
        n_sub = 20
        for k in range(n_sub + 1):
            theta = ts + (te - ts) * k / n_sub
            ef = 0.85
            x3 = (R + r_minor * ef * math.cos(ep)) * math.cos(theta)
            y3 = (R + r_minor * ef * math.cos(ep)) * math.sin(theta)
            z3 = r_minor * ef * math.sin(ep)
            sx, sy = project(x3, y3, z3)
            seg_pts.append((sx, sy))
            if is_front(theta, ep):
                any_visible = True

        if not any_visible:
            d = "M " + " L ".join(f"{F(x)},{F(y)}" for x, y in seg_pts)
            L(f'    <path d="{d}" fill="none" stroke="#e76f51" stroke-width="1.2" opacity="0.15" stroke-dasharray="3,3"/>')
        else:
            d = "M " + " L ".join(f"{F(x)},{F(y)}" for x, y in seg_pts)
            L(f'    <path d="{d}" fill="none" stroke="#e76f51" stroke-width="2" opacity="0.75"/>')
            ex1, ey1 = seg_pts[-2]
            ex2, ey2 = seg_pts[-1]
            angle = math.atan2(ey2 - ey1, ex2 - ex1)
            adx = 7 * math.cos(angle)
            ady = 7 * math.sin(angle)
            px = -ady * 0.5
            py = adx * 0.5
            L(f'    <polygon points="{F(ex2+adx)},{F(ey2+ady)} {F(ex2+px)},{F(ey2+py)} {F(ex2-px)},{F(ey2-py)}" fill="#e76f51" opacity="0.8"/>')

L('  </g>')

# E-field label
elx, ely = torus_2d(math.pi * 0.6, 0)
L(f'  <text x="{F(elx-90)}" y="{F(ely-30)}" font-family="Georgia, serif" font-size="15" fill="#e76f51" font-weight="bold" font-style="italic">E</text>')
L(f'  <text x="{F(elx-90)}" y="{F(ely-13)}" font-family="Georgia, serif" font-size="12" fill="#e76f51">(toroidal)</text>')
L(f'  <line x1="{F(elx-70)}" y1="{F(ely-23)}" x2="{F(elx-30)}" y2="{F(ely-15)}" stroke="#e76f51" stroke-width="0.8" opacity="0.5" stroke-dasharray="3,2"/>')

# ============================================================
# PHOTON HELICAL PATH
# ============================================================
L('  <!-- Photon helical path -->')
L('  <g id="photon-path">')

winding = 3.5
n_photon = 500
all_pts = []

for k in range(n_photon + 1):
    theta = math.pi * k / n_photon
    phi = winding * theta * 2
    x3 = (R + r_minor * math.cos(phi)) * math.cos(theta)
    y3 = (R + r_minor * math.cos(phi)) * math.sin(theta)
    z3 = r_minor * math.sin(phi)
    sx, sy = project(x3, y3, z3)
    front = is_front(theta, phi)
    all_pts.append((sx, sy, front))

seg = []
cur_front = None
for sx, sy, f in all_pts:
    if cur_front is None:
        cur_front = f
        seg = [(sx, sy)]
    elif f == cur_front:
        seg.append((sx, sy))
    else:
        if len(seg) > 1:
            d = "M " + " L ".join(f"{F(x)},{F(y)}" for x, y in seg)
            if cur_front:
                L(f'    <path d="{d}" fill="none" stroke="#e9c46a" stroke-width="2.5" opacity="0.9" stroke-dasharray="8,4" filter="url(#glowGold)"/>')
            else:
                L(f'    <path d="{d}" fill="none" stroke="#e9c46a" stroke-width="1.5" opacity="0.3" stroke-dasharray="5,4"/>')
        cur_front = f
        seg = [(sx, sy)]
if len(seg) > 1:
    d = "M " + " L ".join(f"{F(x)},{F(y)}" for x, y in seg)
    if cur_front:
        L(f'    <path d="{d}" fill="none" stroke="#e9c46a" stroke-width="2.5" opacity="0.9" stroke-dasharray="8,4" filter="url(#glowGold)"/>')
    else:
        L(f'    <path d="{d}" fill="none" stroke="#e9c46a" stroke-width="1.5" opacity="0.3" stroke-dasharray="5,4"/>')

# Arrowhead at end
px1, py1, _ = all_pts[-5]
px2, py2, _ = all_pts[-1]
angle = math.atan2(py2 - py1, px2 - px1)
adx = 10 * math.cos(angle)
ady = 10 * math.sin(angle)
px_p = -ady * 0.5
py_p = adx * 0.5
L(f'    <polygon points="{F(px2+adx)},{F(py2+ady)} {F(px2+px_p)},{F(py2+py_p)} {F(px2-px_p)},{F(py2-py_p)}" fill="#e9c46a" opacity="0.9"/>')

L('  </g>')

# Photon label
mi = n_photon // 3
pmx, pmy, _ = all_pts[mi]
L(f'  <text x="{F(pmx-20)}" y="{F(pmy-60)}" font-family="Georgia, serif" font-size="11" fill="#e9c46a" opacity="0.9">Photon path</text>')
L(f'  <text x="{F(pmx-20)}" y="{F(pmy-47)}" font-family="Georgia, serif" font-size="10" fill="#e9c46a" opacity="0.7">(helical on torus surface)</text>')
L(f'  <line x1="{F(pmx)}" y1="{F(pmy-43)}" x2="{F(pmx+5)}" y2="{F(pmy-18)}" stroke="#e9c46a" stroke-width="0.7" opacity="0.4" stroke-dasharray="3,2"/>')

# ============================================================
# CROSS-SECTION OVERLAY
# ============================================================
L('  <!-- Cross-section overlay -->')
L('  <g id="cross-section">')

cs_theta = 0
cs_cx2, cs_cy2 = project(R, 0, 0)
cs_offset = 195
cs_draw_cx = cs_cx2 + cs_offset
cs_draw_cy = cs_cy2

# Cut indicator on torus
cut_pts = []
for k in range(49):
    phi = 2 * math.pi * k / 48
    x3 = (R + r_minor * math.cos(phi))
    z3 = r_minor * math.sin(phi)
    sx, sy = project(x3, 0, z3)
    cut_pts.append((sx, sy))
d = "M " + " L ".join(f"{F(x)},{F(y)}" for x, y in cut_pts) + " Z"
L(f'    <path d="{d}" fill="#e0e0e0" fill-opacity="0.06" stroke="#e0e0e0" stroke-width="1.2" opacity="0.4"/>')

# Connecting line
L(f'    <line x1="{F(cs_cx2 + r_minor + 8)}" y1="{F(cs_cy2)}" x2="{F(cs_draw_cx - 50)}" y2="{F(cs_draw_cy)}" stroke="#e0e0e0" stroke-width="0.8" opacity="0.3" stroke-dasharray="5,5"/>')

# Cross-section circle
cs_r = 45
L(f'    <circle cx="{F(cs_draw_cx)}" cy="{F(cs_draw_cy)}" r="{F(cs_r)}" fill="#22242c" stroke="#a855f7" stroke-width="1.5" opacity="0.8"/>')

# B-field arrows circulating (poloidal)
n_arr = 8
arr_r = cs_r * 0.7
for k in range(n_arr):
    phi = 2 * math.pi * k / n_arr
    phi_next = phi + 0.3
    ax1 = cs_draw_cx + arr_r * math.cos(phi)
    ay1 = cs_draw_cy - arr_r * math.sin(phi)
    ax2 = cs_draw_cx + arr_r * math.cos(phi_next)
    ay2 = cs_draw_cy - arr_r * math.sin(phi_next)
    angle = math.atan2(ay2 - ay1, ax2 - ax1)
    L(f'    <line x1="{F(ax1)}" y1="{F(ay1)}" x2="{F(ax2)}" y2="{F(ay2)}" stroke="#2a9d8f" stroke-width="2" opacity="0.8"/>')
    adx = 5 * math.cos(angle)
    ady = 5 * math.sin(angle)
    px = -ady * 0.5
    py = adx * 0.5
    L(f'    <polygon points="{F(ax2+adx)},{F(ay2+ady)} {F(ax2+px)},{F(ay2+py)} {F(ax2-px)},{F(ay2-py)}" fill="#2a9d8f" opacity="0.8"/>')

# E-field dot (out of page)
L(f'    <circle cx="{F(cs_draw_cx)}" cy="{F(cs_draw_cy)}" r="7" fill="none" stroke="#e76f51" stroke-width="1.5" opacity="0.9"/>')
L(f'    <circle cx="{F(cs_draw_cx)}" cy="{F(cs_draw_cy)}" r="2.5" fill="#e76f51" opacity="0.9"/>')

# Labels
L(f'    <text x="{F(cs_draw_cx)}" y="{F(cs_draw_cy - cs_r - 12)}" font-family="Georgia, serif" font-size="12" fill="#e0e0e0" text-anchor="middle" font-weight="bold">Cross-section</text>')
L(f'    <text x="{F(cs_draw_cx + cs_r + 10)}" y="{F(cs_draw_cy - 12)}" font-family="Georgia, serif" font-size="11" fill="#2a9d8f" font-weight="bold" font-style="italic">B</text>')
L(f'    <text x="{F(cs_draw_cx + 16)}" y="{F(cs_draw_cy + 5)}" font-family="Georgia, serif" font-size="11" fill="#e76f51" font-style="italic">E &#x2299;</text>')

L('  </g>')

# ============================================================
# DIMENSION LINES
# ============================================================
L('  <!-- Dimension lines -->')

center_sx, center_sy = CX, CY
tube_cx, tube_cy = project(R, 0, 0)
dim_y = center_sy + 145

L(f'  <line x1="{F(center_sx)}" y1="{F(dim_y)}" x2="{F(tube_cx)}" y2="{F(dim_y)}" stroke="#e9c46a" stroke-width="1.5" marker-start="url(#arrowDimRev)" marker-end="url(#arrowDim)"/>')
L(f'  <line x1="{F(center_sx)}" y1="{F(dim_y-5)}" x2="{F(center_sx)}" y2="{F(dim_y+5)}" stroke="#e9c46a" stroke-width="1.2"/>')
L(f'  <line x1="{F(tube_cx)}" y1="{F(dim_y-5)}" x2="{F(tube_cx)}" y2="{F(dim_y+5)}" stroke="#e9c46a" stroke-width="1.2"/>')
L(f'  <line x1="{F(center_sx)}" y1="{F(center_sy+5)}" x2="{F(center_sx)}" y2="{F(dim_y-5)}" stroke="#e9c46a" stroke-width="0.7" opacity="0.3" stroke-dasharray="4,3"/>')
L(f'  <line x1="{F(tube_cx)}" y1="{F(tube_cy + r_minor * COS_T + 5)}" x2="{F(tube_cx)}" y2="{F(dim_y-5)}" stroke="#e9c46a" stroke-width="0.7" opacity="0.3" stroke-dasharray="4,3"/>')

mid_rx = (center_sx + tube_cx) / 2
L(f'  <text x="{F(mid_rx)}" y="{F(dim_y-8)}" font-family="Georgia, serif" font-size="14" fill="#e9c46a" text-anchor="middle" font-style="italic">R</text>')
L(f'  <text x="{F(mid_rx)}" y="{F(dim_y+18)}" font-family="Georgia, serif" font-size="11" fill="#e9c46a" text-anchor="middle" opacity="0.85">R &#x2248; &#x03BB;C</text>')

# Minor radius r
tube_top_sx, tube_top_sy = project(R, 0, r_minor)
rd_x = 32
L(f'  <line x1="{F(tube_cx + rd_x)}" y1="{F(tube_cy)}" x2="{F(tube_top_sx + rd_x)}" y2="{F(tube_top_sy)}" stroke="#e9c46a" stroke-width="1.5" marker-start="url(#arrowDimRev)" marker-end="url(#arrowDim)"/>')
L(f'  <line x1="{F(tube_cx + rd_x - 4)}" y1="{F(tube_cy)}" x2="{F(tube_cx + rd_x + 4)}" y2="{F(tube_cy)}" stroke="#e9c46a" stroke-width="1.2"/>')
L(f'  <line x1="{F(tube_top_sx + rd_x - 4)}" y1="{F(tube_top_sy)}" x2="{F(tube_top_sx + rd_x + 4)}" y2="{F(tube_top_sy)}" stroke="#e9c46a" stroke-width="1.2"/>')

rl_x = tube_cx + rd_x + 14
rl_y = (tube_cy + tube_top_sy) / 2
L(f'  <text x="{F(rl_x)}" y="{F(rl_y)}" font-family="Georgia, serif" font-size="14" fill="#e9c46a" font-style="italic">r</text>')
L(f'  <text x="{F(rl_x+14)}" y="{F(rl_y-8)}" font-family="Georgia, serif" font-size="11" fill="#e9c46a" opacity="0.85">&#x2248; r<tspan dy="4" font-size="8">e</tspan></text>')

# Center dot
L(f'  <circle cx="{F(center_sx)}" cy="{F(center_sy)}" r="3" fill="#e9c46a" opacity="0.8"/>')

# ============================================================
# RATIO NOTE
# ============================================================
L('  <g transform="translate(450, 540)">')
L('    <rect x="-120" y="-18" width="240" height="36" rx="6" fill="#22242c" stroke="#e9c46a" stroke-width="1" opacity="0.7"/>')
L('    <text x="0" y="5" font-family="Georgia, serif" font-size="14" fill="#e9c46a" text-anchor="middle">r / R = &#x03B1; &#x2248; 1/137</text>')
L('  </g>')

# ============================================================
# HOPF LINKING INSET
# ============================================================
L('  <!-- Hopf Linking Inset -->')
L('  <g transform="translate(130, 660)">')
L('    <rect x="-100" y="-55" width="200" height="110" rx="8" fill="#22242c" stroke="#444" stroke-width="1"/>')
L('    <text x="0" y="-35" font-family="Georgia, serif" font-size="12" fill="#bbb" text-anchor="middle" font-weight="bold">Hopf Linking</text>')
L('    <ellipse cx="-8" cy="8" rx="35" ry="18" fill="none" stroke="#e76f51" stroke-width="2" opacity="0.3" stroke-dasharray="3,2"/>')
L('    <ellipse cx="8" cy="8" rx="18" ry="32" fill="none" stroke="#2a9d8f" stroke-width="2.2" opacity="1"/>')
L('    <path d="M -43,8 A 35,18 0 0,0 27,8" fill="none" stroke="#e76f51" stroke-width="2.5"/>')
L('    <text x="-50" y="-3" font-family="Georgia, serif" font-size="10" fill="#e76f51" font-weight="bold">E</text>')
L('    <text x="32" y="-20" font-family="Georgia, serif" font-size="10" fill="#2a9d8f" font-weight="bold">B</text>')
L('    <text x="0" y="52" font-family="Georgia, serif" font-size="9" fill="#999" text-anchor="middle" font-style="italic">Every E-line links every B-line</text>')
L('  </g>')

# ============================================================
# HOPF INVARIANT BOX
# ============================================================
L('  <g transform="translate(730, 655)">')
L('    <rect x="-150" y="-35" width="300" height="75" rx="8" fill="#22242c" stroke="#a855f7" stroke-width="1.5" opacity="0.8"/>')
L('    <text x="0" y="-8" font-family="Georgia, serif" font-size="18" fill="#e0e0e0" text-anchor="middle" font-weight="bold">H = &#xB1;1  &#x2192;  Q = &#xB1;e</text>')
L('    <text x="0" y="16" font-family="Georgia, serif" font-size="12" fill="#bbb" text-anchor="middle" font-style="italic">Topological invariant &#x2192; quantized charge</text>')
L('  </g>')

# ============================================================
# CAPTION
# ============================================================
L('  <text x="500" y="775" font-family="Georgia, serif" font-size="13" fill="#999" text-anchor="middle">')
L('    Figure 2: Toroidal electron &#x2014; B-field (poloidal) and E-field (toroidal) with Hopf linking H = &#xB1;1, showing photon helical path and cross-section')
L('  </text>')

L('</svg>')

svg = '\n'.join(lines)

script_dir = os.path.dirname(os.path.abspath(__file__))
out_path = os.path.join(script_dir, "DarkMatter-Fig2-ToroidalElectron.svg")
with open(out_path, 'w', encoding='utf-8') as f:
    f.write(svg)
print(f"Generated {len(lines)} lines, {len(svg)} chars -> {out_path}")
