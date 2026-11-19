#!/usr/bin/env python3
"""Generate a watertight ASCII STL of a simple sloop for OpenFOAM.

All dimensions in metres. The boat is built from closed manifold bodies;
each body is checked for winding (signed volume) and edge-pairing
(watertightness) before the file is written.

Bodies:
    hull, keel fin, keel bulb, rudder          (waterborne)
    mast, mainsail, headsail                   (aero rig)

v3 additions:
    --roll R        heel the whole boat about its x-axis (deg, + = starboard)
    --yaw Y         rotate heading about vertical (deg, + = bow to starboard)
    --main-sheet A  ease the mainsail about its luff (deg from fully trimmed)
    --jib-sheet A   ease the jib about the forestay (deg from fully trimmed)

Sheet angles are rigid rotations of each sail about its luff axis: 0 deg is
the fully-trimmed (100%) geometry; positive values open the leech outboard.
Roll and yaw rotate the whole boat (hull + rig), applied after sheeting, so
the wind stays fixed along +x while you vary heel / heading / trim.

Usage:
    python3 sailboat_stl_v3.py                           # full boat, flat, trimmed
    python3 sailboat_stl_v3.py --sails-only              # rig only (airflow case)
    python3 sailboat_stl_v3.py --roll 20 --yaw -15       # heeled to starboard
    python3 sailboat_stl_v3.py --sails-only --main-sheet 12 --jib-sheet 6

Standard library only, Python 3.8+.
"""

import argparse
import math
import os
from collections import Counter


# --------------------------------------------------------------- mesh bits

def loft(rings):
    """Triangles between consecutive closed rings (equal point counts)."""
    tris = []
    n = len(rings[0])
    for i in range(len(rings) - 1):
        a, b = rings[i], rings[i + 1]
        for j in range(n):
            k = (j + 1) % n
            tris.append((a[j], b[j], b[k]))
            tris.append((a[j], b[k], a[k]))
    return tris


def cap(ring, center, start=True):
    """Fan-cap a ring from an interior centre point (never on the ring)."""
    n = len(ring)
    tris = []
    for j in range(n):
        k = (j + 1) % n
        if start:
            tris.append((center, ring[j], ring[k]))
        else:
            tris.append((center, ring[k], ring[j]))
    return tris


def sphere(cx, cy, cz, r, nlat=12, nlon=24):
    rings = []
    for i in range(1, nlat):
        th = math.pi * i / nlat
        st, ct = math.sin(th), math.cos(th)
        rings.append([(cx + r * st * math.cos(2 * math.pi * j / nlon),
                       cy + r * st * math.sin(2 * math.pi * j / nlon),
                       cz + r * ct) for j in range(nlon)])
    tris = cap(rings[0], (cx, cy, cz + r), start=True)
    tris += loft(rings)
    tris += cap(rings[-1], (cx, cy, cz - r), start=False)
    return tris


# -------------------------------------------------------- rigid transforms

def rot_x(p, c, s):
    """Heel about the x-axis; positive angle tips +z toward +y (starboard)."""
    x, y, z = p
    return (x, y * c + z * s, -y * s + z * c)


def rot_z(p, c, s):
    """Heading rotation about vertical; positive swings +x toward +y."""
    x, y, z = p
    return (x * c - y * s, x * s + y * c, z)


def orient_bodies(bodies, roll_deg, yaw_deg):
    """Heel the whole boat about its x-axis, then rotate heading about
    vertical. Rigid: watertightness, winding and volumes are preserved."""
    cr = math.cos(math.radians(roll_deg))
    sr = math.sin(math.radians(roll_deg))
    cy = math.cos(math.radians(yaw_deg))
    sy = math.sin(math.radians(yaw_deg))
    out = []
    for name, tris in bodies:
        new = []
        for (a, b, c3) in tris:
            q = (rot_x(a, cr, sr), rot_x(b, cr, sr), rot_x(c3, cr, sr))
            new.append((rot_z(q[0], cy, sy), rot_z(q[1], cy, sy),
                        rot_z(q[2], cy, sy)))
        out.append((name, new))
    return out


def rotate_body(tris, p0, u, ang_deg):
    """Rigidly rotate a closed body about the line through p0 along unit
    vector u (Rodrigues' formula). Used for sheeting a sail about its luff.
    Rigid: watertightness, winding and volume are preserved."""
    c = math.cos(math.radians(ang_deg))
    s = math.sin(math.radians(ang_deg))
    ux, uy, uz = u
    out = []
    for (a, b, cc) in tris:
        def rp(p):
            vx, vy, vz = p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]
            dot = ux * vx + uy * vy + uz * vz
            cx = uy * vz - uz * vy
            cyy = uz * vx - ux * vz
            cz = ux * vy - uy * vx
            return (p0[0] + vx * c + cx * s + ux * dot * (1.0 - c),
                    p0[1] + vy * c + cyy * s + uy * dot * (1.0 - c),
                    p0[2] + vz * c + cz * s + uz * dot * (1.0 - c))
        out.append((rp(a), rp(b), rp(cc)))
    return out


# ------------------------------------------------------------- components

def hull(L, B, D, ns=30, nr=28):
    """Lofted hull. x runs 0 (bow) -> L (stern); deck line sits near z = 0.
    End rings stop short of the tips so bow/stern caps are real cone fans."""
    sheer = 0.02 * L

    def half_beam(s):
        return 0.5 * B * math.sin(math.pi * s) ** 0.85

    def draft(s):
        return D * math.sin(math.pi * s) ** 0.6

    rings = []
    for i in range(1, ns):                     # skip s=0 and s=1 (point rings)
        s = i / ns
        x = s * L
        b, d = half_beam(s), draft(s)
        z_top = sheer * (1.0 - s)              # deck rises to the bow
        ring = []
        for j in range(nr):
            th = 2 * math.pi * j / nr
            flare = 1.0 + 0.35 * max(0.0, math.cos(th)) ** 2
            y = b * math.sin(th) * flare
            z = z_top - 0.5 * d + 0.5 * d * math.cos(th)
            ring.append((x, y, z))
        rings.append(ring)

    tris = cap(rings[0], (0.0, 0.0, sheer), start=True)
    tris += loft(rings)
    tris += cap(rings[-1], (L, 0.0, 0.0), start=False)
    return tris, draft


def fin(x0, x1, hmax, z_top, z_bot, ns=16, nr=14):
    """Vertical flat-plate fin (keel / rudder) with rounded LE and TE."""
    hc = 0.5 * (z_top - z_bot)          # positive half-chord
    zc = 0.5 * (z_bot + z_top)
    rings = []
    for i in range(ns + 1):
        t = i / ns
        h = hmax * max(math.sin(math.pi * t) ** 0.7, 0.12)  # no zero-thickness caps
        x = x0 + t * (x1 - x0)
        rings.append([(x, h * math.sin(2 * math.pi * j / nr),
                       zc + hc * math.cos(2 * math.pi * j / nr))
                      for j in range(nr)])
    tris = cap(rings[0], (x0, 0.0, zc), start=True)
    tris += loft(rings)
    tris += cap(rings[-1], (x1, 0.0, zc), start=False)
    return tris


def mast(xm, r, z0, z1, nz=24, nr=16):
    rings = []
    for i in range(nz + 1):
        z = z0 + (z1 - z0) * i / nz
        rings.append([(xm + r * math.sin(2 * math.pi * j / nr),
                       r * math.cos(2 * math.pi * j / nr), z)
                      for j in range(nr)])
    tris = cap(rings[0], (xm, 0.0, z0), start=True)
    tris += loft(rings)
    tris += cap(rings[-1], (xm, 0.0, z1), start=False)
    return tris


def sail(x_luff, z0, z1, chord_max, th_max, camber, nz=24, nr=20):
    """Mainsail: thin closed cambered lens, luff vertical at the mast.
    The top ring stops one step short of z1 so the head cap is a real cone."""
    dz = (z1 - z0) / nz
    rings = []
    for i in range(nz):
        uz = i / (nz - 1)
        z = z0 + (z1 - z0 - dz) * uz           # last ring at z1 - dz
        chord = chord_max * (1.0 - 0.92 * uz) ** 0.85
        ring = []
        for j in range(nr):
            th = 2 * math.pi * j / nr
            u = 0.5 * (1.0 - math.cos(th))     # 0 luff -> 1 leech
            x = x_luff + chord * u
            y = (camber * 4.0 * u * (1.0 - u)
                 + th_max * math.sin(th) * math.sqrt(u * (1.0 - u)))
            ring.append((x, y, z))
        rings.append(ring)

    tris = cap(rings[0], (x_luff + 0.45 * chord_max,
                          camber * 4 * 0.45 * 0.55, z0), start=True)
    tris += loft(rings)
    tris += cap(rings[-1], (x_luff, 0.0, z1), start=False)
    return tris


def headsail(x_tack, z_tack, x_head, z_head, x_clew, z_foot,
             th_max, camber, roach=0.0, nz=24, nr=18):
    """Jib: luff along the forestay (tack -> head), foot to the clew,
    leech from clew to head with optional roach. Thin closed lens in z."""
    dz = (z_head - z_foot) / nz
    t_foot = (z_foot - z_tack) / (z_head - z_tack)
    x_l_foot = x_tack + (x_head - x_tack) * t_foot
    chord_foot = x_clew - x_l_foot
    min_chord = 0.04 * chord_foot

    rings = []
    for i in range(nz):
        uz = i / (nz - 1)
        z = z_foot + (z_head - z_foot - dz) * uz   # last ring one step short of head
        t = (z - z_tack) / (z_head - z_tack)
        x_l = x_tack + (x_head - x_tack) * t       # luff on the forestay line
        x_le = x_clew + (x_head - x_clew) * uz + roach * math.sin(math.pi * uz)
        chord = max(x_le - x_l, min_chord)         # keep head ring non-degenerate
        ring = []
        for j in range(nr):
            th = 2 * math.pi * j / nr
            u = 0.5 * (1.0 - math.cos(th))
            x = x_l + chord * u
            y = (camber * 4.0 * u * (1.0 - u)
                 + th_max * math.sin(th) * math.sqrt(u * (1.0 - u)))
            ring.append((x, y, z))
        rings.append(ring)

    tris = cap(rings[0], (x_l_foot + 0.45 * chord_foot,
                          camber * 4 * 0.45 * 0.55, z_foot), start=True)
    tris += loft(rings)
    tris += cap(rings[-1], (x_head, 0.0, z_head), start=False)
    return tris


# ---------------------------------------------------------------- assembly

def build(L, B, D, mode, main_sheet=0.0, jib_sheet=0.0):
    """mode: 'full' | 'hull' | 'sails'. Sheet angles in degrees (0 = 100% trimmed)."""
    bodies = []

    if mode in ('full', 'hull'):
        hull_tris, draft_fn = hull(L, B, D)
        bodies.append(('hull', hull_tris))

        kx0, kx1 = 0.40 * L, 0.62 * L
        kz_top, kz_bot = -(D - 0.15), -(D + 0.85 * L / 6.0)
        bodies.append(('keel fin', fin(kx0, kx1, 0.009 * L, kz_top, kz_bot)))

        rb = 0.045 * L
        bodies.append(('keel bulb',
                       sphere(0.5 * (kx0 + kx1), 0.0, kz_bot - 0.35 * rb, rb)))

        rx0, rx1 = 0.86 * L, 0.93 * L
        s_r = 0.5 * (rx0 + rx1) / L
        z_b = 0.02 * L * (1.0 - s_r) - draft_fn(s_r)
        bodies.append(('rudder', fin(rx0, rx1, 0.005 * L,
                                     z_b + 0.08, z_b - 0.45 * L / 6.0)))

    if mode in ('full', 'sails'):
        xm, mr, mz = 0.40 * L, 0.0075 * L, 0.62 * L
        bodies.append(('mast', mast(xm, mr, -0.1, mz)))

        # mainsail; sheeting = rotation about the vertical luff line at the mast
        x_luff = xm + mr + 0.01
        ms = sail(x_luff, 0.1, mz - 0.2, 0.28 * L, 0.035, 0.05)
        if main_sheet:
            ms = rotate_body(ms, (x_luff, 0.0, 0.0), (0.0, 0.0, 1.0), main_sheet)
        bodies.append(('mainsail', ms))

        # headsail (jib); sheeting = rotation about the forestay axis
        x_tack, z_tack = 0.08 * L, 0.12
        x_head, z_head = xm - mr - 0.02, mz - 0.05
        x_clew, z_foot = 0.36 * L, 0.15
        hs = headsail(x_tack, z_tack, x_head, z_head, x_clew, z_foot,
                      0.03, 0.04, roach=0.05)
        if jib_sheet:
            dx, dzv = x_head - x_tack, z_head - z_tack
            nlen = math.hypot(dx, dzv)
            hs = rotate_body(hs, (x_tack, 0.0, z_tack),
                             (dx / nlen, 0.0, dzv / nlen), jib_sheet)
        bodies.append(('headsail', hs))
    return bodies


# --------------------------------------------------------------- checking

def signed_volume(tris):
    v = 0.0
    for (p0, p1, p2) in tris:
        v += (p0[0] * (p1[1] * p2[2] - p1[2] * p2[1])
              + p0[1] * (p1[2] * p2[0] - p1[0] * p2[2])
              + p0[2] * (p1[0] * p2[1] - p1[1] * p2[0]))
    return v / 6.0


def is_watertight(tris):
    """Closed manifold: every undirected edge shared by exactly 2 triangles."""
    c = Counter()
    for (a, b, cc) in tris:
        for e in ((a, b), (b, cc), (cc, a)):
            c[tuple(sorted(e))] += 1
    return all(n == 2 for n in c.values())


def shift(tris, ox, oy, oz):
    out = []
    for (a, b, c) in tris:
        out.append(((a[0] + ox, a[1] + oy, a[2] + oz),
                    (b[0] + ox, b[1] + oy, b[2] + oz),
                    (c[0] + ox, c[1] + oy, c[2] + oz)))
    return out


def write_stl(path, bodies):
    with open(path, 'w') as f:
        f.write('solid sailboat\n')
        for _, tris in bodies:
            for (p0, p1, p2) in tris:
                ux, uy, uz = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
                vx, vy, vz = p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]
                f.write('  facet normal %.6f %.6f %.6f\n'
                        % (uy * vz - uz * vy, uz * vx - ux * vz,
                           ux * vy - uy * vx))
                f.write('    outer loop\n')
                for p in (p0, p1, p2):
                    f.write('      vertex %.6f %.6f %.6f\n' % p)
                f.write('    endloop\n  endfacet\n')
        f.write('endsolid sailboat\n')


def main():
    ap = argparse.ArgumentParser(description='Watertight sloop STL for OpenFOAM')
    ap.add_argument('-o', '--output', default='sailboat.stl')
    ap.add_argument('--split', action='store_true',
                    help="write one STL per body ('<stem>_<body>.stl') so "
                         "snappyHexMesh can assign a named patch per sail; "
                         "skips the combined file")
    grp = ap.add_mutually_exclusive_group()
    grp.add_argument('--hull-only', action='store_true',
                     help='omit the rig (hydro case)')
    grp.add_argument('--sails-only', action='store_true',
                     help='mast + mainsail + headsail only (airflow case)')
    ap.add_argument('--length', type=float, default=6.0, help='LOA in metres')
    ap.add_argument('--beam', type=float, default=1.9, help='max beam, m')
    ap.add_argument('--draft', type=float, default=1.1,
                    help='hull draft excluding keel, m')
    ap.add_argument('--roll', type=float, default=0.0,
                    help='heel in degrees; positive = to starboard')
    ap.add_argument('--yaw', type=float, default=0.0,
                    help='heading rotation about vertical, degrees; '
                         'positive = bow toward starboard (wind stays along +x)')
    ap.add_argument('--main-sheet', type=float, default=0.0,
                    help='mainsail opening about its luff, degrees from the '
                         'fully-trimmed (100%%) position; positive = eased outboard')
    ap.add_argument('--jib-sheet', type=float, default=0.0,
                    help='jib opening about the forestay, degrees from fully '
                         'trimmed; positive = eased outboard')
    ap.add_argument('--offset', nargs=3, type=float, default=(0, 0, 0),
                    metavar=('X', 'Y', 'Z'), help='shift the whole boat, m')
    args = ap.parse_args()

    mode = 'hull' if args.hull_only else ('sails' if args.sails_only else 'full')
    bodies = build(args.length, args.beam, args.draft, mode,
                   args.main_sheet, args.jib_sheet)
    if args.roll or args.yaw:
        bodies = orient_bodies(bodies, args.roll, args.yaw)
    ox, oy, oz = args.offset
    bodies = [(name, shift(tris, ox, oy, oz)) for name, tris in bodies]

    # safety net: flip any body whose winding came out inward
    checked = []
    for name, tris in bodies:
        if signed_volume(tris) < 0:
            print('note: flipped winding of "%s"' % name)
            tris = [(a, c, b) for (a, b, c) in tris]
        checked.append((name, tris))

    print('%-12s %10s %12s  %s' % ('body', 'triangles', 'volume m^3', 'watertight'))
    ok_all, total = True, 0
    for name, tris in checked:
        ok = is_watertight(tris)
        ok_all &= ok
        total += len(tris)
        print('%-12s %10d %12.3f  %s'
              % (name, len(tris), signed_volume(tris), 'yes' if ok else 'NO - do not use'))

    if args.split:
        stem = os.path.splitext(args.output)[0]
        for name, tris in checked:
            p = '%s_%s.stl' % (stem, name.replace(' ', '_'))
            write_stl(p, [(name, tris)])
            print('wrote %s (%d triangles)' % (p, len(tris)))
    else:
        write_stl(args.output, checked)
        print('\n%d triangles -> %s' % (total, args.output))
    print('axes: x bow(0)->stern, y port->starboard, z up; deck ~ z=0')
    if not ok_all:
        raise SystemExit('watertight check failed')


if __name__ == '__main__':
    main()
