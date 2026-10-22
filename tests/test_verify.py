#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright 2026 Zhizheng Ju
# SPDX-License-Identifier: Apache-2.0
"""e213ips.verify -- shipped tests (public code-only candidate).

Every control is CONSTRUCTED in this file, either by an exact formula applied to an elementary
triangle or by hand-solved small integer coordinates. No coordinate dataset is read, so nothing here
depends on data whose rights or provenance are unresolved. The hand-solved sets are commented with the
arithmetic that makes them work, which is what lets a reader check the fixture without trusting it.

Run: python3 tests/test_verify.py
"""

from __future__ import division

import itertools
from fractions import Fraction

from _harness import check, expect_raises, finish

from e213ips import verify as V


def circumcentre(a, b, c):
    """Exact perpendicular-bisector solve for the centre of the circle through three non-collinear
    points; None when they are collinear. A different algorithm from the determinant under test."""
    ax, ay = Fraction(a[0]), Fraction(a[1])
    bx, by = Fraction(b[0]), Fraction(b[1])
    cx, cy = Fraction(c[0]), Fraction(c[1])
    d = 2 * ((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))
    if d == 0:
        return None
    # Solved for the OFFSET (h - ax, k - ay), whose right-hand sides are |b-a|^2 and |c-a|^2 -- not
    # |b|^2 - |a|^2. Writing the un-shifted differences here was this helper's own bug, and the grid
    # agreement leg below is what caught it: 157 quadruples disagreed, all of them with a[0] != (0,0).
    e1 = (bx - ax) ** 2 + (by - ay) ** 2
    e2 = (cx - ax) ** 2 + (cy - ay) ** 2
    return (ax + (e1 * (cy - ay) - e2 * (by - ay)) / d,
            ay + ((bx - ax) * e2 - (cx - ax) * e1) / d)


def hub_from_scratch():
    """LEGAL integral 4-set in the s=1 plane, derived from Heron geometry.

    A 40-40-48 triangle (a 5-5-6 triangle scaled by 8) has 16K^2 = 128*900, K = 768, and
    R = abc/4K = 40*40*48/(4*768) = 25, an integer. Put the base on the x-axis:
    A=(0,0), B=(48,0), C=(24,32), so |AB|=48 and |AC|=|BC|=sqrt(24^2+32^2)=40. A triangle's vertices
    together with its circumcentre form a 4-point set with six integral distances and a degree-3 hub.
    """
    tri = [(0, 0), (48, 0), (24, 32)]
    for a, b in itertools.combinations(tri, 2):
        assert V.integral_distance(a, b, 1) is not None, "fixture triangle must be integral"
    o = circumcentre(*tri)
    r2 = V.dist2(o, tri[0], 1)
    assert V.is_square_rational(r2), "the circumradius must be integral for the hub to be usable"
    return tri + [(o[0], o[1])], int(V._isqrt_fraction(r2))


def tri_s2():
    """LEGAL integral triangle in the s=2 model plane: (0,0),(2,0),(1,2).
    |(2,0)|=2; |(1,2)|^2 = 1 + 2*4 = 9; and the same from (2,0): 1 + 8 = 9. Sides 2,3,3."""
    return [(0, 0), (2, 0), (1, 2)]


def trapezoid_s2():
    """INTEGRAL, non-collinear, and ILLEGAL on concyclicity -- in the s=2 plane, built by hand.

    Vectors of square norm in Z[sqrt(-2)]: (9,0)->9, (7,4)->49+32=81, (2,4)->4+32=36. So
    p0=(0,0), p1=(9,0), p2=(7,4) have lengths 9,9 and |p1p2| = |(2,-4)| = 6. Add q=(2,4): |q|=6,
    |q-p1| = |(7,-4)| = 9, |q-p2| = 5. Six distances {9,9,9,6,6,5}, no three collinear -- and the two
    points at t=4 sit on a line parallel to the two at t=0 with equal legs 6 and 6, which is an
    ISOSCELES TRAPEZOID: parallel chords plus equal legs force one circle. Every distance is integral
    and general position fails on exactly the concyclicity condition -- the obstruction the eight-point
    problem actually hits, reproduced from three lines of arithmetic."""
    return [(0, 0), (9, 0), (7, 4), (2, 4)]


def grid_points(lo, hi):
    return [(Fraction(x), Fraction(y)) for x in range(lo, hi + 1) for y in range(lo, hi + 1)]


def all_heron(P, s):
    rep = V.verify(P, s)
    return [V.heron_16k2(rep.distances[t], rep.distances[u], rep.distances[v])
            for t, u, v in ((a, b, c) for a, b, c in
                            (((i, j), (i, k), (j, k)) for i, j, k in itertools.combinations(range(len(P)), 3)))
            if (t in rep.distances and u in rep.distances and v in rep.distances)]


def main():
    # 1. POSITIVE, constructed above: the hub 4-set, radius derived from the Heron formula.
    hub, radius = hub_from_scratch()
    rep = V.verify(hub, 1)
    check("R-VERIFY the hub 4-set built from a 40-40-48 triangle is legal, its circumradius is %d and "
          "its six distances are exactly 25,25,25,40,40,48" % radius,
          rep.is_legal and sorted(rep.distances.values()) == [25, 25, 25, 40, 40, 48] and radius == 25,
          "reasons=%s dist=%s" % (rep.reasons, sorted(rep.distances.values())))

    # 2. POSITIVE in a characteristic-2 plane, so the s-weighted predicates run on a legal set at s>1.
    t2 = tri_s2()
    rep = V.verify(t2, 2)
    check("R-VERIFY the constructed characteristic-2 triangle is legal with sides 2,3,3",
          rep.is_legal and sorted(rep.distances.values()) == [2, 3, 3], "reasons=%s" % rep.reasons)

    # 3. NEGATIVE on CONCYCLICITY in a characteristic-2 plane: the hand-built isosceles trapezoid.
    trap2 = trapezoid_s2()
    rep = V.verify(trap2, 2)
    check("R-VERIFY the hand-built s=2 trapezoid has 6 integral distances, no collinear triple and "
          "fails on exactly 1 concyclic quad -- the isosceles-trapezoid mechanism, no data needed",
          rep.integral_distance_count == 6 and not rep.collinear_triples
          and len(rep.concyclic_quads) == 1 and rep.reasons == ["1 concyclic quad(s)"]
          and sorted(rep.distances.values()) == [5, 6, 6, 9, 9, 9],
          "reasons=%s" % rep.reasons)

    # 4. NEGATIVE on concyclicity at s=1 (the classic: a rectangle).
    rect = [(0, 0), (6, 0), (6, 8), (0, 8)]
    rep = V.verify(rect, 1)
    check("R-VERIFY the concyclicity detector FIRES on a 6x8 rectangle: integral, non-collinear, "
          "illegal on exactly 1 concyclic quad",
          rep.integral_distance_count == 6 and not rep.collinear_triples
          and len(rep.concyclic_quads) == 1 and rep.reasons == ["1 concyclic quad(s)"],
          "reasons=%s" % rep.reasons)

    # 5. NEGATIVE on collinearity, reported AS collinearity.
    rep = V.verify([(0, 0), (3, 0), (6, 0)], 1)
    check("R-VERIFY a collinear triple is caught as collinearity, not as concyclicity or non-integrality",
          rep.integral_distance_count == 3 and len(rep.collinear_triples) == 1
          and not rep.concyclic_quads and rep.reasons == ["1 collinear triple(s)"],
          "reasons=%s" % rep.reasons)

    # 6. NEGATIVE on integrality: sqrt(2) is refused, not rounded.
    rep = V.verify([(0, 0), (1, 1)], 1)
    check("R-VERIFY sqrt(2) is refused rather than rounded",
          not rep.is_legal and rep.non_integral_pairs == [(0, 1)], "reasons=%s" % rep.reasons)

    # 7. CHARACTERISTIC both ways: each constructed set must FAIL at the other s.
    t2_at1 = V.verify(t2, 1)
    hub_at2 = V.verify(hub, 2)
    check("R-VERIFY the s=2 triangle fails integrality read at s=1 (%d of 3 pairs) and the s=1 hub "
          "fails at s=2 (%d of 6): s is load-bearing, not decoration"
          % (len(t2_at1.non_integral_pairs), len(hub_at2.non_integral_pairs)),
          not t2_at1.is_legal and not hub_at2.is_legal
          and t2_at1.non_integral_pairs and hub_at2.non_integral_pairs,
          "tri@1=%d hub@2=%d" % (len(t2_at1.non_integral_pairs), len(hub_at2.non_integral_pairs)))

    # 8. CONCYCLICITY is s-sensitive in BOTH directions, on coordinates checked by hand:
    #    the diamond {(0,0),(1,1),(0,2),(-1,1)} is a square at s=1 (cyclic) and at s=2 it is a rhombus
    #    with diagonals 2 and 2*sqrt2, which is cyclic only when it is a square -- so not cyclic;
    #    {(3,0),(-3,0),(1,2),(1,-2)} lies on x^2 + 2*t^2 = 9, i.e. the true-plane circle of radius 3
    #    (both x-values give 9, both points at t=+-2 give 1+8), while the circle through its first
    #    three points at s=1 has centre (0,-1) and r^2 = 10, which misses (1,-2) at 2.
    diamond = [(0, 0), (1, 1), (0, 2), (-1, 1)]
    kite2 = [(3, 0), (-3, 0), (1, 2), (1, -2)]
    check("R-VERIFY the circle predicate moves with the characteristic: cyclic at s=1 and not at s=2 "
          "for the diamond; not cyclic at s=1 and cyclic at s=2 for the second fixture",
          V.concyclic(diamond, 1) and not V.concyclic(diamond, 2)
          and not V.concyclic(kite2, 1) and V.concyclic(kite2, 2),
          "diamond=%s/%s kite2=%s/%s" % (V.concyclic(diamond, 1), V.concyclic(diamond, 2),
                                         V.concyclic(kite2, 1), V.concyclic(kite2, 2)))

    # 9. isosceles trapezoid cyclic, scalene trapezoid not (s=1).
    trap = [(0, 0), (6, 0), (4, 4), (2, 4)]
    scat = [(0, 0), (6, 0), (5, 4), (2, 4)]
    check("R-VERIFY an isosceles trapezoid is detected concyclic and a scalene one is not",
          V.concyclic(trap, 1) and not V.concyclic(scat, 1),
          "trap=%s scalene=%s" % (V.concyclic(trap, 1), V.concyclic(scat, 1)))

    # 10. TWO IMPLEMENTATIONS over a dense grid: the fraction-free determinant and the circumcentre
    #     construction must agree on every 4-subset of a 4x4 lattice, with BOTH verdicts represented.
    pts = grid_points(0, 3)
    agree = conc = 0
    disagree, skipped = [], 0
    for q in itertools.combinations(pts, 4):
        d = V.det4(list(q), 1)
        # det4 vanishes for a cyclic quad AND for any quad holding a collinear triple (the circle of
        # infinite radius), so the circle-fit comparison only applies when no triple is collinear.
        if any(V.collinear(*tr) for tr in itertools.combinations(q, 3)):
            skipped += 1
            continue
        c0 = circumcentre(q[0], q[1], q[2])
        if c0 is None:
            skipped += 1
            continue
        r2 = V.dist2(c0, q[0], 1)
        fit = (q[3][0] - c0[0]) ** 2 + (q[3][1] - c0[1]) ** 2 == r2
        if (d == 0) == fit:
            agree += 1
        else:
            disagree.append(q)
        conc += (d == 0)
    check("R-VERIFY det4 (fraction-free) and the circumcentre construction agree on %d grid "
          "quadruples: %d concyclic, %d not, %d skipped because three of the four are collinear"
          % (agree, conc, agree - conc, skipped),
          agree + skipped == 1820 and agree >= 1000 and 0 < conc < agree and not disagree,
          "agree=%d concyclic=%d disagree=%d" % (agree, conc, len(disagree)))

    # 11. Bareiss vs Laplace as linear algebra, including a zero first pivot (forces row pivoting), a
    #     singular matrix, and the determinant that must vanish on a concyclic quadruple.
    cyclic4 = [[x * x + y * y, x, y, 1] for x, y in [(0, 0), (6, 0), (6, 8), (0, 8)]]
    zero_pivot = [[0, 2, 3, 4], [1, 5, 6, 7], [2, 9, 10, 11], [3, 13, 14, 16]]
    ok = (V.det_bareiss(cyclic4) == V.det_laplace(cyclic4) == 0
          and V.det_bareiss(zero_pivot) == V.det_laplace(zero_pivot)
          and V.det_bareiss([[1, 2], [2, 4]]) == 0)
    check("R-VERIFY det_bareiss == det_laplace == 0 on a cyclic quadruple, agrees on a zero-pivot "
          "4x4, and returns 0 on a singular 2x2", ok, "cyclic=%s" % (V.det_bareiss(cyclic4) == 0))

    # 12. MUTATION: nudge one vertex and integrality must break -- the verifier reads coordinates.
    mutated = [hub[0], hub[1], (hub[2][0] + 1, hub[2][1]), hub[3]]
    check("R-VERIFY nudging one vertex by +1 makes the hub set non-integral (sensitivity, not a "
          "stored flag)", not V.verify(mutated, 1).is_legal,
          "reasons=%s" % V.verify(mutated, 1).reasons[:1])

    # 13. exactness at the door.
    expect_raises("R-VERIFY a float coordinate is REFUSED with TypeError", TypeError,
                  V.verify, [(0.0, 0.0), (3.0, 4.0)], 1)

    # 14. the characteristic is a property of the DISTANCES: the hub's four triangles sit in square
    #     class 1, and the s=2 triangle's single Heron value sits in class 2.
    hh = all_heron(hub, 1)
    th = all_heron(t2, 2)
    check("R-VERIFY all %d hub triangles have 16K^2 in square class 1 and the %d characteristic-2 "
          "triangle's in class 2 (pair test, no factoring)" % (len(hh), len(th)),
          len(hh) == 4 and len(th) == 1 and all(V.same_square_class(h, 1) for h in hh)
          and all(V.same_square_class(h, hh[0]) for h in hh)
          and all(V.same_square_class(h, 2) for h in th),
          "hub=%d/%d tri=%d/%d" % (sum(V.same_square_class(h, 1) for h in hh), len(hh),
                                   sum(V.same_square_class(h, 2) for h in th), len(th)))

    return finish("public package tests / verify")


if __name__ == "__main__":
    raise SystemExit(main())
