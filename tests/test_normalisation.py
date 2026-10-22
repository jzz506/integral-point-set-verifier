#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright 2026 Zhizheng Ju
# SPDX-License-Identifier: Apache-2.0
"""Scale normalisation -- shipped tests (public code-only candidate). Run:
python3 tests/test_normalisation.py

Two "primitives" circulate in this literature under one word:

    distance normalisation  divide the C(n,2) distances by their gcd   -> the key convention used by
                            multiset-storing pipelines; the representative generally needs RATIONAL
                            coordinates, because clearing a denominator L makes every distance a
                            multiple of L
    lattice normalisation   divide the coordinate DIFFERENCES by their gcd -> the integer-cluster
                            primitive; the largest d for which the frame is d times another integer
                            frame, and the only one that is translation invariant

These legs test the arithmetic that distinguishes them. The divergence itself is a property of
join-derived frames and is NOT reproducible with the constructed controls shipped here -- that is
recorded as a coverage item, not papered over (README, "coverage reduction").
"""

from __future__ import division

import itertools
import math
from fractions import Fraction

from _harness import check, finish
from test_verify import hub_from_scratch, trapezoid_s2, tri_s2

from e213ips import similarity as S
from e213ips import verify as V

HUB, _R = hub_from_scratch()
CTRL = [(HUB, 1), (tri_s2(), 2), (trapezoid_s2(), 2), ([(0, 0), (3, 0), (0, 4)], 1)]


def gcd_all(vs):
    g = 0
    for v in vs:
        g = math.gcd(g, abs(int(v)))
    return g


def edge_gcd(P, s):
    return gcd_all(V.verify([(Fraction(p[0]), Fraction(p[1])) for p in P], s).distances.values())


def diff_gcd(P, s):
    F = [(Fraction(p[0]), Fraction(p[1])) for p in P]
    d = []
    for i in range(len(F)):
        for j in range(i + 1, len(F)):
            d += [F[i][0] - F[j][0], F[i][1] - F[j][1]]
    return gcd_all(d)


def main():
    # 1. the raw-coordinate gcd is not a property of a set: translating changes it. diff_gcd does not.
    rect = [(0, 0), (6, 0), (6, 8), (0, 8)]
    shifted = [(p[0] + 1, p[1]) for p in rect]
    raw_a, raw_b = gcd_all([v for p in rect for v in p]), gcd_all([v for p in shifted for v in p])
    check("R-NORM translation changes the raw-coordinate gcd (%d -> %d) but not diff_gcd (%d -> %d): "
          "only the difference form is a property of the set"
          % (raw_a, raw_b, diff_gcd(rect, 1), diff_gcd(shifted, 1)),
          raw_a == 2 and raw_b == 1 and diff_gcd(rect, 1) == diff_gcd(shifted, 1) == 2,
          "raw=%d/%d" % (raw_a, raw_b))

    # 2. the two diameters coincide exactly when the two gcds coincide -- recomputed, not asserted.
    rows = []
    bad = []
    for P, s in CTRL:
        eg, dg = edge_gcd(P, s), diff_gcd(P, s)
        mx = max(V.verify([(Fraction(p[0]), Fraction(p[1])) for p in P], s).distances.values())
        rows.append((eg, dg))
        if (mx // eg == mx // dg) != (eg == dg):
            bad.append((eg, dg))
    check("R-NORM the distance-normalised and lattice-normalised diameters agree exactly when "
          "edge_gcd == diff_gcd, on all %d controls (all of them have edge_gcd == diff_gcd == 1, so "
          "the divergence the internal candidate measured in 32 of 36 classes is NOT reproducible "
          "from constructed controls -- disclosed coverage gap)" % len(CTRL),
          not bad and rows == [(1, 1)] * len(CTRL), "rows=%s" % rows)

    # 3. the lattice normalisation machinery: a 7x frame really is 7x an integral set, and dividing it
    #    back out does not change the object.
    F = [(p[0] * 7, p[1] * 7) for p in HUB]
    rep = V.verify([(Fraction(p[0]), Fraction(p[1])) for p in F], 1)
    eg, dg = edge_gcd(F, 1), diff_gcd(F, 1)
    back = [(Fraction(p[0]) / 7, Fraction(p[1]) / 7) for p in F]
    check("R-NORM the hub dilated by 7 has edge_gcd=%d and diff_gcd=%d, and dividing the frame by "
          "diff_gcd returns a legal integral set with the same class_id" % (eg, dg),
          eg == 7 and dg == 7 and V.verify(back, 1).is_legal and rep.is_legal
          and S.class_id(back, 1) == S.class_id(HUB, 1),
          "eg=%d dg=%d" % (eg, dg))

    # 4. clearing denominators is exact and never rounds.
    L, ints = S.clear_denominators([Fraction(1, 3), Fraction(2, 9), Fraction(5)])
    check("R-NORM clear_denominators([1/3, 2/9, 5]) returns L=9 and [3, 2, 45]",
          L == 9 and ints == [3, 2, 45], "L=%d ints=%s" % (L, ints))

    print("-" * 72)
    print("RULE for any text that quotes a diameter: name the normalisation. On a lattice set the two "
          "differ by a factor of edge_gcd/diff_gcd, which is 1 exactly when the frame is both "
          "distance-primitive and lattice-primitive.")
    return finish("public package tests / normalisation")


if __name__ == "__main__":
    raise SystemExit(main())
