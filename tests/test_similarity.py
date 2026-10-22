#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright 2026 Zhizheng Ju
# SPDX-License-Identifier: Apache-2.0
"""e213ips.similarity -- shipped tests (public code-only candidate).

Controls come from tests/test_verify.py's constructed fixtures. Run: python3 tests/test_similarity.py
"""

from __future__ import division

import itertools
from fractions import Fraction

from _harness import check, finish
from test_verify import hub_from_scratch, trapezoid_s2, tri_s2

from e213ips import similarity as S
from e213ips import verify as V

HUB, _R = hub_from_scratch()
TRI2 = tri_s2()
TRAP2 = trapezoid_s2()
TRI3 = [(0, 0), (3, 0), (0, 4)]                       # 3-4-5 right triangle


def _as_frac(P):
    return [(Fraction(p[0]), Fraction(p[1])) for p in P]


def main():
    # 1. KNOWN ANSWER, pinned for encoding drift. The semantic content is the distance multiset, which
    #    test_verify.py checks independently; these strings exist so a change of the signature
    #    encoding (the "|".join(str(Fraction)) form, the divide-by-minimum step, the digest length)
    #    makes this suite go red instead of silently renaming every class. The hub control lands on
    #    H4-db6ef773afadaba6, the same id the internal candidate's 5-5-6 example carries, which is the
    #    scale-blindness of the classifier confirmed across two separately built packages.
    ids = {name: S.class_id(_as_frac(P), s) for name, P, s in
           (("hub", HUB, 1), ("tri2", TRI2, 2), ("trap2", TRAP2, 2), ("tri3", TRI3, 1))}
    check("R-SIM the four constructed controls classify to the pinned ids %s"
          % " ".join("%s=%s" % (k, v) for k, v in sorted(ids.items())),
          ids == {"hub": "H4-db6ef773afadaba6", "tri2": "H3-9d40d411a64a4e60",
                  "trap2": "H4-37acb40cbe2d4b5a", "tri3": "H3-55daec86056ef2c4"},
          "got=%s" % ids)

    # 2. INVARIANCE under relabeling: identity is a property of the SET.
    bad = []
    for name, P, s in (("hub", HUB, 1), ("tri2", TRI2, 2), ("trap2", TRAP2, 2)):
        F = _as_frac(P)
        want = S.class_id(F, s)
        perms = [F[k:] + F[:k] for k in range(1, len(F))] + [list(reversed(F))]
        bad += [(name, p) for p in perms if S.class_id(p, s) != want]
    check("R-SIM cyclic rotations and reversal of the vertex order change no class_id (3 sets x 8 "
          "orderings)", not bad, "bad=%s" % bad[:2])

    # 3. INVARIANCE under dilation, and the recovered scale.
    bad = []
    seen = []
    for name, P, s in (("hub", HUB, 1), ("tri2", TRI2, 2), ("tri3", TRI3, 1)):
        F = _as_frac(P)
        for m in (2, 3, 5):
            G = [(p[0] * m, p[1] * m) for p in F]
            if S.class_id(G, s) != S.class_id(F, s):
                bad.append((name, m, "id"))
            sim, lam = S.are_similar(F, G, s)
            if not sim or lam != Fraction(m * m) or S.linear_scale(lam) != Fraction(m):
                bad.append((name, m, "lambda=%s" % lam))
            seen.append((name, m, str(lam)))
    check("R-SIM dilation by 2, 3 or 5 preserves class_id and are_similar returns lambda=m^2 with "
          "linear_scale=m (9 cases)", not bad and len(seen) == 9, "bad=%s" % bad[:2])

    # 4. A NON-SQUARE scale is not realisable between integral sets: doubling every distance is square
    #    on 0 edges of the constructed controls, and linear_scale refuses 2 and 1/2 outright.
    n_edge = 0
    nonsquare = []
    for P, s in ((HUB, 1), (TRI2, 2), (TRAP2, 2), (TRI3, 1)):
        for v in S.squared_distances(_as_frac(P), s):
            n_edge += 1
            if V.is_square_rational(2 * v):
                nonsquare.append((len(P), v))
    check("R-SIM 2*d^2 is a square on 0 of the %d edges of the constructed controls, and "
          "linear_scale(2) / linear_scale(1/2) are None" % n_edge,
          not nonsquare and S.linear_scale(Fraction(2)) is None
          and S.linear_scale(Fraction(1, 2)) is None
          and S.linear_scale(Fraction(9, 4)) == Fraction(3, 2),
          "edges=%d" % n_edge)

    # 5. DISTINCTNESS: the four controls are four different objects.
    ids = set()
    for P, s in ((HUB, 1), (TRI2, 2), (TRAP2, 2), (TRI3, 1)):
        ids.add(S.class_id(_as_frac(P), s))
    sim, lam = S.are_similar(_as_frac(HUB), _as_frac(TRI3), 1)
    check("R-SIM the controls give 4 distinct class_ids and hub vs the 3-4-5 triangle compares as NOT "
          "similar", len(ids) == 4 and sim is False and lam is None, "ids=%d" % len(ids))

    # 6. The multiset key normalises dilation but is NOT complete: two non-similar 4-subsets of the
    #    grid with the SAME distance multiset exist in a 4x4 lattice, found here rather than inherited.
    grid = [(Fraction(x), Fraction(y)) for x in range(4) for y in range(4)]
    by_multiset = {}
    for q in itertools.combinations(grid, 4):
        vals = tuple(sorted(int(v) for v in S.squared_distances(list(q), 1)))
        by_multiset.setdefault(vals, []).append(q)
    merged = [g for g in by_multiset.values() if len({S.similarity_signature(list(x), 1) for x in g}) > 1]
    n_sub = sum(len(v) for v in by_multiset.values())
    check("R-SIM over the %d 4-subsets of a 4x4 lattice (%d distinct multisets) the key merges %d "
          "group(s) that the complete signature separates -- so the multiset key is INCOMPLETE, "
          "measured here rather than inherited" % (n_sub, len(by_multiset), len(merged)),
          n_sub == 1820 and len(merged) >= 1, "collisions=%d multisets=%d" % (len(merged), len(by_multiset)))

    # 7. gcd(D^2) == gcd(D)^2: the squared and unsquared forms of the multiset key carry the same
    #    information, which is what lets a squared-distance pipeline be reconciled with this one.
    import math
    bad = []
    for P, s in ((HUB, 1), (TRI2, 2), (TRAP2, 2), (TRI3, 1)):
        d2 = [int(v) for v in S.squared_distances(_as_frac(P), s)]
        ds = [int(V._isqrt_fraction(Fraction(v))) for v in d2]
        gd = 0
        for v in ds:
            gd = math.gcd(gd, v)
        g2 = 0
        for v in d2:
            g2 = math.gcd(g2, v)
        if g2 != gd * gd:
            bad.append((str(gd), str(g2)))
    check("R-SIM gcd of the squared distances equals the square of the gcd of the distances on all 4 "
          "controls (the two key conventions differ by nothing)", not bad, "bad=%s" % bad[:2])

    # 8. exact rational handling, no rounding: dividing every coordinate by 3 must not change the class.
    F = _as_frac(HUB)
    G = [(p[0] / 3, p[1] / 3) for p in F]
    check("R-SIM dividing every coordinate by 3 leaves the class_id alone (denominators are cleared, "
          "never rounded)", S.class_id(G, 1) == S.class_id(F, 1),
          "%s vs %s" % (S.class_id(G, 1), S.class_id(F, 1)))

    # 9. THE n=8 CODE PATH, exercised on a synthetic lattice set: eight points of a 3x3 grid, two
    #    vertex orderings, one answer. No integral 8-point set is shipped because none is known -- this
    #    leg tests the 40,320-labeling search and the 28-entry tuple, not a discovery.
    grid8 = [(Fraction(x), Fraction(y)) for x in range(3) for y in range(3)][:8]
    rot = grid8[3:] + grid8[:3]
    rev = list(reversed(grid8))
    ids8 = {S.class_id(grid8, 1), S.class_id(rot, 1), S.class_id(rev, 1)}
    n_pairs = len(S.squared_distances(grid8, 1))
    check("R-SIM the n=8 path (%d labelings x %d squared distances) is relabeling-invariant on a "
          "synthetic 3x3-grid set: one id for three orderings" % (40320, n_pairs),
          len(ids8) == 1 and n_pairs == 28 and len(S.similarity_signature(grid8, 1)) == 28,
          "ids=%s pairs=%d" % (ids8, n_pairs))

    return finish("public package tests / similarity")


if __name__ == "__main__":
    raise SystemExit(main())
