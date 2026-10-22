# SPDX-FileCopyrightText: Copyright 2026 Zhizheng Ju
# SPDX-License-Identifier: Apache-2.0
"""Exact verification of integral point sets in the square-root coordinate model.

    True point          P = (x, t * sqrt(s))        with x, t rational and s a squarefree integer >= 1
    Squared distance    d(P,Q)^2 = (x1-x2)^2 + s * (t1-t2)^2
    Stored pair         (x, t)                      -- s is shared by the whole set, never per point

The stored coordinate is the SECOND STORAGE coordinate t, not the true y. Keeping t rational is what
makes every predicate below an integer / Fraction computation: the true y-coordinate is irrational
whenever s > 1, so there is no way to write these sets down in the plane they actually live in.

MODEL_NOTES -- three facts a reader must not have to rediscover:

1. Collinearity in the model plane equals collinearity in the true plane. The map
   (x, t) -> (x, t*sqrt(s)) is an invertible linear map (a scaling of one axis), and invertible
   linear maps preserve lines. So the 2x2 cross-product test can be run on the stored pairs.

2. Concyclicity does NOT transfer that way, and the naive circle equation is wrong here. Circles in
   the true plane pull back to the conics x^2 + s*t^2 + D*x + E*t + F = 0, so the determinant to test
   is the weighted one, with first column x^2 + s*t^2. `concyclic` implements exactly that -- the same
   weighted determinant used across this problem's literature.

3. The characteristic s is a property of the distance matrix, not of the chosen coordinates
   (Kemnitz's shared-characteristic theorem): every non-degenerate triangle of the set has the same
   square class for 16*area^2. `heron_16k2` and `same_square_class` recompute it from distances alone,
   so a reader can check a stated s against the coordinates that produced it. A set verified at the
   wrong s is not the same geometric object -- tests/test_verify.py reads the characteristic-2 control
   set at s=1 and back, and requires the verifier to refuse the mismatched pairing.

Public API
    dist2(p, q, s=1)              -> Fraction
    integral_distance(p, q, s=1)  -> int | None
    collinear(a, b, c)            -> bool
    concyclic(points4, s=1)       -> bool
    det4(points4, s=1)            -> int              (the exact determinant; 0 <=> concyclic)
    heron_16k2(a, b, c)           -> int | Fraction   (16*area^2 from the three distances)
    is_square_rational(v)         -> bool            (v is a square in Q)
    same_square_class(u, v)         -> bool            (u, v differ by a rational square)
    verify(points, s=1)           -> Report
"""

from __future__ import division

import itertools
import math
from fractions import Fraction

__all__ = ["dist2", "integral_distance", "collinear", "concyclic", "det4", "det_bareiss",
           "det_laplace", "heron_16k2", "is_square_rational", "same_square_class",
           "verify", "Report", "MODEL_NOTES"]

Point = tuple  # (x, t) with x, t in int | Fraction | str-parseable


def _num(v):
    """Accept ints, Fractions and the strings the JSON data uses; nothing else."""
    if isinstance(v, int):
        return Fraction(v)
    if isinstance(v, Fraction):
        return v
    if isinstance(v, str):
        return Fraction(v)
    if isinstance(v, float):
        raise TypeError("float coordinate %r is refused: this package decides with exact "
                        "arithmetic only (Fraction(int(v)) would silently change the answer)" % v)
    raise TypeError("unsupported coordinate type %s" % type(v).__name__)


def dist2(p, q, s=1):
    """Squared distance between two stored points, as an exact Fraction."""
    x1, t1 = _num(p[0]), _num(p[1])
    x2, t2 = _num(q[0]), _num(q[1])
    return (x1 - x2) ** 2 + s * (t1 - t2) ** 2


def _isqrt_fraction(v):
    """sqrt(v) as an int when v is a perfect square of a rational and the root is an integer, else None."""
    if v <= 0:
        return None
    n, d = v.numerator, v.denominator
    rn, rd = math.isqrt(n), math.isqrt(d)
    if rn * rn != n or rd * rd != d:
        return None
    root = Fraction(rn, rd)
    return int(root) if root.denominator == 1 else None


def integral_distance(p, q, s=1):
    """The distance itself if it is a positive integer, else None. Never a float."""
    return _isqrt_fraction(dist2(p, q, s))


def collinear(a, b, c):
    """True iff the three stored points lie on one line (see MODEL_NOTES 1)."""
    ax, at = _num(a[0]), _num(a[1])
    bx, bt = _num(b[0]), _num(b[1])
    cx, ct = _num(c[0]), _num(c[1])
    return (bx - ax) * (ct - at) == (bt - at) * (cx - ax)


def det_bareiss(matrix):
    """Exact determinant of a square integer matrix, fraction-free (Bareiss), with row pivoting.

    Bareiss keeps every intermediate an integer and divides exactly by the previous pivot, so the
    result is the true determinant with no gcd or coefficient blow-up on the 4x4 case that matters.
    A 4x4 determinant is translation-invariant under column operations and scales by lambda^4 under a
    uniform coordinate scaling, which is why clearing denominators one row at a time (a per-row factor)
    can change the magnitude but never the zero set.

    Two notes from this package's build. (1) Bareiss divides exactly by the previous pivot; that is
    asserted here rather than trusted, because a floored ``//`` on a non-exact division would print a
    plausible wrong determinant instead of failing. (2) A first version of this loop transcribed the
    second term as A[i][k]*A[i][k] instead of A[i][k]*A[k][j] and returned a 53-digit nonzero value for a
    quadruple that is concyclic; comparing this function against ``det_laplace`` -- a second
    implementation, by definition -- caught it, and tests/test_verify.py keeps that comparison as a
    shipped leg.
    """
    n = len(matrix)
    A = [list(r) for r in matrix]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            pivot = None
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    pivot = i
                    break
            if pivot is None:
                return 0          # trailing block has a zero column => singular
            A[k], A[pivot] = A[pivot], A[k]
            sign = -sign
        akk = A[k][k]
        for i in range(k + 1, n):
            aik = A[i][k]
            for j in range(k + 1, n):
                num = A[i][j] * akk - aik * A[k][j]
                q, r = divmod(num, prev)
                if r:
                    raise ArithmeticError("Bareiss step %d: %d / %d is not exact (remainder %d)"
                                          % (k, num, prev, r))
                A[i][j] = q
            A[i][k] = 0
        prev = A[k][k]
    return sign * A[n - 1][n - 1] if n else sign


def det_laplace(matrix):
    """Second implementation, by definition (Laplace along row 0) -- the cross-check for det_bareiss."""
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    out = 0
    for j in range(n):
        if matrix[0][j] == 0:
            continue
        minor = [[matrix[i][k] for k in range(n) if k != j] for i in range(1, n)]
        out += ((-1) ** j) * matrix[0][j] * det_laplace(minor)
    return out


def _row_clear(row):
    """Multiply a rational row by the lcm of its denominators -> integer row (nonzero factor)."""
    L = 1
    for v in row:
        L = L * v.denominator // math.gcd(L, v.denominator)
    return [int(v * L) for v in row]


def det4(points4, s=1):
    """The weighted concyclicity determinant det[[x^2+s t^2, x, t, 1]] as an exact integer.

    Zero <=> the four points are concyclic (or degenerate: three of them collinear makes the same
    determinant vanish, which `verify` reports separately so a caller can tell the two apart).
    """
    if len(points4) != 4:
        raise ValueError("det4 takes exactly 4 points")
    rows = []
    for p in points4:
        x, t = _num(p[0]), _num(p[1])
        rows.append(_row_clear([x * x + s * t * t, x, t, Fraction(1)]))
    return det_bareiss(rows)


def concyclic(points4, s=1):
    return det4(points4, s) == 0


def heron_16k2(a, b, c):
    """16 * area^2 of the triangle with sides a, b, c, from the distances alone: (a+b+c)(-a+b+c)(a-b+c)(a+b-c)."""
    a, b, c = _num(a), _num(b), _num(c)
    return (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c)


def is_square_rational(v):
    """True iff the positive rational v is a square in Q: both numerator and denominator perfect squares."""
    v = _num(v)
    if v <= 0:
        return False
    rn, rd = math.isqrt(v.numerator), math.isqrt(v.denominator)
    return rn * rn == v.numerator and rd * rd == v.denominator


def same_square_class(u, v):
    """True iff u and v differ by a rational square, i.e. lie in ONE square class.

    This is the test the characteristic claim needs, and it is deliberately not implemented as
    "factor both and compare squarefree parts": the Heron values here run to 40+ digits, and a
    squarefree representative would require factoring them. Membership in a square class is a pair test
    on a single division, exact and fast.
    """
    u, v = _num(u), _num(v)
    if u <= 0 or v <= 0:
        return False
    return is_square_rational(Fraction(u, v))


class Report(object):
    """What `verify` decided, with the witnesses, not just the verdict."""

    def __init__(self, **kw):
        self.__dict__.update(kw)

    @property
    def is_legal(self):
        return not self.reasons

    def as_dict(self):
        d = dict(self.__dict__)
        d["is_legal"] = self.is_legal
        d["collinear_triples"] = [list(t) for t in d["collinear_triples"]]
        d["concyclic_quads"] = [list(q) for q in d["concyclic_quads"]]
        d["non_integral_pairs"] = [list(p) for p in d["non_integral_pairs"]]
        return d


def verify(points, s=1, require_integral=True, require_general_position=True):
    """Decide a point set against the definitions, exactly.

    Legal means: n distinct points; every pairwise distance a positive integer (when
    require_integral); no three collinear and no four concyclic (when require_general_position).
    The Report carries the offending index tuples, so a caller can see WHY.
    """
    pts = [(p[0], p[1]) for p in points]
    n = len(pts)
    seen = set()
    dup = []
    for i, p in enumerate(pts):
        key = (_num(p[0]), _num(p[1]))
        if key in seen:
            dup.append(i)
        seen.add(key)

    pairs = list(itertools.combinations(range(n), 2))
    distances = {}
    non_integral = []
    for a, b in pairs:
        root = integral_distance(pts[a], pts[b], s)
        if root is None:
            non_integral.append((a, b))
        else:
            distances[(a, b)] = root

    triples, quads = [], []
    if n >= 3:
        for t in itertools.combinations(range(n), 3):
            if collinear(pts[t[0]], pts[t[1]], pts[t[2]]):
                triples.append(t)
    if n >= 4:
        for q in itertools.combinations(range(n), 4):
            if det4([pts[i] for i in q], s) == 0:
                quads.append(q)

    reasons = []
    if dup:
        reasons.append("repeated point at index %s" % dup)
    if require_integral and non_integral:
        reasons.append("%d of %d pairs at non-integral distance" % (len(non_integral), len(pairs)))
    if require_general_position and triples:
        reasons.append("%d collinear triple(s)" % len(triples))
    if require_general_position and quads:
        reasons.append("%d concyclic quad(s)" % len(quads))

    return Report(n=n, s=s, total_pairs=len(pairs), integral_distance_count=len(distances),
                  distances={k: v for k, v in distances.items()},
                  non_integral_pairs=non_integral, collinear_triples=triples,
                  concyclic_quads=quads, reasons=reasons)
