# SPDX-FileCopyrightText: Copyright 2026 Zhizheng Ju
# SPDX-License-Identifier: Apache-2.0
"""Exact similarity classification of integral point sets.

Two point sets are SIMILAR iff a vertex permutation and one positive scale carry one distance matrix
onto the other. For sets with integral distances the scale is necessarily rational: if the squared
distance matrices are proportional then lambda = d'^2 / d^2 = (d'/d)^2 for any pair, a rational square,
and the true ratio d'/d is rational. (Conversely a non-square rational lambda, e.g. 2, maps an integral
distance matrix out of the domain: 2*d^2 is never a square for d >= 1 -- tests/test_similarity.py
checks that on every edge of the planted control sets.) So no algebraic scale handling is needed, and dilations are the same
class -- which is the point of this module, because the frame generators emit dilated copies and those
must not be counted as discoveries.

PUBLIC API
    similarity_signature(points, s=1)  -> tuple of Fraction   complete invariant, see below
    class_id(points_or_sig, s=1)       -> str                 "H<n>-" + sha256(signature)[:16]
    distance_gcd_key(points, s=1)      -> (s, tuple)          dilation-normalised multiset, NOT complete
    are_similar(A, B, s=1)             -> (bool, lambda)      exact squared scale factor
    linear_scale(lam)                  -> Fraction | None     sqrt of a rational square, else None

COMPLETENESS, and why both keys ship.
    `similarity_signature` is the complete invariant: the lexicographic minimum over ALL n! vertex
    labelings of the C(n,2) squared distances in upper-triangular order, divided by its own minimum.
    Minimum-then-divide and divide-then-minimise agree because the minimum is permutation-invariant and
    multiplying every entry by one positive constant preserves lexicographic order -- so the search runs
    on plain integers and only the winner is normalised.
    `distance_gcd_key` is the SORTED distance multiset divided by its gcd. It normalises dilation,
    which is why multiset-storing pipelines use it, but it is NOT a similarity proof: distinct
    similarity classes can share a distance multiset. The internal candidate measured 2 such collisions
    among 1,260 four-point subsets of its 36-class benchmark; that data is not shipped here, so in this
    package the caveat is DOCUMENTED, not re-measured (README, "coverage reduction").
    Never key a population on the multiset.

Cost: n! labelings -- 5,040 for n=7, 40,320 for n=8, each a tuple comparison on ints.
"""

from __future__ import division

import hashlib
import itertools
import math
from fractions import Fraction

from .verify import _num, dist2

__all__ = ["squared_distances", "clear_denominators", "similarity_signature", "class_id",
           "distance_gcd_key", "are_similar", "linear_scale"]


def _pair_order(n):
    return list(itertools.combinations(range(n), 2))


def squared_distances(points, s=1):
    """The C(n,2) squared distances in upper-triangular (i<j) order, exact Fractions."""
    n = len(points)
    return [dist2(points[a], points[b], s)
            for a, b in _pair_order(n)]


def clear_denominators(values):
    """(L, [ints]): the lcm of the denominators and the scaled integers. Never a float multiply."""
    L = 1
    for v in values:
        v = _num(v)
        L = L * v.denominator // math.gcd(L, v.denominator)
    return L, [int(v * L) for v in values]


def distance_matrix(points, s=1):
    """Symmetric integer matrix of CLEARED squared distances, plus the clearing denominator L."""
    n = len(points)
    pairs = _pair_order(n)
    L, ints = clear_denominators(squared_distances(points, s))
    D2 = [[0] * n for _ in range(n)]
    for (a, b), v in zip(pairs, ints):
        D2[a][b] = D2[b][a] = v
    return D2, L


def _canon_winner(D2, n):
    """The lex-min tuple over all n! labelings, still in integers."""
    pairs = _pair_order(n)
    best = None
    for perm in itertools.permutations(range(n)):
        cand = tuple(D2[perm[a]][perm[b]] for a, b in pairs)
        if best is None or cand < best:
            best = cand
    return best


def _canon_from_integer_matrix(D2, n):
    pairs = _pair_order(n)
    m = min(D2[a][b] for a, b in pairs)
    if m <= 0:
        raise ValueError("a similarity signature needs strictly positive distances (repeated point?)")
    return tuple(Fraction(v, m) for v in _canon_winner(D2, n))


def similarity_signature(points, s=1):
    """Complete similarity signature: ratios of squared distances, invariant under relabeling and
    under ONE positive rational scale. Sorted distances are NOT enough -- module docstring."""
    n = len(points)
    if n < 3:
        raise ValueError("need at least 3 points for a signature")
    D2, _L = distance_matrix(points, s)
    return _canon_from_integer_matrix(D2, n)


def signature_digest(sig):
    return hashlib.sha256(("|".join(str(x) for x in sig)).encode()).hexdigest()


def _n_from_pair_count(m):
    """Invert m = C(n,2) -> n; used only when a bare signature is handed to class_id."""
    n = 2
    while (n * (n - 1)) // 2 < m:
        n += 1
    if (n * (n - 1)) // 2 != m:
        raise ValueError("%d is not C(n,2) for any n" % m)
    return n


def class_id(points_or_sig, s=1):
    """'H<n>-' + the first 16 hex digits of sha256 of the signature.

    Accepts a point list or an already-computed signature tuple. The encoding is
    '|'.join(str(Fraction)) -- so "1" for the minimum entry and "25/4" style for ratios -- which makes
    the digest a function of the signature alone, reproducible across runs and hash seeds.
    tests/test_similarity.py pins the exact ids of the planted control sets as known answers.
    """
    seq = list(points_or_sig)
    if seq and isinstance(seq[0], (tuple, list)):
        pts = seq
        sig = similarity_signature(pts, s)
        n = len(pts)
    else:
        sig = tuple(Fraction(x) for x in seq)
        n = _n_from_pair_count(len(sig))
    return "H%d-%s" % (n, signature_digest(sig)[:16])


def distance_gcd_key(points, s=1):
    """(s, sorted DISTANCES after clearing denominators and dividing by their gcd).

    A dilation-normalised multiset key. It discards which distance belongs to which pair, so it
    identifies a set only up to its multiset -- see the module docstring. Use `class_id` for identity
    and this when reconciling with a pipeline that stores multisets.
    """
    D2, _L = distance_matrix(points, s)
    n = len(points)
    roots = []
    for a, b in _pair_order(n):
        r = math.isqrt(D2[a][b])
        if r * r != D2[a][b]:
            return (s, None)                   # not integral under this characteristic
        roots.append(r)
    g = 0
    for v in roots:
        g = math.gcd(g, v)
    return (s, tuple(sorted(v // g for v in roots)))


def linear_scale(lam):
    """sqrt of a rational if it is a rational square, else None. lambda=9/4 -> 3/2; lambda=2 -> None."""
    lam = _num(lam)
    if lam <= 0:
        return None
    rn, rd = math.isqrt(lam.numerator), math.isqrt(lam.denominator)
    if rn * rn == lam.numerator and rd * rd == lam.denominator:
        return Fraction(rn, rd)
    return None


def are_similar(A, B, s=1):
    """(bool, lambda): similar, and the exact ratio of the cleared squared-distance matrices.

    Similarity is decided by the complete signatures, which quantify over all labelings. lambda is then
    read off the canonical winners themselves -- both winner tuples are indexed by the same pair order
    and both realise the same canonical form, so their entrywise ratio IS the scale correspondence.
    Comparing the raw matrices pair-by-pair instead would need A and B to arrive already aligned, which
    is exactly what a caller cannot assume.
    """
    if len(A) != len(B):
        return False, None
    n = len(A)
    Da, _ = distance_matrix(A, s)
    Db, _ = distance_matrix(B, s)
    wa, wb = _canon_winner(Da, n), _canon_winner(Db, n)
    if min(wa) <= 0 or min(wb) <= 0:
        return False, None
    sig_a = tuple(Fraction(v, min(wa)) for v in wa)
    sig_b = tuple(Fraction(v, min(wb)) for v in wb)
    if sig_a != sig_b:
        return False, None
    lam = Fraction(wb[0], wa[0])
    for x, y in zip(wa, wb):
        if _num(x) * lam != y:
            return False, None                 # unreachable given equal signatures; kept as a self-check
    return True, lam
