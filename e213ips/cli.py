#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright 2026 Zhizheng Ju
# SPDX-License-Identifier: Apache-2.0
"""Command line: python3 -m e213ips.cli <input.json> [--s N] [--quiet]

<input.json> is either a list of [x, t] pairs, or an object with "points" and optionally "s".
Coordinates may be integers, strings, or rationals written as "123/45" -- all parsed with Fraction, so
nothing rounds through a float.

Prints the verifier's report (including WHICH pairs are non-integral, which triples are collinear and
which quads are concyclic), then the similarity class_id. Exit code 0 if the set is legal, 1 otherwise,
3 if the input could not be read.
"""

from __future__ import division

import itertools
import json
import sys
from fractions import Fraction

from . import similarity as S
from . import verify as V


def load(path):
    with open(path) as f:
        doc = json.load(f)
    if isinstance(doc, list):
        raw, s = doc, 1
    else:
        raw, s = doc["points"], int(doc.get("s", 1))
    return [(Fraction(str(x)), Fraction(str(t))) for x, t in raw], s


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    quiet = "--quiet" in argv
    if quiet:
        argv.remove("--quiet")
    s_override = None
    if "--s" in argv:
        i = argv.index("--s")
        s_override = int(argv[i + 1])
        del argv[i:i + 2]
    if not argv:
        print(__doc__)
        return 3
    path = argv[0]
    try:
        pts, s = load(path)
    except (OSError, ValueError, KeyError) as e:
        print("cannot read %s: %s" % (path, e))
        return 3
    if s_override is not None:
        s = s_override
    rep = V.verify(pts, s)
    print("input      : %s" % path)
    print("n, s       : %d, %d" % (rep.n, rep.s))
    print("distances  : %d of %d pairs integral" % (rep.integral_distance_count, rep.total_pairs))
    if rep.non_integral_pairs:
        print("             non-integral at %s" % rep.non_integral_pairs[:8])
    print("collinear  : %d of %d triples" % (len(rep.collinear_triples),
                                             len(list(itertools.combinations(range(rep.n), 3)))))
    print("concyclic  : %d of %d quads%s"
          % (len(rep.concyclic_quads),
             len(list(itertools.combinations(range(rep.n), 4))),
             (" at %s" % rep.concyclic_quads[:8]) if rep.concyclic_quads else ""))
    print("verdict    : %s" % ("LEGAL" if rep.is_legal else "ILLEGAL -- " + "; ".join(rep.reasons)))
    if rep.integral_distance_count != rep.total_pairs:
        print("class_id   : not computed (some pair is not at an integral distance)")
    elif len(pts) < 3:
        print("class_id   : not computed (a similarity signature needs at least 3 points)")
    else:
        print("class_id   : %s" % S.class_id(pts, s))
    return 0 if rep.is_legal else 1


if __name__ == "__main__":
    raise SystemExit(main())
