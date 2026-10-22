# SPDX-FileCopyrightText: Copyright 2026 Zhizheng Ju
# SPDX-License-Identifier: Apache-2.0
"""A deliberately small check harness: named legs, both ways, honest tally.

Every leg states what it can indict. A leg that only shows the shipped code passing proves nothing --
the rectangle and collinear-triple legs exist to prove the detectors can say NO, and the
Bareiss-vs-Laplace leg exists because it is the one that caught a real bug in this package during the
build (a transcribed matrix factor, see e213ips/verify.py:det_bareiss).

This is this package's own test suite. It is not a project-wide regression gate, and its tally must
never be quoted as one: it covers the code in this directory and nothing else.
"""

import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_state = {"pass": 0, "fail": 0, "failures": [], "mark": 0}


def check(name, condition, detail=""):
    if condition:
        _state["pass"] += 1
        print("PASS %s%s" % (name, ("  [%s]" % detail) if detail else ""))
        return 1
    _state["fail"] += 1
    _state["failures"].append(name)
    print("FAIL %s  [%s]" % (name, detail))
    return 0


def expect_raises(name, exc, fn, *a, **kw):
    try:
        fn(*a, **kw)
    except exc:
        return check(name, True, exc.__name__)
    except Exception as e:                                    # noqa: BLE001 - the point is to NAME it
        return check(name, False, "raised %s, expected %s" % (type(e).__name__, exc.__name__))
    return check(name, False, "no exception raised")


def finish(label):
    """Report THIS file's legs and the running total, so a per-file number can never be misread as the
    suite's (the first version counted cumulatively and printed '31/31' under a 7-leg file)."""
    total = _state["pass"] + _state["fail"]
    section = total - _state["mark"]
    passed_here = _state["pass"] - _state["mark"]
    _state["mark"] = total
    print("-" * 72)
    print("%s: %d/%d checks passed  (suite so far %d/%d)%s"
          % (label, passed_here, section, _state["pass"], total,
             "" if not _state["failures"] else " -- FAILED: " + "; ".join(_state["failures"])))
    return 1 if _state["fail"] else 0


def load_data(name):
    import json
    with open(os.path.join(_ROOT, "data", name)) as f:
        return json.load(f)


def pts(class_row, coord_field="witness_coordinates"):
    """JSON stores coordinates as strings so nothing can round through a float on the way in."""
    from fractions import Fraction
    return [(Fraction(x), Fraction(y)) for x, y in class_row[coord_field]]
