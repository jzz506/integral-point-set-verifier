#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright 2026 Zhizheng Ju
# SPDX-License-Identifier: Apache-2.0
"""Run this package's test suite: python3 tests/run_all.py

Exit code 0 only if every leg passes. Totalled from the shared harness, so the number below is the sum
of the three files' legs and not a re-count.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import _harness                      # noqa: E402
from _harness import finish          # noqa: E402

import test_verify                   # noqa: E402
import test_similarity               # noqa: E402
import test_normalisation            # noqa: E402

rc = 0
for mod in (test_verify, test_similarity, test_normalisation):
    rc |= mod.main()
    print()

print("package tests TOTAL: %d/%d legs passed"
      % (_harness._state["pass"], _harness._state["pass"] + _harness._state["fail"]))
if _harness._state["fail"]:
    print("FAILURES: " + "; ".join(_harness._state["failures"]))
raise SystemExit(rc)
