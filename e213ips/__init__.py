# SPDX-FileCopyrightText: Copyright 2026 Zhizheng Ju
# SPDX-License-Identifier: Apache-2.0
"""e213ips — exact verification and similarity classification of integral point sets.

Public code-only release candidate. Standard library only, Python 3.9+. Every decision path uses
exact integer / Fraction arithmetic; no float appears on one, and no geometric predicate is
approximated. This package ships code, tests and constructed controls -- no coordinate benchmark.

Modules
    verify      -- integral distances, no-three-collinear, no-four-concyclic, from the definitions.
    similarity  -- complete similarity signature (class_id) and a dilation-normalised multiset key.

The coordinate model is documented in verify.MODEL_NOTES; read it before using the results.
"""

__version__ = "0.1.0"
__all__ = ["verify", "similarity", "__version__"]
