# PROVENANCE and rights — public code-only candidate

Written 2026-10-22. This file exists so that a recipient can answer, for every byte in the archive,
"who wrote it and on whose authority may I use it".

## 1. Code: original to this campaign, and provably unchanged from the verified sibling

Three modules implement the mathematics. Each is **identical to the corresponding file of this
project's internal candidate once module/class/function docstrings are removed** — the public copies
differ only in comments, which were rewritten so no internal path, dataset name or round identifier
appears in shipped code. That claim is machine-checked, not asserted:

| file | executable AST identical to the internal sibling | size of the normalised AST |
|---|---|---|
| `e213ips/verify.py` | **True** | 33,280 chars |
| `e213ips/similarity.py` | **True** | 20,867 chars |
| `e213ips/cli.py` | True except one disclosed message fix, below | 10,224 chars |
| `e213ips/__init__.py` | differs only in `__version__` (`0.1.0-rc1` → `0.1.0`) | — |

**One public-only divergence, found by the clean-environment battery.** Given a set with fewer than
three points, the internal CLI printed `class_id : not computed (needs every distance integral)` even
when every distance *was* integral -- a false reason for a true verdict. The public copy now prints the
actual cause (`a similarity signature needs at least 3 points`, or the non-integral case). It touches a
printed string only, no decision path; `verify.py` and `similarity.py` remain byte-for-behaviour
identical to the internal siblings.

The check: parse both files, pop the first `Expr(Constant(str))` of every `Module`/`FunctionDef`/
`ClassDef`, compare `ast.dump`. Reproduce with:

```python
import ast
def strip(t):
    for n in ast.walk(t):
        if isinstance(n, (ast.Module, ast.FunctionDef, ast.ClassDef)):
            b = n.body
            if b and isinstance(b[0], ast.Expr) and isinstance(b[0].value, ast.Constant) \
               and isinstance(b[0].value.value, str):
                b.pop(0)
    return t
a = strip(ast.parse(open('fileA.py').read())); b = strip(ast.parse(open('fileB.py').read()))
print(ast.dump(a) == ast.dump(b))
```

**Third-party code: none.** No package is imported beyond the standard library (`__future__`, `fractions`, `math`,
`itertools`, `hashlib`, `json`, `os`, `sys`). No source was extracted from any paper, repository or
dataset. The algorithms implemented are the standard ones for this problem — the weighted
concyclicity determinant and the lexicographically minimal squared-distance signature — each written
here from the definitions, with a second independent implementation of the determinant and of the
circle test inside the test suite so that either can indict the other.

One bug found this way is documented in `e213ips/verify.py`'s `det_bareiss` docstring and in
`tests/test_verify.py`'s `circumcentre` helper: an early draft of the helper solved for the centre
using un-shifted differences (`|b|² − |a|²`) where the offset form needs `|b − a|²`, and the
grid-agreement leg caught it as 157 disagreements. The shipped verifier was not affected.

## 2. Tests and examples: constructed, with the arithmetic shown

| file | content | how it was obtained | rights |
|---|---|---|---|
| `tests/test_verify.py` | 14 legs | derived in-file: the hub 4-set comes from `R = abc/4K` on a 40-40-48 triangle; the characteristic-2 triangle and trapezoid from hand-solved `Δx² + 2Δt² = square`; the rectangle, the collinear triple and `√2` are toys | constructed here |
| `tests/test_similarity.py` | 9 legs | same fixtures, plus 4×4-lattice enumerations | constructed here |
| `tests/test_normalisation.py` | 4 legs | same fixtures | constructed here |
| `tests/_harness.py`, `tests/run_all.py` | test reporting | written here | original |
| `examples/hub_4set_legal.json` | legal `s=1` set | the circumradius construction above | constructed here |
| `examples/trapezoid_s2_illegal.json` | integral, illegal-on-one-circle `s=2` set | hand-solved norm-9 / norm-36 vectors | constructed here |

Nothing in this candidate is a transcription of someone else's coordinates. The published literature is
*referenced* in prose (Kreisel–Kurz arXiv:0804.1303; Kurz–Noll–Rathbun–Simmons arXiv:1312.2318;
Kurz–Wassermann) for definitions and for the one constant quoted in the README, and no figure or
dataset is copied from them.

## 3. Excluded on purpose

* the 36-class seven-point coordinate benchmark, and the two eight-point near-miss classes;
* any dataset whose attribution this project has not verified, and anything read from a published
  coordinate collection whose licence terms were never checked;
* enumeration engines, sieves, certifiers, cost models, coverage certificates, ledgers and caches;
* every private path, round identifier, username, absolute filesystem path and credential — verified
  by a scan over the package for the repository directory name, for a home-directory prefix, for
  the internal dataset and module names, and for credential-shaped tokens: no match in shipped code,
  tests, examples or documents;
* the licence is NOT excluded: Apache-2.0 is applied and the official text ships as `LICENSE` (see §4).

## 4. Licence applied, and the one attribution field still unverified

`LICENSE` is the complete official Apache License 2.0 text, fetched from
`https://www.apache.org/licenses/LICENSE-2.0.txt` on 2026-10-22: **11,358 bytes, 202 lines, sections 1
through 9 plus the appendix, sha256
`cfc7749b96f63bd31c3c42b5c471bf756814053e847c10f3eb003417bc523d30`**, left byte-for-byte unedited — the
bracketed boilerplate in its appendix is intentionally unfilled. Nine code and test files carry the
two-line SPDX header (`SPDX-FileCopyrightText: Copyright 2026 Zhizheng Ju`,
`SPDX-License-Identifier: Apache-2.0`), and the two JSON fixtures carry the same two identifiers as
leading object fields, since JSON has no comment syntax — so every shipped file is self-identifying
rather than relying on prose.

Adding those headers changed no behaviour: the three functional modules remain AST-identical to the
internal siblings with docstrings stripped, re-checked after the edit, and the suite still passes 27/27.

**Attribution.** The copyright line reads `Copyright 2026 Zhizheng Ju`, confirmed by the project owner
on 2026-10-22 along with the authorisation to license this original code under Apache-2.0. Until that
confirmation the header carried a placeholder, and no name was invented in its absence. The grant covers
this campaign's own material; it is not a grant over anyone else's, and the audit found no third-party
material for one to be needed.

**No `NOTICE` file is required**, because no third-party licensed material is included (see the import
audit above). No third-party copyright notices exist to preserve.

## 5. Rights position

Everything shipped here was authored by this campaign and describes elementary mathematics. On
2026-10-22 the owner authorised releasing this independently owned verification code under Apache-2.0,
which is now applied to the whole archive as delivered -- code, tests, constructed fixtures and the two
documents. The **data** exclusion is not
conceded necessity — it is unresolved provenance: the benchmark coordinate sets are held with published
attribution that this project has never evidenced locally, and their upstream licence terms are unread.
Adding that benchmark to a public release is therefore blocked until either the attribution is resolved
or the benchmark is dropped permanently. If it is ever included, it needs its own licence decision
(CC0-1.0 is the usual choice for dedicating coordinate tables) and its own provenance file.
