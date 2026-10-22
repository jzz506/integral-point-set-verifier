# e213ips — exact verification and similarity classification of integral point sets

**Public code release `e213ips` version 0.1.0**, prepared 2026-10-22. Standard library only; declared
Python 3.9 or later. **Built and tested on CPython 3.9.6 only** — 3.10 through 3.13 have not been
executed here, so CI on those versions is the first thing that should confirm the "or later".

**Licence: Apache License 2.0.** The complete official text is in [`LICENSE`](LICENSE); every code,
test and example file carries an SPDX header. See *Licence and attribution* below for what is covered,
whose it is, and what a distributor must not assume.

This candidate contains **code, tests and constructed controls only**. It ships no coordinate
benchmark, no research dataset and no path into any research repository: every point set in here is
either derived from an elementary formula in the test file or is a hand-solved small integer
configuration, with the arithmetic written down next to it so you can check the fixture yourself.

## What this is, and what it is not

An **integral point set** is a set of distinct plane points whose pairwise distances are all positive
integers. Erdős Problem #213 asks whether eight such points can exist with **no three collinear and no
four concyclic**. **Nobody has found such a set, and this package does not contain one.** It ships:

1. `e213ips.verify` — a verifier deciding, from the definitions and in exact arithmetic, whether a
   point set qualifies: integral distances, no collinear triple, no concyclic quadruple;
2. `e213ips.similarity` — a classifier deciding when two such sets are the *same object up to scaling*,
   which is the unit the literature counts in;
3. `e213ips.cli` — a command line around both.

Nothing here enumerates, generates or searches. It decides and classifies. The problem stays open.

## Quick start

```sh
shasum -a 256 -c SHA256SUMS.txt                              # manifest verifies
python3 tests/run_all.py                                     # 27 legs, <1 s, exit 0 only if all pass
python3 -m e213ips.cli examples/hub_4set_legal.json           # exit 0
python3 -m e213ips.cli examples/trapezoid_s2_illegal.json     # exit 1, names the concyclic quad
```

```python
from fractions import Fraction
from e213ips import verify as V, similarity as S

# a 40-40-48 triangle (a 5-5-6 triangle scaled by 8) has circumradius abc/4K = 25,
# so its three vertices plus its circumcentre are a legal integral 4-set.
pts = [(0, 0), (48, 0), (24, 32), (Fraction(24), Fraction(7))]
rep = V.verify(pts, s=1)
rep.is_legal                      # True
sorted(rep.distances.values())    # [25, 25, 25, 40, 40, 48]
S.class_id(pts, 1)                # 'H4-db6ef773afadaba6' -- invariant under relabeling and scaling
```

## The coordinate model

A stored pair `(x, t)` means the true point **`(x, t·√s)`**, with `x, t ∈ Q` and one squarefree `s ≥ 1`
shared by the whole set (its **characteristic**). Squared distance is `Δx² + s·Δt²`.

* **Collinearity** transfers directly: `(x,t) → (x, t√s)` is an invertible linear map, so the cross
  product test runs on the stored pairs.
* **Concyclicity** does not. Circles pull back to `x² + s·t² + Dx + Et + F = 0`, so the predicate is the
  weighted determinant `det[[x²+s·t², x, t, 1]] = 0`. Using the unweighted circle equation on stored
  pairs is wrong whenever `s > 1`. `tests/test_verify.py` shows the effect in **both** directions on
  hand-checked coordinates: one quadruple is cyclic at `s=1` and not at `s=2`, another the other way.
* `s` is **load-bearing**: the same coordinate rows give different verdicts at different `s`, and the
  verifier refuses to pretend otherwise.

The characteristic is also a property of the *distances*: all of a set's triangles have `16K²` in one
square class, and that class equals the characteristic (Kemnitz's shared-characteristic theorem).
`heron_16k2` plus `same_square_class` recompute it from distances alone — as a pair test, because
getting a squarefree *representative* of a 40-digit Heron value would mean factoring it.

## API

`verify` — `dist2`, `integral_distance`, `collinear`, `concyclic`, `det4`, `det_bareiss`,
`det_laplace`, `heron_16k2`, `is_square_rational`, `same_square_class`, `verify(points, s=1)` →
`Report(n, s, total_pairs, integral_distance_count, distances, non_integral_pairs, collinear_triples,
concyclic_quads, reasons, is_legal)`.

`similarity` — `squared_distances`, `clear_denominators`, `distance_matrix`, `similarity_signature`,
`class_id`, `distance_gcd_key`, `are_similar`, `linear_scale`.

Coordinates may be `int` or `Fraction`; **floats are refused with a `TypeError`**. No decision path in
this package touches a float and no geometric predicate is approximated. Determinism: the classifier is
exhaustive over `n!` labelings with a `min`, so ids do not depend on input order, dict ordering or hash
seed. Cost is `n!` labelings — 24 at `n=4`, 5,040 at `n=7`, 40,320 at `n=8`.

## Tests: 27 legs, every detector answering both ways

| group | legs | what it can indict |
|---|---|---|
| `tests/test_verify.py` | 14 | three integrality/concyclicity/collinearity negatives, `√2`, a +1 vertex mutation, floats, a wrong-characteristic read in both directions, and a 1,278-quadruple agreement between two independent concyclicity implementations |
| `tests/test_similarity.py` | 9 | relabeling and dilation invariance, a non-square scale refused on 18 edges, distinctness, the `n=8` path, exact denominator clearing, and the multiset key's incompleteness measured on a lattice |
| `tests/test_normalisation.py` | 4 | raw-coordinate gcd is not translation invariant while `diff_gcd` is; the two primitives coincide exactly when `edge_gcd == diff_gcd`; a 7× frame divided back is the same object |

Controls constructed in the suite: a **hub 4-set** derived from the circumradius formula `R = abc/4K`
of a 40-40-48 triangle (legal, distances 25,25,25,40,40,48); a **characteristic-2 triangle** `(0,0),(2,0),(1,2)`
with sides 2,3,3; and a **characteristic-2 isosceles trapezoid** `(0,0),(9,0),(7,4),(2,4)` with six
integral distances 5,6,6,9,9,9, no three collinear and **exactly one concyclic quadruple** — the
mechanism the open problem keeps hitting, reproduced from three lines of arithmetic instead of a
dataset.

## Coverage reduction relative to the internal candidate

The internal candidate (15 files, 31 legs) verifies against a 36-class seven-point benchmark and two
eight-point near-miss classes. Those datasets are excluded here, so the public suite is **weaker in
seven specific ways**, listed rather than glossed:

1. **No external known-answer for the classifier.** The internal suite reproduces 36 `class_id` strings
   written by a different program in an earlier round. Here the four pinned ids are regression constants
   for this package only, not a cross-instrument check. (One coincidence is worth noting: the
   constructed hub lands on `H4-db6ef773afadaba6`, the same id the internal example of the same
   classical configuration carries.)
2. **No seven-point set is verified.** All controls here have `n ≤ 4` except the synthetic `n=8`
   lattice set, which is not integral. The 21-distance / 35-triple / 35-quad geometry is unexercised.
3. **No real eight-point near-miss.** Internally, two classes with 28 integral distances and exactly 3
   concyclic quads of 70 are re-decided. Here the same *mechanism* is reproduced at `n=4` with 1 quad of
   1, which is a far weaker test of the quad loop.
4. **The scale-normalisation divergence is not reproducible.** Internally `edge_gcd` and `diff_gcd`
   disagree for 32 of 36 classes; every constructed control here has both equal to 1, so the public
   suite tests the relation and the machinery, not the phenomenon.
5. **No reproduction of a published constant.** Internally the lattice normalisation of one class yields
   exactly `2,262,000`, the smallest primitive 7-cluster diameter printed in arXiv:1312.2318. That
   positive control cannot be rebuilt without the data.
6. **The square-class check is thin.** Internally 35 Heron values × 36 classes; here 4 triangles at
   `s=1` and 1 at `s=2`.
7. **No reconciliation against a stored key convention.** Internally the distance key, squared
   elementwise, equals a stored `campaign_squared_gcd_key` on 36/36. Here only the identity
   `gcd(D²) = gcd(D)²` is checked, on four controls.

The incompleteness of the sorted-distance-multiset key, which the internal suite measured on 1,260 real
four-subsets, is measured here **more strongly** on independently constructed data: among the 1,820
four-subsets of a 4×4 lattice the multiset key merges 15 groups that the complete signature separates.

## Limitations

* Eight points have not been found. Nothing in this package is a non-existence claim, and the shipped
  near-miss fixture is a 4-point toy, not evidence about `n=8`.
* "Legal" means the definitions hold for the coordinates given. It says nothing about whether a set is
  new, published, or whose it is.
* Only characteristics 1 and 2 are exercised by the tests. The code is written for any squarefree `s`;
  that is a claim about the code, not a measurement.
* `verify` is `O(n²)` determinants and `class_id` is `O(n!·n²)`; for `n ≥ 9` the classifier will be slow
  by design.

## Licence and attribution

**Apache License 2.0**, applied to this archive as delivered: the four modules in `e213ips/`, the five
files in `tests/`, the two fixtures in `examples/`, and the two documents (`README.md`,
`PROVENANCE.md`). The owner of the project authorised the Apache-2.0 release of this campaign's
independently authored verification code on 2026-10-22; that authorisation covers material this
campaign wrote, and nothing else.

Copyright attribution currently reads:

```
SPDX-FileCopyrightText: Copyright 2026 Zhizheng Ju
SPDX-License-Identifier: Apache-2.0
```

This attribution was confirmed by the project owner on 2026-10-22, together with the authorisation to
license the campaign's independently authored original verification code under Apache-2.0. The grant
covers that material and nothing else: no third-party code is included, so there is nothing here whose
ownership would need to be asserted. `LICENSE` is the verbatim official text — 11,358 bytes, 202 lines,
sections 1-9 plus the appendix, sha256
`cfc7749b96f63bd31c3c42b5c471bf756814053e847c10f3eb003417bc523d30` — so its bracketed boilerplate is
deliberately left unfilled.

**Third-party material: none, so no `NOTICE` file is required.** The audit behind that sentence is in
`PROVENANCE.md` §1 and is re-runnable: the package imports only standard-library modules (`__future__`,
`fractions`, `math`, `itertools`, `hashlib`, `json`, `os`, `sys`), bundles no source, data file or binary obtained
from elsewhere, and reproduces no third party's expression. Using CPython's standard library at runtime
is not redistribution of third-party code, so no Python licence text is embedded; a distributor who
packages this into a larger artifact may still want to note CPython itself.

The published literature is **cited, never copied**: arXiv:0804.1303, arXiv:1312.2318 and
Kurz-Wassermann appear as references, and the one constant quoted in the README (`2,262,000`, the
smallest primitive 7-cluster diameter printed in arXiv:1312.2318) is a numerical fact attributed to its
source, not reproduced content. Mathematical facts, definitions and coordinates derived independently
from the same formulas are not copyrightable expression; where a fixture happens to coincide with a
published example, it was constructed here from the formula shown next to it in the test file.

If the owner later prefers a different grant for the non-code material (for example CC0-1.0 for the
fixtures, or CC-BY-4.0 for the documentation), that is a one-line change to the SPDX headers plus this
section — it has not been assumed here.

