"""Ma-Zheng dihedral-angle libraries (block-pasting Stage 4).

Loads the precomputed spherical / Euclidean / Lanner dihedral-angle libraries that
ship with the Ma-Zheng HCPdm repository (cloned under ``scratchpad/HCPdm``). These are
the constraint sets that make per-combinatorial-type enumeration tractable: instead of
brute-forcing every integer label assignment, the block-paste enumerator (see
``pipeline/stage4_blockpaste.py``) prunes any partial Gram matrix whose sub-angle tuples
violate these libraries.

Encoding facts (verified against the files):
  * Each line of a ``*lis.txt`` file is a tuple of integer dihedral labels ``m`` (the
    Gram entry is ``-cos(pi/m)``).  ``S{k}lis`` lists the angle tuples of the ``C(k,2)``
    pairs among ``k`` facets that form a *spherical* (finite Coxeter, i.e. elliptic)
    rank-``k`` configuration; ``E{k}lis`` the *Euclidean* (parabolic) ones; ``L4lis`` the
    Lanner tetrahedra (compact hyperbolic 4-simplices = size-4 missing faces).
  * The files are PERMUTATION-CLOSED (every ordering of each diagram appears), so
    membership is a plain ``tuple in frozenset`` lookup -- no per-row sorting needed,
    provided the query tuple is built in the same column order (``itertools.combinations``
    of the sorted facet subset; see :func:`facet_pair_columns`).
  * The alphabet is ``{2,3,4,5,6,7}`` where ``7`` is a WILDCARD for "label >= 7".  The
    join therefore operates entirely in this alphabet; concrete labels >= 7 are only
    introduced at the solver hand-off (expand 7 -> {7,8,9,10,12}) and resolved by the
    exact signature solve.  :func:`collapse7` maps a concrete tuple back into the wildcard
    alphabet for membership tests / cross-checks.

Facet indices here are 0-based (matching the rest of the pipeline).  The HCPdm files are
pure label values (index-free), so there is no 0-vs-1 conversion at the value level.
"""
from __future__ import annotations

import itertools
from pathlib import Path

# Repo-root-relative location of the cloned Ma-Zheng libraries.
_TOOL = Path(__file__).resolve().parents[2] / "scratchpad" / "HCPdm" / "ToolPolytope"

# Expected sizes (rows incl. all orderings) -- sanity-checked on load.
_EXPECTED = {
    ("Slis", 3): 31, ("Slis", 4): 242, ("Slis", 5): 1946,
    ("Slis", 6): 20206, ("Slis", 7): 227676,
    ("Elis", 3): 10, ("Elis", 4): 27, ("Elis", 5): 257,
    ("Elis", 6): 870, ("Elis", 7): 6870,
    ("Lanner", 4): 392, ("Lanner", 5): 420,
}

_cache: dict[tuple[str, int], frozenset[tuple[int, ...]]] = {}


def _path(kind: str, k: int) -> Path:
    prefix = {"Slis": "S", "Elis": "E", "Lanner": "L"}[kind]
    return _TOOL / kind / f"{prefix}{k}lis.txt"


def load_library(kind: str, k: int) -> frozenset[tuple[int, ...]]:
    """Return the library ``kind`` ('Slis'|'Elis'|'Lanner') of rank ``k`` as a frozenset
    of integer tuples.  Cached.  Validates the row count against the known size."""
    key = (kind, k)
    if key in _cache:
        return _cache[key]
    path = _path(kind, k)
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rows.append(tuple(int(x) for x in line.split()))
    lib = frozenset(rows)
    expected = _EXPECTED.get(key)
    if expected is not None and len(rows) != expected:
        raise ValueError(f"{path}: expected {expected} rows, got {len(rows)}")
    # every tuple must have length C(k,2)
    width = k * (k - 1) // 2 if kind != "Slis" or k != 3 else 3
    width = k * (k - 1) // 2
    bad = next((r for r in lib if len(r) != width), None)
    if bad is not None:
        raise ValueError(f"{path}: row {bad} has width {len(bad)}, expected {width}")
    _cache[key] = lib
    return lib


# Convenience accessors -----------------------------------------------------

def S(k: int) -> frozenset[tuple[int, ...]]:
    return load_library("Slis", k)


def E(k: int) -> frozenset[tuple[int, ...]]:
    return load_library("Elis", k)


def L(k: int) -> frozenset[tuple[int, ...]]:
    """Lannér library of rank k (compact hyperbolic k-simplices = size-k missing faces).
    Defined for k = 4, 5 (Lannér diagrams exist only up to order 5)."""
    return load_library("Lanner", k)


def L4() -> frozenset[tuple[int, ...]]:
    return load_library("Lanner", 4)


def collapse7(tup) -> tuple[int, ...]:
    """Map a concrete label tuple into the wildcard alphabet: every label >= 7 -> 7.

    Used for membership tests / cross-checks when a tuple may carry concrete labels
    (8,9,10,12).  The block-paste join itself stays in the {2..7} alphabet and does not
    need this."""
    return tuple(7 if m >= 7 else int(m) for m in tup)


def facet_pair_columns(subset) -> list[tuple[int, int]]:
    """The ordered list of facet pairs (i, j), i < j, for a facet ``subset``, in
    ``itertools.combinations`` order of the SORTED subset.  This fixes the column order of
    a sub-config's angle tuple, matching the HCPdm convention (their pair key 10*i+j is
    enumerated in the same combinations order).

    >>> facet_pair_columns([0, 1, 3, 6])
    [(0, 1), (0, 3), (0, 6), (1, 3), (1, 6), (3, 6)]
    """
    s = sorted(subset)
    return [(a, b) for a, b in itertools.combinations(s, 2)]


if __name__ == "__main__":
    # quick self-test
    assert facet_pair_columns([0, 1, 3, 6]) == [(0, 1), (0, 3), (0, 6), (1, 3), (1, 6), (3, 6)]
    assert collapse7((2, 8, 12, 3, 7, 5)) == (2, 7, 7, 3, 7, 5)
    for (kind, k) in _EXPECTED:
        lib = load_library(kind, k)
        print(f"{kind}{k}: {len(lib)} rows, width {len(next(iter(lib)))}")
    print("mazheng_lib self-test OK")
