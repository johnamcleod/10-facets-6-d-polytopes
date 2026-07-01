"""CoxIter adapter: verify a solved Coxeter diagram is a COMPACT polytope of the right
dimension.  CoxIter (github.com/rgugliel/CoxIter, R. Guglielmetti) is the standard verifier
(CLAUDE.md 6); signature+isolation is only a proxy.  Built binary lives under
scratchpad/CoxIter/build/coxiter and reads the graph from STDIN (it uses stdin whenever
stdin is not a tty).

Graph format:
    n d                # facets, dimension
    i j m              # facets i,j (1-indexed) meet at angle pi/m  (m>=3; m=2 omitted,
                       #   since CoxIter defaults unlisted pairs to perpendicular)
    i j 1              # facets i,j are ultraparallel / non-intersecting (dotted)
"""
from __future__ import annotations

import subprocess
from pathlib import Path

_BIN = Path(__file__).resolve().parents[2] / "scratchpad" / "CoxIter" / "build" / "coxiter"


def diagram_text(label_assignment, dotted_pairs, n, d) -> str:
    """Build the CoxIter graph text from a solved config.  ``label_assignment`` maps
    ordinary pairs (i,j) (0-indexed) -> integer label m; ``dotted_pairs`` is the list of
    ultraparallel pairs.  m=2 (perpendicular) edges are omitted per Coxeter convention."""
    lines = [f"{n} {d}"]
    for (i, j), m in sorted(label_assignment.items()):
        if int(m) >= 3:
            lines.append(f"{i + 1} {j + 1} {int(m)}")
    for (i, j) in sorted(tuple(sorted(p)) for p in dotted_pairs):
        lines.append(f"{i + 1} {j + 1} 1")   # weight 1 = dotted (ultraparallel)
    return "\n".join(lines) + "\n"


def check(label_assignment, dotted_pairs, n, d, timeout=60):
    """Return (cocompact: bool, dimension: int|None) for the solved diagram, via CoxIter."""
    txt = diagram_text(label_assignment, dotted_pairs, n, d)
    out = subprocess.run([str(_BIN), "-c"], input=txt, capture_output=True,
                         text=True, timeout=timeout).stdout
    cocompact, dim = False, None
    for line in out.splitlines():
        s = line.strip()
        if s.startswith("Cocompact:"):
            cocompact = s.split(":", 1)[1].strip().lower().startswith("y")
        elif s.startswith("Dimension:"):
            try:
                dim = int(s.split(":", 1)[1].strip())
            except ValueError:
                dim = None
    return cocompact, dim


def is_compact_polytope(label_assignment, dotted_pairs, n, d, timeout=60) -> bool:
    """True iff CoxIter certifies a COMPACT polytope of dimension d."""
    cocompact, dim = check(label_assignment, dotted_pairs, n, d, timeout=timeout)
    return cocompact and dim == d


def available() -> bool:
    return _BIN.exists()
